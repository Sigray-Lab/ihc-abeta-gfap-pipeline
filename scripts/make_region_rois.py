#!/usr/bin/env python3
"""Turn the delivered hippocampus/isocortex outlines into measurement ROIs.

The 2026-08-19 delivery needs five repairs before it can be measured against, all of them
mechanical and none needing the annotator:

1. Region identity lives in the GeoJSON ``name`` property, not ``classification``. QuPath
   displays it either way; our importer reads ``classification``, so it would have seen 246
   unclassified objects.
2. Spelling varies: TISSUE, TISSSUE, Cortex, Isocortes.
3. Eight objects carry no name. Six are 18-25 mm2 in files that already have a hippocampus
   and an isocortex, so they are the tissue outline by elimination. Two are small and
   genuinely ambiguous and are left out rather than guessed.
4. Two hippocampi are degenerate (0.305 and 0.000 mm2).
5. Outlines were drawn on the ORIGINAL scan for five animals whose measurement uses the
   rescan. Those are reported and skipped here; recovering them needs registration between
   the two acquisitions (they are the same tissue -- mask Dice 0.989-0.996).

Each region is intersected with the DAPI tissue mask so numerator and denominator share a
support, exactly as ADR-0029 requires for whole tissue.

    python3 scripts/make_region_rois.py --src <unzipped dir> --out <roi dir>
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from shapely.geometry import shape, mapping, MultiPolygon, Polygon
from shapely.ops import unary_union

REPO = Path(__file__).resolve().parent.parent
CANON = {"tissue": "Tissue", "tisssue": "Tissue", "hippocampus": "Hippocampus",
         "isocortex": "Isocortex", "isocortes": "Isocortex", "cortex": "Isocortex"}
MIN_REGION_MM2 = 0.5
REGION_COLOR = {"Hippocampus": -65281, "Isocortex": -16776961}


def load_regions(path: Path, px_um: float) -> dict:
    d = json.loads(path.read_text())
    feats = d.get("features", d) if isinstance(d, dict) else d
    named, unnamed = {}, []
    for ft in feats:
        raw = (ft.get("properties", {}).get("name") or "").strip()
        g = shape(ft["geometry"])
        if not g.is_valid:
            g = g.buffer(0)
        key = CANON.get(raw.lower())
        if key:
            named[key] = unary_union([named[key], g]) if key in named else g
        else:
            unnamed.append(g)
    # an unnamed object in a file that already has both cortical regions is the tissue outline
    if unnamed and "Tissue" not in named and {"Hippocampus", "Isocortex"} <= set(named):
        big = max(unnamed, key=lambda g: g.area)
        if big.area * px_um ** 2 / 1e6 > 10:
            named["Tissue"] = big
    return named


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--src", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path.home() / "ihc_work/region_rois")
    ap.add_argument("--tissue-rois", type=Path, default=Path.home() / "ihc_work/tissue_rois")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    man = pd.read_csv(REPO.parent / "BlindingKey" / "provenance_manifest.csv")
    man = man[man.row_kind == "section"].copy()
    man["stem"] = ("Image_" + man.tube_id.astype(str)
                   + np.where(man.scan == "rescan", "r", "")
                   + "_s" + man.section_label.astype(str).str.zfill(2))
    lab = np.where(man.scan == "rescan", man.physical_section_label,
                   man.section_label.astype(str).str.zfill(2))
    man["image"] = (man.code + np.where(man.scan == "rescan", "b", "")
                    + "_s" + pd.Series(lab).astype(str).values)
    by_stem = man.set_index("stem")

    rows, skipped = [], []
    for f in sorted(a.src.rglob("*.geojson")):
        if "MACOSX" in str(f):
            continue
        stem = f.name.replace("_annotations.geojson", "")
        if stem not in by_stem.index:
            skipped.append((stem, "no manifest row")); continue
        m = by_stem.loc[stem]
        if isinstance(m, pd.DataFrame):
            m = m.iloc[0]
        px = float(m.pixel_size_um)
        regions = load_regions(f, px)

        tis_path = a.tissue_rois / f"{m.image}.geojson"
        if not tis_path.exists():
            skipped.append((stem, f"{m.image} not in the measurement set")); continue
        tj = json.loads(tis_path.read_text())
        tissue = shape(tj["features"][0]["geometry"])
        if not tissue.is_valid:
            tissue = tissue.buffer(0)

        for name in ("Hippocampus", "Isocortex"):
            g = regions.get(name)
            if g is None:
                continue
            mm2 = g.area * px ** 2 / 1e6
            if mm2 < MIN_REGION_MM2:
                skipped.append((stem, f"{name} degenerate at {mm2:.3f} mm2")); continue
            inter = g.intersection(tissue)
            # An intersection can return a GeometryCollection -- polygons plus zero-width
            # line or point contacts where the two boundaries touch. QuPath's GeoJSON
            # reader rejects GeometryCollection outright ("Json object does not contain
            # coordinates"), so keep the polygonal parts and drop the degenerate ones.
            # They carry no area, so nothing measurable is lost.
            if inter.geom_type == "GeometryCollection":
                polys = [q for q in inter.geoms if q.geom_type in ("Polygon", "MultiPolygon")]
                inter = MultiPolygon([q for p_ in polys
                                      for q in (p_.geoms if p_.geom_type == "MultiPolygon" else [p_])]) \
                        if polys else Polygon()
            if inter.is_empty:
                skipped.append((stem, f"{name} does not meet the tissue mask")); continue
            sub = a.out / name
            sub.mkdir(parents=True, exist_ok=True)
            out = sub / f"{m.image}.geojson"
            out.write_text(json.dumps({"type": "FeatureCollection", "features": [{
                "type": "Feature", "id": f"{m.image}_{name}",
                "geometry": mapping(inter),
                # QuPath's PathClass deserialiser requires colorRGB; without it
                # PathIO.readObjects throws "Unable to parse PathClass". Match the format
                # make_tissue_rois.py already writes.
                "properties": {"objectType": "annotation",
                               "classification": {"name": name,
                                                  "colorRGB": REGION_COLOR[name]}}}]}))
            rows.append(dict(image=m.image, stem=stem, region=name,
                             drawn_mm2=mm2, roi_mm2=inter.area * px ** 2 / 1e6,
                             kept_frac=inter.area / g.area))

    d = pd.DataFrame(rows)
    d.to_csv(a.out / "region_roi_index.csv", index=False)
    print(f"  layout: {a.out}/<Region>/<image>.geojson  (matches eval_in_tissue.groovy)")
    print(f"  wrote {len(d)} region ROIs to {a.out}")
    for r, g in d.groupby("region"):
        print(f"    {r:<12} n={len(g):>3}  median {g.roi_mm2.median():.2f} mm2  "
              f"kept {g.kept_frac.median():.1%} of the drawn outline")
    print(f"  sections covered: {d.image.nunique()}")
    if skipped:
        print(f"\n  skipped {len(skipped)}:")
        from collections import Counter
        for reason, n in Counter(r for _, r in skipped).most_common():
            print(f"    {n:>3}  {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
