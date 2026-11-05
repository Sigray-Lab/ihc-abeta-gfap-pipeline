#!/usr/bin/env python3
"""Per-image, per-region normalisation anchors for the REGIONNORM arm.

Writes a CSV with columns image, region, glass_p50, r_p50 -- the input that
make_region_classifiers.py needs (docs/normalisation.md, option 2, regional variant):
glass_p50 is the section's glass level from the section-anchor CSV; r_p50 is the median
Cy3 inside the region outline AND inside the DAPI tissue mask, read at a pyramid level of
about 2.6 um/px. The file that produced the August anchors lived only in /tmp
(/tmp/region_anchors.csv) and is gone; this script replaces it in the repository.

    python3 scripts/make_region_anchors.py --out /path/region_anchors.csv
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
from shapely.geometry import shape
from skimage.draw import polygon as draw_polygon

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from ihc.qc.stain_check import read_pyramid_plane, CH_DAPI, CH_CY3   # noqa: E402
from ihc.qc.tissue_mask import tissue_mask                            # noqa: E402
import json

LEVEL = 3   # 2^3 = 8x downsample, ~2.6 um/px


def _img_id(r):
    pl = str(r.physical_section_label)
    return f"{r.code}{'b' if r.scan == 'rescan' else ''}_s{pl.zfill(2) if pl.isdigit() else pl}"


def rasterise(geom, shape_hw, scale):
    m = np.zeros(shape_hw, bool)
    polys = [geom] if geom.geom_type == "Polygon" else list(geom.geoms)
    for pg in polys:
        x, y = np.asarray(pg.exterior.coords).T
        rr, cc = draw_polygon(y / scale, x / scale, shape=shape_hw); m[rr, cc] = True
        for hole in pg.interiors:
            x, y = np.asarray(hole.coords).T
            rr, cc = draw_polygon(y / scale, x / scale, shape=shape_hw); m[rr, cc] = False
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rois", type=Path, default=Path.home() / "ihc_work/region_rois")
    ap.add_argument("--section-anchors", type=Path,
                    default=REPO.parent / "DerivedData/aris_v2_2026-08-18/global_norm_anchors.csv")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    man = pd.read_csv(REPO.parent / "BlindingKey" / "provenance_manifest.csv")
    man = man[man.row_kind == "section"].copy(); man["image"] = man.apply(_img_id, axis=1)
    man = man.drop_duplicates("image").set_index("image")
    glass = pd.read_csv(a.section_anchors).set_index("image").glass_p50
    rows = []
    for region in ("Hippocampus", "Isocortex"):
        for f in sorted((a.rois / region).glob("*.geojson")):
            img = f.stem
            if img not in man.index:
                print(f"  SKIP {img}: not in provenance manifest"); continue
            r = man.loc[img]
            g = shape(json.loads(f.read_text())["features"][0]["geometry"]).buffer(0)
            kw = dict(pixel_size_um=float(r.pixel_size_um), true_size=(int(r.width_px), int(r.height_px)))
            dapi, lv = read_pyramid_plane(r.ets_path, CH_DAPI, level=LEVEL, **kw)
            cy3, _ = read_pyramid_plane(r.ets_path, CH_CY3, level=LEVEL, **kw)
            scale = 2 ** lv
            m = rasterise(g, dapi.shape, scale) & tissue_mask(dapi) & np.isfinite(cy3)
            if m.sum() < 500:
                print(f"  SKIP {img}/{region}: {int(m.sum())} px inside region+tissue"); continue
            rows.append(dict(image=img, region=region, glass_p50=float(glass.get(img, np.nan)),
                             r_p50=float(np.median(cy3[m])), n_px=int(m.sum()), level=lv))
            print(f"  {img:10s} {region:11s} r_p50={rows[-1]['r_p50']:6.1f} glass={rows[-1]['glass_p50']:5.1f} n={m.sum()}", flush=True)
    d = pd.DataFrame(rows); a.out.parent.mkdir(parents=True, exist_ok=True); d.to_csv(a.out, index=False)
    print(f"wrote {len(d)} anchors -> {a.out}")


if __name__ == "__main__":
    main()
