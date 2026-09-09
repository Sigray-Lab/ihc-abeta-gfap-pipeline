#!/usr/bin/env python3
"""Fixed validation boxes for blinded human plaque tracing (plan §4; review 2026-09-09 F4/F5).

Picks stained sections from animals NOT used to train any classifier, spanning dim / mid /
bright Cy3 background in each arm, and places non-overlapping 500 x 500 um boxes inside the
existing Hippocampus and Isocortex outlines, away from the outline edge. Boxes are chosen
by a seeded RNG from a grid of admissible centres -- before anyone has looked at any
classifier mask -- and written as QuPath-importable GeoJSON (class RefBox) plus an index.

    python3 scripts/make_reference_boxes.py --out ~/ihc_work/reference_boxes
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from shapely.geometry import shape, box, mapping

REPO = Path(__file__).resolve().parent.parent
PX = 0.325                       # um per native pixel
BOX_UM, EDGE_UM = 500.0, 100.0   # box side, minimum distance from the outline edge
PER_REGION = {"Hippocampus": 3, "Isocortex": 4}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=Path.home() / "ihc_work/reference_boxes")
    ap.add_argument("--rois", type=Path, default=Path.home() / "ihc_work/region_rois")
    ap.add_argument("--seed", type=int, default=20260909)
    ap.add_argument("--per-arm", type=int, default=3)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    key = REPO.parent / "BlindingKey"
    man = pd.read_csv(key / "provenance_manifest.csv"); man = man[man.row_kind == "section"].copy()
    def _img(r):
        pl = str(r.physical_section_label); s = pl.zfill(2) if pl.isdigit() else pl
        return f"{r.code}{'b' if r.scan == 'rescan' else ''}_s{s}"
    man["image"] = man.apply(_img, axis=1)
    man = man[(man.use_for_measurement.astype(str) == "True") & (man.condition == "positive")]
    # animals used for ANY classifier training are excluded from validation
    trained = set()
    for l in (key / "training_annotation_map.csv").read_text().splitlines():
        if l.startswith("#") or "," not in l: continue
        trained.add(l.split(",")[1].strip()[:3])
    brt = (pd.read_csv(REPO.parent / "DerivedData/aris_v2_2026-08-18/three_arm_section_level.csv")
             .drop_duplicates("image").set_index("image").cy3_med)
    have = {img for img in man.image
            if (a.rois / "Hippocampus" / f"{img}.geojson").exists() and (a.rois / "Isocortex" / f"{img}.geojson").exists()}
    cand = man[man.image.isin(have) & ~man.code.isin(trained) & ~man.image.isin({"Y13_s01", "S17_s04"})].copy()
    cand["cy3"] = cand.image.map(brt); cand = cand.dropna(subset=["cy3"])

    picks = []
    for arm in ("diet", "ip"):
        c = cand[cand.arm == arm].sort_values("cy3")
        tert = [c.iloc[ix] for ix in np.array_split(np.arange(len(c)), a.per_arm)]   # dim / mid / bright
        used = set()
        for t in tert:
            t = t[~t.code.isin(used)]
            if not len(t): continue
            r = t.iloc[len(t) // 2]                   # the median section of the tertile
            picks.append(r); used.add(r.code)
    picks = pd.DataFrame(picks)

    a.out.mkdir(parents=True, exist_ok=True)
    index = []
    side, edge = BOX_UM / PX, EDGE_UM / PX
    for _, r in picks.iterrows():
        feats = []
        for region, n in PER_REGION.items():
            g = shape(json.loads((a.rois / region / f"{r.image}.geojson").read_text())["features"][0]["geometry"]).buffer(0)
            inner = g.buffer(-edge)                  # every box must sit >= EDGE_UM inside the outline
            chosen, side_used = [], None
            for side_um in (BOX_UM, 300.0):         # fall back to 300 um boxes in small regions
                sd = side_um / PX
                if inner.is_empty: break
                minx, miny, maxx, maxy = inner.bounds
                grid = [(x, y) for x in np.arange(minx, maxx - sd, sd / 4) for y in np.arange(miny, maxy - sd, sd / 4)
                        if inner.contains(box(x, y, x + sd, y + sd))]
                rng.shuffle(grid)
                chosen = []
                for (x, y) in grid:
                    b = box(x, y, x + sd, y + sd)
                    if all(b.intersection(o).is_empty for o in chosen):
                        chosen.append(b)
                    if len(chosen) == n: break
                side_used = side_um
                if len(chosen) >= 2: break
            side = (side_used or BOX_UM) / PX
            for i, b in enumerate(chosen, 1):
                bid = f"{region[:3].upper()}{i}"
                feats.append({"type": "Feature", "id": f"{r.image}_{bid}", "geometry": mapping(b),
                              "properties": {"objectType": "annotation", "name": bid,
                                             "classification": {"name": "RefBox", "colorRGB": -256}}})
                bx = b.bounds
                index.append(dict(image=r.image, arm=r.arm, region=region, box=bid,
                                  x_px=round(bx[0]), y_px=round(bx[1]), w_px=round(side), h_px=round(side),
                                  cy3_background=float(r.cy3)))
        (a.out / f"{r.image}.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
        print(f"  {r.image}  arm={r.arm}  cy3 background {r.cy3:.0f}  boxes {len(feats)}")
    pd.DataFrame(index).to_csv(a.out / "reference_boxes_index.csv", index=False)
    print(f"wrote {len(picks)} sections, {len(index)} boxes -> {a.out}")


if __name__ == "__main__":
    main()
