#!/usr/bin/env python3
"""Machine-drawn Negative annotations on secondary-only controls, and overlays to inspect them.

A secondary-only control carries no primary antibody, so every tissue pixel in it is a true
negative by construction. That makes two kinds of Negative annotation legitimate to draw
without a human: (1) HARD negatives -- the objects the current classifier wrongly calls Abeta
on the control, taken as polygons; (2) EASY negatives -- short strokes at seeded random
positions inside the tissue, away from the edge and from the false positives. Positives can
never be drawn this way. Over-labelling bright granules can cost sensitivity to small real
plaques, so the count per section is capped and the effect must be checked on the reference
boxes after retraining.

Inputs: <dir>/hosts.csv (image, cy3_med), <dir>/fp/<image>.geojson (classifier Abeta objects
on the control), ~/ihc_work/tissue_rois/<image>.geojson. Outputs: <dir>/annotations/<image>_auto_annotations.geojson
(+ map.csv), <dir>/png/<image>_overview.png and <image>_zoom.png.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from shapely.geometry import shape, mapping, LineString, Point
from shapely.ops import unary_union
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPoly

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from ihc.qc.stain_check import read_pyramid_plane, CH_CY3     # noqa: E402

PX = 0.325
N_HARD, N_EASY, MIN_FP_UM2, EDGE_UM, STROKE_UM = 25, 12, 20.0, 200.0, 30.0
NEG_RGB = -16776961      # QuPath colour int for the Negative class (blue), as in the delivered files


IGN_RGB = -4144960       # QuPath colour int for Ignore* (grey)
EDGE_IGNORE_UM = 20.0    # a false positive this close to the tissue boundary is section-rim glow -> Ignore*


def feat(geom, name, cls="Negative"):
    return {"type": "Feature", "geometry": mapping(geom),
            "properties": {"objectType": "annotation", "name": name,
                           "classification": {"name": cls, "colorRGB": NEG_RGB if cls == "Negative" else IGN_RGB}}}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, default=Path.home() / "ihc_work/auto_negatives")
    ap.add_argument("--seed", type=int, default=20260909)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    (a.dir / "annotations").mkdir(exist_ok=True); (a.dir / "png").mkdir(exist_ok=True)
    man = pd.read_csv(REPO.parent / "BlindingKey" / "provenance_manifest.csv"); man = man[man.row_kind == "section"].copy()
    def _img(r):
        pl = str(r.physical_section_label); s = pl.zfill(2) if pl.isdigit() else pl
        return f"{r.code}{'b' if r.scan == 'rescan' else ''}_s{s}"
    man["image"] = man.apply(_img, axis=1); man = man.drop_duplicates("image").set_index("image")
    hosts = pd.read_csv(a.dir / "hosts.csv"); maplines = []; summary = []
    for _, h in hosts.iterrows():
        img = h.image
        tissue = shape(json.loads((Path.home() / "ihc_work/tissue_rois" / f"{img}.geojson").read_text())["features"][0]["geometry"]).buffer(0)
        raw = [shape(f["geometry"]).buffer(0) for f in json.loads((a.dir / "fp" / f"{img}.geojson").read_text())["features"]]
        fp = []                                   # QuPath exports one MultiPolygon per class: explode it
        for g in raw:
            fp.extend(list(g.geoms) if g.geom_type == "MultiPolygon" else [g])
        fp = [g for g in fp if g.area * PX * PX >= MIN_FP_UM2 and tissue.contains(g.representative_point())]
        fp.sort(key=lambda g: -g.area)
        hard = fp[:N_HARD]
        inner = tissue.buffer(-EDGE_UM / PX); fp_union = unary_union(fp).buffer(50 / PX) if fp else None
        minx, miny, maxx, maxy = inner.bounds; easy = []
        for _ in range(20000):
            if len(easy) == N_EASY: break
            p = Point(rng.uniform(minx, maxx), rng.uniform(miny, maxy))
            if not inner.contains(p) or (fp_union is not None and fp_union.contains(p)): continue
            ang = rng.uniform(0, np.pi); d = STROKE_UM / PX / 2
            ln = LineString([(p.x - d * np.cos(ang), p.y - d * np.sin(ang)), (p.x + d * np.cos(ang), p.y + d * np.sin(ang))])
            if inner.contains(ln): easy.append(ln)
        boundary = tissue.boundary
        edge = [g for g in hard if g.distance(boundary) * PX < EDGE_IGNORE_UM]      # rim glow -> Ignore*
        hard_in = [g for g in hard if g not in edge]                                # vessels, granules -> Negative
        feats = ([feat(g, f"hardneg{i+1}") for i, g in enumerate(hard_in)] + [feat(g, f"rim{i+1}", "Ignore*") for i, g in enumerate(edge)]
                 + [feat(l, f"easyneg{i+1}") for i, l in enumerate(easy)])
        (a.dir / "annotations" / f"{img}_auto_annotations.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
        maplines.append(f"{img}_auto,{img}")
        summary.append(dict(image=img, cy3_med=h.cy3_med, fp_objects=len(fp), hard_negative=len(hard_in), rim_ignore=len(edge), easy=len(easy),
                            hard_area_um2=round(sum(g.area for g in hard) * PX * PX), largest_fp_um2=round(fp[0].area * PX * PX) if fp else 0))
        # ---------- overlays ----------
        r = man.loc[img]; kw = dict(pixel_size_um=float(r.pixel_size_um), true_size=(int(r.width_px), int(r.height_px)))
        cy3, lv = read_pyramid_plane(r.ets_path, CH_CY3, level=2, **kw); s = 2 ** lv
        lo, hi = np.nanpercentile(cy3, [1, 99.7]); im = np.clip((np.nan_to_num(cy3, nan=lo) - lo) / (hi - lo), 0, 1)
        fig, ax = plt.subplots(figsize=(11, 11 * im.shape[0] / im.shape[1])); ax.imshow(im, cmap="gray", vmin=0, vmax=1); ax.set_axis_off()
        for tp in ([tissue] if tissue.geom_type == "Polygon" else tissue.geoms):
            tx, ty = tp.exterior.xy; ax.plot(np.array(tx) / s, np.array(ty) / s, color="white", lw=.6, alpha=.6)
        for g in fp:
            for poly in ([g] if g.geom_type == "Polygon" else g.geoms):
                x, y = poly.exterior.xy; ax.plot(np.array(x) / s, np.array(y) / s, color="red", lw=.7)
        for g in hard:
            col = "orange" if g in edge else "cyan"
            for poly in ([g] if g.geom_type == "Polygon" else g.geoms):
                ax.add_patch(MplPoly(np.c_[np.array(poly.exterior.xy[0]) / s, np.array(poly.exterior.xy[1]) / s], closed=True, facecolor=col, alpha=.45, edgecolor=col, lw=.8))
        for l in easy:
            x, y = l.xy; ax.plot(np.array(x) / s, np.array(y) / s, color="deepskyblue", lw=2.2)
        ax.set_title(f"{img}  (secondary-only control, Cy3 background {h.cy3_med:.0f})   red = corrected GLOBAL calls 'Abeta' ({len(fp)} objects)   "
                     f"cyan = {len(hard_in)} hard Negatives   orange = {len(edge)} rim objects -> Ignore*   blue = {len(easy)} easy Negative strokes", fontsize=8)
        fig.savefig(a.dir / "png" / f"{img}_overview.png", dpi=220, bbox_inches="tight"); plt.close(fig)
        # zooms on the three largest false positives, at 0.65 um/px
        cy1, lv1 = read_pyramid_plane(r.ets_path, CH_CY3, level=1, **kw); s1 = 2 ** lv1
        lo1, hi1 = np.nanpercentile(cy1, [1, 99.7])
        fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
        for ax, g in zip(axes, hard[:3] if len(hard) >= 3 else hard + [None] * (3 - len(hard))):
            ax.set_axis_off()
            if g is None: continue
            cx, cy = g.centroid.x / s1, g.centroid.y / s1; w = 300 / PX / s1     # 600 um window
            y0, y1 = int(max(cy - w, 0)), int(min(cy + w, cy1.shape[0])); x0, x1 = int(max(cx - w, 0)), int(min(cx + w, cy1.shape[1]))
            crop = np.clip((np.nan_to_num(cy1[y0:y1, x0:x1], nan=lo1) - lo1) / (hi1 - lo1), 0, 1)
            ax.imshow(crop, cmap="gray", vmin=0, vmax=1, extent=(x0, x1, y1, y0))
            for gg in fp:
                if gg.centroid.x / s1 < x0 or gg.centroid.x / s1 > x1 or gg.centroid.y / s1 < y0 or gg.centroid.y / s1 > y1: continue
                for poly in ([gg] if gg.geom_type == "Polygon" else gg.geoms):
                    x, y = poly.exterior.xy; ax.plot(np.array(x) / s1, np.array(y) / s1, color="red", lw=1.0)
            col = "orange" if g in edge else "cyan"
            for poly in ([g] if g.geom_type == "Polygon" else g.geoms):
                x, y = poly.exterior.xy; ax.plot(np.array(x) / s1, np.array(y) / s1, color=col, lw=1.6)
            ax.set_title(f"{img}: false positive {round(g.area*PX*PX)} um2 -> {'Ignore* (rim)' if g in edge else 'Negative'}  (600 um window)", fontsize=9)
        fig.savefig(a.dir / "png" / f"{img}_zoom.png", dpi=200, bbox_inches="tight"); plt.close(fig)
        print(f"  {img}: {len(fp)} false-positive objects, {len(hard)} hard + {len(easy)} easy Negatives drawn", flush=True)
    (a.dir / "annotations" / "map.csv").write_text("# auto negatives: stem -> image (already coded)\n" + "\n".join(maplines) + "\n")
    pd.DataFrame(summary).to_csv(a.dir / "summary.csv", index=False); print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
