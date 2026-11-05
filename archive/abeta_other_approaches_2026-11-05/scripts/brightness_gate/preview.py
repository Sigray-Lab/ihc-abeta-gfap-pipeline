import sys, json
from pathlib import Path
import numpy as np, cv2
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from shapely.geometry import shape
sys.path.insert(0, str(Path(__file__).parent)); from gate import gate_section, read_pyramid_plane, CH_CY3, manifest
IMG, GJ, OUT = sys.argv[1], sys.argv[2], Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
meas, M = gate_section(IMG, GJ, keep_masks=True); s, px = M["s"], M["px"]
print({k: (round(v / 1e6, 3) if isinstance(v, float) and "um2" in k else v) for k, v in meas.items() if "um2" in k and ("Hippocampus" in k or k.startswith(("tissue", "abeta")))})
r = manifest().loc[IMG]; kw = dict(pixel_size_um=px, true_size=(int(r.width_px), int(r.height_px)))
cy1, lv1 = read_pyramid_plane(r.ets_path, CH_CY3, level=1, **kw); s1 = 2 ** lv1
boxes = [(ft["properties"]["name"], shape(ft["geometry"])) for ft in json.loads((Path.home() / f"ihc_work/reference_boxes/{IMG}.geojson").read_text())["features"] if ft["properties"]["name"].startswith("HIP")]
objs = []
for ft in json.loads(Path(GJ).read_text())["features"]:
    g = shape(ft["geometry"]).buffer(0); objs.extend(list(g.geoms) if g.geom_type == "MultiPolygon" else [g])
for bname, b in boxes[:2]:
    x0, y0, x1, y1 = [int(v) for v in b.bounds]
    crop = cy1[y0 // s1:y1 // s1, x0 // s1:x1 // s1]; lo, hi = np.nanpercentile(crop, [1, 99.7]); disp = np.clip((np.nan_to_num(crop, nan=lo) - lo) / (hi - lo), 0, 1)
    g2 = M["gates"][2.0][y0 // s:y1 // s, x0 // s:x1 // s]; a2 = M["abeta"][y0 // s:y1 // s, x0 // s:x1 // s]
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.6))
    for ax, title in zip(axes, ["v4 as it is (red outline)", "v4 gated at 2x local background (yellow = kept, red = v4 outline)"]):
        ax.imshow(disp, cmap="gray", vmin=0, vmax=1, extent=(x0, x1, y1, y0)); ax.set_axis_off(); ax.set_title(title, fontsize=10); ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
        for o in objs:
            if o.intersects(b):
                for pg in (list(o.geoms) if o.geom_type == "MultiPolygon" else [o]): x, y = pg.exterior.xy; ax.plot(x, y, color="red", lw=0.7)
    ov = np.zeros(g2.shape + (4,)); ov[g2] = (1, 1, 0, 0.45); axes[1].imshow(ov, extent=(x0, x1, y1, y0), interpolation="nearest")
    fig.suptitle(f"{IMG} {bname}: v4 Abeta {100*a2.mean():.1f}% of box -> gated {100*g2.mean():.1f}%", fontsize=12)
    fig.savefig(OUT / f"{IMG}_{bname}_gate2x.png", dpi=170, bbox_inches="tight"); plt.close(fig)
# gated result as GeoJSON annotations for QuPath (drag the file onto the viewer)
m = M["gates"][2.0].astype(np.uint8); cs, hier = cv2.findContours(m, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE); feats = []
if hier is not None:
    hier = hier[0]
    for i, c in enumerate(cs):
        if hier[i][3] != -1 or len(c) < 3: continue
        rings = [[(float(x * s), float(y * s)) for x, y in c[:, 0, :]] + [(float(c[0, 0, 0] * s), float(c[0, 0, 1] * s))]]
        j = hier[i][2]
        while j != -1:
            h = cs[j]
            if len(h) >= 3: rings.append([(float(x * s), float(y * s)) for x, y in h[:, 0, :]] + [(float(h[0, 0, 0] * s), float(h[0, 0, 1] * s))])
            j = hier[j][0]
        feats.append({"type": "Feature", "geometry": {"type": "Polygon", "coordinates": rings}, "properties": {"objectType": "detection", "classification": {"name": "Abeta_gated2x", "color": [255, 220, 0]}}})
Path(OUT / f"{IMG}_v4_gated2x.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats})); print(len(feats), "gated polygons written")
