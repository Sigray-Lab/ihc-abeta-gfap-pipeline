#!/usr/bin/env python3
"""Prepare the data the plaque coder GUI shows: full-resolution three-channel crops of each box, the v4 classifier
objects inside it, the 3x-gate and DAPI-veto masks, and a list of 'practice spots' with what the machine says about
them. Reads only tiles that overlap the box (level 0, 0.325 um/px), so it is quick. Output: ~/ihc_work/manual_coding/data/."""
import sys, json, argparse
from pathlib import Path
import numpy as np, pandas as pd
from shapely.geometry import shape, box as sbox
from scipy import ndimage as ndi
from PIL import Image
REPO = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(REPO / "scripts" / "brightness_gate"))
from ihc.qc.stain_check import read_ets_index, _decode_tile, CH_DAPI, CH_FITC, CH_CY3
from gate import gate_section, manifest, geoms, rasterise
from veto import veto_section
M = Path.home() / "ihc_work/manual_coding"; PX = 0.325
def read_region(ets_path, channel, x0, y0, x1, y1, level=0):
    """Assemble one channel over [x0,x1) x [y0,y1) (level-0 pixel coordinates) from the tiles that overlap it."""
    idx = read_ets_index(ets_path); n_dim = idx["n_dim"]; ch_axis = 2 if n_dim >= 4 else None; tw, th = idx["tile"][0], idx["tile"][1]
    s = 2 ** level; X0, Y0, X1, Y1 = x0 // s, y0 // s, -(-x1 // s), -(-y1 // s)
    c0, c1, r0, r1 = X0 // tw, (X1 - 1) // tw, Y0 // th, (Y1 - 1) // th
    out = np.zeros((Y1 - Y0, X1 - X0), np.float32)
    with open(ets_path, "rb") as fh:
        for coords, off, nb in idx["chunks"]:
            if coords[-1] != level or (ch_axis is not None and coords[ch_axis] != channel): continue
            c, r = coords[0], coords[1]
            if not (c0 <= c <= c1 and r0 <= r <= r1): continue
            if off < 0 or nb <= 0 or off + nb > idx["file_bytes"]: continue
            fh.seek(off); t = _decode_tile(fh.read(nb))
            if t is None: continue
            ty, tx = r * th, c * tw                                   # tile origin in level coords
            ya, yb = max(ty, Y0), min(ty + t.shape[0], Y1); xa, xb = max(tx, X0), min(tx + t.shape[1], X1)
            if yb > ya and xb > xa: out[ya - Y0:yb - Y0, xa - X0:xb - X0] = t[ya - ty:yb - ty, xa - tx:xb - tx]
    return out
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--images", nargs="*"); a = ap.parse_args()
    idx = pd.read_csv(M / "boxes_index.csv"); man = manifest(); (M / "data").mkdir(exist_ok=True)
    for img, g in idx.groupby("image", sort=False):
        if a.images and img not in a.images: continue
        r = man.loc[img]; px = float(r.pixel_size_um); kind = "control" if r.condition != "positive" else "stained"
        print(f"{img}: {len(g)} boxes, {kind}", flush=True)
        gj = M / "v4_objects" / f"{img}.geojson"; gj5 = M / "v5_objects" / f"{img}.geojson"
        objs5 = []
        if gj5.exists():
            for ft in json.loads(gj5.read_text())["features"]:
                sh = shape(ft["geometry"]).buffer(0); objs5.extend(list(sh.geoms) if sh.geom_type == "MultiPolygon" else [sh])
        _, GM = gate_section(img, gj, keep_masks=True); _, VM = veto_section(img, gj, keep_masks=True); s2 = GM["s"]
        objs = []
        for ft in json.loads(gj.read_text())["features"]:
            sh = shape(ft["geometry"]).buffer(0); objs.extend(list(sh.geoms) if sh.geom_type == "MultiPolygon" else [sh])
        for _, b in g.iterrows():
            x0, y0, x1, y1 = int(b.x0), int(b.y0), int(b.x1), int(b.y1); W, Hh = x1 - x0, y1 - y0; win = sbox(x0, y0, x1, y1)
            chans = {nm: read_region(r.ets_path, ch, x0, y0, x1, y1) for nm, ch in (("cy3", CH_CY3), ("dapi", CH_DAPI), ("fitc", CH_FITC))}
            crop2 = lambda m: np.repeat(np.repeat(m[y0 // s2:-(-y1 // s2), x0 // s2:-(-x1 // s2)], s2, 0), s2, 1)[: Hh, : W]
            gate3 = crop2(GM["gates"][3.0]); gate2 = crop2(GM["gates"][2.0]); gate15 = crop2(GM["gates"][1.5]); veto = crop2(VM["veto"] & VM["abeta"]); v4m = crop2(GM["abeta"])
            def clip_polys(src):
                out = []
                for o in src:
                    if not o.intersects(win): continue
                    for pg in (list(o.geoms) if o.geom_type == "MultiPolygon" else [o]):
                        out.append([[float(x - x0), float(y - y0)] for x, y in pg.exterior.coords])
                return out
            polys = clip_polys(objs); polys5 = clip_polys(objs5)
            # practice spots from v4 objects inside the box
            dapi_bg = float(np.median(chans["dapi"])); cy3_bg = float(np.median(chans["cy3"][~v4m])) if (~v4m).any() else 1.0
            lab, n = ndi.label(v4m); spots = []
            if n:
                ids = np.arange(1, n + 1); area = ndi.sum(np.ones_like(lab), lab, ids) * px * px
                g3 = ndi.sum(gate3, lab, ids) / np.maximum(ndi.sum(np.ones_like(lab), lab, ids), 1)
                dmed = np.array(ndi.median(chans["dapi"], lab, ids)) / max(dapi_bg, 1); cmed = np.array(ndi.median(chans["cy3"], lab, ids)) / max(cy3_bg, 1)
                cents = ndi.center_of_mass(v4m, lab, ids)
                cats = []
                for i in range(n):
                    A, G3, D, C = area[i], g3[i], dmed[i], cmed[i]
                    if A < 20: continue
                    if kind == "control": cat, txt = "control: false positive", "Antibody-free section: nothing here is amyloid. The classifier called it anyway (rim glow, autofluorescent cell, debris or dim texture)."
                    elif D >= 1.5: cat, txt = "DAPI-bright: autofluorescence suspect", "Bright in the DAPI channel too. Antibody signal is Cy3-only, so broad glow points to an autofluorescent cell or debris, not a deposit."
                    elif A >= 300 and G3 >= 0.4: cat, txt = "typical deposit", "Bright, granular, most of it survives the 3x brightness gate. Every method calls this amyloid."
                    elif A >= 300 and G3 >= 0.1: cat, txt = "deposit with wide halo", "A bright granular core that survives the gate, inside a much larger v4 outline that runs into dim surroundings. Count the deposit; the halo is the over-classification you circled."
                    elif 20 <= A < 300 and G3 >= 0.5: cat, txt = "small bright deposit", "Small but clearly brighter than its surroundings. Kept by the gate; v5 would have dropped it."
                    elif A >= 100 and G3 < 0.1: cat, txt = "dim, v4 only", "Only the v4 classifier calls this amyloid. Dim texture with no bright core; the 3x gate removes it. Look for a real deposit inside before counting."
                    else: continue
                    cats.append((cat, txt, i))
                ORDER = ["typical deposit", "deposit with wide halo", "small bright deposit", "dim, v4 only", "DAPI-bright: autofluorescence suspect", "control: false positive"]
                for cat in [c for c in ORDER if any(cc == c for cc, _, _ in cats)]:
                    members = [(t, i) for c, t, i in cats if c == cat]; chosen = []
                    for _ in range(min(3, len(members))):
                        best = max(members, key=lambda m: (min([np.hypot(cents[m[1]][1] - cents[j][1], cents[m[1]][0] - cents[j][0]) for _, j in chosen] or [1e9]), area[m[1]]))
                        chosen.append(best); members.remove(best)
                    for t, i in chosen:
                        spots.append(dict(x=float(cents[i][1]), y=float(cents[i][0]), r=float(max(15, np.sqrt(area[i] / np.pi) / px * 1.6)), category=cat, text=t,
                                          area_um2=round(float(area[i])), gate3_survival=round(float(g3[i]), 2), dapi_contrast=round(float(dmed[i]), 2), cy3_contrast=round(float(cmed[i]), 2)))
            stem = f"{img}_{b.box}"
            np.savez_compressed(M / "data" / f"{stem}.npz", cy3=chans["cy3"].astype(np.uint16), dapi=chans["dapi"].astype(np.uint16), fitc=chans["fitc"].astype(np.uint16))
            for nm, m in (("gate3", gate3), ("gate2", gate2), ("gate15", gate15), ("veto", veto), ("v4", v4m)): Image.fromarray((m.astype(np.uint8) * 255)).save(M / "data" / f"{stem}_{nm}.png")
            json.dump(dict(image=img, box=b.box, bounds=[x0, y0, x1, y1], px_um=px, side_um=float(b.side_um), kind=kind, label=b.label, set=str(b.get("set", "practice")), v4_polys=polys, v5_polys=polys5, spots=spots,
                           v4_pct_of_box=round(100 * float(v4m.mean()), 2), gate3_pct_of_box=round(100 * float(gate3.mean()), 2)), open(M / "data" / f"{stem}.json", "w"))
            print(f"   {stem}: {W}x{Hh} px, v4 {100*v4m.mean():.1f}% gate3 {100*gate3.mean():.1f}%, {len(polys)} v4 objects, {len(spots)} practice spots", flush=True)
if __name__ == "__main__": main()
