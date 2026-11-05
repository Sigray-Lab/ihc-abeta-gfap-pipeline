#!/usr/bin/env python3
"""Focus score per section, for every fluorescence section in the manifest: width of nuclear edges in um (higher = blurrier).

Two earlier scores were dropped after failing validation: a raw gradient and a power-spectrum ratio both let pixel noise in
dim sections pass as fine detail, so blurred sections scored like sharp ones. Edge width has physical units, grows with
defocus whatever the brightness, and is computed only on tiles with enough signal. Tiles are spread evenly over tissue at
least 150 um from the edge. Usage: focus_scores.py [--images A B C] [--tiles 48] [--procs 4] [--tag _x].
Output: ~/ihc_work/focus_20261105/focus_{sections,tiles}.csv"""
import sys, argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy import ndimage as ndi
from multiprocessing import Pool
REPO = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(REPO / "src")); sys.path.insert(0, str(REPO / "scripts" / "brightness_gate"))
from ihc.qc.stain_check import read_ets_index, read_pyramid_plane, _decode_tile, CH_DAPI, CH_CY3
from ihc.qc.tissue_mask import tissue_mask
from gate import manifest
OUT = Path.home() / "ihc_work/focus_20261105"; OUT.mkdir(exist_ok=True); MIN_SNR = 8.0
def tile_score(t, px):
    """Edge width in um: step height of the brightest structures over their steepest gradient, after light smoothing.
    A blurred edge of Gaussian width sigma gives width = 2.5 sigma, so defocus widens it whatever the brightness; the
    smoothing keeps pixel noise out of the gradient. Returns (edge_um, snr); tiles with low snr carry no usable edges."""
    a = t.astype(np.float32); d = ndi.gaussian_filter(a, 2.0)
    noise = 1.4826 * float(np.median(np.abs(a - d))) + 1e-6
    H = float(np.percentile(d, 99.5) - np.percentile(d, 50))
    g = np.hypot(*np.gradient(d)) / px                                    # counts per um
    G = float(np.percentile(g, 99.5))
    return (H / G if G > 0 else np.nan), H / noise
def score_section(args):
    img, ets, px, w, h, ntiles = args
    try:
        d5, lv = read_pyramid_plane(ets, CH_DAPI, level=5, pixel_size_um=px, true_size=(w, h)); s = 2 ** lv
        tissue = np.asarray(tissue_mask(d5), bool); tissue = ndi.binary_erosion(tissue, iterations=max(1, int(150 / (px * s))))
        idx = read_ets_index(ets); tw, th = idx["tile"][0], idx["tile"][1]; lut = {}
        for coords, off, nb in idx["chunks"]:
            if coords[-1] == 0: lut[(coords[0], coords[1], coords[2])] = (off, nb)
        cand = []
        for (c, r, ch) in lut:
            if ch != CH_DAPI: continue
            y0, y1, x0, x1 = r * th // s, -(-(r + 1) * th // s), c * tw // s, -(-(c + 1) * tw // s)
            blk = tissue[y0:y1, x0:x1]
            if blk.size and blk.mean() > 0.9: cand.append((r, c))
        cand.sort(); pick = [cand[i] for i in np.unique(np.linspace(0, len(cand) - 1, min(ntiles, len(cand))).astype(int))] if cand else []
        rows = []
        with open(ets, "rb") as fh:
            for r, c in pick:
                row = dict(image=img, tile_row=r, tile_col=c)
                for nm, ch in (("dapi", CH_DAPI), ("cy3", CH_CY3)):
                    if (c, r, ch) not in lut: continue
                    off, nb = lut[(c, r, ch)]; fh.seek(off); t = _decode_tile(fh.read(nb))
                    if t is not None and t.shape[0] == th and t.shape[1] == tw:
                        row[nm], row[f"{nm}_snr"] = tile_score(t, px)
                rows.append(row)
        return img, rows, None
    except Exception as e:
        return img, [], repr(e)
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--images", nargs="*"); ap.add_argument("--tiles", type=int, default=48); ap.add_argument("--procs", type=int, default=4); ap.add_argument("--tag", default=""); a = ap.parse_args()
    man = manifest(); man = man[man.index.isin(a.images)] if a.images else man
    jobs = [(i, r.ets_path, float(r.pixel_size_um), int(r.width_px), int(r.height_px), a.tiles) for i, r in man.iterrows()]; tiles = []; sec = []
    with Pool(a.procs) as pool:
        for img, rows, err in pool.imap_unordered(score_section, jobs):
            if err: print(f"{img}: FAILED {err}", flush=True); sec.append(dict(image=img, n_tiles=0, error=err)); continue
            t = pd.DataFrame(rows); tiles += rows
            if t.empty or "dapi_snr" not in t:
                print(f"{img:12s} NO USABLE TILES (tissue mask left none at least 150 um from the edge)", flush=True); sec.append(dict(image=img, n_tiles=0, n_usable=0)); continue
            d = t[t.dapi_snr >= MIN_SNR].dapi if "dapi" in t else pd.Series(dtype=float); c = t[t.cy3_snr >= MIN_SNR].cy3 if "cy3" in t else pd.Series(dtype=float)
            sec.append(dict(image=img, n_tiles=len(t), n_usable=len(d), edge_um=d.median(), edge_p10=d.quantile(.1), edge_p90=d.quantile(.9), cy3_edge_um=c.median(), cy3_usable=len(c), dapi_snr=t.dapi_snr.median()))
            print(f"{img:12s} tiles {len(t):3d} usable {len(d):3d}  DAPI edge {d.median():5.2f} um [{d.quantile(.1):5.2f} to {d.quantile(.9):5.2f}]  snr {t.dapi_snr.median():5.1f}  Cy3 edge {c.median():5.2f} (n {len(c)})", flush=True)
    pd.DataFrame(sec).sort_values("edge_um", ascending=False).to_csv(OUT / f"focus_sections{a.tag}.csv", index=False); pd.DataFrame(tiles).to_csv(OUT / f"focus_tiles{a.tag}.csv", index=False); print("DONE", len(sec), "sections")
if __name__ == "__main__": main()
