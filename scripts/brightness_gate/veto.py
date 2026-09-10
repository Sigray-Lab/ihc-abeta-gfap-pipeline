"""v6 = v4 minus a DAPI-channel autofluorescence veto (Aris's 'second detector for the negatives', minimal deterministic form).
A v4 Abeta pixel is dropped where DAPI >= 1.5 x local DAPI background (400-um window of tissue that is not classifier-Abeta).
Rationale (v4_objects_channel_contrast.csv, 2026-09-10): on antibody-free controls ~30 % of v4's false-positive area is
DAPI-bright, on stained sections only 2-5 % of the Abeta area is; FITC cannot be used (GFAP gliosis lights up 23-63 % of
real plaque area). Also reports v6 combined with the 2x Cy3 gate. Same pipeline as gate.py."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import ndimage as ndi
sys.path.insert(0, str(Path(__file__).parent))
from gate import manifest, geoms, rasterise, read_pyramid_plane, CH_CY3, LEVEL, WIN_UM, HALO_UM
from ihc.qc.stain_check import CH_DAPI
DAPI_K = 1.5; GATE_K = 2.0
def local_bg(plane, valid, px_eff):
    f = 4; H4, W4 = plane.shape[0] // f, plane.shape[1] // f
    blk = lambda a: a[:H4 * f, :W4 * f].reshape(H4, f, W4, f).mean(axis=(1, 3))
    v = blk(valid.astype(np.float32)); num = blk(plane * valid); win = int(round(WIN_UM / (px_eff * f))) | 1
    vs = ndi.uniform_filter(v, win, mode="constant"); ns = ndi.uniform_filter(num, win, mode="constant")
    med = float(np.median(plane[valid])) if valid.any() else 1.0
    bg4 = np.where(vs > 0.05, ns / np.maximum(vs, 1e-6), med); bg = np.repeat(np.repeat(bg4, f, 0), f, 1)
    return np.pad(bg, ((0, plane.shape[0] - bg.shape[0]), (0, plane.shape[1] - bg.shape[1])), mode="edge")
def veto_section(image, geojson, keep_masks=False):
    r = manifest().loc[image]; px = float(r.pixel_size_um); kw = dict(pixel_size_um=px, true_size=(int(r.width_px), int(r.height_px)))
    cy3, lv = read_pyramid_plane(r.ets_path, CH_CY3, level=LEVEL, **kw); dapi, _ = read_pyramid_plane(r.ets_path, CH_DAPI, level=LEVEL, **kw)
    cy3 = np.nan_to_num(cy3.astype(np.float32)); dapi = np.nan_to_num(dapi.astype(np.float32)); s = 2 ** lv; pxA = (px * s) ** 2; shape = cy3.shape
    tissue = rasterise(geoms(Path.home() / f"ihc_work/tissue_rois/{image}.geojson"), shape, s); abeta = rasterise(geoms(geojson), shape, s)
    halo = ndi.binary_dilation(abeta, iterations=int(round(HALO_UM / (px * s)))); valid = tissue & ~halo
    rel_c = cy3 / np.maximum(local_bg(cy3, valid, px * s), 1e-6); rel_d = dapi / np.maximum(local_bg(dapi, valid, px * s), 1e-6)
    veto = rel_d >= DAPI_K; v6 = abeta & ~veto; g2 = abeta & (rel_c >= GATE_K); v6g2 = v6 & (rel_c >= GATE_K)
    out = dict(image=image, tissue_um2=tissue.sum() * pxA, abeta_v4_um2=(abeta & tissue).sum() * pxA, abeta_v6_um2=(v6 & tissue).sum() * pxA,
               abeta_gate2_um2=(g2 & tissue).sum() * pxA, abeta_v6_gate2_um2=(v6g2 & tissue).sum() * pxA, vetoed_share=float((abeta & veto & tissue).sum() / max((abeta & tissue).sum(), 1)))
    for reg in ("Hippocampus", "Isocortex"):
        rf = Path.home() / f"ihc_work/region_rois/{reg}/{image}.geojson"
        if not rf.exists(): continue
        rm = rasterise(geoms(rf), shape, s); out[f"{reg}_um2"] = rm.sum() * pxA
        for nm, m in [("v4", abeta), ("v6", v6), ("gate2", g2), ("v6_gate2", v6g2)]: out[f"{reg}_abeta_{nm}_um2"] = (m & rm).sum() * pxA
    if keep_masks: return out, dict(abeta=abeta, v6=v6, g2=g2, v6g2=v6g2, veto=veto, s=s, px=px)
    return out
if __name__ == "__main__":
    print(json.dumps(veto_section(sys.argv[1], sys.argv[2])))
