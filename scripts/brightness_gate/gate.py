"""Brightness gate on v4 (GLOBAL) Abeta output: keep Abeta pixels whose Cy3 >= k x local background.
Local background = mean Cy3 over a 400-um window of tissue pixels that are NOT classifier-Abeta (10-um dilation excluded),
computed at 5.2 um/px and upsampled; falls back to the section's tissue median where <5% of the window is usable.
Everything at the classifier's working resolution (level 2, 1.3 um/px). Same rule on stained and control sections."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import ndimage as ndi
REPO = Path("/Users/pontusps/Dropbox/longevity_project/projects/Rapamycin_CSF_cell_mice_Per/IHC_analysis_pipeline/pipeline"); sys.path.insert(0, str(REPO / "src"))
from ihc.qc.stain_check import read_pyramid_plane, CH_CY3
THRESH = (1.5, 2.0, 3.0); WIN_UM = 400.0; HALO_UM = 10.0; LEVEL = 2
_man = None
def manifest():
    global _man
    if _man is None:
        m = pd.read_csv(REPO.parent / "BlindingKey/provenance_manifest.csv"); m = m[m.row_kind == "section"].copy()
        def _img(r):
            pl = str(r.physical_section_label); s = pl.zfill(2) if pl.isdigit() else pl
            return f"{r.code}{'b' if r.scan == 'rescan' else ''}_s{s}"
        m["image"] = m.apply(_img, axis=1); _man = m.drop_duplicates("image").set_index("image")
    return _man
def geoms(path):
    out = []
    for ft in json.loads(Path(path).read_text())["features"]:
        g = ft["geometry"]
        if g["type"] == "Polygon": out.append(g["coordinates"])
        elif g["type"] == "MultiPolygon": out.extend(g["coordinates"])
    return out                                  # list of [exterior, hole, hole...]
def rasterise(polys, shape, s):
    """Pixel-centre rasterisation (vertices are pixel corners at the classifier grid, hence the -0.5 shift).
    Reproduces QuPath's own 'Abeta area' measurement within 2% on T18_s03, C18_s01, Z23_s04 (2026-09-09);
    PIL/OpenCV fills include the boundary ring and overstate small objects by up to 35%."""
    from skimage.draw import polygon as skpoly
    m = np.zeros(shape, bool)
    for rings in polys:
        ext = np.asarray(rings[0], float) / s - 0.5
        if len(ext) < 3: continue
        rr, cc = skpoly(ext[:, 1], ext[:, 0], shape); m[rr, cc] = True
        for h in rings[1:]:
            hh = np.asarray(h, float) / s - 0.5
            if len(hh) >= 3: rr, cc = skpoly(hh[:, 1], hh[:, 0], shape); m[rr, cc] = False
    return m
def gate_section(image, geojson, regions=("Hippocampus", "Isocortex"), keep_masks=False):
    r = manifest().loc[image]; px = float(r.pixel_size_um); kw = dict(pixel_size_um=px, true_size=(int(r.width_px), int(r.height_px)))
    cy3, lv = read_pyramid_plane(r.ets_path, CH_CY3, level=LEVEL, **kw); s = 2 ** lv; cy3 = np.nan_to_num(cy3.astype(np.float32), nan=0.0)
    pxA = (px * s) ** 2; shape = cy3.shape
    tissue = rasterise(geoms(Path.home() / f"ihc_work/tissue_rois/{image}.geojson"), shape, s)
    abeta = rasterise(geoms(geojson), shape, s)
    # ---- local background at level+2 (4x coarser) ----
    f = 4; H4, W4 = shape[0] // f, shape[1] // f
    blk = lambda a: a[:H4 * f, :W4 * f].reshape(H4, f, W4, f).mean(axis=(1, 3))
    halo = ndi.binary_dilation(abeta, iterations=int(round(HALO_UM / (px * s))))
    valid = blk((tissue & ~halo).astype(np.float32)); num = blk(cy3 * (tissue & ~halo))
    win = int(round(WIN_UM / (px * s * f))) | 1
    vsum = ndi.uniform_filter(valid, win, mode="constant"); nsum = ndi.uniform_filter(num, win, mode="constant")
    tmed = float(np.median(cy3[tissue & ~halo])) if (tissue & ~halo).any() else float(np.median(cy3[tissue]))
    bg4 = np.where(vsum > 0.05, nsum / np.maximum(vsum, 1e-6), tmed)
    bg = np.repeat(np.repeat(bg4, f, 0), f, 1); bg = np.pad(bg, ((0, shape[0] - bg.shape[0]), (0, shape[1] - bg.shape[1])), mode="edge")
    rel = cy3 / np.maximum(bg, 1e-6)
    out = dict(image=image, tissue_um2=tissue.sum() * pxA, abeta_v4_um2=(abeta & tissue).sum() * pxA, tissue_median_cy3=tmed)
    gates = {k: abeta & (rel >= k) for k in THRESH}
    for k, g in gates.items(): out[f"abeta_gate{k:g}_um2"] = (g & tissue).sum() * pxA
    for reg in regions:
        rf = Path.home() / f"ihc_work/region_rois/{reg}/{image}.geojson"
        if not rf.exists(): continue
        rm = rasterise(geoms(rf), shape, s); out[f"{reg}_um2"] = rm.sum() * pxA; out[f"{reg}_abeta_v4_um2"] = (abeta & rm).sum() * pxA
        for k, g in gates.items(): out[f"{reg}_abeta_gate{k:g}_um2"] = (g & rm).sum() * pxA
    if keep_masks: return out, dict(cy3=cy3, bg=bg, abeta=abeta, tissue=tissue, gates=gates, s=s, px=px)
    return out
if __name__ == "__main__":
    img, gj = sys.argv[1], sys.argv[2]; print(json.dumps(gate_section(img, gj)))
