#!/usr/bin/env python3
"""Blob-recovery ("halo") experiment, stage 1: build the synthetic images.

Design: docs/PLAN_halo_experiment.md. Secondary-only control sections carry no primary
antibody, so anything a classifier reports there is either an injected blob or a known false
positive. Each host is exported at pyramid level 2 -- 1.3 um/px, the classifiers' own working
resolution, so QuPath resamples nothing -- as a 3-channel OME-TIFF, once clean and once with
Gaussian blobs added to Cy3 in raw counts: FWHM 20/40/80/160/320 um x peak 1.25/1.5/2/4 x the
LOCAL tissue background, placed inside the DAPI tissue mask, apart from each other and from
the tissue edge. A per-image JSON records every blob (centre, sigma, amplitude, footprint).
Footprint = disc where the added signal is >= 25 % of the blob peak (radius 1.665 sigma).

    python3 scripts/halo_inject.py --hosts W12_s03 S17_s01 L03_s01 --out ~/ihc_work/halo_20260909
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
import tifffile
from scipy import ndimage as ndi

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from ihc.qc.stain_check import read_pyramid_plane, CH_DAPI, CH_FITC, CH_CY3   # noqa: E402
from ihc.qc.tissue_mask import tissue_mask                                     # noqa: E402

LEVEL = 2                       # 0.325 * 2**2 = 1.3 um/px
UM = 0.325 * 2 ** LEVEL
FWHM_UM = [20, 40, 80, 160, 320]
AMPS = [1.25, 1.5, 2.0, 4.0]
MIN_SEP_UM, EDGE_UM = 400.0, 300.0
CHANNELS = ["DAPI", "FITC", "Cy3"]


MASK_LEVEL = 5                  # tissue_mask() is calibrated for the pipeline's coarse mask level


def load_host(image, man):
    r = man.loc[image]
    kw = dict(pixel_size_um=float(r.pixel_size_um), true_size=(int(r.width_px), int(r.height_px)))
    planes = [read_pyramid_plane(r.ets_path, ch, level=LEVEL, **kw)[0] for ch in (CH_DAPI, CH_FITC, CH_CY3)]
    support = np.isfinite(planes[2])
    # the tissue mask is built where the pipeline builds it (level 5) and upsampled, because
    # tissue_mask()'s smoothing and size rules are in pixels of that level
    d5, lv5 = read_pyramid_plane(r.ets_path, CH_DAPI, level=MASK_LEVEL, **kw)
    m5 = tissue_mask(np.nan_to_num(d5, nan=0.0))
    f = 2 ** (lv5 - LEVEL)
    m = np.repeat(np.repeat(m5, f, axis=0), f, axis=1)
    H, W = planes[0].shape
    m = np.pad(m, ((0, max(0, H - m.shape[0])), (0, max(0, W - m.shape[1]))))[:H, :W]
    return [np.nan_to_num(p, nan=0.0) for p in planes], support, m & support


def write_ome(path, planes):
    arr = np.stack([np.clip(np.rint(p), 0, 65535).astype(np.uint16) for p in planes])   # C,Y,X
    md = {"axes": "CYX", "Channel": {"Name": CHANNELS}, "PhysicalSizeX": UM, "PhysicalSizeXUnit": "µm",
          "PhysicalSizeY": UM, "PhysicalSizeYUnit": "µm"}
    tifffile.imwrite(str(path), arr, ome=True, metadata=md, compression="zlib", tile=(512, 512))


def place(mask, support, rng):
    """Greedy placement of the 20 (width x amplitude) blobs inside the eroded tissue."""
    px = UM
    dist = ndi.distance_transform_edt(mask & support) * px          # um to nearest non-tissue
    blobs, tries = [], 0
    specs = [(f, a) for f in FWHM_UM for a in AMPS]
    rng.shuffle(specs)
    for fwhm, amp in specs:
        sigma_um = fwhm / 2.355; foot_um = 1.665 * sigma_um
        ok = dist >= (EDGE_UM + foot_um)
        ys, xs = np.nonzero(ok)
        if not len(ys): continue
        for _ in range(4000):
            i = rng.integers(len(ys)); y, x = ys[i], xs[i]
            if all(np.hypot(x - b["cx"], y - b["cy"]) * px >= MIN_SEP_UM + foot_um + b["foot_um"] for b in blobs):
                blobs.append(dict(cx=int(x), cy=int(y), fwhm_um=fwhm, sigma_px=sigma_um / px, amp_factor=amp,
                                  foot_um=foot_um, foot_px=foot_um / px)); break
    return blobs


def inject(cy3, mask, blobs):
    out = cy3.astype(np.float64).copy()
    yy, xx = np.mgrid[0:cy3.shape[0], 0:cy3.shape[1]]
    for b in blobs:
        r_ann = 200 / UM
        ann = ((xx - b["cx"]) ** 2 + (yy - b["cy"]) ** 2 <= r_ann ** 2) & mask
        local = float(np.median(cy3[ann]))
        b["local_bg"] = local; b["added_peak"] = (b["amp_factor"] - 1) * local
        g = np.exp(-((xx - b["cx"]) ** 2 + (yy - b["cy"]) ** 2) / (2 * b["sigma_px"] ** 2))
        out += b["added_peak"] * g
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--hosts", nargs="+", required=True, help="coded image names of secondary-only sections")
    ap.add_argument("--out", type=Path, default=Path.home() / "ihc_work/halo_20260909")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(REPO.parent / "BlindingKey" / "provenance_manifest.csv")
    man = man[man.row_kind == "section"].copy()
    def _img(r):
        pl = str(r.physical_section_label); s = pl.zfill(2) if pl.isdigit() else pl
        return f"{r.code}{'b' if r.scan == 'rescan' else ''}_s{s}"
    man["image"] = man.apply(_img, axis=1); man = man.drop_duplicates("image").set_index("image")
    rng = np.random.default_rng(a.seed)
    anchors = []
    for host in a.hosts:
        assert man.loc[host].condition == "negative", f"{host} is not a secondary-only control"
        planes, support, mask = load_host(host, man)
        dapi, fitc, cy3 = planes
        glass = float(np.median(cy3[support & ~ndi.binary_dilation(mask, iterations=int(200 / UM))]))
        t50 = float(np.median(cy3[mask]))
        blobs = place(mask, support, rng)
        cy3_inj = inject(cy3, mask, blobs)
        write_ome(a.out / f"{host}_clean.ome.tif", [dapi, fitc, cy3])
        write_ome(a.out / f"{host}_blobs.ome.tif", [dapi, fitc, cy3_inj])
        (a.out / f"{host}_blobs.json").write_text(json.dumps(dict(host=host, um_per_px=UM, glass_p50=glass, t_p50=t50,
                                                                   tissue_px=int(mask.sum()), blobs=blobs), indent=1))
        anchors += [dict(image=f"{host}_clean", glass_p50=glass, t_p50=t50), dict(image=f"{host}_blobs", glass_p50=glass, t_p50=t50)]
        print(f"  {host}: tissue {mask.sum()*UM*UM/1e6:.1f} mm2, cy3 tissue median {t50:.0f}, glass {glass:.0f}, "
              f"{len(blobs)}/20 blobs placed", flush=True)
    pd.DataFrame(anchors).to_csv(a.out / "anchors.csv", index=False)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
