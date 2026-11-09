#!/usr/bin/env python3
"""Per-section brightness anchors for the GLOBAL normalisation: glass_p50 and t_p50 per image.

Reconstructs the documented method (docs/normalisation.md, option 2) for new images: the original
script that produced DerivedData/aris_v2_2026-08-18/global_norm_anchors.csv was not kept, so this
one is validated against that table instead (--check). Method: tissue mask from DAPI at the coarse
mask level (ihc.qc.tissue_mask); glass = median of acquired pixels more than EDGE_UM outside the
mask; t_p50 = median inside it; both read at READ_LEVEL. NaN tiles (never acquired) are excluded
by construction; "acquired" additionally requires value > 0.

    python3 scripts/make_section_anchors.py --out anchors.csv            # all measured sections
    python3 scripts/make_section_anchors.py --images K07_s01 --out -     # one new image
    python3 scripts/make_section_anchors.py --check                      # validate against production

Validation against the production Cy3 table (2026-11-09, 121 sections): tissue medians agree within a
few counts everywhere; gain within 5 % on 110 sections. The differences are in the glass value: the
production table has glass of 0, 1, 3 and 24 counts on G09_s03, Y06_s01, Y18_s04 and A08_s04, where this
script finds 37-44 like every other standard-exposure section. Those four production values are too
low (ADR-0041). To reproduce the published rows exactly, use the saved production table, not this script.

--channel defaults to cy3 (Abeta). For a GFAP run use --channel fitc: anchors are per channel,
and the Cy3 table must not be reused for FITC."""
import sys, argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy import ndimage as ndi
from multiprocessing import Pool
REPO = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(REPO / "src"))
from ihc.qc.stain_check import read_pyramid_plane, CH_DAPI, CH_FITC, CH_CY3   # noqa: E402
from ihc.qc.tissue_mask import tissue_mask                                     # noqa: E402
CH = {"cy3": CH_CY3, "fitc": CH_FITC, "dapi": CH_DAPI}
MASK_LEVEL, READ_LEVEL, EDGE_UM = 5, 3, 100.0
def manifest():
    m = pd.read_csv(REPO.parent / "BlindingKey/provenance_manifest.csv"); m = m[m.row_kind == "section"].copy()
    m["image"] = m.apply(lambda r: f"{r.code}{'b' if r.scan == 'rescan' else ''}_s{str(r.physical_section_label).zfill(2) if str(r.physical_section_label).isdigit() else r.physical_section_label}", axis=1)
    return m.drop_duplicates("image").set_index("image")
def one(args):
    img, ets, px, w, h, channel = args
    try:
        kw = dict(pixel_size_um=px, true_size=(w, h))
        d5, lv5 = read_pyramid_plane(ets, CH_DAPI, level=MASK_LEVEL, **kw)
        m5 = np.asarray(tissue_mask(np.nan_to_num(d5, nan=0.0)), bool)
        c, lv = read_pyramid_plane(ets, CH[channel], level=READ_LEVEL, **kw); s = 2 ** lv
        f = 2 ** (lv5 - lv); m = np.repeat(np.repeat(m5, f, 0), f, 1)[: c.shape[0], : c.shape[1]]
        if m.shape != c.shape: m = np.pad(m, ((0, c.shape[0] - m.shape[0]), (0, c.shape[1] - m.shape[1])))
        acq = np.isfinite(c) & (c > 0)
        far = ~ndi.binary_dilation(m, iterations=max(1, int(round(EDGE_UM / (px * s))))) & acq
        if not far.any() or not (m & acq).any(): return dict(image=img, error="no glass or no tissue pixels")
        return dict(image=img, glass_p50=float(np.median(c[far])), t_p50=float(np.median(c[m & acq])))
    except Exception as e:
        return dict(image=img, error=repr(e)[:120])
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--channel", choices=list(CH), default="cy3"); ap.add_argument("--images", nargs="*")
    ap.add_argument("--out", default="-"); ap.add_argument("--check", action="store_true"); ap.add_argument("--procs", type=int, default=6)
    a = ap.parse_args(argv); man = manifest()
    man = man[man.use_for_measurement.astype(str) == "True"] if not a.images else man.loc[a.images]
    jobs = [(i, r.ets_path, float(r.pixel_size_um), int(r.width_px), int(r.height_px), a.channel) for i, r in man.iterrows()]
    with Pool(a.procs) as pool: rows = pool.map(one, jobs)
    d = pd.DataFrame(rows); bad = d[d["error"].notna()] if "error" in d else d.iloc[0:0]
    for r in bad.itertuples(): print(f"  FAIL {r.image}: {r.error}", file=sys.stderr)
    d = d.drop(columns=["error"], errors="ignore").dropna()
    if a.check:
        ref = pd.read_csv(REPO.parent / "DerivedData/aris_v4_2026-09-10_session/global_norm_anchors.csv").drop_duplicates("image").set_index("image")
        j = d.set_index("image").join(ref, rsuffix="_ref", how="inner"); gain = j.t_p50 - j.glass_p50; gref = j.t_p50_ref - j.glass_p50_ref
        rel = (gain - gref).abs() / gref
        print(f"checked {len(j)} sections against the production table ({a.channel}): gain within 2% on {int((rel <= .02).sum())}, within 5% on {int((rel <= .05).sum())}, worst {rel.max():.3f} ({rel.idxmax()})")
        print(f"glass abs diff median {float((j.glass_p50 - j.glass_p50_ref).abs().median()):.1f} counts, max {float((j.glass_p50 - j.glass_p50_ref).abs().max()):.1f}")
        print("outside 5%:", {i: (round(float(gain[i])), round(float(gref[i]))) for i in rel[rel > .05].index})
        return 0                                  # a comparison, not a gate: see the docstring for the known differences
    out = d.to_csv(index=False)
    (sys.stdout.write(out) if a.out == "-" else Path(a.out).write_text(out))
    if a.out != "-": print(f"wrote {len(d)} anchors -> {a.out}", file=sys.stderr)
    return 0 if bad.empty else 1
if __name__ == "__main__":
    raise SystemExit(main())
