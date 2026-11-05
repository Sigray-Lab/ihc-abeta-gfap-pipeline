#!/usr/bin/env python3
"""Blob-recovery experiment, stage 3: recovery versus width per classifier.

Reads the CSVs written by halo_eval.groovy. Recovery = Abeta area inside the footprint disc of
an injected blob divided by the disc area, after subtracting the host's own false-positive
rate in that disc measured on the clean image. Writes halo_recovery.csv and halo_recovery.png.
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--dir", type=Path, default=Path.home() / "ihc_work/halo_20260909")
    a = ap.parse_args()
    d = pd.concat([pd.read_csv(f) for f in sorted(a.dir.glob("halo_*.csv")) if f.name != "halo_recovery.csv"], ignore_index=True)
    d = d.drop_duplicates()          # a run launched twice appends identical rows
    d["host"] = d.image.str.replace(r"_(clean|blobs)$", "", regex=True); d["kind"] = d.image.str.extract(r"_(clean|blobs)$")
    d["arm"] = d.classifier.str.replace(r"^gx_.*", "GLOBAL", regex=True)
    d["frac"] = d.abeta_um2 / d.foot_area_um2
    key = ["host", "arm", "blob", "scale"]
    clean = d[d.kind == "clean"].set_index(key).frac
    b = d[d.kind == "blobs"].copy()
    b["baseline"] = b.set_index(key).index.map(clean).values
    b["recovery"] = (b.frac - b.baseline.fillna(0)).clip(lower=0)
    b.to_csv(a.dir / "halo_recovery.csv", index=False)
    s = b[b.scale == 1].groupby(["arm", "fwhm_um", "amp_factor"]).recovery.median().unstack("amp_factor")
    print("median recovery inside the footprint (scale 1), rows=fwhm um, cols=peak/background:\n", s.round(2).to_string())
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.2), sharey=True)
    for ax, amp in zip(axes, sorted(b.amp_factor.unique())):
        for arm, g in b[(b.scale == 1) & (b.amp_factor == amp)].groupby("arm"):
            m = g.groupby("fwhm_um").recovery.median()
            ax.plot(m.index, m.values, marker="o", label=arm)
        ax.set_xscale("log"); ax.set_title(f"peak = {amp}x background", fontsize=9); ax.set_xlabel("blob FWHM (um)")
        ax.axvline(130, color="grey", lw=.8, ls="--")
    axes[0].set_ylabel("recovered fraction of footprint"); axes[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(a.dir / "halo_recovery.png", dpi=200); print("wrote", a.dir / "halo_recovery.png")


if __name__ == "__main__":
    main()
