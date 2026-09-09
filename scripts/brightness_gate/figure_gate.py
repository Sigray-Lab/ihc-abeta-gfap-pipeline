"""Cohort figure: v4 vs gated (2x) per section, stained vs controls; and the four gate levels as a strip."""
import sys
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
IN, OUT = Path(sys.argv[1]), Path(sys.argv[2])
d = pd.read_csv(IN / "section_level_gate.csv", index_col=0)
cols = [("v4", "pct_abeta_v4_um2"), ("gate 1.5x", "pct_abeta_gate1.5_um2"), ("gate 2x", "pct_abeta_gate2_um2"), ("gate 3x", "pct_abeta_gate3_um2")]
st, ct = d[d.condition == "positive"], d[d.condition != "positive"]
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
ax = axes[0]; ax.scatter(st.pct_abeta_v4_um2, st.pct_abeta_gate2_um2, s=22, c="tab:red", label=f"stained (n={len(st)})"); ax.scatter(ct.pct_abeta_v4_um2, ct.pct_abeta_gate2_um2, s=22, c="tab:blue", label=f"antibody-free controls (n={len(ct)})")
lim = max(d.pct_abeta_v4_um2.max(), 1) * 1.05; ax.plot([0, lim], [0, lim], "k:", lw=.8, label="no change"); ax.plot([0, lim], [0, 0.38 * lim], "k--", lw=.8, label="x0.38 (median stained ratio)")
ax.set_xlabel("v4 GLOBAL, Abeta % of tissue"); ax.set_ylabel("v4 gated at 2x local background, %"); ax.legend(fontsize=8); ax.set_title("Per section: gate removes ~60% of v4 area, ranking kept", fontsize=10)
ax = axes[1]; xs = np.arange(len(cols))
for i, (name, c) in enumerate(cols):
    ax.scatter(np.full(len(st), i - 0.15) + np.random.default_rng(1).uniform(-.06, .06, len(st)), st[c], s=12, c="tab:red", alpha=.7)
    ax.scatter(np.full(len(ct), i + 0.15) + np.random.default_rng(2).uniform(-.06, .06, len(ct)), ct[c], s=12, c="tab:blue", alpha=.7)
    ax.hlines(st[c].median(), i - 0.3, i, colors="tab:red", lw=2); ax.hlines(ct[c].max(), i, i + 0.3, colors="tab:blue", lw=2)
ax.set_yscale("log"); ax.set_xticks(xs); ax.set_xticklabels([n for n, _ in cols]); ax.set_ylabel("Abeta % of tissue (log)")
ax.set_title("Red bar = stained median, blue bar = worst control", fontsize=10)
fig.suptitle("Brightness gate on v4 GLOBAL output, all sections, 2026-09-09", fontsize=12); fig.savefig(OUT / "gate_cohort.png", dpi=160, bbox_inches="tight"); print("figure written")
