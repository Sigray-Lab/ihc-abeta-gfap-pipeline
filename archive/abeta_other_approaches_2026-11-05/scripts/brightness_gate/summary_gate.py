"""Cohort table: v4 GLOBAL vs brightness-gated v4 (1.5x / 2x / 3x local background)."""
import sys, json, glob
from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import spearmanr, mannwhitneyu
P = Path("/Users/pontusps/Dropbox/longevity_project/projects/Rapamycin_CSF_cell_mice_Per/IHC_analysis_pipeline"); R = Path(sys.argv[1]); OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
BAD = {"Y13_s01", "S17_s04"}
man = pd.read_csv(P / "BlindingKey/provenance_manifest.csv"); man = man[man.row_kind == "section"].copy()
def _img(r):
    pl = str(r.physical_section_label); s = pl.zfill(2) if pl.isdigit() else pl
    return f"{r.code}{'b' if r.scan == 'rescan' else ''}_s{s}"
man["image"] = man.apply(_img, axis=1)
meta = man[man.use_for_measurement.astype(str) == "True"].drop_duplicates("image").set_index("image")[["code", "group", "arm", "condition"]]
brt = pd.read_csv(P / "DerivedData/aris_v2_2026-08-18/three_arm_section_level.csv").drop_duplicates("image").set_index("image")["cy3_med"]
d = pd.DataFrame([json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(R / "*.json")))]).set_index("image")
d = d.join(meta).join(brt); d = d[d.code.notna() & ~d.index.isin(BAD)]
VARS = [("v4", "abeta_v4_um2"), ("gate 1.5x", "abeta_gate1.5_um2"), ("gate 2x", "abeta_gate2_um2"), ("gate 3x", "abeta_gate3_um2")]
out = []; W = lambda *a: out.append(" ".join(str(x) for x in a))
W(f"BRIGHTNESS GATE ON v4 GLOBAL, 2026-09-09  (n sections = {len(d)}; gate = Cy3 >= k x local background, 400-um window of non-Abeta tissue)\n")
W("WHOLE SECTION inside tissue ROI")
for name, col in VARS:
    pct = 100 * d[col] / d.tissue_um2; st, ct = pct[d.condition == "positive"], pct[d.condition != "positive"]
    W(f"  {name:9s} stained n={len(st)} median {st.median():6.2f}% [{st.min():.2f}-{st.max():.2f}]   controls n={len(ct)} median {ct.median():.3f}% max {ct.max():.3f}%   "
      f"stained median / worst control {st.median()/ct.max():4.1f}x   brightness Spearman stained {spearmanr(d.cy3_med[st.index], st).correlation:+.2f} controls {spearmanr(d.cy3_med[ct.index], ct).correlation:+.2f}")
    d[f"pct_{col}"] = pct
v4 = d["pct_abeta_v4_um2"]; g2 = d["pct_abeta_gate2_um2"]; st = d.condition == "positive"
W(f"  v4 -> gate 2x, stained sections: rank Spearman {spearmanr(v4[st], g2[st]).correlation:+.3f}, median ratio {(g2[st]/v4[st]).median():.2f} [{(g2[st]/v4[st]).min():.2f}-{(g2[st]/v4[st]).max():.2f}]")
W("\nHIPPOCAMPUS (manual regions), positives only, area-weighted per animal, medians and MWU p (descriptive, uncorrected)")
h = d[d.Hippocampus_um2.notna() & (d.condition == "positive")]
rows = []
for name, col in VARS:
    hc = col.replace("abeta_", "Hippocampus_abeta_")
    pa = h.groupby(["code", "group", "arm"]).apply(lambda g: 100 * g[hc].sum() / g.Hippocampus_um2.sum(), include_groups=False).rename("pct").reset_index()
    for arm in ["diet", "ip"]:
        g = pa[pa.arm == arm]; r = g[g.group.str.contains("Rapamycin")].pct; c = g[~g.group.str.contains("Rapamycin")].pct
        W(f"  {arm:4s} {name:9s} ctrl {c.median():6.2f} (n={len(c)})  rapa {r.median():6.2f} (n={len(r)})  ratio {r.median()/c.median():.3f}  p {mannwhitneyu(r, c).pvalue:.3f}")
    rows.append(pa.assign(method=name))
pd.concat(rows).drop(columns=["group"]).to_csv(OUT / "hippocampus_per_animal_gate.csv", index=False)
d.drop(columns=["group"]).to_csv(OUT / "section_level_gate.csv")
(OUT / "summary.txt").write_text("\n".join(out) + "\n"); print("\n".join(out))
