import sys, json, glob
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
P = Path("/Users/pontusps/Dropbox/longevity_project/projects/Rapamycin_CSF_cell_mice_Per/IHC_analysis_pipeline"); R = Path(sys.argv[1]); OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
BAD = {"Y13_s01", "S17_s04"}
man = pd.read_csv(P/"BlindingKey/provenance_manifest.csv"); man = man[man.row_kind=="section"].copy()
man["image"] = man.apply(lambda r: f"{r.code}{'b' if r.scan=='rescan' else ''}_s{str(r.physical_section_label).zfill(2) if str(r.physical_section_label).isdigit() else r.physical_section_label}", axis=1)
meta = man.drop_duplicates("image").set_index("image")[["code","group","arm","condition","use_for_measurement","exposure_is_standard"]]
d = pd.DataFrame([json.loads(Path(f).read_text()) for f in sorted(glob.glob(str(R/"*.json")))]).set_index("image").join(meta)
d = d[d.code.notna() & ~d.index.isin(BAD)]   # note: evaluated on the 121-section set as exported (pre-swap)
VARS = [("v4","abeta_v4_um2"),("v6 (DAPI veto)","abeta_v6_um2"),("gate 2x","abeta_gate2_um2"),("v6 + gate 2x","abeta_v6_gate2_um2")]
out = []; W = lambda *a: out.append(" ".join(str(x) for x in a))
W(f"v6 = v4 minus DAPI veto (DAPI >= 1.5x local background), 2026-09-10; n sections = {len(d)} (section set as exported before the session swap)\n\nWHOLE SECTION inside tissue ROI")
for name, col in VARS:
    pct = 100*d[col]/d.tissue_um2; st = pct[d.condition=="positive"]; ct = pct[d.condition!="positive"]; cts = pct[(d.condition!="positive") & (d.exposure_is_standard.astype(str)=="True")]
    W(f"  {name:15s} stained n={len(st)} median {st.median():5.2f}% [{st.min():.2f}-{st.max():.2f}]   controls n={len(ct)} median {ct.median():.3f}% max {ct.max():.3f}% (standard exposure n={len(cts)} max {cts.max():.3f}%)   stained median / worst control {st.median()/ct.max():4.1f}x")
    d[f"pct_{col}"] = pct
W(f"  vetoed share of v4 area: stained median {100*d[d.condition=='positive'].vetoed_share.median():.1f}%, controls median {100*d[d.condition!='positive'].vetoed_share.median():.1f}%")
W("\nHIPPOCAMPUS, positives only, area-weighted per animal, medians and MWU p")
h = d[d.Hippocampus_um2.notna() & (d.condition=="positive")]
for name, col in VARS:
    hc = col.replace("abeta_", "Hippocampus_abeta_")
    pa = h.groupby(["code","group","arm"]).apply(lambda g: 100*g[hc].sum()/g.Hippocampus_um2.sum(), include_groups=False).rename("pct").reset_index()
    for arm in ["diet","ip"]:
        g = pa[pa.arm==arm]; r = g[g.group.str.contains("Rapamycin")].pct; c = g[~g.group.str.contains("Rapamycin")].pct
        W(f"  {arm:4s} {name:15s} ctrl {c.median():6.2f} (n={len(c)})  rapa {r.median():6.2f} (n={len(r)})  ratio {r.median()/c.median():.3f}  p {mannwhitneyu(r,c).pvalue:.3f}")
d.drop(columns=["group"]).to_csv(OUT/"section_level_v6.csv"); (OUT/"summary.txt").write_text("\n".join(out)+"\n"); print("\n".join(out))
