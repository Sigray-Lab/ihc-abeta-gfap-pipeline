#!/usr/bin/env python3
"""Coverage QC: for every fluorescence section in the manifest, what has actually been done to it. One row per section,
built only from files on disk (ROIs, classifiers, measurement rows), never from what the docs say should exist.
Output: DerivedData/coverage_qc_2026-11-05/coverage_by_section.csv (no group or tube columns)."""
import json, sys
from pathlib import Path
import pandas as pd
R = Path(__file__).resolve().parents[3]; H = Path.home(); DD = R / "DerivedData"; EV = H / "ihc_work/eval_20260909"; SW = DD / "session_swap_2026-09-10"; OUT = DD / "coverage_qc_2026-11-05"; OUT.mkdir(exist_ok=True)
man = pd.read_csv(R / "BlindingKey/provenance_manifest.csv"); man = man[man.row_kind == "section"].copy()
man["image"] = man.apply(lambda r: f"{r.code}{'b' if r.scan == 'rescan' else ''}_s{str(r.physical_section_label).zfill(2) if str(r.physical_section_label).isdigit() else r.physical_section_label}", axis=1)
assert man.image.is_unique, "image name is not a unique key"
d = man.set_index("image")[["scan", "condition", "exposure_Cy3_ms", "exposure_is_standard", "use_for_measurement", "bench_action", "flagged_lower_quality", "has_flagged_artefact"]].copy()
d["in_use"] = d.use_for_measurement.astype(str) == "True"; d = d.drop(columns="use_for_measurement")
ids = lambda f, col="image": set(pd.read_csv(f)[col]) if Path(f).exists() else set()
def ids_for(f, clf):
    t = pd.read_csv(f); return set(t[t.classifier.astype(str).str.fullmatch(clf)].image)
d["in_qupath_project"] = d.index.isin([e["image_name"] for e in json.load(open(H / "ihc_work/qupath/project_spec.json"))["images"]])
d["tissue_roi"] = [(H / f"ihc_work/tissue_rois/{i}.geojson").exists() for i in d.index]
d["hip_outline"] = [(H / f"ihc_work/region_rois/Hippocampus/{i}.geojson").exists() for i in d.index]
d["iso_outline"] = [(H / f"ihc_work/region_rois/Isocortex/{i}.geojson").exists() for i in d.index]
d["anchors"] = d.index.isin(ids(DD / "aris_v4_2026-09-10_session/global_norm_anchors.csv"))
CL = H / "ihc_work/qupath/qupath/classifiers/pixel_classifiers"; d["gx_classifier"] = [(CL / f"gx_{i}.json").exists() for i in d.index]
aug = DD / "aris_v2_2026-08-18/full_cohort_in_tissue_RAW.csv"
d["nonorm_whole"] = d.index.isin(ids_for(aug, "v2_nonorm_s1") | ids(SW / "tissue_nonorm.csv"))
d["z130_whole"] = d.index.isin(ids_for(aug, "v2_z130_s1") | ids(SW / "tissue_z130.csv"))
d["v4_whole"] = d.index.isin(ids(EV / "tissue_GLOBAL_v4.csv") | ids(SW / "tissue_GLOBAL_v4_july_originals.csv"))
d["v5_whole"] = d.index.isin(ids(EV / "tissue_GLOBAL_v5.csv"))
d["gate_whole"] = d.index.isin(ids(DD / "brightness_gate_2026-09-09/section_level_gate.csv") | ids(H / "ihc_work/gate_20260910_session/section_level_gate.csv"))
d["v6_whole"] = d.index.isin(ids(DD / "v6_dapi_veto_2026-09-10/section_level_v6.csv"))
d["v4_hip"] = d.index.isin(ids(EV / "region_Hippocampus_GLOBAL_v4.csv") | ids(SW / "region_Hippocampus_GLOBAL_v4_july_originals.csv"))
d["v4_iso"] = d.index.isin(ids(EV / "region_Isocortex_GLOBAL_v4.csv") | ids(SW / "region_Isocortex_GLOBAL_v4_july_originals.csv"))
sec = pd.read_csv(DD / "aris_v4_2026-09-10_session/section_level_in_tissue.csv"); g = sec[sec.arm_method == "GLOBAL"].set_index("image")
d["in_working_table"] = d.index.isin(g.index); d["valid_in_table"] = d.index.map(g.valid.to_dict()).fillna(False).astype(bool)
reg = pd.read_csv(DD / "aris_v4_2026-09-10_session/regional_section_level.csv"); rg = reg[reg.arm_method == "GLOBAL"]
d["in_hip_table"] = d.index.isin(rg[rg.region == "Hippocampus"].image); d["in_iso_table"] = d.index.isin(rg[rg.region == "Isocortex"].image)
d["v4_pct"] = d.index.map((g.pct).to_dict()); d["v4_coverage"] = d.index.map(g.coverage.to_dict())
tr = set()
for f in ("training_annotation_map.csv", "training_annotation_map_v5.csv"):
    for l in (R / "BlindingKey" / f).read_text().splitlines():
        if "," in l and not l.startswith("#"): tr.add(l.split(",")[1].strip())
d["hand_annotated_for_training"] = d.index.isin(tr) & ~d.index.isin(["W12_s03", "J21_s02", "J01_s01", "G22_s02", "D19_s02", "G09_s04"])
f = H / "ihc_work/focus_20261105/focus_sections.csv"
if f.exists():
    fo = pd.read_csv(f).set_index("image"); d = d.join(fo[["n_usable", "edge_um", "edge_p10", "edge_p90", "cy3_edge_um"]])
    tl = pd.read_csv(f.with_name("focus_tiles.csv")); tl = tl[tl.dapi_snr >= 8]; d["blurred_tile_share"] = tl.groupby("image").dapi.apply(lambda x: float((x > 2.9).mean()))
    rf = pd.read_csv(f.with_name("focus_regions.csv")).set_index("image"); d = d.join(rf[["Hippocampus_edge_um", "Isocortex_edge_um"]])
d.sort_values(["in_use", "condition", "scan"], ascending=[False, False, True]).to_csv(OUT / "coverage_by_section.csv")
u = d[d.in_use]; n = d[~d.in_use]; st = u[u.condition == "positive"]; ct = u[u.condition != "positive"]
print(f"fluorescence sections {len(d)}: in use {len(u)} (stained {len(st)}, controls {len(ct)}), not in use {len(n)}\n")
print("not in use:"); print(n[["scan", "condition", "exposure_Cy3_ms", "bench_action", "nonorm_whole", "z130_whole", "v4_whole", "v5_whole", "gate_whole", "v6_whole", "in_working_table"]].to_string())
cols = ["in_qupath_project", "tissue_roi", "anchors", "gx_classifier", "nonorm_whole", "z130_whole", "v4_whole", "v5_whole", "gate_whole", "v6_whole", "in_working_table", "valid_in_table", "hip_outline", "iso_outline", "v4_hip", "v4_iso", "in_hip_table", "in_iso_table"]
print("\ncoverage among sections in use:"); print(pd.DataFrame({"stained": st[cols].sum(), "controls": ct[cols].sum(), "all": u[cols].sum()}).astype(int).to_string())
print("\nin use but missing:")
for c in cols[:12]:
    m = u.index[~u[c]].tolist()
    if m: print(f"  {c}: {len(m)} -> {m[:12]}")
