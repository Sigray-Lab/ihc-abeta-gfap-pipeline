# Plan: v6 — a second detector for the negatives (Aris's proposal, 2026-09-09), and how every method is compared

**Aris's idea (paraphrased).** Keep the plaque classifier trained on stained slides only. Train a second classifier on
the antibody-free control slides to recognise what lights up there (autofluorescence, rims, debris), and let it veto
the first. Do not fold control annotations into the plaque classifier itself; ADR-0036 (v5) showed that this removes
small deposits along with the false positives.

**Why the veto must not look at Cy3 alone.** On a single Cy3 channel an autofluorescent cell body and a small deposit
are the same object; a veto trained on Cy3 texture from control slides re-creates v5. The information that separates
them is spectral: lipofuscin-type autofluorescence glows in every channel, 82E1 signal only in Cy3. Measured on v4's
objects (`~/ihc_work/manual_qc_20260909/v4_objects_channel_contrast.csv`, 2026-09-10):

| v4 objects on | DAPI-bright (>1.5× local) share of area | FITC-bright share of area |
|---|---|---|
| 3 antibody-free controls | 29–33 % | 21–38 % |
| 3 stained sections | 2–5 % | 23–63 % |

FITC is unusable as a veto: reactive astrocytes (GFAP) surround real plaques, so a FITC veto would remove up to 63 % of
genuine plaque area. DAPI separates cleanly: ~30 % of control false-positive area, 2–5 % of stained area.

**v6 as implemented first (deterministic, no training).** v6 = v4 Aβ pixels minus pixels where DAPI ≥ 1.5 × its local
background (400-µm window of non-Aβ tissue, same machinery as the brightness gate). Script:
`scripts/brightness_gate/veto.py` + `run_veto.py`. Also reported: v6 combined with the 2× Cy3 gate. A trained DAPI-only
veto classifier (RTrees on DAPI texture features, trained on the six machine-drawn control files) is the next step if
the rule shows promise; it is not built yet.

**Expected ceiling.** Because only ~30 % of control false-positive area is DAPI-bright, v6 alone can at best cut the
control ceiling by about a third. The remaining control false positives are Cy3-only (dim texture, gain-amplified
noise, rims that are not DAPI-bright) and need the gate, the gain/exposure rule, or rim erosion.

## The reference set every method is compared against

| Method | What it is | Where |
|---|---|---|
| v4 GLOBAL | working classifier, one offset/gain per section | ADR-0032, `DerivedData/aris_v3_2026-09-09_globalfix/` |
| v5 | v4 annotations + six machine-drawn control files, capped 700 samples each | ADR-0036, `DerivedData/v5_capped_negatives_2026-09-09/` |
| v4 + gate k× | v4 pixels kept where Cy3 ≥ k × local background, k = 1.5/2/3 | ADR-0037, `DerivedData/brightness_gate_2026-09-09/` |
| v6 (DAPI veto) | v4 pixels minus DAPI-bright pixels; also v6 + gate 2× | this plan, `DerivedData/v6_dapi_veto_2026-09-10/` |

Compared on the same three things, always: (1) whole-section table — stained median, worst control (all exposures and
standard exposure), stained-to-control ratio, brightness Spearman; (2) hippocampus per-animal diet and IP ratios;
(3) the same reference boxes rendered side by side (raw | v4 | v5 | gate | v6), starting with T18_s03 HIP1.

**How the choice will be made (PI, 2026-09-10).** Not on the treatment p-value. The PI and the wet-lab scientist look at
a random sample of sections and boxes with all methods side by side and agree on which outline matches what they call
a plaque; that agreed definition (D-3) is then frozen and, where possible, checked against traced reference boxes.
The random sample should include bright rescans, dim originals, antibody-free controls and tissue rims, which the
current 38-box panel does not cover.
