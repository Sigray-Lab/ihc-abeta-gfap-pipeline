# Hand-off to the GFAP session (written 2026-11-09)

Aβ is finished and locked. This note is everything the next session needs to start GFAP quantification without
rediscovering what the Aβ work already established. Paths: `<project>` is the Dropbox folder
`…/Rapamycin_CSF_cell_mice_Per/IHC_analysis_pipeline/`; the git repo is `<project>/pipeline/`; working files are in `~/ihc_work/`.

## 1. State

- **Aβ: done.** Two approaches are presented side by side, unnormalised and v4 GLOBAL, with neither called the better one
  (ADR-0040). Results: `<project>/DerivedData/abeta_final_2026-11-05/` (read its README before quoting a number).
  Figures: `<project>/Results/Figures/abeta_final_2026-11-05/`. Runbook: `docs/RUNBOOK_abeta.md`.
- **GFAP: nothing measured.** No classifier, no training annotations, no evaluation script. The plan and the training guide exist.
- **Everything else that was tried for Aβ is archived**, not deleted: `<project>/Archive/abeta_other_approaches_2026-11-05/README.md`
  is the index; git tag `abeta-all-approaches-2026-11-05` restores the repo with all of it in place.

## 2. Read these first, in this order

1. `docs/decisions.md`: ADR-0040 and ADR-0041 (the lock), ADR-0038 (session rule), ADR-0017 and ADR-0020 (GFAP-specific), ADR-0025/0029/0030 (denominator and QC rules).
2. `docs/classifier_gfap.md`: the step-by-step GFAP training guide written for the wet-lab scientist. Section 4 (where faint stops counting) and section 8 (negative-control check) carry the method.
3. `config/config.yaml`: block 6 `gfap`, block 7 `intensity` (still PENDING_PI_DECISION), block 8 `gfap_enrichment` (bands, exploratory, D-8 not approved), `classifier.resolution_downsample` (GFAP candidates 1 and 2, not locked), `classifier.input_channels` (GFAP sees FITC only).
4. `<project>/DerivedData/coverage_qc_2026-11-05/README.md` and `coverage_by_section.csv`: what every section is, and its focus flags.
5. `docs/EXTERNAL_REVIEW_ROUND2_ASSESSMENT_2026-09-09.md`: what two rounds of external review found; section 5 lists actions, several never applied.
6. `docs/RUNBOOK_abeta.md`: the commands; most transfer to GFAP unchanged.

## 3. How the PI wants to work

- Brief answers, plain language, simple tables; "explain like I'm 12" when asked. Show PNGs when a picture settles it.
- This run is exploratory and will be repeated. **Do not raise** blinding, pre-registration, backups or data-integrity items as things to act on. Give unblinded names (`Image_NN_sXX`) when asked, without treatment groups.
- Scope is image processing and quantification, not inferential statistics or study design.
- Region outlines are manual until an atlas step exists; do not propose atlas work.
- He decides and locks; present options with a recommendation, then wait. Do not let a method choice follow the treatment p-value.
- Aris (the wet-lab scientist) trains classifiers and judges biology; neither he nor the PI claims plaque or astrocyte expertise, so show images for anything that needs their eye.

## 4. Where everything is

| What | Where |
|---|---|
| Raw images, read-only | `<project>/RawData/` (37 VSI files; 130 fluorescence sections; never write here) |
| Key and manifest | `<project>/BlindingKey/provenance_manifest.csv` (one row per section; `use_for_measurement` is the only correct filter) |
| QuPath project, rebuilt 2026-11-09 | `~/ihc_work/qupath/qupath/project.qpproj`: the 121 measured sections; 123 Aβ classifiers; reference boxes on 6 sections |
| Tissue masks | `~/ihc_work/tissue_rois/<image>.geojson` (DAPI-derived; `tissue_rois_pairs/` for original+rescan pairs) |
| Region outlines | `~/ihc_work/region_rois/{Hippocampus,Isocortex}/<image>.geojson` (60 stained sections + 1 control) |
| Training annotations, Aβ | `~/ihc_work/annotations/annotations_training_v2/` (12 files, named `Image_NN_sXX`); hand-drawn region outlines in `Annotations_Hippo_Cortex_Tissue_v1/` |
| Saved Aβ measurement rows | `~/ihc_work/eval_20260909/` (copy in `DerivedData/aris_v4_2026-09-10_session/eval_inputs_backup_20261109/`) |
| Focus scores | `~/ihc_work/focus_20261105/` and `DerivedData/coverage_qc_2026-11-05/` |
| Frozen Aβ classifiers | `classifiers/abeta/` in the repo, with SHA-256 |
| Active scripts | `scripts/` (`run_eval_batch.sh`, `make_section_anchors.py`, `make_global_classifiers.py`, `make_tissue_rois.py`, `make_region_rois.py`, `register_rescan_regions.py`, `retrain*.groovy`, `qc/`) and `qupath/scripts/` |
| Archived tools that still run | Plaque Coder viewer: `archive/abeta_other_approaches_2026-11-05/scripts/manual_coding/plaque_coder.py` (data in `~/ihc_work/manual_coding/`); single-image comparison project `~/ihc_work/qupath_J09_s04/` |
| Archived working outputs | `~/ihc_work/archive_abeta_other_approaches_20261105/` (README.txt inside) |

Image names: `K07_s01` is animal code + section; a `b` after the code (`A09b_s03`) is an August rescan; `snew` marks a section first imaged in August.
Channels: 0 DAPI, 1 FITC = GFAP, 2 Cy3 = Aβ. 16-bit, 0.325 µm/px; pyramid level 2 = 1.3 µm/px (the Aβ working resolution), level 5 = tissue-mask level.

## 5. The Aβ numbers, for orientation

Hippocampus, diet arm, per mouse (control n 8, rapamycin n 11): unnormalised ratio of medians 0.39, p 0.051; v4 GLOBAL 0.54, p 0.016
(Mann-Whitney, uncorrected, exploratory). IP arm: no difference under either. Isocortex: not significant. The two approaches rank mice almost identically (Spearman 0.95).
Neither is validated against traced reference masks; on one section about half of v4's area in hippocampal boxes was near background.

## 6. Cohort facts GFAP will need

- 31 mice, four groups (control diet 8, rapamycin diet 12, control IP 6, rapamycin IP 5). 121 measured sections: 69 stained, 52 secondary-only controls. The controls lack the GFAP primary too, so they are valid GFAP negatives (ADR-0020: GFAP brightness is what separates stained from control sections).
- **FITC exposure:** 397.9 ms on 119 sections; 240.8 ms on two controls, A09_s01 and A09_s02 (their DAPI is also short). The Cy3 exposure problem additionally covered R22's controls, but R22's FITC is standard.
- **Sessions:** July originals for 114 sections; 7 stained sections in use come from the 5-August rescan session (A09b ×2, R22b ×2, V22b_snew01, L03b_snew02, B16b_snew01). Mice A09 and R22 have stained data only from August. For Cy3 the August scans read very differently (ADR-0038). For FITC the audit found a milder change (background ×1.2–1.35, tails close to proportional), but measure the GFAP session effect on the original+rescan pairs before mixing sessions.
- **Failed tissue masks:** Y13_s01 (DAPI blown out) and S17_s04. Excluded everywhere.
- **Focus:** blur is common and patchy. 12 of 58 measured hippocampi are blurred; worst sections Z23_s02, F15 (all four), J09_s04, S17_s03, V07_s01, F02_s01. Blurred sections read lower for Aβ. GFAP processes are thin, so expect blur to matter more. The focus score is nuclear edge width from DAPI (`scripts/qc/focus_scores.py`); crisp is about 2 µm, blurred above 2.9 µm.
- **Regions:** hand-drawn, not atlas-matched; hippocampal outline area varies about tenfold between sections and about 30 % between diet groups. Only one control has outlines.
- **Age:** diet controls are 24.4 days older than diet rapamycin mice. The IP arm's ages were not found in the project files; ask the PI.

## 7. What exists for GFAP

- **Endpoint (config):** GFAP-immunoreactive percent area, explicitly not astrocyte number or activation state. Report absolute GFAP area per region alongside any enrichment contrast.
- **Local enrichment around Aβ (config block 8, exploratory, not PI-approved):** bands 0–25, 25–50 and >50 µm from the Aβ mask plus the Aβ core as its own compartment; primary contrast near minus distant. It needs an Aβ mask, and there are now two (unnormalised and v4): which mask defines the bands is an open decision.
- **Accepted fallback:** whole-region FITC intensity, agreed in advance as a legitimate standard method if a classifier proves unworkable (`<project>/EXECUTION_PLAN_v4.md`).
- **Known signal facts:** GFAP is strongly elevated around plaques (up to 2.85× near deposits) and this is biology, not bleed-through (ADR-0017). Far from deposits, Cy3 in GFAP-bright tissue is only 1.07–1.32× GFAP-dim tissue. On antibody-free sections, autofluorescent cells glow in FITC and DAPI as well as Cy3.
- **Reusable as is:** the QuPath project, tissue masks, region outlines, `run_eval_batch.sh` (now fails loudly), `make_section_anchors.py --channel fitc`, `make_global_classifiers.py`, focus and coverage QC, the table and figure builders as templates.
- **Needs writing:** a GFAP trainer (adapt `scripts/retrain.groovy`: FITC channel, GFAP class names) and a GFAP evaluation script (`qupath/scripts/eval_in_tissue.groovy` looks up measurements named `Abeta area` and `Negative area`; it needs the GFAP class name, ideally as an argument).

## 8. Lessons from Aβ that apply directly

1. **Define the boundary before training.** For Aβ, where the outline stops decided the number, and the classifier was generous at the edges. GFAP is a gradient with no edge, so this is the whole problem. Follow section 4 of the guide; train with annotations that mark edges and matched nearby negatives, not short strokes through bright centres (the Aβ positive class was only about 950 pixels from 87 strokes on 8 mice, which is the root of most of what followed).
2. **Judge specificity on the controls, at standard exposure.** Report the false-positive floor from the antibody-free sections next to every stained number. A very dim or short-exposure section breaks any divide-by-brightness scheme.
3. **Brightness normalisation: do not argue from stained sections.** Correlating the marker channel's brightness with the marker's area is circular. Use antibody-free sections and same-tissue rescan pairs. For Aβ the evidence was mixed and the decision was to present normalised and unnormalised together; expect the same question for GFAP and consider the same answer. If normalising, compute FITC anchors (never reuse the Cy3 table) and check that every glass value is plausible: four values in the Cy3 production table were wrong (ADR-0041).
4. **Do not train against control false positives in the same channel.** Doing so for Aβ (v5) removed the false positives only by removing small real deposits. A second channel separates them better: autofluorescence is broad-spectrum, antibody signal is not.
5. **Post-hoc gates have a free parameter.** The Aβ brightness gate moved the treatment ratio between 0.55 and 0.70 depending on its setting. Any threshold must be fixed before looking at group results.
6. **Check what the classifier was actually applied to.** The project once held the wrong scans for three mice by a default nobody chose; six sections were never classified; four failed controls were counted as zeros. `scripts/qc/coverage_qc.py` builds the truth from files on disk; extend it for GFAP.
7. **Validate every new metric against cases with a known answer** before using it. Two focus metrics and one rasteriser looked reasonable and were wrong; each was caught only by a direct test (known blur added to sharp tiles; areas compared with QuPath's own).
8. **QuPath operations.** `QuPath script -p <project>` without `-i <image>` runs once per image. Never run headless QuPath on a project that is open in the GUI, and never several QuPath processes beside an open window: they crash each other and the GUI overwrites headless saves. The import script refuses to add to an existing project: move the old one aside and rebuild. Pixel classifiers cannot compute per-image constants, hence one `gx_<image>` file per section. `createTrainingData` per image numbers classes differently per image: use the canonical label map (ADR-0032). Set `JAVA_TOOL_OPTIONS=-Djava.awt.headless=true` for scripted runs on macOS.
9. **Reading the images.** `read_pyramid_plane` returns NaN for never-acquired tiles; use nan-aware statistics and pass `true_size`. Rasterise polygons at pixel centres (PIL and OpenCV fills overstate small objects by up to a third). The slide label and overview images inside each VSI are brightfield, not data. Dropbox can leave files as empty placeholders; `./ihc doctor` checks.
10. **Keep results reproducible while iterating.** Pin input file lists (a glob once swallowed later files), assert headline values in the table builder, tag before restructuring, and keep superseded outputs with a pointer instead of overwriting.

## 9. Decisions the PI still has to make for GFAP

- The GFAP positive-class boundary: how faint still counts (guide section 4), and whether thin processes are in.
- Classifier or the whole-region intensity fallback, or both side by side.
- Normalised, unnormalised, or both, mirroring Aβ.
- Working resolution for GFAP (candidates: native and 0.65 µm/px).
- Which Aβ mask defines the enrichment bands, and whether the enrichment endpoint is in scope now.
- Handling of blurred sections and of the two August-only mice for GFAP.
- The `intensity` block in the config (population, background adjustment, saturation, empty-mask value).

## 10. Known issues carried forward

- Neither Aβ approach is validated against traced reference masks; 38 reference boxes exist on 6 sections (class `RefBox` in the project) and were never traced.
- `make_global_classifiers.py` has no minimum-gain check although `docs/normalisation.md` describes one.
- Table builders in `DerivedData/` are pinned to this cohort and sit outside the repo.
- Text corrections listed in the round-2 assessment (older ADRs, including the withdrawn "no photobleaching" statement) were not all applied.
- The IP arm has no age record in the project files.

## 11. Suggested first steps

1. Run `./ihc doctor`, open the project once, and look at FITC on three stained sections and two controls at full resolution, including one blurred section.
2. Measure the cheap baselines before training anything: whole-region FITC intensity per section (the agreed fallback), on stained and control sections, with the session and focus flags beside it.
3. Agree the boundary rule with the PI and Aris on a few boxes, then have Aris annotate following `docs/classifier_gfap.md`.
4. Train, then run the negative-control check (guide section 8) and the rescan-pair check before any group comparison.
5. Record each decision as an ADR as it is made.
