# Aβ runbook: reproduce the results, or quantify new images

Two approaches are run and presented side by side (ADR-0040): **unnormalised** (`v2_nonorm_s1`) and **v4 GLOBAL**
(`gx_<image>`, the frozen v4 forest with one offset and gain per section). Both classifiers are frozen in `classifiers/abeta/`.
Paths below assume `config/paths.yaml` as shipped: raw data and the custodian folder in the project folder, working files in `~/ihc_work`.

## A. Reproduce the published tables and figures

No image processing needed; everything is rebuilt from saved measurement rows.

```bash
cd <project>/DerivedData/aris_v4_2026-09-10_session && python3 build_tables.py      # all four arms, from the saved QuPath rows
cd ../abeta_final_2026-11-05 && python3 build_final_tables.py                        # the two presented approaches + flags
cd <project>/Results/Figures/abeta_final_2026-11-05 && python3 make_group_difference.py
```

`build_final_tables.py` asserts the headline values (hippocampus diet p 0.051 unnormalised, 0.016 v4 GLOBAL), so a drift fails loudly.
The saved QuPath rows are in `~/ihc_work/eval_20260909/`, with a copy in `aris_v4_2026-09-10_session/eval_inputs_backup_20261109/`.

## B. Re-measure from the raw images

```bash
./ihc doctor && ./ihc verify                                   # environment and raw containers
./ihc manifest && ./ihc blind --seed <recorded seed> --force   # manifest; same seed keeps the same codes
python3 scripts/build_qupath_project.py                        # project spec + coded symlinks (move any old project aside first)
QuPath script --args ~/ihc_work/qupath/project_spec.json qupath/scripts/import_blinded_project.groovy
python3 scripts/make_tissue_rois.py                            # DAPI tissue masks as GeoJSON, ~/ihc_work/tissue_rois
python3 scripts/make_section_anchors.py --out anchors.csv      # glass and tissue medians per section (Cy3); see the note below
cp classifiers/abeta/*.json ~/ihc_work/qupath/qupath/classifiers/pixel_classifiers/
python3 scripts/make_global_classifiers.py --classifier v4_global --anchors anchors.csv --prefix gx
scripts/run_eval_batch.sh 4 ~/ihc_work/tissue_rois 'gx_{image}'  out_v4.csv       # v4 GLOBAL, whole section
scripts/run_eval_batch.sh 4 ~/ihc_work/tissue_rois 'v2_nonorm_s1' out_nonorm.csv  # unnormalised, whole section
scripts/run_eval_batch.sh 4 ~/ihc_work/region_rois/Hippocampus 'gx_{image}' out_hip_v4.csv Hippocampus   # per region
```

**Anchors.** To reproduce the published v4 rows exactly, use the saved table `DerivedData/abeta_final_2026-11-05/global_norm_anchors.csv`
in place of a fresh `anchors.csv`. The script reproduces it within 5 % of gain on 110 of 121 sections; four production glass values
are too low (ADR-0041). Checked 2026-11-09: with the saved anchors the rebuilt project returns values identical to the saved rows.

`run_eval_batch.sh` exits non-zero and keeps its work folder if any section fails or a row is missing. Close the QuPath window
first: several QuPath processes on one machine crash each other. `QuPath script -p <project>` without `-i <image>` runs the
script once per image.

## C. Quantify new images

1. Put the new `.vsi` files and their `_Image_NN_` folders in `RawData/`, and add the slides to `config/slides.csv` and the wet-lab slide record.
2. Run section B from the top. New sections get new coded names; existing codes are unchanged with the same seed.
3. Region measurements need outlines: `scripts/make_region_rois.py` turns delivered hippocampus and isocortex outlines into ROIs
   (`register_rescan_regions.py` carries them to a rescan). Outlines are drawn by hand; there is no atlas step.
4. Run `scripts/qc/focus_scores.py` and `scripts/qc/coverage_qc.py` and read the flags before trusting a number: blurred
   sections read lower, and a section imaged in a different session or at a different exposure is not comparable (ADR-0038).
5. The table builders in `DerivedData/` are pinned to this cohort's files. Copy them for a new cohort and point them at the new rows.

**Limits to keep in mind.** The classifiers were trained on this cohort's staining and imaging run; a new staining batch is a new
validation question, not a rerun. Neither approach is validated against traced reference masks. Standard exposure here is
Cy3 1840 ms; `make_global_classifiers.py` has no minimum-gain check, so a short-exposure or near-black section gets a
meaningless v4 result (the four short-exposure controls of this cohort are the example).
