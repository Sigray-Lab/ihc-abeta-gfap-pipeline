# Independent audit of the reviewer's same-tissue rescan experiment (REVIEW_RESPONSE_ASSESSMENT §2.3 and §4)

Auditor: image-analysis audit, 2026-09-09. Read-only on all project trees; QuPath not run.
Working files (this directory): task1_pair_table.csv, task2_anchor_comparison.csv, task3_project_pairs_fullframe.csv,
task3_commonroi_vs_wholesection.csv, task4_pixel_qq_compact.csv, qq_pairs_summary.csv, qq_pairs_quantiles.csv,
qq_pairs.py (my own pixel re-read of the 18 acquisitions at 2.6 um/px inside the reviewer's common ROIs).

## 0. Bottom line

The reviewer's numbers are correct to the last digit, the geometry and classifiers are what they say they are, and the
conclusion ("the corrected model still changes across sessions; mixing acquisitions is not automatically harmless") is
sound and, if anything, understated. My own pixel-level re-read shows why an affine-invariant normaliser fails here:
the July->August change is not affine. Cy3 tissue background rises 2.0-2.2x at byte-identical exposure while the plaque
band (p99-p99.9) rises only 1.0-1.3x; in normalised units the plaque quantiles fall to 0.54-0.60x and the fraction of
pixels above a fixed normalised threshold falls to 0.30-0.56x. The forest's 0.37-0.59 simply reports that. The project's
own August notes already contain the same observation ("compressive background pedestal ... not a lamp or camera gain
change", docs/PLAN_after_vacation_2026-08-19.md) and the same design fact ("5 animals sit wholly on a session that reads
1.7-2.1x brighter ... zero mixed ... no model carries a session term"), yet ADR-0034 and decisions.md still carry
"No photobleaching between the original scans and the rescans" and no session rule. That is the real gap.

## 1. Task 1 - reproducing the pair table

Source: common_roi/measurements.csv (27 rows; matches common_roi_eval_success.log line for line), registration.csv,
pair_results.csv; exposures from BlindingKey/provenance_manifest.csv.

| Pair | Original -> rescan | Cy3 exp orig/resc (ms) | Same exp | Orig GLOBAL % | Rescan GLOBAL % | Ratio | Delta pp | Common mm2 | ROI area ratio resc/orig | Fit Dice | Coverage orig/resc | Usable |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| 1 | V22_s01 -> V22b_s01 | 1840/1840 | yes | 5.065 | 3.004 | 0.593 | -2.06 | 19.589 | 0.99998 | 0.9976 | 0.977/0.980 | yes |
| 2 | L03_s03 -> L03b_s03 | 1840/1840 | yes | 6.025 | 2.631 | 0.437 | -3.39 | 20.150 | 1.00002 | 0.9930 | 0.954/0.955 | yes |
| 3 | Y13_s01 -> Y13b_s01 | 1840/1840 | yes | 1.009 | 0.000 | 0.000 | -1.01 | 0.119 | 0.99998 | 0.0239 | 0.992/1.000 | NO |
| 4 | Y13_s04 -> Y13b_s04 | 1840/1840 | yes | 15.303 | 7.929 | 0.518 | -7.37 | 22.124 | 1.00000 | 0.9974 | 0.973/0.975 | yes |
| 5 | A09_s03 -> A09b_s03 | 146/1840 | no | 32.344 | 12.948 | 0.400 | -19.40 | 20.455 | 1.00002 | 0.9960 | 0.977/0.989 | yes |
| 6 | A09_s04 -> A09b_s04 | 146/1840 | no | 18.329 | 9.050 | 0.494 | -9.28 | 20.923 | 1.00000 | 0.9964 | 0.937/0.960 | yes |
| 7 | B16_s01 -> B16b_s01 | 1840/1840 | yes | 16.330 | 6.006 | 0.368 | -10.32 | 22.699 | 0.99999 | 0.9956 | 0.986/0.988 | yes |
| 8 | R22_s03 -> R22b_s03 | 398/1840 | no | 7.658 | 1.353 | 0.177 | -6.31 | 22.228 | 1.00002 | 0.9433 | 0.981/0.985 | yes |
| 9 | R22_s04 -> R22b_s04 | 398/1840 | no | 14.266 | 2.831 | 0.198 | -11.43 | 16.486 | 1.00005 | 0.9867 | 0.962/0.964 | yes |

Checks: percentages recomputed from abeta_um2/tissue_um2 agree with pct_tissue to 1e-6 on all 27 rows; ratios and
deltas agree with pair_results.csv to 4 decimals. Same-exposure usable pairs (1,2,4,7): 0.3678-0.5931, median 0.4774,
delta -2.06..-10.32 pp, i.e. "41-63% lower" - exactly as reported. All eight usable: 0.1767-0.5931 (median 0.42, geometric
mean 0.37). Exposure flags agree with the provenance manifest (A09 originals also have DAPI 60.5 vs 128.5 ms). Both
members of every pair are measured on the same physical domain: ROI areas agree to <=5e-5, and the rigid transform
(angle, dy, dx, Dice) is identical to the project's own ~/ihc_work/region_rois/rescan_registration.csv for all seven
pairs the project registered. The 50-um erosion removes 4.6-5.7% of the smaller DAPI silhouette; 95% of the section
remains. Pair 3 is rightly excluded: the July Y13_s01 DAPI channel is blown out (DAPI tissue median 12,026 vs ~75
elsewhere; auto mask 0.22 mm2), so the "common domain" is a 0.119 mm2 sliver. (Side finding: that same broken mask feeds
the production anchor for Y13_s01 - global_norm_anchors.csv tissue 0.23 mm2, t_p50 2233.5 - and tissue_GLOBAL_v4.csv
carries Y13_s01 with a 0.22 mm2 ROI and qc_status=OK; it is only removed by a hard-coded BAD list in build_tables.py.)
Classifier check: all 18 reviewer classifiers are byte-identical to v4_global.json once the two constants are masked, the
constants equal anchors.csv, and the nine production gx_*.json are the same forest with global_norm_anchors.csv
constants. "The forest is the exact corrected v4 forest" is verified.

## 2. Task 2 - uniform vs production anchors

Production anchors exist only for the 121 selected images (111 originals + 10 rescans); nine of the 18 pair images have
one. No repository script produces whole-section glass_p50/t_p50 (only consumers: make_global_classifiers.py;
make_region_anchors.py reads the August CSV), so the reviewer's statement is correct.

| Image | Uniform glass/gain | Production glass/gain | dGlass | gain ratio uni/prod | GLOBAL% uniform | GLOBAL% production | dpp | rel% |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| V22b_s01 | 62/467 | 59/471 | +3 | 0.992 | 3.004 | 2.972 | -0.032 | -1.1 |
| L03b_s03 | 70/505 | 63/512 | +7 | 0.986 | 2.631 | 2.583 | -0.048 | -1.8 |
| Y13_s01 (0.23 mm2 mask) | 42/2257 | 60/2173.5 | -18 | 1.038 | 1.009 | 1.159 | +0.150 | +14.9 |
| Y13_s04 | 44/218 | 42/223 | +2 | 0.978 | 15.303 | 14.982 | -0.321 | -2.1 |
| A09b_s03 | 54/412 | 53/417 | +1 | 0.988 | 12.948 | 12.760 | -0.188 | -1.5 |
| A09b_s04 | 59/491 | 56/498 | +3 | 0.986 | 9.050 | 8.875 | -0.174 | -1.9 |
| B16b_s01 | 58/462 | 57/465 | +1 | 0.994 | 6.006 | 5.945 | -0.060 | -1.0 |
| R22b_s03 | 54/242 | 51/246 | +3 | 0.984 | 1.353 | 1.301 | -0.052 | -3.9 |
| R22b_s04 | 71/289 | 69/292 | +2 | 0.990 | 2.831 | 2.764 | -0.068 | -2.4 |

The uniform protocol reproduces production anchors within 1-7 counts of glass and 0.6-2.2% of gain (Y13_s01 excepted,
where both are garbage). "At most 0.32 percentage points" is exactly right (Y13_s04); relative effect <=3.9% on the eight
real sections - an order of magnitude below the session effect. Supported. Cross-check: the reviewer's common-ROI
values with production anchors are within 1-8% of the project's own whole-section v4 values for the same images
(B16b_s01 5.95 vs 6.20; R22b_s03 1.30 vs 1.42; Y13_s04 14.98 vs 15.02), so the eroded domain is representative.
Caveat the reviewer under-states: the four changed-exposure originals have uniform gains of 10, 12, 22, 26 counts
against a cohort range of 201-512; the normalised image is quantised in 0.04-0.10 unit steps and noise is amplified
20-50x. docs/normalisation.md promises a "configurable floor" on gain; make_global_classifiers.py implements none
(only gain<=0 is refused). Pairs 5,6,8,9 test GLOBAL far outside its operating domain and confound a 4.6-12.6x exposure
change with the session change; the headline must rest on pairs 1,2,4,7 and "0.18-0.59 over all eight" should be
labelled a stress test.

## 3. Task 3 - what the project's own data say

Where pairs exist: nowhere for GLOBAL/REGIONNORM. The seven replaced July originals have no production anchors, no
gx_ classifiers and no region outlines measured on them (outlines drawn on originals were transferred to and measured
on the rescans only), so the regional *_rescan_nonorm/z130 rows (6 hippocampus, 7 isocortex rescans) have no original
counterpart. The only paired evidence in the project is the August full-frame table
(DerivedData/aris_v2_2026-08-18/pairs_results_nonorm.csv, pairs_results_z130.csv; README: not the ADR-0029 endpoint).
Ratios use abeta_pct_of_nonignored; frame-fraction ratios are within +-5% (task3_project_pairs_fullframe.csv).

| Arm | Same-exposure pairs 1,2,4,7 rescan/original | median | Changed-exposure pairs 5,6,8,9 | Direction |
|---|---|---:|---|---|
| unnormalised production v2_nonorm_s1 | 1.19-1.59 | 1.29 | 0.24-8.2 (original coverage 0.000-0.003: abstained on dark images) | UP |
| unnormalised variants noIgOnly / noCtrl / v1files / abeta_rt_nonorm | 1.22-1.46 / 1.06-1.46 / 2.5-3.4 / 6.6-10.4 | - | meaningless | UP |
| local z130 production v2_z130_s1 | 0.37-0.74 | 0.64 | 0.31-1.27 | DOWN |
| GLOBAL v4 (reviewer, common ROI) | 0.37-0.59 | 0.48 | 0.18-0.49 | DOWN |
| GLOBAL/REGIONNORM (project) | no pairs | - | - | - |

Reading: the project's own August numbers already showed no arm is session-stable - unnormalised up 1.2-1.6x, z130
down 0.4-0.7x, corrected GLOBAL down 0.4-0.6x (same direction and similar size as z130). The 19-August plan computed an
in-tissue "z130 session ratio 0.813", called the session "a genuine design defect", and judged it "real but inert" on
the strength of the then-defective GLOBAL. ADR-0034 is silent on it; decisions.md still records "No photobleaching
between the original scans and the rescans" - a finding based on tissue-median counts per millisecond (background) on
the two slides re-acquired for exposure, contradicted by every pair in section 4.
Unpaired cohort view (tissue_GLOBAL_v4.csv, positives): 59 July originals median 6.34% (IQR 3.24-10.45); 10 August
rescans median 2.93% (IQR 2.67-5.85). Per animal the five rescan-only animals rank 2, 3, 6, 15, 24 of 31 (three of five
among the six lowest). Confounded with biology, but exactly what a ~0.5x technical depression produces.

## 4. Task 4 - design and cause

Common domain: fair. Same physical pixels (Dice 0.993-0.998 on the four clean pairs; rigid fit identical to the
project's), 95% of the section, symmetric erosion (cannot manufacture a ratio; at most removes pial/edge deposits from
both), coverage 0.95-0.99 in both members. Coarse (~10 um/px) registration is irrelevant for an area fraction over 20 mm2.
Uniform anchors: do not undermine the conclusion (section 2).

Is the change affine? No. docs/normalisation.md: "if the two differ by an approximately affine intensity map, a method
that is affine-invariant by construction cannot fail that test. Establish what the acquisitions actually differ by."
I re-read all three channels at 2.6 um/px inside the reviewer's common ROIs, fitted rescan = a + b*original on the
bulk quantiles (p5-p75) and compared the bright tail with the affine prediction:

| Pair | Cy3 median ratio | Cy3 affine gain | Cy3 p99/p99.5/p99.9 observed / affine prediction | normalised p99.5 orig -> resc | frac pixels I'>3 resc/orig | forest ratio |
|---|---:|---:|---|---|---:|---:|
| 1 V22 | 2.03 | 2.10 | 0.64/0.62/0.59 | 8.2 -> 5.0 (x0.60) | 0.56 | 0.59 |
| 2 L03 | 1.73 | 1.72 | 0.64/0.57/0.51 | 6.7 -> 3.7 (x0.55) | 0.46 | 0.44 |
| 4 Y13_s04 | 2.04 | 2.05 | 0.60/0.58/0.55 | 11.9 -> 6.5 (x0.54) | 0.48 | 0.52 |
| 7 B16 | 2.19 | 2.14 | 0.55/0.53/0.50 | 8.2 -> 4.0 (x0.49) | 0.30 | 0.37 |

The changed-exposure pairs show the same tail deficit (p99.5 observed/affine 0.53-0.81) once exposure is absorbed in
the gain. No saturation in the four clean pairs (max 15-35k of 65,535). Background doubled, plaque contrast relative to
background halved, and a one-offset-one-gain normaliser (or a fixed-window z-score) cannot undo that. A crude threshold
on the normalised image (0.30-0.56) accounts for the forest's 0.37-0.59; anchor choice is not the mechanism.

Scanner or specimen? Pure illumination/gain change (affine) is rejected by the tail deficit. Glass rose at identical
exposure in every channel (Cy3 x1.3-1.5, DAPI x2.6-6.5, FITC x1.1-1.4) and channel factors differ (DAPI median
x2.2-3.1, Cy3 x1.7-2.2, FITC/GFAP x1.2-1.35 with FITC tails nearly affine, 0.85-1.02): glass carries no specific
fluorophore, so the instrument state (per-channel light source, stray light, dirt/haze, focus) differed in August.
ADR-0006 says lamp intensities were identical across the July files; nobody has reported the light-source/camera fields
of the August VSI files - the cheap decisive check still undone. A specimen-side loss of specific Cy3/DAPI signal
between 19 July and 5 August fits the Cy3/DAPI deficits and near-affine FITC, but wavelength-dependent haze from "dirt on
the slide and the microscope slightly out of focus" (the wet-lab note for that session, config/section_notes.csv) would
also raise medians more than peaks; defocus alone cannot raise glass 1.5-6x. Verdict: instrument state changed AND
plaque-to-background contrast fell; the data cannot apportion the contrast loss between fading and haze. For the study
it does not matter: August is on a different scale in every arm, and rescanning again would not restore July.
The decisions.md "no photobleaching" finding measured background per millisecond; in the plaque band the same-exposure
rescans carry 0.50-0.64 of the affine expectation. It should be withdrawn or restated.

## 5. Task 5 - macro rescan-range ruling

Verified from codex_audit/macro_rescan_pairs.csv: same-exposure pairs 1,2,3,4,7 -> 7.98, 14.64, 5.30, 5.80, 6.52x (range
5.303-14.638); changed-exposure pairs 5,6,8,9 -> 1577, 1896, 622, 3270x (">600x", max 3,269.7); original macro endpoints
0.00733, 0.01338, 0.00087, 0.00044% - the near-zero quantity is the original Cy3-positive area, not DAPI area; DAPI-area
ratios 1.089-1.426 and Cy3-area ratios 6.96-15.95 on same-exposure pairs as stated. The stratified restatement is
arithmetically right, concedes the project's like-for-like range, and corrects two imprecisions in ADR-0034 ("two
originals" -> four sections from two slides; "denominator" -> original Cy3-positive area). Fair, and more careful than
the response it rules on. Residuals: pair 3 is among the "five same-exposure pairs" despite its blown-out July DAPI
channel (macro's DAPI denominator unreliable) - over four clean pairs the range is 5.8-14.6x, nothing changes; and the
5-15x is what a fixed absolute threshold must do on a 2x brighter background (affine gain 2.05-2.14; unnormalised
pipeline arm up 1.2-1.6x on the same pairs) - a property of the threshold rule on a brighter session, not of the tissue.

## 6. Study-design implications and recommendation (auditor's judgment)

1. 121 selected = 111 July originals (52 secondary-only controls, 59 positives) + 10 August rescans, all positive, none
   negative. The 10 come from 5 animals (A09, B16, L03, R22, V22): for each, BOTH positives are August and BOTH controls
   are July - five animals have their entire endpoint on the other session. Seven are re-acquisitions; three
   (V22b_snew01, L03b_snew02, B16b_snew01) were first imaged in August, flagged "not super nice" by the wet lab, who asked
   for a sensitivity analysis without them that was never run.
2. The five animals span all four groups (2/1/1/1; group sizes 12/8/6/5): the ~0.5x depression is diluted, not balanced.
   With 5-12 animals per group, one or two animals at half their true burden can move a group median; the diet-arm
   hippocampal GLOBAL ratio (0.579, MWU p 0.020, aris_v3 summary) is not protected. The project's 19-August estimate of
   the session propagation (28% of between-animal spread) used z130 and the defective GLOBAL; under v4 GLOBAL at 0.48x it
   will be larger.
3. For V22, L03, B16 the standard-exposure July originals exist, are in focus, have manual outlines and Dice 0.993-0.998
   to their rescans; they were displaced only by a "cohort default of preferring a rescan" that the wet lab says "exists
   only for the two slides re-acquired to fix exposure (the two exposure-fix slides)". Three animals sit on the wrong session for no reason.
4. Rule: PRIMARY ANALYSIS ON THE JULY SESSION ONLY. Substitute V22_s01, L03_s03, B16_s01 for their rescans (anchors are two
   numbers per image; tissue ROIs exist in ~/ihc_work/tissue_rois_pairs/), drop the three August-only sections, treat
   A09 and R22 (July positives unusable at 146/398 ms) as missing in the primary (n=29 animals) and carry them in a
   labelled sensitivity analysis. Apply to every arm: nonorm (x1.2-1.6), z130 (x0.4-0.7), GLOBAL (x0.4-0.6) are all
   session-dependent; "GLOBAL fixes it" is not available.
5. Do not correct August by 1/0.48: four pairs spanning 0.37-0.59 give +-25% on the factor, and the changed-exposure
   animals would need a different one. Use the pairs as a bound, not a calibration.
6. Before any further rescanning: dump light-source/LED intensity and camera settings from RawData/Rescan/Image_*.vsi and
   compare with July. Different -> instrumental, July is the reference; identical -> specimens changed, and any future
   rescan is invalid as a replacement.
7. Record it: a new ADR with the session effect per arm and the selection rule; withdraw the "no photobleaching" finding;
   fix docs/normalisation.md (the promised gain floor does not exist in make_global_classifiers.py); reword ADR-0012's
   "clean repeatability pairs" to what they measure - session compatibility, which fails.
8. The reviewer's section 2.3 conclusion stands; it should say "demonstrably harmful to the endpoint for five of 31
   animals", not merely "not automatically harmless".
