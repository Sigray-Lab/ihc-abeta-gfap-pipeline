# Independent audit: training annotations, reference-box panel, normalisation argument
**Scope:** claims 1–6 from the caller's brief, checked against the reviewer's round-2 assessment
(`/Users/pontusps/Documents/Codex/2026-09-09/re/outputs/REVIEW_RESPONSE_ASSESSMENT.md`, §2.1, §2.7, §2.9, §3.2).
**Method:** everything below was recomputed from the GeoJSON/CSV files with Python (shapely 2.1.2) in this
scratch directory; QuPath was not run; nothing under Dropbox, `~/ihc_work` or `~/Documents/Codex` was modified.
Pixel size 0.325 µm; classifier working resolution 1.30 µm/px (from `retrain_global.groovy`).
**Blinding:** counts are keyed by blinded image name only; no annotation-map rows, no tube/group links are reproduced.
The only QuPath-dependent number I could not regenerate is the sample count per file; for that I rely on the reviewer's
reproduction log (byte-identical `v4_global.json`, 4,782 samples) and ADR-0036, which quotes the same 4,782.

Supporting files written here: `stroke_counts_by_blinded_image.csv`, `box_geometry_recomputed.csv`.

---

## Claim 1 — the two Ignore*-only files hold 4 strokes (3 + 1), not 30 — **Confirmed**

Source set: `/Users/pontusps/ihc_work/annotations/combined_v5/` has exactly 12 human files (`Image_*_annotations.geojson`)
plus 6 machine files (`*_auto_annotations.geojson`, excluded here). The 12 human files are byte-identical to the 12 files
in `annotations/annotations_training_v2/` (12 identical, 0 differing), and the old and v5 blinding maps list the same 12
human images, so the reviewer and I assessed the same annotation state.

All 202 features are `LineString` pen strokes. Counts per class per file (blinded image name):

| image | Abeta | Negative | Ignore* | total | classes present |
|---|---|---|---|---|---|
| A09b_s03 | 8 | 9 | 4 | 21 | 3 |
| F02_s01 | 11 | 4 | 1 | 16 | 3 |
| G09_s03 | 0 | 7 | 3 | 10 | Negative + Ignore* (secondary-only control) |
| L03_s01 | 0 | 6 | 2 | 8 | Negative + Ignore* (secondary-only control) |
| L03b_s03 | 15 | 10 | 4 | 29 | 3 |
| Q07_s01 | 13 | 14 | 5 | 32 | 3 |
| R22b_s03 | 11 | 8 | 2 | 21 | 3 |
| V07_s01 | 5 | 4 | 2 | 11 | 3 |
| **V07_s02** | 0 | 0 | **3** | 3 | **Ignore* only** |
| W12_s02 | 12 | 13 | 1 | 26 | 3 |
| Y06_s02 | 12 | 10 | 2 | 24 | 3 |
| **Y13_s01** | 0 | 0 | **1** | 1 | **Ignore* only** |
| **Total** | **87** | **85** | **30** | **202** | |

- Single-class files: two, both Ignore*-only, **3 + 1 = 4 strokes**. **30 is the Ignore* total over all 12 files.**
- Stroke length (µm, polyline length × 0.325): Abeta 1,158 total, median 10.2 (IQR 4.5–18.1, max 54.6);
  Negative 3,004 total, median 14.1 (IQR 6.5–39.0, max 201.8); Ignore* 2,025 total, median 53.4 (IQR 24.6–93.0, max 216).
  All 202 strokes: 6,186 µm, median 13.9 µm. **Abeta + Negative (172 strokes): median 12.2 µm** — this is the
  "172 short lines (median 12 µm)" of ADR-0036 and the reviewer's "≈170 strokes, ~12 µm". 119/202 strokes are < 20 µm.
- The 4 dropped strokes carry 449 µm (22 % of all Ignore* length). Calibrating stroke length against the reviewer's
  per-file sample counts (0.70–1.07 samples/µm, median 0.875) gives an estimated pool composition of
  ≈ 950 Abeta / 2,480 Negative / 1,350 Ignore* samples (of 4,782); adding the 416 dropped samples would raise Ignore* from
  ≈28 % to ≈34 % of the pool.
- Abeta strokes come from only **8 images = 8 animals**; the entire positive class is ≈ 950 pixel samples.

## Claim 2 — ADR-0032 "losing 30 strokes"; the label-map API; is the fix trivial — **Confirmed (ADR sentence is wrong)**

Exact sentence, `pipeline/docs/decisions.md` lines 1506–1508 (ADR-0032, "Fix" paragraph):

> The two `Ignore*`-only files are still dropped: QuPath yields no training data from a single class, and inventing a
> second class in them would be worse than losing 30 strokes.

Judgement: **wrong on the number** (4 strokes, ≈416 samples, not 30) and **wrong on the premise**. "QuPath yields no
training data from a single class" is true only of the method the script calls. `retrain_global.groovy` line 128 calls
`training.createTrainingData([imageData])`; that method emits "Unlocked annotations for at least two classes are required
to train a classifier!" and returns null, and line 135 turns that into `SKIP … single-class file, no contrast to learn`.
The same class in the installed QuPath (`~/Applications/QuPath.app/Contents/app/qupath-extension-processing-0.7.0.jar`,
`qupath/process/gui/commands/ml/PixelClassifierTraining.class`) also exports **`createTrainingDataForLabelMap`**
(confirmed from the class-file string table without running QuPath). The reviewer's probe
(`…/re/work/round2/model_tables/single_class_probe.groovy` line 121, log `single_class_raw.log`) calls exactly that with the
canonical map on the two files and gets **316 + 100 = 416 samples, all label 1 (Ignore\*)**, no invented annotations.
So the reviewer is right: the exclusion is an implementation choice.

Why the script skips them: it trains per image (one offset/gain per section, lines 97–102, 126–128) so each image is a
separate `createTrainingData` call with its own auto-numbered label map — the very cause of the ADR-0032 scramble. The
**fix is a one-line change**: replace line 128 with `createTrainingDataForLabelMap([imageData], CANON-as-PathClass-map)`.
That both admits single-class files and makes the remap block (lines 137–176) a no-op safety check, because every image
then comes back already in canonical order (the probe log shows `remapped 0`).

Two caveats the ADR should carry instead of the "30 strokes" sentence: (1) one of the two files maps to a section whose
tissue anchor is broken — `global_norm_anchors.csv`: Y13_s01 glass 60, **t_p50 2,233.5, tissue area 0.2 mm²** (cohort
t_p50 141–575, tissue ≈ 20 mm²); its Ignore* features would be divided by a gain of 2,173 and are garbage under the GLOBAL
transform. It is already excluded from the 67-section measurement set for that reason. So only the 3-stroke file (V07_s02,
gain 196, 316 samples) should be added. (2) The trainer comment (lines 129–134) cites a "2026-08-18 ablation" showing the
drop changes repeat error by < 1 SE; I found no ADR or DerivedData record of that ablation (the word appears only in
review documents), so it is an unrecorded claim.

## Claim 3 — reference panel geometry — **Confirmed in full**

Inputs: `/Users/pontusps/ihc_work/reference_boxes/{F15_s03,G07_s01,J09_s04,T18_s03,V02_s01,Z23_s01}.geojson`,
region outlines `/Users/pontusps/ihc_work/region_rois/{Hippocampus,Isocortex}/<image>.geojson`, tissue outlines
`/Users/pontusps/ihc_work/tissue_rois/<image>.geojson`, generator `pipeline/scripts/make_reference_boxes.py`.

- 6 sections, **6 distinct animals**; none in the old map (10 training animals) nor the v5 map (14 animals incl. auto files): overlap 0 / 0.
- **38 features, all class `RefBox`: 18 HIP + 20 ISO.** Per section (HIP/ISO): F15_s03 3/4, G07_s01 3/2, J09_s04 3/4,
  T18_s03 3/4, V02_s01 3/2, Z23_s01 3/4. The 42-box ceiling (3+4 per section) is missed because the script accepts a
  500-µm layout as soon as ≥ 2 boxes fit (line 86), giving two sections only 2 ISO boxes.
- Side lengths: exactly 500.0 or 300.0 µm at 0.325 µm/px, square: **17 × 500 µm** (9 HIP, 8 ISO) and **21 × 300 µm** (9 HIP, 12 ISO).
- Containment: every box `contains`-tested inside its own region polygon (max area outside = 0.0 µm²) and inside the tissue
  outline. **Pairwise overlap: 0** (all pairs within each image, cross-region included).
- **Minimum edge clearance to the region outline: 100.65 µm** (max 285.8); minimum clearance to the tissue outline 108.6 µm.
  The reviewer's 100.658 differs from my 100.651 only because they used QuPath's per-image calibration instead of 0.325.
- Total box area 6.14 mm²; boxes cover 6–27 % of each region outline (median ≈ 13 %).

## Claim 4 — brightness stratification of the panel — **Confirmed; the limitation is real and larger than §2.7 states**

"67 valid stained sections" = rows with `condition == positive` and `valid == True` in
`DerivedData/aris_v2_2026-08-18/three_arm_section_level.csv` (= the 69 manifest rows with `use_for_measurement == True`
minus S17_s04 and Y13_s01). Their `cy3_med`: **143–582 (4.07×)**, median 259.

- Panel `cy3_med`: **204, 236, 245, 293, 314, 325** → cohort percentile ranks 16–76. **15/67 above 325; 11/67 below 204.**
- Panel provenance (manifest): 6/6 `scan == original`, 6/6 `condition == positive`, 3 diet + 3 IP animals, Cy3 exposure 1840 ms.
- **No rescans in the panel, yet 10/67 valid sections are rescans and they occupy the bright end (cy3_med 298–582)** —
  i.e. most of the 15 sections above the panel maximum are rescans, which is the regime the halo/rescan-drift findings
  concern. No secondary-only sections: the 52 valid controls (cy3_med 11–381) — the population on which the false-positive
  ceiling (up to 2.9 % under v4) is defined — are unrepresented.
- Edges/folds: all boxes are ≥ 100.65 µm inside the region outline and ≥ 108.6 µm inside the tissue outline, so tissue
  rims — the leading false-positive source named in ADR-0036 — cannot be scored. One panel section carries
  `has_flagged_artefact == True` in the manifest; whether any box overlaps that artefact cannot be decided from geometry.
- Why the panel is narrow: the generator draws the median section of each brightness tertile *within arm* after excluding
  14 training animals and images lacking both outlines; taking tertile medians guarantees the extremes are never sampled.
  For contrast, the stained human-annotated **training** images span 143–582, the full cohort range.
- What the panel can validate: pixel-level agreement on interior parenchyma of mid-brightness original scans of stained
  tissue. What it cannot: the false-positive rate on antibody-free tissue, rim/fold behaviour, the bright-rescan regime,
  and whether box-level accuracy transfers to the whole-region percent-area endpoint (boxes cover 6–27 % of each region).

## Claim 5 — the median-shift argument and the group gain ratios — **Confirmed arithmetic; practical effect small but systematic; ratios reproduced**

Arithmetic: background U[100,200], fraction p replaced by pixels above the median. The median m solves
(1 − p)(m − 100)/100 = 0.5 → **m = 100 + 50/(1 − p)**: 150.00 (p = 0), 155.56 (p = 0.1), 162.50 (p = 0.2). Correct.
So `normalisation.md` lines 46–47 and the trainer header ("the median … is plaque-free by construction and barely moves")
are false as stated — the median pixel is background, but its value moves.

Actual sensitivity here (67 sections, `global_norm_anchors.csv`; new median = background quantile at 0.5/(1 − p),
interpolated linearly between t_p50 and t_p75; assumes deposit pixels lie above the median):

| plaque fraction p | median shift (counts), median across sections [range] | as % of gain (t_p50 − glass) |
|---|---|---|
| 1 % | 0.8 [0.4–2.1] | 0.4 % [0.3–0.8] |
| 2 % | 1.6 [0.8–4.3] | 0.8 % [0.5–1.7] |
| 5 % | 4.2 [2.1–11.1] | 2.0 % [1.4–4.3] |
| 10 % | 8.9 [4.4–23.3] | 4.1 % [2.9–9.1] |
| 25 % | 26.7 [13.3–70.0] | 12.4 % [8.8–27.3] |

Observed between-section spread of t_p50: **141–575, 4.08× (p90/p10 = 2.75×)**, IQR 208–313. Observed burden estimates:
z130 median 2.1 % (max 3.1 %), v3_global median 5.2 % (p95 12 %), unnormalised median 6.1 % (p95 27 %, inflated by rim
false positives). So at realistic burdens the anchor moves by ≈ 1–4 % of gain against a 4× between-section spread: the
shift is 1–2 orders of magnitude below the variation the normalisation is meant to remove. It is, however, systematic
and burden-correlated: a 10 %-burden section gets a ≈4 % higher gain, its normalised intensities fall ≈4 %, and its
measured area is compressed toward the null — a bias, small, conservative in direction. The calculation cannot address the
reviewer's second point (diffuse immunoreactivity raising the whole neuropil distribution), which no anchor statistic can
detect; that is a validation question, not an arithmetic one.

Group gain ratios (gain = t_p50 − glass_p50 per section, averaged per animal over its valid sections, median across
animals, larger group median / smaller): **diet 1.250 (8 vs 12 animals); IP 1.211 (5 vs 6 animals)** — reproduced
exactly. (Using the 69-row manifest set instead of the 67 gives IP 1.285; using t_p50 instead of gain gives 1.20/1.16.)
Which group is brighter is not reported here.

## Claim 6 — the reference instruction excludes vessels — **Confirmed; the scope is indeed unstated at the endpoint level**

`/Users/pontusps/ihc_work/reference_boxes/README.md` lines 17–19:

> Inside each yellow `RefBox`, outline every Aβ deposit with the brush or polygon tool, class `Abeta_ref`. Include diffuse
> deposits; exclude vessels, folds, edge glow and autofluorescent granules. Do not trace outside the box.

That defines the reference target as parenchymal deposits. The endpoint, by the project's own words
(`pipeline/README.md` lines 197–199), is "Aβ-immunoreactive area, not plaque load. A pixel classifier cannot separate
parenchymal deposits from vascular or intracellular signal without object rules". Decision D-3 (positive-class definition:
parenchymal only / + vascular / all immunoreactivity) is still 🟡 open in `pipeline/docs/decisions.md` line 270 and
`docs/DECISION_REGISTER.md` lines 62–71. Consequence: every vessel wall the classifier calls positive will be scored as a
false positive against this reference even though it is inside the declared endpoint. The reviewer is right that the
choice must be made explicit — and it must be made *before* tracing, because it changes what the tracer draws.

---

## Auditor's judgement

1. **Training set — not adequate for the generalisation being asked of it.** 4,782 samples is enough to fit a 12-feature
   random forest, but the positive class is ≈ 950 pixels from 87 lines of median 10 µm (≈ 8 px at 1.3 µm/px) on 8 animals;
   a line through a plaque core teaches "bright textured centre", not where a deposit ends, and percent-area is decided
   at the edge. The brightness range is covered (training images 143–582 = cohort range) but by one image per stratum.
   The RTrees fit is texture-driven (halo experiment: no arm recovers smooth blobs), and ADR-0036 shows the model cannot
   drop control false positives without dropping small deposits — the signature of an under-specified boundary.
   **Single most valuable addition:** polygon (area) annotations of whole deposits, core plus diffuse margin, with a
   matched Negative polygon on the adjacent 20–30 µm of neuropil and on the nearest vessel wall / cell body / rim, on
   ~10 sections chosen at the brightness extremes and including ≥ 3 rescans. That converts ≈ 950 core pixels into tens of
   thousands of boundary-informative ones and directly targets both failure modes.
2. **Reference panel — a sound but insufficient accuracy test.** Geometry is clean and blind-by-construction; but 6
   original mid-brightness stained sections (percentile 16–76) with interior-only boxes cannot certify the endpoint.
   Before trusting it, add: 12 boxes on 2 bright rescans (cy3 ≥ 400), 6 on 2 dim originals (< 180), 8 on 2 secondary-only
   controls (everything traced there is a false positive by definition), 6 rim/fold boxes deliberately straddling the
   tissue edge, whole-region tracing of at least 2 sections to tie box Dice to section percent-area, and a second tracer on
   ≥ 12 boxes for the human ceiling — roughly 70 boxes total. Lock D-3 first; otherwise the score is uninterpretable.
3. **Median-shift argument — right in principle, minor in practice.** At 2–10 % burden the anchor moves 1–4 % of gain
   versus a 4× between-section spread, so it is not the mechanism that could invalidate GLOBAL. Reword the "by construction"
   claim to the measured bound and move on. The 1.21–1.25× between-group gain differences are the point that matters:
   if treatment changes diffuse neuropil immunoreactivity, GLOBAL divides part of the effect out, and only an
   independent reference (including diffuse deposits) or a burden-independent anchor can settle that.
