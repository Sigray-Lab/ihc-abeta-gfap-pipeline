# Assessment of the external round-2 review (2026-09-09)

**What was assessed.** `/Users/pontusps/Documents/Codex/2026-09-09/re/outputs/REVIEW_RESPONSE_ASSESSMENT.md`
(Codex, round 2) and its companion `OUTSIDE_TISSUE_BRIGHTNESS.md`. The reviewer saw the repository at
commit e0fd0cc; it did not see ADR-0036 (machine-drawn negatives, v5), the reference-box QC of T18_s03, or the
brightness-gate experiment (ADR-0037), all of which happened later the same day.

**How.** Five independent audits, each given a disjoint set of the reviewer's claims and told to recompute
from the raw rows, GeoJSONs and code rather than trust either the reviewer or our own ADRs. Their reports are in
`docs/external_review_2026-09-09/round2_audits/` (tube identifiers scrubbed). What follows is my own reading of
those audits plus the checks I did myself. No blinding-key content appears here.

## 1. Verdict on the review

The review is careful, almost entirely correct on numbers, and fair in its rulings. Every one of the ~40 numeric
claims the audits re-derived reproduced, most to four decimals. Its two retractions (one continuous original
session; macro ratios must be stratified by exposure) are right. Its central scientific verdict, that F1 and F6
are fixed but the measurement's accuracy is not established, is right and is now reinforced by evidence the
reviewer did not have (section 4).

Where it is wrong or understated:

- **"Current real batches were complete."** False. `region_Hippocampus_REGIONNORM_v4.csv.log` already contains
  `WORKER1 FAILED rx_Hippocampus_S17_s04` followed by `DONE` with 60 of 61 rows. It was harmless only because
  S17_s04 is in the excluded set. The silent-failure path the reviewer describes as prospective has fired.
- **The session effect is understated.** The reviewer calls mixing acquisitions "not automatically harmless".
  The rescan audit shows it is harmful to the endpoint as saved: 10 of the 121 measured sections are August
  rescans, all antibody-positive, from five animals whose entire endpoint sits on the August scale while their
  controls sit on July; same-exposure rescans read 0.37–0.59× under GLOBAL v4, 0.37–0.74× under z130 and
  1.2–1.6× under the unnormalised arm. Three of those five animals have good July originals that were
  displaced by a "prefer the rescan" default the wet-lab notes restrict to the two exposure-fix slides.
- **Brightness dependence of GLOBAL is understated.** The reviewer reports the stained Spearman (−0.06) and a
  control Spearman of −0.21 within standard exposure. The control audit applied the project's own documented
  diagnostic (`docs/normalisation.md` §19–22: correlation must be judged on controls, log scale) and gets
  r = −0.48 (p 0.0003) overall and −0.29 (p 0.048) within standard exposure. The README's "brightness-
  independent" uses the correlation the docs say not to use. The sign is reversed relative to the unnormalised
  arm: GLOBAL rescales detector noise with the gain, so dim sections gain texture the forest reads as plaque.
- **"No photobleaching" is not addressed.** `decisions.md` still records that finding; it was based on
  background counts per millisecond and is contradicted by every pair in the reviewer's own table.

Where the reviewer is right and it changes what we say:

| Item | Status | Consequence |
|---|---|---|
| Worst GLOBAL control (2.86 %) is the dimmest section (Cy3 median 11, gain 6, exposure 146 ms), not the brightest | confirmed | README:31 wording wrong; the number is a 146-ms artefact of the missing gain floor |
| Standard-exposure control ceiling is 1.16 % | confirmed | quote 1.16 % (n=48) as the ceiling, 2.86 % only with the exposure caveat |
| `make_global_classifiers.py` has no gain floor although `normalisation.md` promises one; regional factory has MIN_GAIN = 50 | confirmed | a floor of 50 removes exactly the four short-exposure controls and no stained section (stained gains 105–512) |
| Four collapsed unnormalised controls counted as valid zeros by `build_tables.py` | confirmed, regression against August's own `valid=False` flags | control median 0.036 % (n=48), not 0.030 % (n=52); treatment endpoints unaffected |
| ADR-0032 "losing 30 strokes" | wrong: 4 strokes, ≈416 Ignore* samples; `createTrainingDataForLabelMap` would take them | one-line trainer change; one of the two files sits on a broken anchor and should stay out anyway |
| Halo prose: footprint 0.857 mm² not 5.1; baselines not ~0 everywhere; weak-amplitude responses not "nothing"; "<1 % of Aβ totals" wrong denominator; no 80-µm hard cut-off | all confirmed | ADR-0035 sentences listed in the halo audit must be rewritten; the scientific conclusion (smooth plateaus poorly detected; z130 best on compact blobs) stands |
| Reference panel: 6 mid-brightness original stained sections, interior boxes only; 15/67 cohort sections brighter than the panel, all 10 rescans and all 52 controls unrepresented | confirmed | the panel cannot certify the control ceiling, rims or the rescan regime; add boxes there before trusting it |
| Vessel exclusion in the tracing instructions vs "all immunoreactive area" endpoint (D-3 open) | confirmed | lock D-3 before Aris traces |
| Median-shift argument (median = 100 + 50/(1−p)) | arithmetic confirmed; at 2–10 % burden the anchor moves 1–4 % of gain against a 4× between-section spread | reword "by construction"; minor in practice; the 1.21–1.25× between-group gain difference is the real normalisation question |
| BH across the 16 displayed comparisons: min q 0.27 | confirmed | as already stated: exploratory |
| ADR-0033 "≤31 µm" refers to a different comparison; archived→current max 51.7 µm, +2.05 % area | confirmed | fix the sentence |
| README "82E1 stains 2.5–4.5× more than the age-curve antibody" cites a 6E10 comparison | confirmed | remove or re-source; the sentence is regenerated by `make_age_adjusted.py:141` |

## 2. Findings beyond the review

Each of these came from an audit going past its checklist; none is in the reviewer's text.

1. **The halo plateaus are relabelled `Ignore*`, not `Negative`.** At 4×, the host-median `Ignore*` fraction of the
   injected footprint is 59–82 % for z130 (40–160 µm) and 55–80 % for GLOBAL (20–40 µm), ≈0 on clean images.
   The forests treat a bright smooth region as fold/debris and abstain. ADR-0035's mechanism sentence ("called
   Negative whatever its height") is wrong, and the recovery metric counts abstention as zero, which
   `normalisation.md` §5 forbids. A real broad, smooth deposit would be dropped from the measurement by every arm.
2. **Mechanism of the low-gain false positives, measured on raw pixels.** Background texture (LoG σ1 at 1.3 µm) is
   gain-invariant in standard sections; the model n² = T² + (0.86 ADU/gain)² fits within ±33 %. At gain 6–7 the
   texture is 2.5–2.7× the learned level and is classified as plaque. ≤5 % added noise corresponds to gain ≈ 50,
   which is the regional factory's existing floor.
3. **The session effect is a design defect with five animals on the wrong scale** (see above). The project's own
   19-August plan had computed a z130 session ratio of 0.81 and called the session "a genuine design defect";
   ADR-0034 is silent on it. The rescan audit's pixel re-read shows the July→August change is non-affine:
   bulk Cy3 background ×2.0–2.2 at identical exposure but the plaque band only ×1.0–1.3, so plaque-to-background
   contrast halved; a plain threshold on normalised pixels reproduces the forest's 0.4–0.6× without any forest.
   Glass rose in all channels with channel-specific factors, so the instrument state changed as well. The
   decisive cheap test, never run, is to dump the light-source and camera fields from the August VSI headers.
4. **Blinding leak in derived CSVs.** `DerivedData/aris_v2_2026-08-18/global_norm_anchors.csv` (`tube` column),
   `DerivedData/brightness_investigation_2026-08-15/abeta_sections_all_resolutions.csv` (`tube_id`),
   `DerivedData/aris_v2_2026-08-18/pair_tissue.csv`, and `~/ihc_work/region_rois/region_roi_index.csv`
   (`stem` = raw scan name) carry the code→tube half of the key next to blinded codes, in Dropbox-synced folders.
   This is a larger hole than the symlink issue the reviewer flagged. Custodian's decision.
5. **The training set's positive class is ≈950 pixels** from 87 pen strokes (median 10 µm) on 8 animals. That is
   the common root of the halo over-classification seen in the T18_s03 boxes (about half of v4's Aβ area at
   <1.5× local background), of the v5 collapse (ADR-0036), and of the reviewer's "accuracy not established".
6. **`build_tables.py` no longer reproduces its own outputs**: its `region_*.csv` glob now also matches the v5
   files written today, so a rerun would add 61 `g5_*` rows to the regional CSVs. Pin the inputs.

## 3. What the reviewer did not know

- **Reference-box QC of T18_s03 (this afternoon).** Inside three 500-µm hippocampal boxes, 39–52 % of v4's Aβ
  area is dimmer than 1.5× the box background; only 17–22 % is ≥3×. The extra area is faint textured tissue
  next to deposits and dim somata, not plaque. v4 is over-loose at the periphery; v5 (ADR-0036) was over-tight.
- **Brightness gate (ADR-0037, in progress).** Keeping v4's Aβ pixels only where Cy3 ≥ 2× a 400-µm local
  background of non-Aβ tissue: cohort table to be inserted below when the run completes.

| whole section, 119 sections | v4 | gate 2× |
|---|---|---|
| stained median (n=67) | 5.46 % | 2.22 % |
| worst control, all exposures (n=52) | 2.87 % | 0.61 % |
| worst control, standard exposure (n=48) | 1.16 % | 0.35 % |
| stained median / worst control | 1.9× | 3.6× |
| section ranking, stained (Spearman v4 vs gated) | | +0.979 |
| hippocampus diet ratio, MWU p (n=8 vs 11) | 0.578, 0.020 | 0.713, 0.109 |
| hippocampus IP ratio, p (n=6 vs 5) | 1.28, 0.66 | 1.32, 0.66 |

  The gate removes ~56 % of v4's Aβ area on stained sections (median ratio 0.44, range 0.29–0.61) and ~80 % on
  controls, keeps section ranking, leaves the brightness correlations unchanged, and weakens the diet difference
  (gate 1.5×: 0.75, p 0.09; gate 3×: 0.57, p 0.09). Part of the diet difference under v4 therefore lives in the
  dim halo. Full tables: `DerivedData/brightness_gate_2026-09-09/`.

## 4. What this means for the numbers as saved

- The whole-section and regional tables are arithmetically reproducible (reviewer and audits agree). They are
  not yet interpretable as amyloid burden: half the measured area in dense hippocampus is dim halo, five of 31
  animals sit on a session that reads ~0.5×, and the control ceiling depends on an exposure/gain rule that the
  code does not implement.
- The hippocampal diet ratio (0.58, MWU p 0.02 under v4) is not robust to any of these: the age scenario moves it
  to 0.72; the gate at 2× moves it to 0.71 (p 0.11); the
  session rule (July-only primary, n=29) has not been run. It should be presented as exploratory and
  method-dependent, which is what ADR-0034 already says.
- Nothing in either round changes the direction of the diet difference in hippocampus; everything changes its
  size and its p-value.

## 5. Actions, ranked (PI decides; nothing below has been applied unless marked)

1. **Session rule.** Primary analysis on the July session only: substitute the three July originals that were
   displaced without cause, treat the two exposure-fix animals as missing (n=29), carry all-31 as a labelled
   sensitivity analysis, all arms. Do not rescale August by 1/0.48. Dump the August VSI light-source/camera
   fields first (cheap, decisive). New ADR; withdraw "no photobleaching".
2. **Gain/exposure rule for GLOBAL.** MIN_GAIN = 50 (or `exposure_is_standard`) in `make_global_classifiers.py`,
   matching the regional factory; removes exactly the four short-exposure controls. Quote the control ceiling
   as 1.16 % (n=48). Fix `normalisation.md`.
3. **Boundary definition, not more classifier tuning.** Lock D-3 (parenchymal vs all immunoreactivity), then get
   whole-deposit polygon annotations with matched adjacent-negative polygons on ~10 sections at the brightness
   extremes including rescans; extend the reference panel with boxes on bright rescans, dim originals,
   secondary-only controls and tissue rims (~70 boxes total), second tracer on ≥12. This is the only route to
   "accuracy established".
4. **Table hygiene.** `build_tables.py`: pin inputs, assert one row per key, gate every arm on coverage
   (restores August's `valid=False` for the four collapsed controls). `run_eval_batch.sh`: compare rows with
   the image count and exit non-zero on any worker failure. Both are ~10 lines.
5. **Blinding.** Remove or move the tube columns from the four derived files listed in §2.4; regenerate
   `region_roi_index.csv` without raw stems. Custodian's call.
6. **ADR text corrections.** ADR-0032 (4 strokes, label-map API), ADR-0033 (51.7 µm), ADR-0035 (six sentences
   listed in the halo audit; add signed responses and the `Ignore*` fraction to `halo_analyse.py`), README:31
   (dimmest, not brightest), README:70 (6E10 comparator), the "brightness-independent" sentence.
7. **Not worth doing now:** the reviewer's full factorial stress test (texture × width × contrast × location).
   Signed responses and the `Ignore*` fraction are free and should be added; a one-day textured-patch run is
   reasonable; the rest waits for real reference masks.

## 6. Deferrals re-examined

Pre-registration, atlas regions, backups, the 129-row manifest, `verified_ok`, hash comparison and config wiring
remain deferred as decided. Two deferrals are unsafe as framed and are re-flagged: "atomic writes" (the batch
driver's silent failure is its live instance) and "symlink leak" (the key is already in plain derived CSVs).
