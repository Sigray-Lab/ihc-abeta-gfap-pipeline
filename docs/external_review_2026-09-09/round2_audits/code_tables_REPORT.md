# Independent check of the round-2 reviewer: code robustness and derived tables

Date 2026-09-09. Method: read the code; small read-only checks on scratch copies of CSVs and
polygons (`check.py`, output in `check_output.txt`, inputs under `data/` in this directory).
No QuPath runs; nothing under Dropbox, `~/ihc_work` or `Codex` was modified. Blinding-key content
is not reproduced; coded image names (e.g. `A09_s02`) are used exactly as they appear in the
derived tables and ADRs.

Abbreviations: `P` = `.../IHC_analysis_pipeline/pipeline`, `DD` = `.../IHC_analysis_pipeline/DerivedData`,
`RF` = `.../IHC_analysis_pipeline/Results/Figures/20260909_globalfix`, `EV` = `~/ihc_work/eval_20260909`.

## Verdicts

| # | Claim | Verdict |
|---|---|---|
| 1 | Batch worker failure looks like success | **Confirmed — and it already fired in the production run** (S17_s04) |
| 2 | build_tables.py ignores QC; no key/uniqueness checks | **Confirmed**; whole-section path is weaker than the reviewer tested |
| 3 | Factories accept NaN/Inf/duplicates/missing names; no global gain floor | **Confirmed**; the missing floor sits under the "2.9 % worst control" headline |
| 4 | Missing halo baseline becomes zero | **Confirmed** (and the clip hides 38/360 negative responses) |
| 5 | Rejected rerun leaves stale polygon; dry-run overwrites CSV; index not rebuilt | **Confirmed**, with the nuance that the index is not on the measurement path |
| 6 | verified_ok never enforced; digests never compared; 129-row manifest; hardcoded thresholds | **Confirmed from code** |
| 7 | ADR-0033 "≤2 %, ≤31 µm"; README "2.5–4.5×" citing ADR-0018 | ADR-0033: **wrong for the on-disk change** (fair on area); README: **misattributed** |
| 8 | Minimum BH q ≈ 0.27 over 16 comparisons | **Confirmed**: 0.269976 |

## 1. `run_eval_batch.sh` — Confirmed, stronger than stated

`P/scripts/run_eval_batch.sh:13` uses `set -u` only (no `-e`, no `pipefail`). `:35-36` consumes
QuPath's exit status with `|| say "WORKER$w FAILED ..."` while sending stdout/stderr to `/dev/null`,
so the error text is lost. `:38` reports "worker done (n rows)" for whatever accumulated. `:40-47`
`wait`s, then unconditionally concatenates the worker CSVs and writes `$OUT`. `:48` logs `DONE` and
the script's exit status is that of `rm -rf`, i.e. 0. The merged row count is never compared with
`TOTAL` (`:26`). If a QuPath script exception does not yield a non-zero process exit (not verifiable
without running QuPath), even the log line is absent; the only universal signal — rows vs images — is
never checked.

Production evidence: `EV/region_Hippocampus_REGIONNORM_v4.csv.log` contains
`16:51:53 WORKER1 FAILED rx_Hippocampus_S17_s04 S17_s04` followed by `DONE`; the CSV holds 60 rows for
61 images. Cause: `P/scripts/make_region_anchors.py:69-70` skipped S17_s04/Hippocampus (< 500 px
inside region ∩ DAPI mask — this is one of the two failed-DAPI sections), so
`~/ihc_work/region_anchors_20260909.csv` has 60 Hippocampus / 61 Isocortex rows, no
`rx_Hippocampus_S17_s04.json` exists (`gx_S17_s04.json` and `g5_S17_s04.json` do), and QuPath could
not load the classifier. Harmless only because S17_s04 ∈ `BAD` (`build_tables.py:26,71`). The
reviewer's "current real batches were complete" is not accurate: the failure path fired, and nothing
but one log line recorded it.

## 2. `build_tables.py` — Confirmed

**Validity.** `DD/aris_v3_2026-09-09_globalfix/build_tables.py:26` defines `BAD` as two images;
`:48` `valid = ~image.isin(BAD)`; `:43` `cols` drops `qc_status`/`classified_coverage` from the new
GLOBAL file; `:47` computes `coverage` and never gates on it; `:59` and `:85` use `valid` only. The old
input (`DD/aris_v2_2026-08-18/full_cohort_in_tissue_RAW.csv`) has no `qc_status` column at all
(10 columns), so the gate has to be coverage-based — and coverage is already computed at `:47`.
Numbers (check §1, from the saved `section_level_in_tissue.csv`): `v2_nonorm_s1` controls n = 52,
median 0.029953 %; four rows with coverage ≤ 0.0002 (A09_s01, A09_s02, R22_s01, R22_s02; pct
0–0.0006 %) carry `valid = True`; excluding them n = 48, median 0.035989 %. GLOBAL and z130 have zero
collapsed rows; no stained section collapses under any arm. This contradicts ADR-0030
(`P/docs/decisions.md:1401-1408`: "Coverage below 1 % is a failed measurement, never a zero") and the
evaluator's own comment (`P/qupath/scripts/eval_in_tissue.groovy:70-73`). `summary.txt:4`
("controls n=52 median 0.030%") is therefore a saved number that the project's own rule says is wrong.

**Keys.** `:44` concatenates whole-section rows with no de-duplication and no completeness check: a
duplicated row double-counts in the area-weighted sum (`:57`), a missing row silently shrinks
`n_sections`. `:68` de-duplicates regional rows with `keep="last"`, which resolves conflicts by the glob
sort order of `:64` — an override, not a check. Live demonstration (check §4): since 20:35 the glob
`region_*.csv` also matches `EV/region_Hippocampus_GLOBAL_v5.csv` (classifier prefix `g5_`), which
`arm_method()` (`:28`) passes through unchanged, producing 61 per-image `arm_method` values. A rerun
today would add 61 junk rows to `regional_section_level.csv` / `regional_per_animal.csv` (saved:
118 × 4 rows). `summary.txt` and the figures would not change because they select methods by name, but
the CSVs would no longer be byte-reproducible. The whole-section path is weaker than the reviewer's
test, which only duplicated a regional row.

**Current inputs are sound** (check §2-3): `tissue_GLOBAL_v4.csv` 121/121 unique images, every
classifier equals `gx_<image>`; regional GLOBAL 61/61 in both regions; REGIONNORM 60/61 Hippocampus
(S17_s04, see §1) and 61/61 Isocortex; the 26 rescan regional rows removed at `:63` are matched exactly
by the 26 replacement rows. The single `LOW_COVERAGE` row in `tissue_GLOBAL_v4.csv` is S17_s04
(coverage 0.49, 33.9 %), excluded by `BAD`. No current treatment endpoint is affected.

## 3. Classifier factories — Confirmed, with a consequence the reviewer under-stated

`P/scripts/make_global_classifiers.py:59-63` rejects only `gain <= 0`; NaN and +Inf compare False
and pass. `:72-73` `float()` of NaN reaches `json.dumps`, which writes the non-JSON literal `NaN`.
`:74` `f"{a.prefix}_{r.image}.json"` with a missing image name writes `gx_nan.json`. `:66-76` iterate
every row, so a duplicated image overwrites one file while `written` counts two (`:78`). There is no
gain floor, yet `P/docs/normalisation.md:53-55` states "`make_global_classifiers.py` refuses below a
configurable floor" — nothing is configurable and nothing refuses. The regional factory does have a
floor: `P/scripts/make_region_classifiers.py:28` `MIN_GAIN = 50.0`, applied at `:59-62` — but its test
`gain < MIN_GAIN` also lets NaN and Inf through, and duplicates behave the same way.
(`make_region_anchors.py:69-70` has a separate 500-px floor on the anchor side.)

Current anchors (check §6): `DD/aris_v2_2026-08-18/global_norm_anchors.csv` has 121 unique, finite
rows, gains 6.0–2173.5; the four gains below 50 are A09_s02 (6), A09_s01 (7), R22_s02 (21), R22_s01
(30) — the same four sections that collapse under the unnormalised arm. `region_anchors_20260909.csv`:
121 rows, all finite, none below 50, no duplicate keys. So the reviewer's "gain 6 is present in the
worst current control" is correct, and more pointed than they made it: the worst GLOBAL control,
2.857 % (`RF/README.md:21` "up to 2.9 %", and the premise of ADR-0036), is A09_s02 at gain 6 (tissue
Cy3 median 11 counts over glass 4); the third worst, 0.947 %, is A09_s01 at gain 7. Spearman(gain,
GLOBAL control %) = −0.44 over 52 controls. Dividing by 6 where the cohort median gain is ~250
amplifies detector noise roughly 40×, and the documented-but-unimplemented floor would have refused
exactly these sections; with a ≥ 50 floor the worst control would be C18_s01 at 1.16 %. Not proof of
cause, but the headline specificity number and the design gap coincide on the same section.

## 4. `halo_analyse.py` — Confirmed

`P/scripts/halo_analyse.py:25` maps baselines by key (NaN when absent); `:26`
`(b.frac - b.baseline.fillna(0)).clip(lower=0)` turns a missing baseline into a zero baseline, and
the clip discards negative responses. Saved `DD/aris_v3_2026-09-09_globalfix/halo_recovery.csv`
(check §7): 360 rows, 0 NaN baselines — the current grid is complete — and 38 rows whose signed net
response is negative (min −0.0735) are reported as 0. `:17-18` de-duplicates exact rows only.

## 5. `register_rescan_regions.py` — Confirmed, with a nuance

`P/scripts/register_rescan_regions.py:177` `if a.dry_run or after <= before: continue` — a rejected
pair deletes nothing, so a polygon written by an earlier run stays in place. `:209-211` write
`rescan_registration.csv` unconditionally, dry-run included. `region_roi_index.csv` is written only by
`P/scripts/make_region_rois.py:140`, from that script's own drawn-outline rows, so a rerun of
*that* script would drop the 13 transferred entries again; the registration script never touches
the index. Nuance: nothing on the measurement path reads the index — `run_eval_batch.sh:23`
enumerates `$ROI/*.geojson` — so the stale-polygon path is the real risk and the index is provenance
only. ADR-0033's "can no longer be absent" describes a one-off repair, not an invariant. Current
state (check §9-10): 7 pairs, all accepted, 3 with non-zero angle; index 122 rows = 122 files (old
index 109).

## 6. Verification and ingest — Confirmed from code

- **`verified_ok`** is computed at `P/src/ihc/ingest/manifest.py:472` inside `_verified_hashes()`.
  The only other occurrence in the repository is `P/docs/decisions.md:506` (ADR-0016, 2026-08-07,
  "written but never read"). It is not in `MANIFEST_COLUMNS` (`:205-235`) and is not consulted at
  `:739-803` where digests are copied into rows. Selection is `finalise_scan_preference`
  (`:926-971`): `analysis_include & scan_is_preferred & condition-resolved`. `validate_manifest`
  (`:1243-1498`) contains no digest or verification check (grep of that range: none). Poisoned
  flags and digests therefore cannot change selection or validation; the reviewer's 121 selected /
  0 problems follows from the code without needing their run.
- **Digests are recorded, never compared.** `P/src/ihc/ingest/verify.py:759-762` says hashing
  "records the digest, so a later re-verification can prove the bytes have not changed"; no code
  performs that comparison. `inputs_fingerprint` (`manifest.py:1006-1046`) hashes file-name lists
  plus the bytes of `slides.csv`/`section_notes.csv`, and `cmd_doctor` (`P/ihc:206`) compares only
  that — a changed `.ets` under an unchanged name passes "manifest is current".
  `tests/test_verify.py:642-651` tests change detection at function level only.
- **129-row manifest.** `write_manifest` (`manifest.py:1564`) filters on `analysis_include` only.
  `~/ihc_work/manifest/manifest_analysis.csv`: 129 rows, `use_for_measurement` True 121 / False 8;
  the blinded manifest has 121. `P/README.md:127` documents the filter and `build_tables.py:35`
  applies it. A documented contract, and a footgun for any new consumer.
- **Thresholds.** `verify.py:143-156` (pixel tolerance 0.002 µm, dehydration ratio 0.10),
  `manifest.py:148-172` (pixel tolerance 0.001, exposure tolerance 0.5 ms, rescan match 0.25),
  `:1082-1083`, `:1640`; `P/ihc:64` (ratio 0.5). `P/config/config.yaml:31-32` (pixel tolerance
  0.0005) and `:410` (ratio 0.5) are read by nothing (grep over `ihc`, `src`, `scripts`: no consumer).
  Three different pixel-size tolerances coexist.

## 7. Documentation claims

- **ADR-0033** (`P/docs/decisions.md:1533-1536`): "Effect on the six non-zero-angle outlines is ≤ 2 %
  area and ≤ 31 µm vertex displacement (reviewer's measurement)". Archived tarball vs current files
  (check §8, clipped polygons): six outlines changed (A09b_s03, A09b_s04, R22b_s03 × two regions),
  seven are identical; maximum area change +2.0526 % (Isocortex A09b_s04) — reproduces the reviewer's
  figure exactly; Hausdorff distance up to 50.3 µm (Hippocampus A09b_s03), consistent with the
  reviewer's 51.7 µm on unclipped source vertices. Hippocampal areas change by < 0.003 % while the
  outlines move 30–50 µm: area is the wrong metric for a rigid displacement. Verdict: the "≤ 31 µm"
  refers to a different comparison (round-1: centre-corrected-but-direction-wrong formula vs true
  inverse) and is wrong as a description of what changed on disk; the ADR itself says the archived
  files predated the offset fix, so the sentence is internally inconsistent. "≤ 2 %" is fair rounding.
- **README** (`RF/README.md:69-70`): "anti-Aβ42 vs our 82E1, which stains 2.5–4.5× more area in this
  line, ADR-0018". ADR-0018 (`decisions.md:558-560`) compares 82E1 with **6E10** (Araki et al.), not
  with an anti-Aβ42 antibody. Misattributed. The sentence is generated by
  `RF/make_age_adjusted.py:141-142` and rewritten into the README between the AGEADJ markers
  (`:150-152`), so editing the README alone will be overwritten on the next run.

## 8. BH correction — Confirmed

Exact MWU from `regional_per_animal.csv` (check §5): p = 0.020349 / 0.040883 / 0.050621 for
hippocampus-diet GLOBAL / REGIONNORM / unnormalised; Benjamini–Hochberg over the 16 displayed
comparisons gives min q = 0.269976 (the three tie at 0.050621 × 16 / 3). From the rounded p-values in
`summary.txt`, min q = 0.272. scipy's default `mannwhitneyu` (two-sided, exact at these n without
ties) is what `build_tables.py:96` uses, so the p-values match the summary.

## Beyond the checklist

- **Blinding.** `DD/aris_v2_2026-08-18/global_norm_anchors.csv` carries a `tube` column (integers
  29–60, 31 unique) beside coded `image` names; so do `DD/aris_v2_2026-08-18/pair_tissue.csv`,
  `DD/brightness_investigation_2026-08-15/abeta_sections_all_resolutions.csv` (`tube_id`) and
  `~/ihc_work/region_rois/region_roi_index.csv` (`stem` = `Image_<tube>_sNN`). That is the code→tube
  half of the key, sitting in Dropbox-synced derived folders. The reviewer's blinding audit covered
  the blinded manifest only. Whoever can read `DerivedData` can un-blind, independently of the
  symlink question.
- **Tests.** None exist for `run_eval_batch.sh`, `build_tables.py`, `make_global_classifiers.py`,
  `make_region_classifiers.py`, `halo_analyse.py` or `register_rescan_regions.main()`;
  `tests/test_rescan_transform.py` covers only the matrix functions.

## Judgment

1. **Saved numbers, ranked.** (a) `summary.txt:4` unnormalised-control line (n = 52, median
   0.030 %) violates ADR-0030; the rule-consistent value is n = 48, 0.036 %. Real, small, control-only.
   (b) The "worst GLOBAL control 2.9 %" (README:21, ADR-0036 premise) is A09_s02 at gain 6, which the
   documented floor would have excluded; any specificity argument resting on it is fragile.
   (c) Everything else: inputs complete, unique and finite, the 26-row rescan swap exact, arithmetic
   reproducible — the reviewer's "none" holds for every treatment endpoint but not for (a), and (b)
   deserves a footnote wherever 2.9 % is quoted.
2. **Future reruns.** Batch driver (has already failed silently once), glob ingestion (already
   matches v5 files), no key/duplicate checks, factories, dry-run overwrite. All are cheap fixes.
3. **Three changes first.** (i) `run_eval_batch.sh`: record each QuPath exit code, compare merged
   rows with `TOTAL`, write `$OUT.partial` and rename only when complete, exit non-zero otherwise.
   (ii) `build_tables.py`: pin the input file list (no glob) and assert classifier prefixes ∈
   {gx, rx, v2}; assert exactly one row per (image[, region], arm_method) and raise on duplicates
   instead of `keep="last"`; set `valid = ~BAD & coverage ≥ 0.01` for every arm and report the
   failures in `summary.txt`. (iii) `make_global_classifiers.py`: reject non-finite, blank or
   duplicate anchors and enforce `MIN_GAIN` (50, as the regional factory) — or delete the sentence in
   `normalisation.md`; mirror the NaN/duplicate rejection in `make_region_classifiers.py`.
4. **Deferrals.** Pre-registration, atlas regions, backups: fine to defer (no code-safety
   consequence). `verified_ok`, hash baseline, config wiring, 129-row manifest: acceptable given the
   reviewer's fresh 167-file hash match, with one caveat — the raw tree lives in a sync folder, so a
   20-line comparison of `verify_raw.json` against a frozen copy is worth running before any rerun.
   "Atomic writes" should not be deferred as framed: the batch driver's failure-as-success is its
   concrete form and it has already happened. The "symlink leak" deferral is unsafe as framed: the
   code→tube mapping is already in plain CSVs under `DerivedData`, a larger hole than the symlinks.
