# Independent verification — secondary-only controls, exposure and normalisation gain under corrected GLOBAL (v4)

Auditor: independent re-computation, 2026-09-09. Read-only on all project paths; QuPath not run.
Everything below was recomputed from the raw evaluation rows, the August raw rows, the anchor file and
the manifest (allowed columns only: `condition`, `use_for_measurement`, `scan`, exposure fields). The
reviewer's own tables were read only after my numbers existed. No blinding-key content appears here.

Root `D` = `/Users/pontusps/Dropbox/longevity_project/projects/Rapamycin_CSF_cell_mice_Per/IHC_analysis_pipeline`.
Scratch scripts and logs: this directory (`pass1.py` … `pass6_fit.py`, `*.log`, `global_v4_joined.csv`,
`all_arms_joined.csv`, `noise_by_gain.csv`).

## 0. Data reconciliation (so the verdicts rest on the same rows as the project tables)

| item | value | source |
|---|---|---|
| manifest section rows / `use_for_measurement == True` | 130 / 121 (69 positive, 52 negative) | `D/BlindingKey/provenance_manifest.csv` |
| GLOBAL v4 rows | 121, all matched to the manifest; `qc_status` OK ×120, LOW_COVERAGE ×1 (S17_s04) | `~/ihc_work/eval_20260909/tissue_GLOBAL_v4.csv` |
| August nonorm / z130 rows | 121 + 121 | `D/DerivedData/aris_v2_2026-08-18/full_cohort_in_tissue_RAW.csv` |
| valid = not in {Y13_s01, S17_s04} | 119 = 67 stained + 52 controls | `D/DerivedData/aris_v3_2026-09-09_globalfix/build_tables.py:26,48` |
| my `pct` vs project `section_level_in_tissue.csv` (363 rows) | max abs diff 1.8e-15 pp; `valid` identical | pass2.log |
| `gain` in `three_arm_section_level.csv` vs `t_p50 − glass_p50` in `global_norm_anchors.csv` | identical on all 121 | pass1.log |
| `cy3_med` (three_arm, from the 2026-08-15 `allchan.csv`) vs anchors `t_p50` | different measurements, agree only on 3 %; both used exactly as the reviewer used them | pass1.log |

## 1. Claim 1 — the worst GLOBAL control is the dimmest, not the brightest — **Confirmed**

52 valid controls under GLOBAL v4, sorted by percent of tissue (full list in `pass1.log`):

| rank | image | % tissue | coverage | cy3_med | glass_p50 | t_p50 | gain | Cy3 exposure (ms) |
|---|---|---|---|---|---|---|---|---|
| 1 | **A09_s02** | **2.856537** | 0.9515 | **11** | 4 | 10 | **6** | **145.677** |
| 2 | C18_s01 | 1.160959 | 0.9668 | 167 | 34 | 164 | 130 | 1839.999 |
| 3 | A09_s01 | 0.947342 | 0.9812 | 11 | 4 | 11 | 7 | 145.677 |
| … | | | | | | | | |
| 23 | **L03_s01** (brightest control) | **0.267636** | 0.9561 | **381** | 40 | 366 | 326 | 1839.999 |

* A09_s02 and A09_s01 share `cy3_med = 11`, the minimum among the 52 controls **and** among all 121
  sections ("tied for the dimmest" is exact). Gain 6 = 10 − 4 from the anchor file, and the shipped
  classifier `gx_A09_s02.json` carries exactly `subtract 4.0 / divide 6.0`.
* The raw row (`abeta_in_um2 625452.1`, `tissue_roi_um2 21895465.2`) reproduces 2.856537 %.
* Faintest valid stained section: F15_s03, 0.707073 %.
* `README.md` line 31 (`D/DerivedData/aris_v3_2026-09-09_globalfix/README.md`) reads: "2.86 % on the
  brightest control exceeds the faintest genuinely stained sections (0.71 %)". The two numbers are right;
  the word "brightest" is wrong — it is the dimmest, shortest-exposure control. The reviewer's citation
  (`README.md:31`) is accurate.

## 2. Claim 2 — standard-exposure ceiling 1.161 %, still above 0.707 % — **Confirmed**

Cy3 exposure across the 121 measured sections (manifest `exposure_Cy3_ms`, `exposure_is_standard`):

| condition | scan | exposure (ms) | standard? | n |
|---|---|---|---|---|
| negative | original | 145.677 | False | 2 (A09_s01, A09_s02) |
| negative | original | 397.927 | False | 2 (R22_s01, R22_s02) |
| negative | original | 1839.999 | True | 48 |
| positive | original / rescan | 1839.999 | True | 59 / 10 |

* Four short-exposure controls, from two animals; every stained section is standard exposure.
* GLOBAL controls by exposure: 145.677 ms n=2, max 2.857 %, median 1.902 %; 397.927 ms n=2, max 0.278 %,
  median 0.250 %; 1839.999 ms n=48, **max 1.160959 % (C18_s01)**, median 0.2226 %.
* 1.161 % > 0.707 % (F15_s03): restricting to standard exposure removes the A09 outlier but not the
  overlap between controls and the faintest stained sections. Note that C18_s01 (gain 130), Q07_s04
  (0.764 %, gain 215) and S17_s02 (0.643 %, gain 169) are all standard-exposure, mid-gain sections.

## 3. Claim 3 — no configurable minimum gain in the global factory — **Confirmed**

`D/pipeline/scripts/make_global_classifiers.py` lines 58–63 (the only gain check in the file):
```
58    anchors = pd.read_csv(a.anchors)
59    anchors["gain"] = anchors.t_p50 - anchors.glass_p50
60    bad = anchors[anchors.gain <= 0]
61    if len(bad):
62        print(f"  FAIL non-positive gain on {list(bad.image)} -- tissue not above glass")
63        return 1
```
Its argparse options are `--classifier`, `--anchors`, `--dir`, `--prefix` only; nothing sets a floor.

`D/pipeline/scripts/make_region_classifiers.py`:
```
28  MIN_GAIN = 50.0     # dividing by a smaller gain amplifies noise beyond usefulness
57      d = pd.read_csv(a.anchors)
58      d["gain"] = d.r_p50 - d.glass_p50
59      skip = d[(d.gain < MIN_GAIN) | ~d.glass_p50.notna()]
60      for _, r in skip.iterrows():
61          print(f"  SKIP {r.image}/{r.region}: gain {r.gain:.0f} below {MIN_GAIN:.0f}")
62      d = d.drop(skip.index)
```
`D/pipeline/docs/normalisation.md` lines 53–55: "*The gain must be large enough to divide by.* A section
whose tissue barely exceeds glass produces a tiny divisor and amplifies noise. `make_global_classifiers.py`
refuses below a configurable floor."

Additional facts (read-only `git log` in `D/pipeline`): both scripts and the doc were introduced in one
commit (`074af90`, 2026-08-19) and never touched since; no revision of the repository has ever contained
`min_gain`/`min-gain`; no test references either factory. So the doc described an intention that was never
coded, and even the regional floor is a hard-coded constant, not "configurable". The mismatch is real but
has no numeric consequence for the shipped tables: the factory produced all 121 `gx_*` files, including
gain 6 and 7.

## 4. Claim 4 — brightness correlations — **Confirmed** (all four numbers reproduced exactly)

Section level, 67 valid stained, GLOBAL % vs `cy3_med`: Spearman **−0.0580**, p = 0.641.

Animal level: the reviewer's `stratify.py` (line 18) aggregates by the per-animal **mean of section
percentages** and mean `cy3_med`. With that definition I reproduce every figure to four decimals:

| quantity | reviewer | recomputed (mean of sections) | recomputed (area-pooled, the project's own per-animal definition) |
|---|---|---|---|
| GLOBAL stained, all 31 animals | −0.0865 | −0.0865 (p 0.644, n 31) | −0.078 (mean cy3) / −0.067 (median cy3) |
| GLOBAL stained, original + standard, 26 animals | +0.1262 | +0.1262 (p 0.539) | +0.141 |
| unnormalised controls, standard exposure, 25 animals | +0.6485 (p≈0.00046) | +0.6485 (p 0.000455) | +0.658 (p 0.00035) |
| GLOBAL controls, standard exposure, 25 animals | −0.2138 (p≈0.305) | −0.2138 (p 0.305) | −0.202 (p 0.33) |

Remarks: (i) "standard exposure" is a no-op for stained sections (all 67 are 1 840 ms); "original" drops
the 10 rescans. (ii) For controls, "standard", "original+standard" and "non-collapsed" are the same 48
sections / 25 animals. (iii) Section-level rho with p-values on repeated sections is descriptive only, as
the reviewer says.

**Finding of my own, not in the reviewer's text.** The project's own diagnostic (`normalisation.md`
lines 19–22: "r(log background intensity, log measured area) in controls should be ~0 … Do not use stained
sections for this") gives, on the 52 valid controls: unnormalised r = **+0.819**; GLOBAL r = **−0.482**
(p = 0.0003); within standard exposure GLOBAL r = **−0.287** (p = 0.048, n = 48); z130 −0.206 / −0.288.
Spearman: GLOBAL −0.425 (p 0.0017) / −0.359 (p 0.012). GLOBAL has not removed the brightness dependence
of the false-positive rate; it has reversed its sign. `README.md` lines 26 and 30 report only the stained
correlation (−0.06), exactly the quantity the documentation says not to use for this purpose.

## 5. Claim 5 — per-animal stained vs secondary-only readings uncorrelated — **Confirmed**

Area-pooled per animal (Σ Aβ area / Σ tissue area) over valid sections, animals with both conditions:

| method | rows used | n animals | Spearman | p |
|---|---|---|---|---|
| GLOBAL | all valid | 27 | **+0.0006** | 0.998 |
| z130 | all valid | 27 | **+0.1227** | 0.542 |
| unnormalised | all valid | 27 | +0.1545 | 0.442 |
| unnormalised | collapsed controls (coverage < 0.01) removed | **25** | **+0.2177** | 0.296 |

The reviewer's 0.218 / n = 25 is the collapsed-excluded variant (A09 and R22 lose all their controls),
which is consistent with their own §2.5 position, though the JSON does not state the exclusion. Either way
the conclusion holds: no animal-level coupling between stained and control readings under any method.

## 6. Claim 6 — build_tables.py counts four collapsed unnormalised controls as valid zeros — **Confirmed**

`build_tables.py` line 26 `BAD = {"Y13_s01", "S17_s04"}` and line 48 `sec["valid"] = ~sec.image.isin(BAD)`
are the whole-section validity rule; `coverage` is computed on line 47 and never used for QC; no `qc_status`
is read. (The reviewer cites line 63, which is the regional carry-forward filter; the whole-section rule
is line 48 — immaterial.)

The four collapsed unnormalised controls (coverage = (Aβ + Negative)/tissue ROI, from
`full_cohort_in_tissue_RAW.csv`), all `valid = True` in the September `section_level_in_tissue.csv`:

| image | coverage | Aβ µm² | Negative µm² | tissue ROI µm² | % tissue | Cy3 exposure |
|---|---|---|---|---|---|---|
| A09_s01 | 0.000000 | 0.0 | 0.0 | 21 039 944 | 0.000000 | 145.677 ms |
| A09_s02 | 0.000004 | 52.4 | 43.9 | 21 895 465 | 0.000239 | 145.677 ms |
| R22_s01 | 0.000012 | 103.1 | 104.8 | 17 328 215 | 0.000595 | 397.927 ms |
| R22_s02 | 0.000195 | 5.1 | 3 496.6 | 17 923 925 | 0.000028 | 397.927 ms |

* They are exactly the four short-exposure controls: the unnormalised forest sees a near-black image and
  assigns > 99.98 % of it to `Ignore*`. GLOBAL and z130 do not collapse on them (coverage 0.95–0.98).
* On `qc_status`: the August raw file `full_cohort_in_tissue_RAW.csv` has **no** `qc_status` column
  (10 columns; `eval_v4.groovy` did not export one) — the premise that it does is wrong. The status field
  exists only in the later exporter `D/pipeline/qupath/scripts/eval_in_tissue.groovy` lines 74–76
  (`COLLAPSED` if coverage < 0.01, `LOW_COVERAGE` if < 0.50) and in the August GLOBAL/regional files.
  However, the August derived table `aris_v2_2026-08-18/section_level_in_tissue.csv` (from
  `report_in_tissue.py`, `MIN_COVERAGE = 0.01`) flags all four `fail_cov = True, valid = False`. The
  September `build_tables.py` silently reversed that — a genuine regression, as the reviewer says.
* Control median unnormalised: **0.029953 %** among 52 vs **0.035989 %** among 48; maximum unchanged
  (1.988224 %). No stained unnormalised section collapses (minimum coverage 0.857), so no endpoint moves.

## 7. Gain distribution across the 121 sections (anchors: gain = t_p50 − glass_p50)

| condition | n | min | 2 % | 5 % | 10 % | 25 % | median | 75 % | 90 % | max |
|---|---|---|---|---|---|---|---|---|---|---|
| controls | 52 | **6** | 7.3 | 26 | 117 | 141.5 | 163.5 | 193.8 | 227.9 | 339 |
| stained | 69 | 63* | 110 | 129 | 147.8 | 170 | 211 | 274 | 448 | 2173.5* |

\* S17_s04 (63) and Y13_s01 (2173.5; its `t_p50` of 2233.5 was computed on a failed mask) are the two
failed-DAPI sections already excluded. **Lowest valid stained gain: 105 (Y06_s02).** Standard-exposure
valid sections span 105–512 (4.9×); the whole cohort spans 6–512 (85×) plus the failed-mask outlier.

Lowest gains: A09_s02 6, A09_s01 7, R22_s02 21, R22_s01 30 — then nothing until S17_s04 63 (excluded) and
Y06_s02 105 (stained). Consequences for any floor `g_min`:

| floor | sections removed | GLOBAL controls n / max / median | stained affected |
|---|---|---|---|
| none | — | 52 / 2.857 % / 0.251 % | — |
| 8–20 | A09_s01, A09_s02 | 50 / 1.161 % / 0.236 % | none |
| 22–30 | + one of R22 | 49 / 1.161 % | none |
| 31–105 | + both R22 | 48 / 1.161 % / 0.223 % | none |
| > 105 | | | starts removing valid stained sections |

The control-brightness Spearman under GLOBAL moves only from −0.425 to −0.359 when the four are removed:
the residual (negative) brightness dependence lives in the standard-exposure range.

## 8. What the gain-6 section does to GLOBAL — measured on the raw pixels

To avoid arguing from quantisation arithmetic, I read the Cy3 plane of 13 controls (gains 6–339) and 2
stained references at the classifier's working level (1.3 µm/px, pyramid level 2) with the pipeline's own
`read_pyramid_plane`, built the pipeline's DAPI tissue mask (level 5, upsampled), eroded it 52 µm from the
rim, applied each section's exact (glass, gain), and computed the forest's texture features (Laplacian,
weighted SD, gradient magnitude at σ = 1–2 px) inside the interior — in an antibody-free section these
are, by construction, background texture plus detector noise (`pass5_noise.py`, `noise_by_gain.csv`).

* Normalised texture in standard-exposure controls is essentially **gain-invariant**: LoG σ1 robust SD
  0.037–0.065 across gains 115–339 (median 0.056); weighted SD σ1 0.065–0.115. The normalisation does its
  job for real tissue structure. Overall the dependence is only n ∝ gain^−0.26, not 1/gain.
* A two-component model, n² = T² + (F/gain)², fits all 13 controls (median |residual| 9 %, max 33 %):
  LoG σ1 **T = 0.0553 ± 0.0031** (intrinsic, normalised), **F = 0.86 ± 0.05 ADU** (detector +
  quantisation floor at the 1.3 µm level); weighted SD σ1 T = 0.097, F = 1.88 ADU; gradient σ1
  T = 0.042, F = 0.59 ADU.
* Observed/predicted inflation of background texture over the intrinsic level: **A09_s02 2.7× (gain 6),
  A09_s01 2.5× (7), R22_s02 1.3× (21), R22_s01 1.4× (30)**, standard sections 0.7–1.2×. The forest was
  trained on standard sections; a control whose flat background carries 3× the texture it learned as
  background is read as plaque texture. That is the mechanism behind 2.86 %.
* Consistency checks: (i) regressing GLOBAL control % on the section's normalised p99.9 tail (48 standard
  controls, r = 0.57) predicts 0.90 % for A09_s02 and 0.44 % for A09_s01 from their (genuinely heavy)
  tails; the observed values are 3.2× and 2.2× higher — the same factor as the measured noise inflation.
  (ii) Under z130, which divides by the local SD and therefore normalises noise itself, A09_s02 is 2.2×
  the control median (rank 2 of 52), not 11×. A09's tissue is intrinsically busy; GLOBAL turns busy into
  2.9 %.
* Where the model puts a floor: added noise ≤ 5 % of intrinsic texture at gain ≈ 48 (LoG) / 60 (WSD);
  ≤ 10 % at 34 / 42; ≤ 20 % at 23 / 29; ≤ 50 % at 14 / 17; +100 % at 9 / 11. The two R22 sections
  (+13–36 %) read normally; the two A09 sections (+140–240 %) do not. The breakdown of the classifier
  therefore lies somewhere between +40 % and +140 % added noise, i.e. between gain ≈ 8 and ≈ 20; with two
  sections on each side it cannot be located more precisely from this cohort.

## 9. Senior-analyst judgment

1. None of the six findings moves a treatment endpoint: all 67 valid stained sections are standard exposure with gain 105–512, so any gain floor up to 105, or an exposure rule, changes no stained number. What they change is the honesty of the control story.
2. The consequential problem is specificity in the standard-exposure range, and no gain rule touches it: with the four short-exposure controls gone the ceiling is still 1.161 % (C18_s01, gain 130), with 0.764 % (Q07_s04) and 0.643 % (S17_s02) behind it, all above the faintest stained section (0.707 %).
3. Under GLOBAL the control false-positive rate is brightness-dependent with the sign flipped (log r −0.48 overall, −0.29 within standard exposure) and tracks the section's normalised upper tail (rho +0.65): rim, vessel and cell-body autofluorescence does not scale with the neuropil median, which is the ratio model's intrinsic limit. The README's "brightness-independent" rests on the stained correlation the project's own docs forbid for this test; the controls' correlation belongs in that table.
4. Cosmetic: the word "brightest" (README:31), the doc/code mismatch on the floor (no output changed), the reviewer's line-number slips, and mean-of-sections vs area-pooled aggregation (same conclusions either way).
5. Not cosmetic, not endpoint-affecting: `build_tables.py` re-admitted four sections that the exporter labels COLLAPSED and that August's own table had excluded. Honour the status field; report them as failures (control median 0.036 %, n = 48).
6. The gain-6 control shows that GLOBAL rescales detector noise together with signal: measured background texture at gains 6–7 is 2.5–2.7× the gain-invariant level (T ≈ 0.055) of every standard section, and the forest reads that as plaque texture; after allowing for A09's heavy tail the excess is still 2–3×, and z130, which self-normalises noise, sees only 2.2× the median.
7. A minimum-gain rule is justified — the doc promises it, the regional factory has it, and the physics is now measured. Reuse `MIN_GAIN = 50` for the whole-section factory: it is the ≤ 5 % added-noise point, it is the existing constant, and every value in 31–105 partitions this cohort identically (exactly the four non-standard-exposure sections out, every stained section in), so it cannot be accused of tuning on results. A floor of 8–20 that keeps R22 is data-consistent but would look drawn around A09.
8. The cleanest pre-specifiable form is the acquisition rule already in the manifest — GLOBAL is defined only for `exposure_is_standard == True` — with the gain floor as the code-level guard behind it; here they are the same four sections.
9. With either rule the GLOBAL control row becomes n = 48, max 1.161 %, median 0.223 %, brightness Spearman −0.36 — and the README's headline sentence remains true in substance: the corrected model's false-positive ceiling still exceeds its faintest stained sections. That is the sentence to keep; the word "brightest" and the gain-6 outlier are distractions from it.
10. Specificity has to be attacked where ADR-0036 already points (rim erosion applied identically to both conditions, an independent reference set), not in the normalisation constant.
