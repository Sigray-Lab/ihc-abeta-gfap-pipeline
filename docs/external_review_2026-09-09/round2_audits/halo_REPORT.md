# Independent audit of the reviewer's critique of the halo/blob-injection experiment (ADR-0035)

**Auditor role:** senior image-analysis auditor, independent of both the project and the external reviewer.
**Date:** 2026-09-09. **Mode:** read-only on the project, `~/ihc_work` and the reviewer's tree; no QuPath runs.
**Everything below was recomputed from the raw QuPath CSVs and blob JSONs in `/Users/pontusps/ihc_work/halo_20260909/`**, not from
`halo_recovery.csv`, the reviewer's derived files, or the ADR. Scripts and outputs are in this directory:
`recompute.py`, `paired_recompute.csv` (360 paired discs), `frame_totals_recompute.csv`, `median_clipped_all_amps.csv`,
`median_signed_all_amps.csv`, `zscore_diagnostic.csv`, `footprint_composition.csv`.

Sources read: `pipeline/scripts/halo_inject.py`, `pipeline/scripts/halo_analyse.py`, `pipeline/qupath/scripts/halo_eval.groovy`,
`docs/PLAN_halo_experiment.md` (note: it lives in `IHC_analysis_pipeline/docs/`, not `pipeline/docs/` as the ADR implies),
`pipeline/docs/decisions.md` lines 1562-1608 (ADR-0035), `Results/Figures/20260909_globalfix/README.md` lines 30-39,
`pipeline/docs/normalisation.md`, `pipeline/docs/classifier_abeta.md:91`, `pipeline/scripts/retrain.groovy`,
classifier JSONs `~/ihc_work/qupath_pairs/qupath/classifiers/pixel_classifiers/v2_{z130,nonorm}_s1.json` and
`~/ihc_work/halo_20260909/gx_*.json`; the reviewer's `REVIEW_RESPONSE_ASSESSMENT.md` §2.4/§3.4 and `HALO_ROUND2_REVIEW.md`.

Data integrity checks passed before anything else: 18 raw CSVs, 720 rows, 720 after exact de-duplication, 360 unique
(host, arm, blob, scale) keys per condition, 360/360 clean-injected pairs, no missing baseline; every raw CSV carries
non-trivial `Abeta`/`Negative` measurements (so the `meas()` fallback to 0 in `halo_eval.groovy:37` did not fire silently);
`DerivedData/aris_v3_2026-09-09_globalfix/halo_recovery.csv` and `Results/Figures/20260909_globalfix/halo_recovery.png`
are byte-identical to the `~/ihc_work` originals; the reviewer's six fresh QuPath CSVs for W12_s03
(`.../halo_audit/fresh_eval/`) are byte-identical (`cmp`) to the saved ones, so the saved measurements are reproducible
for at least one host under all three arms and both conditions.

---

## 1. Claim-by-claim verdicts

Notation: "net" = (Aβ area in the footprint disc on the injected image − same disc on the clean image) / disc area,
signed. "recovery" = max(0, net), which is what `halo_analyse.py:26` computes and what the ADR table and figure show.
Scale 1 = the footprint disc (radius 1.665 σ); scale 2 = the disc of twice that radius. 180 scale-1 cases = 3 hosts × 3 arms × 20 blobs.

### Claim 1 — the metric is a clipped net area change, not pixel recall; 16/180 scale-1 cases lose area (max 7.35 pp) and are reported as 0
**Verdict: Confirmed.**
`halo_analyse.py:25-26`: `b["baseline"] = ...map(clean)`; `b["recovery"] = (b.frac - b.baseline.fillna(0)).clip(lower=0)`.
Nothing in the pipeline identifies which pixels became positive, so it is not recall. My recompute: 16 of 180 scale-1
cases have net < 0; the largest is **−7.3526 pp** (L03_s01, v2_nonorm_s1, blob 7, 160 µm / 4×: clean 15.93 % → injected
8.58 %). By arm: unnormalised 6 (three of them ≥ 0.85 pp), GLOBAL 2 (max 0.50 pp), z130 8 (all ≤ 0.34 pp). On the 2×
discs 22/180 are negative (max −1.85 pp). The maximum "recovery" is 1.0037 (pixel count vs geometric disc area, harmless).
The `fillna(0)` means a missing clean row would silently become a zero baseline; the saved data has no missing row, so
this is a latent defect, not a present error.

### Claim 2 — scale-1 injected footprint is 0.857 mm² per host, not 5.1 mm²
**Verdict: Confirmed.**
From the blob JSONs: Σ π·foot_um² over 20 discs = **0.856783 mm²**, identical in all three hosts (four blobs at each of
the five widths; footprint radii 14.14 / 28.28 / 56.56 / 113.12 / 226.24 µm; `foot_um / sigma_um = 1.665` as in
`halo_inject.py:68`). The CSV column `foot_area_um2` sums to the same value. Over three hosts: 2.570 mm². The footprint is
3.5–4.1 % of each host's tissue mask (20.8 / 24.7 / 20.8 mm²). The only simple arithmetic on these files that yields 5.1
is 6 × 0.856783 = **5.1407 mm²** = the sum of `foot_area_um2` over the scale-1 rows of one arm's six CSVs (three hosts ×
clean + injected). ADR line 1589 therefore counts the clean images as injected area, sums across hosts, and calls the
result "per host". Off by 6×.

### Claim 3 — clean baselines are non-zero in 39/180, > 1 % in 11, up to 15.93 % (unnormalised) and 3.41 % (GLOBAL)
**Verdict: Confirmed.**
Baseline > 0: **39/180**; > 0.1 %: 21; > 1 %: **11**; maxima **15.93 %** (v2_nonorm_s1, L03 blob 7), **3.41 %** (GLOBAL,
S17 blob 5), 0.59 % (z130). Nine of the eleven > 1 % cases are the unnormalised arm on the bright host L03 (tissue
median 363 counts). The median baseline is 0 in all arms. So "baselines were ~0 everywhere" (ADR:1576) is true for z130
(≤ 0.6 %) and for the dim host, and false for the unnormalised arm on the bright host, whose false-positive floor inside
individual discs reaches 16 % — larger than every recovery value in the unnormalised column of the ADR table.

### Claim 4 — at 1.25×/1.5×, 31/90 cases have positive net response, 11 exceed 1 %; z130 host-median at 40 µm is 1.14 % / 3.30 %
**Verdict: Confirmed, with a qualification the reviewer under-states.**
Weak-amplitude scale-1 cases: 90; net > 0: **31**; > 1 %: **11**; > 5 %: 5; net < 0: 8. Host-median z130 at 40 µm:
**1.142 %** (1.25×) and **3.295 %** (1.5×) — the only cells anywhere below 2× whose host-median exceeds 1 %
(per host: 4.10 / 1.14 / 0.00 and 0.88 / 3.77 / 3.30 %). Qualification: 20 of the 31 "positive" cases are below 1 %
(several are single-pixel-level, e.g. 0.001 %), and the largest weak responses are one blob — L03 blob 4 (20 µm, 1.5×),
which sits on a bright patch (local background 580 counts vs host median 363) and is seen by all three arms
(GLOBAL 28.2 %, z130 15.3 %, unnormalised 7.5 %). "Nothing is recovered" (ADR:1588) is false; the accurate sentence is
"host-median recovery ≤ 3.3 % of footprint below 2×, single blobs up to 28 %".

### Claim 5 — "whole-frame Aβ totals change by < 1 %" is wrong relative to the clean Aβ total (GLOBAL +3.4 / +16.4 / +32.6 %) though < 0.09 pp of frame area
**Verdict: Confirmed.**
`frame` is the full exported image incl. glass (`halo_eval.groovy:32`, `createFullImageAnnotation`). Relative to the clean
frame total: GLOBAL **+32.60 %** (L03), **+3.43 %** (S17), **+16.41 %** (W12); unnormalised −0.14 / +0.82 / **+1.40 %**;
z130 +0.16 / +0.22 / +0.25 %. In percentage points of frame area (dims 3766×5039, 4345×4895, 4319×4853 px at 1.3 µm):
max **0.0900 pp** (L03 GLOBAL, 28,858 µm² on 32.07 mm²). A third, more informative denominator — the injected footprint
(0.857 mm²) — gives GLOBAL 3.37 / 0.28 / 1.93 %, z130 0.14 / 0.12 / 0.43 %, unnormalised −0.07 / 0.04 / 0.005 %.
Under the only reading in which "Aβ totals" means Aβ totals, "< 1 %" is false for GLOBAL in every host and for the
unnormalised arm in W12. The sentence needs a denominator; see §3 for the replacement.

### Claim 6 — GLOBAL median response at 320 µm / 4× is 8.19 %, so "universal intensity blindness" is overstated
**Verdict: Confirmed.**
GLOBAL 320 µm / 4×: 8.19 / 0.65 / 8.27 % per host, median **8.186 %** (no clipping involved — all positive). GLOBAL also
responds at 80 µm / 4× in all three hosts (4.93 / 6.17 / 5.72 %) and at 40 µm / 4× (2.69 / 8.54 / 9.28 %); at 160 µm / 4×
it is 8.64 / 0.00 / 0.82 %. GLOBAL's annulus (2× disc minus footprint) gain in L03 is 7.3 k µm² against 21.5 k in the
footprints, i.e. the response extends into the Gaussian tails, as an intensity-driven response should. A plateau at 4×
background is called Aβ over ~5–8 % of its area by GLOBAL; that is weak, not invisible. The wording "intensity elevation
alone is invisible to every current arm" (ADR:1601) is overstated. See also Finding F1: most of the remaining plateau
pixels are not called `Negative` either.

### Claim 7 — FWHM vs σ: a FWHM-80 µm Gaussian retains ~93.6 % of its peak after subtracting a σ = 130 µm Gaussian local mean
**Verdict: Confirmed, and extended.**
The z130 classifier JSON contains `op.normalize.local` with `sigmaMean = 100.0`, `sigmaStdDev = 100.0` (pixels at
1.3 µm/px = 130 µm); the plot's x-axis is blob FWHM and the dashed line is at 130 (`halo_analyse.py:36`). For a 2-D
Gaussian bump (σ_s = FWHM/2.3548) minus its Gaussian-σ_n local mean, the centre retains σ_n²/(σ_n²+σ_s²):
**0.9957 / 0.9832 / 0.9361 / 0.7854 / 0.4779** at FWHM 20 / 40 / 80 / 160 / 320 µm. Numerically on a 1.3 µm grid
(scipy, truncate 6): 0.9361 and 0.4779. So mean subtraction alone predicts a fall-off between 160 and 320 µm, not at 80 µm.
Extension (F4 below): I applied the same mean-and-SD normalisation (σ = 100 px) to the real clean and injected Cy3
planes. The variance term matters for bright broad blobs — the local SD at the centre inflates from ~20–110 to
115–430 counts for 160–320 µm / 4× blobs because the blob's own residual enters the variance — so the centre z at
160–320 µm / 4× is only 1.3–2.4. But at 80 µm / 4× the centre z is 4.8–5.3 (footprint mean 2.2–2.8), comparable to
40 µm / 2× (3.8–6.8; recovered 4–22 %) and not far below 20 µm / 2× (6.3–6.6; recovered 44–56 %), yet recovery is
0.3–1.2 %. z130's collapse at 80 µm is therefore a property of the forest/feature stage (gradient, Laplacian and
weighted-SD features at σ 1.3–5.2 µm; a Gaussian's peak gradient scales as 1/σ_s, its Laplacian as 1/σ_s²), not of the
130-µm neighbourhood. The ADR's "falls off above its own neighbourhood exactly as its construction predicts" (ADR:1592)
is wrong on both scale and mechanism. (Caveat: my normalisation approximates QuPath 0.7.0's `LocalNormalization`;
all three classifier JSONs carry `inputPadding: 0`, so whether a σ = 100 px normalisation is effectively tile-bounded at
512 px was not verified and could add position dependence.)

### Claim 8 — all 30 cells of the 2×/4× table reproduce
**Verdict: Confirmed — all 30 checked, not a sample.**
Every host-median clipped recovery rounds to the ADR value (measured: unnormalised 0.1264 / 0.0726 / 0 / 0.0088 / 0 /
0 / 0 / 0 / 0.00005 / 0; GLOBAL 0.0108 / 0.2044 / 0 / 0.0854 / 0 / 0.0572 / 0 / 0.0082 / 0.0003 / 0.0819;
z130 0.5490 / 1.0037 / 0.1863 / 0.4116 / 0.0355 / 0.0045 / 0.0078 / 0 / 0 / 0). My recompute matches the saved
`halo_recovery.csv` to < 1e-16 in `recovery`, `baseline` and `frac`. Note the per-host spread behind the "clean" cells:
unnormalised 20 µm / 4× = 35.5 / 2.2 / 7.3 %; GLOBAL 160 µm / 4× = 8.6 / 0.0 / 0.8 %; z130 40 µm / 4× = 27.0 / 47.1 / 41.2 %.

---

## 2. Findings that neither the ADR nor the reviewer reports

**F1. The plateaus are not "called `Negative` whatever their height" — they are mostly re-labelled `Ignore*`.**
`Ignore*` fraction of a disc = 1 − (Aβ + Negative)/disc area. On clean images it is ≈ 0 (−0.4 … +1.5 %, rasterisation).
On injected images at 4×, host-median `Ignore*` fraction of the footprint:

| arm | 20 µm | 40 µm | 80 µm | 160 µm | 320 µm |
|---|---|---|---|---|---|
| z130 | 0 % | **59 %** | **82 %** | **64 %** | 19 % |
| GLOBAL | **80 %** | **55 %** | 3 % | 2 % | 19 % |
| unnormalised | **62 %** | 21 % | 31 % | 1 % | 1 % (but 82–84 % in the bright host) |

At 2×, z130: 44 / 68 / 55 / 22 / 2 %. Composition at 4× (host-median Aβ / Negative / Ignore*): z130 80 µm = 0.5 / 18 / 82 %;
z130 160 µm = 0 / 36 / 64 %; GLOBAL 20 µm = 20 / 0 / 80 %; GLOBAL 40 µm = 8.5 / 36 / 55 %. `Ignore*` is trained on "folds,
tears, edges, debris, out-of-focus patches" (`classifier_abeta.md:91`). So the forests do respond to a bright smooth
plateau: they classify it as artefact and abstain. The ADR's central mechanistic sentence (ADR:1597-1598) is factually
wrong, and the recovery metric silently counts abstention as non-detection — precisely what `normalisation.md` §5
forbids ("a near-total abstention must be recorded as a failure, never as a zero"). Biologically this is a different
failure mode from the one the ADR describes: a bright, smooth immunoreactive region in real tissue would be dropped as a
fold/debris by every arm (and, since QuPath excludes `Ignore*` from measurements — decisions.md:995 — from the tissue
denominator as well), rather than measured as `Negative`.

**F2. z130 is non-additive at frame level (negative halo).** In-footprint gain per host 3.3–3.5 k µm², frame-level gain
1.2 k (L03), 1.0 k (S17), 3.6 k (W12): losses of −389 / −653 µm² in the annuli and −1742 / −1835 / −333 µm² beyond the
2× discs, concentrated around the 320-µm blobs (e.g. S17 blob 16, 320 µm / 4×: −1031 µm² in the annulus). A bright broad
addition switches off existing positives several hundred µm away, as a σ = 130 µm mean/SD normalisation must. GLOBAL is
additive (frame gain ≈ 2×-disc gain; outside ≈ 0). The unnormalised arm on the bright host: footprints −2.4 k, annuli
+1.8 k, frame −0.6 k µm² — adding 20 blobs of signal *reduced* its total Aβ area.

**F3. Host-brightness dependence — the plan's third diagnostic ("recovery versus host brightness should be flat",
PLAN §"What is measured") was never reported.** z130 20 µm / 4×: 100.4 / 100.4 / 100.4 %; 20 µm / 2×: 54.9 / 43.6 / 55.7 %
— flat. GLOBAL 20 µm / 4×: 25.8 / 20.2 / 20.4 % — flat; 160–320 µm / 4× erratic (8.6 / 0.0 / 0.8; 8.2 / 0.7 / 8.3).
Unnormalised 20 µm / 4×: 35.5 / 2.2 / 7.3 %; 40 µm / 4×: 17.0 / 0.9 / 0.0 % — the bright host dominates, and all its
large baselines and losses are in that host. The unnormalised arm fails the design's own flatness criterion; z130 passes.

**F4. z-score diagnostic** (see Claim 7): host-median local z gained at the blob centre under z130-style normalisation,
rows FWHM 20 / 40 / 80 / 160 / 320 µm, columns 1.25 / 1.5 / 2 / 4×: 1.3 / 2.3 / 6.6 / 13.3; 1.7 / 2.8 / 5.4 / 8.7;
1.3 / 2.5 / 4.0 / 5.1; 1.2 / 2.3 / 2.4 / 2.4; 0.7 / 0.9 / 1.4 / 1.3 (`zscore_diagnostic.csv`). In GLOBAL's own units the
injected height is the same at every width by design (median 0.3 / 0.7 / 1.2 / 3.7–4.9 tissue units), so GLOBAL's
fall-off with width is entirely a feature/forest effect.

**F5. Design fragility.** Each table cell is the median of three blobs at three arbitrary placements; local background
differs (e.g. 580 vs 363 counts), so nominally equal "1.5×" blobs differ 3.5-fold in raw added counts (36–290). Single-host
responses up to 28 % vanish in the medians. Amplitude and location are confounded, as the reviewer says.

---

## 3. Sentences in ADR-0035 / the figure README that must be rewritten

1. ADR:1576 "Baselines were ~0 everywhere." → "Baselines were 0 in 141/180 discs and ≤ 0.6 % for z130; on the clean bright host the unnormalised arm reads up to 15.9 % and GLOBAL up to 3.4 % of a disc."
2. ADR:1574-1576 recovery definition → add "clipped at zero; 16/180 discs lost positive area (max 7.4 pp) and appear as 0; signed values are in `halo_recovery.csv`" (after adding them).
3. ADR:1588 "Nothing is recovered by any arm at 1.25× or 1.5×." → "Host-median recovery is ≤ 3.3 % of footprint below 2× (z130 at 40 µm: 1.1 % / 3.3 %); single blobs reach 28 %."
4. ADR:1588-1589 "Whole-frame Aβ totals change by < 1 % after injecting 5.1 mm² of footprint per host." → "0.857 mm² of footprint per host (3.5–4.1 % of tissue). Whole-frame Aβ area rises by ≤ 0.09 percentage points of frame area (≤ 3.4 % of the injected footprint); for GLOBAL that is +3 % to +33 % of its small clean false-positive total."
5. ADR:1595-1598 "All three forests use the same texture features … a smooth intensity plateau has no texture at those scales and is called `Negative` whatever its height." → "The Gaussian feature is an intensity feature; at 4× most plateau pixels leave `Negative` for `Ignore*` (z130 59–82 % of the footprint at 40–160 µm, GLOBAL 80 % at 20 µm), i.e. the forests treat a bright smooth region as a fold/debris artefact and abstain."
6. ADR:1592 "z130 falls off above its own neighbourhood exactly as its construction predicts" → "z130's recovery collapses at 80 µm FWHM (σ 34 µm), well inside its 130-µm normalisation σ; mean subtraction keeps 94 % of such a blob, and the collapse is a feature/forest effect. The normalisation's variance term does suppress 160–320 µm blobs (local z ≈ 1–2.5) and produces negative halos that remove positives up to several hundred µm away."
7. ADR:1594-1595 "the limitation is not the normalisation" → keep for GLOBAL/unnormalised; qualify for z130 as in 6.
8. ADR:1601-1602 "It settles that intensity elevation alone is invisible to every current arm, and that z130 is additionally blind to anything broader than ~80 µm." → "Smooth Gaussian plateaus at ≤ 4× background are scored as Aβ over at most ~8 % of their area (GLOBAL) and ≤ 4 % (other arms) at ≥ 80 µm, and are largely labelled `Ignore*`; z130's recovery falls from 100 % (20 µm) to 41 % (40 µm) to ≤ 4 % (≥ 80 µm)." No hard size cut-off was tested.
9. Figure README lines 33-36: "No arm recovers blobs below 2x background" → "≤ 3.3 % median"; "nothing from 80 um up" → "≤ 4 % median (single blobs to 5.5 %)"; "fully" → "median 100 % in three of three hosts". The dashed 130-µm line sits on an FWHM axis: label it as σ, move it to the equivalent FWHM (306 µm), or drop it.
10. The 30-cell table stands. Add the `Ignore*` fraction and the per-host values alongside it.

---

## 4. Judgment

1. The experiment is sound as a paired perturbation test: images, placements, measurements and the 30-cell table all reproduce, and the reviewer's numerical corrections are all correct. What the ADR gets wrong is the prose around the table, and in one place the mechanism.
2. What it establishes about the arms: (a) all three forests respond to compact (20–40 µm) high-contrast Gaussian bumps far better than to broad ones — z130 best (100 % / 41 % at 4×, 55 % / 19 % at 2×), GLOBAL weakly (20 % / 9 %), the unnormalised arm weakly and mainly in the bright host; (b) at ≥ 80 µm every arm scores ≤ 8 % of a plateau as Aβ, and most of the plateau is re-labelled `Ignore*`, i.e. artefact — the forests are not blind to the intensity, they reject it; (c) z130's response is brightness-invariant across hosts and the unnormalised arm's is not; (d) z130 removes positives around bright additions (negative halo), so it is non-additive.
3. It establishes nothing about textured diffuse deposits, nothing about which arm is right on real tissue, and no hard object-size limit for z130; the observed 80-µm collapse is a property of the trained forest and its 1.3–5.2 µm feature scales, not of the 130-µm normalisation.
4. Do the reviewer's corrections change the scientific conclusion? No. "Smooth plateaus are poorly detected as Aβ by all arms" and "z130 detects compact blobs best" both survive every one of the eight points; those are wording and number fixes. My F1 finding does change the *mechanism* sentence (abstention as artefact, not "Negative whatever its height") and the "not the normalisation" sentence for z130 (F2/F4), but it strengthens rather than weakens the practical caveat: a bright smooth immunoreactive region would be dropped as a fold by every arm, so no arm's absence of confluent Aβ can be quoted as biology.
5. The reviewer's "better stress test": do the free part now — report signed net, raw fractions, baselines, `Ignore*` fraction, per-host values, and fail on missing baselines (≈ 10 lines in `halo_analyse.py`, no QuPath). A small textured-patch run with matched locations and an out-of-footprint readout (the ADR's own follow-up) is worth a day. The full factorial (texture × width × contrast × location × anchor mode) is not worth doing before human masks exist: it can only characterise responses to synthetic stimuli and cannot decide which arm is right about real diffuse deposits, while the reference boxes for six non-training animals already exist and need only a tracer. Trace first; run the factorial only if the tracing shows diffuse deposits on which the arms disagree.
