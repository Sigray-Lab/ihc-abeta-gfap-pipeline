# Decision log (ADRs)

Two parts. **Part 1** records decisions already taken and why — one entry per decision,
append-only, never edited in place. If a decision is reversed, add a new entry that
supersedes the old one and mark the old one Superseded. **Part 2** is the register of
decisions reserved for the PI, mirrored from `docs/DECISION_REGISTER.md` in the project
folder, so that the repository states its own open questions.

Spec `CLAUDE_v1.2.md` §12 draws the line: scientific choices belong to humans,
technical choices belong to whoever writes the code but must be recorded here.

---

# Part 1 — Decisions taken

## ADR-0001 — Storage split: raw in Dropbox, git and scratch outside it

**Status** Accepted · **Date** 2026-08-06 · **Type** technical

Raw data stays in Dropbox (shared, backed up, read-only original). The git repository,
all derived artefacts (`~/ihc_work`), all results (`~/ihc_work/results`) and the
blinding key (`~/ihc_custodian`, mode 700) live on local disk outside Dropbox.

**Why.** Dropbox syncs file-by-file with no notion of a commit, so a repository inside
a synced folder eventually acquires "conflicted copy" files inside `.git`, corrupting
the object store. A full run writes tens to hundreds of gigabytes of masks, overlays,
caches and tiles, which would saturate the sync queue and can dehydrate the raw data
the run is currently reading. The custodian tree must additionally be un-pushable and
un-shareable by accident.

**Consequence.** Every path is resolved through `config/paths.yaml`, so a move is a
one-file edit. `./ihc doctor` fails if `work_root`, `results_root` or `custodian_root`
resolves to a path containing "Dropbox". Derived artefacts are never backed up in
preference to re-running.

**Accepted cost.** Raw data remains subject to Dropbox dehydration, which is why
`./ihc doctor` checks `st_blocks` against `st_size` on every raw file before any run.

---

## ADR-0002 — Box membership is derived from stage coordinates, not section number

**Status** Accepted · **Date** 2026-08-06 · **Type** scientific, taken in `CLAUDE_v1.2.md` §2

Sections are assigned to a PAP-pen box by sorting them on stage X (VSI tag 2018,
micrometres) and splitting at the largest gap. The low-stage-X box is `near_label`.
The split must yield 2+2, 2+1 or 1+2, and the between-box gap is 1.7–2.4× the
within-box gaps in every animal checked.

**Why.** Spec v1.1 asserted that `_01`+`_02` and `_03`+`_04` always pair. **That is
false.** In tube 49 the sections sharing a condition are `_01` and `_04`: sorted by
stage position the layout is `[_02 _03]` — 14.65 mm gap — `[_01 _04]`, and the measured
GFAP/Aβ signal groups accordingly. The section number is **acquisition order**; the
operator scanned the physically-third section first.

**Consequence.** Series-to-section mapping is verified, never assumed. Stage
coordinates are recorded in the manifest. A split that is not 2+2, 2+1 or 1+2 is a hard
failure, not a warning.

**Supersedes** the pairing rule stated in `CLAUDE_v1.1.md`.

---

## ADR-0003 — Condition comes from `slides.csv`; pixels audit it and never relabel

**Status** Accepted · **Date** 2026-08-06 · **Type** scientific, `CLAUDE_v1.2.md` §2 non-negotiable

Positive/negative status is read from `config/slides.csv`. The image-based check
(GFAP 99.9th percentile within tissue: positives 3358–5736, negatives 196–388 — 8.7×
with no overlap) runs as an **auditor only**. A contradiction blocks and is escalated
to a human; it never relabels a section.

**Why.** Deriving condition from pixels would silently convert a failed positive stain
into a fake negative control, and would make negative-control verification circular —
the controls would have been selected for the very property being tested.

**Consequence.** `slides.csv` is version-controlled and validated by
`config/validate_config.py`. Rows flagged `needs_confirmation` are reported loudly and
are not usable until confirmed at the bench. Tube 49 currently carries such a flag.

---

## ADR-0004 — The near GFAP band starts at 10 µm, not 0

**Status** Accepted · **Date** 2026-08-06 · **Type** scientific, `EXECUTION_PLAN_v3.md` §3

Distance bands are: Aβ-positive area (its own compartment), >0–10 µm (reported,
flagged bleed-through-suspect, excluded from the primary contrast), **10–25 µm (the
primary near band)**, 25–50 µm, and >50 µm (distant baseline). Primary contrast is
near minus distant in percentage points; the ratio is secondary and is undefined —
never pseudocounted — when the distant denominator is too small.

**Why.** Profiling GFAP against distance from Aβ deposits on six positive sections
from three animals gave a monotonic decline **peaking inside the deposit**: in-plaque
1.26–2.85× the far-field baseline, 0–10 µm 1.11–1.97×, 10–25 µm 1.02–1.64×, 25–50 µm
1.00–1.42×. A peak inside the deposit is the signature of optical bleed-through; real
astrogliosis would be expected to dip at the amyloid core and peak in a ring outside
it. The test could not separate the two readings — against bleed-through, the ratio of
GFAP excess to Cy3 brightness varies threefold between sections (5.3–15.7 %), whereas
bleed through a fixed filter should be a constant fraction, and plaque-infiltrating
astrocyte processes are well documented.

**What is settled** is the part that matters: **the elevation extending to 25–50 µm
cannot be optical**, because there is no Cy3 signal at that distance to bleed from.
That component is genuinely biological and the measure is viable.

**Consequence.** The primary enrichment contrast must not sit on the in-plaque or
0–10 µm compartment. This supersedes the spec's proposed >0–25 µm near band and
D-8's "keep >0–25 / >25–50". The whole enrichment measure stays **exploratory** until
bleed-through is characterised further.

**Supersedes** the band edges in `CLAUDE_v1.2.md` §8 and in D-8 of the decision register.

---

## ADR-0005 — Lossy JPEG2000 is recorded as a limitation, not corrected

**Status** Accepted · **Date** 2026-08-06 · **Type** technical + methods

The raw tiles are lossy and stay lossy. No attempt is made to compensate.

**Why.** Codestream inspection across ~47,000 tiles gives COD transform = 0 (9/7
irreversible wavelet), scalar-derived quantisation, 2 quality layers, quality 98. The
ETS compression code 3 means "JPEG 2000" and carries no lossless/lossy information —
reading code 3 as lossless would have been an easy and invisible mistake. The loss is
irreversible; there is nothing to correct.

**Consequence.** It goes in the methods. It is uniform across the cohort, which removes
variation in codec *settings* — but compression artefacts are **signal-dependent**, so
ringing at sharp bright edges may still interact with exposure and intensity.
Validation therefore stratifies by exposure and by burden rather than assuming the
codec is neutral.

---

## ADR-0006 — Rescan tubes 51 and 60 rather than rely on numerical correction

**Status** Accepted · **Date** 2026-08-06 · **Type** resource, PI (D-1)

Tubes 51 and 60 were re-acquired at standard exposure (`RawData/Rescan/`, 2 tissue
series each). The originals are retained and corrected numerically as a cross-check,
as a pre-specified sensitivity analysis.

**Why.** Exposure was swept across all 31 index files and validated against
Bio-Formats on the 8 animals with payloads (8/8 exact). Exactly two deviate: tube 51
(DAPI 60.52 / FITC 240.82 / Cy3 145.68 ms) and tube 60 (Cy3 397.93 instead of 1840).
Autoexposure was off and lamp intensities identical, so these were manual settings, and
the metadata is real rather than a stale preview value — glass background outside the
tissue scales with the recorded exposure almost exactly.

`I* = I_raw / t_exposure` rescales but does not restore photons: on slide 51 the tissue
sits ~7 grey levels above background where on slide 29 it sits ~208, so quantisation
and SNR remain worse whatever the arithmetic. Two slides is ~1 hour of scope time,
comfortably within what the wet-lab scientist said they would rescan.

**Consequence.** Rescanned animals are flagged in the manifest. Exposure is recorded
per animal per channel. The exposure check is assignment-free: dump every exposure
record in file order and require N identical triplets, where N is the number of tissue
series — never collapse to a set, which cannot distinguish "a value is missing" from
"two channels share a value". That flaw sank the first sweep.

---

## ADR-0007 — Undecided config values are a sentinel string, not a default

**Status** Accepted · **Date** 2026-08-06 · **Type** technical

Every parameter nobody has decided yet is the string `PENDING_PI_DECISION` in
`config/config.yaml`, with a comment naming the decision it waits on.

**Why.** A plausible default is the dangerous option: a stage would run, produce a
number, and nothing would indicate that the number rests on an unconsidered choice. A
string where a number belongs raises `TypeError` at first arithmetic use. Loud and
immediate beats silent and wrong.

**Consequence.** `./ihc check-config` lists every sentinel and every block still marked
`pi_approved: false`. A stage that depends on one must refuse to run. Blocks with
`pi_approved: false` carry a concrete recommended value and may be used for blinded
pilots and engineering work, but not for any number that reaches the manuscript.

---

## ADR-0008 — No packaging step; the entry command manipulates `sys.path`

**Status** Accepted · **Date** 2026-08-06 · **Type** technical

`./ihc` inserts `src/` on `sys.path` and imports `ihc.*` directly. There is no
`setup.py`, no `pyproject.toml`, no `pip install -e .`.

**Why.** The pipeline is handed to a wet-lab scientist. "Clone it, create the conda
env, run `./ihc doctor`" is a shorter and more reliable instruction than anything
involving an editable install, and one fewer thing to go wrong on the cluster. The
repository is an application, not a library — nothing else imports it.

**Consequence.** Scripts under `scripts/` and tests must do the same insert, or run
from the repository root. If this ever becomes a library, revisit.

---

## ADR-0009 — `config/slides.csv` is version-controlled, including its `group` column

**Status** Accepted · **Date** 2026-08-06 · **Type** technical, with a blinding caveat

`slides.csv` lives in git. It carries `group` (treatment arm), which is unblinding
information.

**Why.** The spec makes `config/slides.csv` the authoritative source of condition, so
it must be versioned and diffable alongside the code. The `group` column leaks nothing
that `tube_id` does not already leak — tube IDs run in contiguous treatment blocks, so
anyone holding the file can infer group regardless.

**Consequence.** The blinded analysis manifest must strip `tube_id`, `group`, `arm`,
original path and acquisition order (`blinding.strip_from_blinded_manifest` in
`config/config.yaml`). Coded IDs come from a random permutation with a recorded seed
held by the custodian — never sequential in file-iteration order, never a hash or
arithmetic transform of `tube_id`, because any order-preserving scheme reproduces the
group structure exactly. Anyone training classifiers, drawing regions or annotating
validation data works from the blinded manifest and should not open this file.

---

## ADR-0010 — Never-acquired tiles are missing support, not background

**Status** Accepted · **Date** 2026-08-06 · **Type** technical, forced by the data

The percent-area denominator carries an `acquired_support` term. Tile positions inside
the bounding box that were never acquired are excluded from the denominator rather than
counted as negative tissue.

**Why.** 5–11 % of tile positions are absent because the scanner used a sample mask.
Counting them as background would inflate the denominator by up to a tenth and deflate
every percent-area number, non-uniformly across sections.

**Consequence.** This is not a PI decision and has no sentinel; it is forced. The
verifier reports the sparsity fraction per series, and a fraction above 0.25 — far
outside the observed 0.05–0.11 band — is a failure.

---

## ADR-0011 — Four subcommands, and the surface stays that size

**Status** Accepted · **Date** 2026-08-06 · **Type** technical

`./ihc doctor | check-config | verify | meta`. New functionality goes into `src/ihc/`
and, if it needs a human-facing runner, into `scripts/` — not into a fifth verb.

**Why.** Spec §1: "a single documented entry command". A CLI that grows a verb per
function becomes a second, undocumented API, and the person who inherits this reads the
README once.

**Consequence.** Each subcommand prints its plan and its output directory before acting.
Exit codes are fixed: 0 ok, 1 problems, 2 usage, 3 a required module is not written yet.

---

# Part 2 — Open decisions reserved for the PI

Mirrored from `docs/DECISION_REGISTER.md` in the project folder. That file is the
source of truth for the reasoning; this table is the version the repository ships so
that the code states its own open questions. Each row names the config key that holds
`PENDING_PI_DECISION` until it is answered.

**Legend** 🔴 only the PI · 🟡 needs real thought · 🔵 recommendation exists, needs one line · 🟢 mine to decide · ⚪ deferred

| ID | Decision | Owner | Config key | Recommendation on the table |
|---|---|---|---|---|
| **D-11** | **Δ — the smallest treatment effect that would matter** | 🔴 | `validation.delta_smallest_meaningful_effect_relative`, `validation.acceptance_thresholds_percentage_points` | Pre-register **20–30 % relative reduction**, converted to percentage points from our own blinded pilot baseline rather than from the literature. Everything keys off this: validation acceptance, the classifier release gate, the negative-control threshold. |
| **D-14** | **Blinding key custodian** | 🔴 | `blinding.custodian`, `blinding.seed` | Name a person. Must not be the wet-lab scientist (they stained and imaged, so they can infer group) and must not be whoever trains classifiers, draws regions or annotates. |
| **D-1** | **Rescan slides 51 and 60?** | 🔴 | — | **Answered — see ADR-0006.** Rescans are in `RawData/Rescan/`. |
| **D-3** | **Aβ positive-class definition** | 🟡 | `abeta.positive_class`, `abeta.objects.vascular_amyloid_policy` | Annotate four classes (parenchymal extracellular, vascular, intracellular/ambiguous, background/artefact); use **parenchymal extracellular** for the primary endpoint and the peri-plaque GFAP mask; fall back to total immunoreactive area if experts cannot separate the classes reliably. Turns on whether CAA is scientifically interesting here or purely a nuisance. |
| **D-4** | **Region definitions** | 🟡 | `regions.*` | Whole hippocampal formation with subiculum excluded; full-thickness isocortex, white matter excluded; hemispheres pooled by summing numerator and denominator. Start with whole structures, not subfields, unless there is a subfield hypothesis. (Partly resolved already: **there is no frontal cortex** — the endpoints are hippocampal formation and isocortex at the sampled level.) |
| **D-2** | **Percent-area denominator** | 🔵 | `denominator.*`, `masks.tissue.method`, `masks.artefact.method` | `ROI ∩ acquired-support ∩ tissue − artefact`. The acquired-support term is forced by the sparse tile grid (ADR-0010). Export every component separately so any alternative denominator is recomputable without re-segmenting. |
| **D-7** | **Section weighting / hemisphere pooling** | 🔵 | `aggregation.*` | Area-weighted (Σ positive area / Σ valid area) as primary; mean-of-section-percentages as sensitivity. Supported by the wet-lab scientist: the two sections in a box are a **depth check**, not independent replicates. Pool hemispheres by summing. |
| **D-8** | **GFAP distance bands and boundary rules** | 🔵 | `gfap_enrichment.*` | Band edges **revised by ADR-0004** — near band starts at 10 µm. 25 µm has direct support in this model (astrocyte increases extend to ~25 µm from plaques, Tomlin et al. 2025, *Brain Commun*). Still open: minimum band area, whether deposits outside the ROI may generate bands inside it, and how bands crossing region boundaries are handled. Measure stays exploratory. |
| **D-9** | **Plaque object rules** | 🔵 | `abeta.objects.*` | Connected components, no aggressive splitting, provisional min 20 µm² with sensitivity at 10 and 50 µm², centroid-in-ROI counting, edge-touching objects excluded from size summaries but their pixels retained in percent area. Secondary endpoint; should not block anything. |
| **D-10** | **Negative-control failure rule** | 🔵 | `negative_control.*` | Threshold at `min(0.10 pp, Δ/10)` — so it cannot be set until D-11 is. Failure triggers visual review and a flag, never automatic exclusion. **Sub-rule needed:** tubes 35, 38, 45 and 53 are double-positive with **no negative control at all**, so the gate is defined per animal-where-negatives-exist and needs an explicit rule for animals that have none. |
| **D-13** | **If the global classifier fails the stratified gate** | 🔵 | — | Fixed hierarchy: improve on development animals → retest on the untouched reserve → if the failure is specifically *region*-dependent, allow separately locked per-region models (accepting loss of absolute cortex-vs-hippocampus comparability) → if it is exposure-, group- or burden-dependent, do **not** introduce stratum-specific classifiers; restrict the endpoint instead. |
| **D-5** | Classifier resolution | 🟢 | `classifier.resolution_downsample.*` | Blinded pilot, two candidates: Aβ at ds2, GFAP at ds1. Different resolution per marker is fine — "frozen" applies *within* a marker. |
| **D-6** | Animal split mechanics | 🟢 | `validation.split_animals` | ~6 development / 9 validation / 4 reserve / rest production-only, custodian-balanced. **All animals still enter the biological analysis** — training costs validation capacity, not sample size. |
| **D-15** | QuPath / environment pin | 🟢 | `env/tool_versions.yaml` | QuPath 0.7.0 + ABBA 0.5.0 + Bio-Formats 8.5.0 on Java 8. arm64 first. |
| — | Missing-section handling | 🟢 | `sections.missing_section_rule` | Use whatever valid positive sections exist, flag animals contributing only one, no imputation. Tube 42 (one positive section) is the live case. |
| **D-12** | Statistical model and multiplicity | ⚪ | — | Out of scope for the image pipeline; inference comes later. The one thing worth capturing **now**, because it decays: the original treatment-allocation procedure — was diet assigned per animal or per cage? That determines the experimental unit. On the list for the wet-lab records. |

**Summary from the register:** genuinely need the PI — 3 (Δ, a custodian name, rescan
authorisation, the last now answered). One line each — 6. Think properly about — 2.
Mine to handle — 4.

**Still the only hard blocker for stage 3:** `config/slides.csv` from the wet-lab bench log.
A copy is now in the repository and passes structural validation, with **tube 49
flagged `needs_confirmation`** — the L/R cell says `near_label` but the annotation says
`01+04`, which is `far_label`, and the image agrees with the annotation. That row is
not usable until the wet-lab scientist confirms it.

---

## Open technical questions (not PI decisions)

| Question | Where it bites |
|---|---|
| The FITC emission filter is unrecorded in all 31 files. | Methods only, no computation depends on it. Ask the imaging core. |
| Excitation lines conflict between sources: the instrument documentation gives an Excelitas X-Cite NOVEM with lines at 405, 485, 551, 584, 639, 730 nm; `CLAUDE_v1.2.md` §5 quotes the KI BIC facility document as 385, 430, 475, 545, 635, 735 nm. | Methods only. Both are recorded in `config/channels.yaml` with the conflict flagged. Do not quote either set until the core confirms which describes this unit. |
| A wet-lab channel table (DAPI 440±40, FITC 521±21, Cy3 607±34) does not match the in-file values. | Rejected as a generic VS200 spec, recorded in `config/channels.yaml` so nobody rediscovers it and mistakes it for data. |
| Payloads for 23 of 31 animals have not been transferred. | `./ihc verify` and any pixel work can only cover the 8 animals present plus the 2 rescans. `./ihc meta` works on all 31, because the index carries the metadata. |

---

## ADR-0009 — Blinding key custodian: the PI

**Date:** 2026-08-06 · **Status:** accepted

The PI holds the code key and the random seed. the wet-lab scientist performs region delineation and
must therefore stay blinded to group.

Consequence for storage: the custodian directory is `~/ihc_custodian` (mode 700),
**outside Dropbox and outside the git repo**, so it is not in the folder the wet-lab scientist works
from. Blinded QuPath projects may live in shared storage; the key may not.

Note the residual limit, stated rather than glossed: the wet-lab scientist performed the staining and
imaging, so he can in principle infer group from tube identity if he ever sees it. The
blinded projects therefore carry coded IDs only, exclude the slide-label image (which
shows the tube ID as text *and* a DataMatrix barcode), and rename internal series names
(tube 60's are prefixed `60_`).

## ADR-0010 — Region delineation: try atlas registration first, hand-drawing as fallback

**Date:** 2026-08-06 · **Status:** accepted

Attempt ABBA atlas registration on a small number of sections first. If it needs heavy
manual correction, fall back to hand-drawn QuPath regions. This is a half-day test, not
a commitment — the decision is made on the result, and recorded per section so that
fallback frequency can later be checked against treatment group.

the wet-lab scientist performs the delineation, blinded. The marker channels (FITC/Cy3) must be hidden
while boundaries are drawn, so that visible plaque burden cannot influence where the
boundary goes.

## ADR-0011 — Tube 49 positive box: CONFIRMED far_label (sections 01 + 04)

**Date:** 2026-08-06 · **Status:** accepted · **Source:** the wet-lab scientist, in writing

Tube 49 genuinely deviates from the usual 01+02 / 03+04 pairing: sections **01 and 04
are the same PAP-pen box**, and that box is the **right / far-from-label** one. The `L`
tick on scan sheet 1 was only meant to draw attention to 01 being the first positive;
sheet 2's `R` is correct.

This independently confirms the stage-coordinate rule: box membership follows physical
position, never section number. Record, geometry and pixel signal now agree 8/8.

## ADR-0012 — Rescanned slides double as an imaging-repeatability check

**Date:** 2026-08-06 · **Status:** accepted · **Source:** the wet-lab scientist's suggestion

Tubes 33, 42, 49, 51, 54 and 60 exist as two scans of the *same physical slide* on
different days. That is a free estimate of **technical repeatability**: how much a
percent-area number moves for imaging reasons alone, with the biology held constant.

Worth having, because it bounds the measurement noise independently of any classifier
validation, and it is the natural denominator for judging whether a treatment
difference is large enough to matter. Report it as a QC output.

Caveat: the two scans are not interchangeable. Fluorophores photobleach between
sessions, so a systematic drop on the second scan is expected and is not measurement
error. Tubes 51 and 60 additionally changed exposure between scans by design, so only
33, 42, 49 and 54 are clean repeatability pairs.

## ADR-0013 — `manifest` is a subcommand, and `condition` has three states

**Date:** 2026-08-07 · **Status:** accepted · **Type:** technical
**Supersedes** the four-verb rule in ADR-0011.

`./ihc manifest` builds the section manifest. It is a verb rather than a script under
`scripts/` because it is a named build stage in the spec (§6 step 2) that every later
stage reads, and because it is where `condition` is derived — the one column whose
corruption produces numbers that are wrong and completely plausible. A step nobody can
find is a step somebody skips.

**Three states, not two.** `condition` is `positive`, `negative` or **`unresolved`**.
The third exists because the alternative to representing "we do not know" is guessing,
and a guess here is indistinguishable from knowledge downstream. A slide is unresolved
when `slides.csv` carries a non-empty `needs_confirmation`, when its `positive_box` is
unreadable, when box membership cannot be derived from stage X, or when the tube has no
`slides.csv` row at all. Unresolved rows are kept in `manifest.csv` (they are provenance)
and are absent from `manifest_analysis.csv` by construction, so a later stage cannot
quantify them by forgetting to filter. Tube 37 is the live case.

**`has_negative_control` is a property of the animal, not of a scan.** The rescans carry
only the positive box, so a per-scan flag would report "no negative control" for tube 51
while two negative sections of tube 51 sit on the original slide. D-10 defines the gate
per animal-where-negatives-exist, so animal-level is the semantics every consumer wants.

**Rescan sections are matched to the original by stage X, never by section number.** The
rescan renumbers its sections `01`/`02` whatever they were on the original slide, so
tube 51's rescan `01` is the original's `03`. Matching by label would bind the wrong
physical section and, because both rescanned slides are `far_label`, the opposite
condition. Nearest-neighbour on stage X with a tolerance of a quarter of the slide's
smallest inter-section spacing; observed offsets 82–349 µm against a 6.7 mm spacing. A
match that fails the tolerance yields `unresolved` rather than a guess.

**Consequence for every consumer.** `(tube_id, section_label)` is *not* a key — tubes 51
and 60 appear twice. Key on `(tube_id, scan, section_label)`; pair the two scans of one
physical section on `(tube_id, physical_section_label)` for the ADR-0012 repeatability
check; and filter on `scan_is_preferred` before aggregating, blinding, or building a
QuPath project, or those two animals contribute twice.

## ADR-0014 — The image cross-check is a separate function that returns a report

**Date:** 2026-08-07 · **Status:** accepted · **Type:** technical, implements ADR-0003

`crosscheck_condition_against_pixels()` returns a DataFrame. It has no write path to
`condition`, and the manifest builder never calls it. This is the code-level expression
of ADR-0003 and is deliberate structure rather than convention: a function that *could*
write the column would eventually be asked to.

Statistic: at pyramid level 3 (downsample 8), tissue is thresholded on DAPI at the
geometric midpoint between the 25th and 99.9th percentile over the acquired support —
scale-free, so it survives tube 51's 2.1× lower DAPI exposure. The 99.9th percentile of
GFAP (FITC) within tissue is then divided by the 99th percentile of DAPI over the same
support. Dividing by DAPI cancels the per-slide gain that would otherwise let an
exposure difference look like a staining difference; the 99th rather than the 99.9th
percentile is used in the denominator because DAPI carries occasional saturated debris
that would swamp the extreme tail. Never-acquired tile positions are excluded throughout
(ADR-0010).

Measured over all 34 sections with payloads: negatives 0.57–1.41, positives 4.86–20.07,
smallest within-slide separation 8.0×. Thresholds are set at ≤ 2.0 negative and ≥ 3.0
positive, leaving an explicit `inconclusive` band that nothing in this cohort occupies.
34 of 34 sections agree with `slides.csv`.

Cy3 is reported as `abeta_index` but does **not** feed the verdict. Amyloid burden is
biology and varies legitimately between animals, whereas GFAP immunoreactivity in a
no-primary section is near-zero by construction; using Cy3 would confuse a low-burden
animal with a failed stain.

---

## ADR-0015 — Blinding scope: keep it simple, accept the residual

**Date:** 2026-08-06 · **Status:** accepted

**Decision: blinding here is one mechanism, not a programme.** Every animal gets a random
code; the key lives with the PI; nobody else sees a tube ID. That is the whole thing.

The mechanism was checked once and works: no column or column pair in the blinded
manifest recovers the group partition, and nothing orders the animals by tube ID. The
check is available as `./ihc blind --audit` if reassurance is wanted; it is not a gate
and it is off by default.

**Residual risk, accepted and not mitigated further.** Roughly 4 of 30 animals could in
principle be identified by someone who both holds the staining record and goes looking,
because the slide structure (three sections, or a double-positive) is itself distinctive.
This is not worth engineering around. The person drawing regions is not tracking which
slide is which, has a great deal else to do, and would gain nothing from knowing. More
obfuscation would add friction to a task that is already the project's bottleneck, in
exchange for closing a route nobody is walking down.

What actually protects the measurement is cheaper and already in place: **boundaries are
drawn on DAPI alone with the marker channels hidden**, so even a fully de-blinded animal
cannot have its region nudged toward a burden the drawer cannot see.

Note it as a one-line limitation in the methods and move on. This is an exploratory study;
the risk of over-engineering the process exceeds the risk of the leak.

---

## ADR-0016 — Adversarial review 2026-08-07: what was accepted, what was not

An eight-reviewer adversarial pass raised 17 red flags. Verdict: *do not proceed to
measurement yet, but the engineering is sound*. What changed:

**Accepted and fixed.**

- **The manifest had no correct filter.** `analysis_include` alone double-counted tubes
  49, 51 and 60 (importing the bad-exposure scans the rescans exist to replace);
  `scan_is_preferred` alone deleted tubes 33, 42 and 54 outright — three animals from
  three different groups. Root cause: preference was decided from the *existence* of a
  file in `Rescan/`, before anything knew whether that rescan resolved. Now computed
  after conditions are known, **per physical section rather than per tube**, exposed as
  a single `use_for_measurement` column, with an invariant that fails if any animal
  vanishes, loses its positives, or double-counts.
- **`Rescan/` was never six rescans.** For tubes 33, 42 and 54 it holds a section the
  first scan never captured — confirmed at the bench: *"because I had some extra time on
  the booking slot I put it also the ones I had one positive."* Those sections are now
  placed by stage position relative to the box split, given their own physical identity
  so they cannot merge with an original, and carried as positives. **The three animals
  that had one positive section now have two.** That is the single largest improvement
  the dataset had available.
- **The bleed-through direction was backwards, and the 10 µm band is reverted.** See
  ADR-0017.
- **The real blinding leak was `group` in committed `slides.csv`** — plain text in every
  clone and reflog, while `.gitignore` carefully protected the key. Moved to the
  custodian tree; the manifest now reads it from there and yields an empty group, i.e. a
  blinded manifest, for anyone without it. Re-seeded with `secrets.randbits(128)`; the
  previous seed was date-shaped, which is ~365 guesses given a public cohort roster.
  **Correction (ADR-0019): the guard claimed here was not actually implemented until
  2026-08-08.**
- **The provenance manifest was not the master record** it was described as — it read
  the analysis subset and silently dropped bench-excluded rows. It now reads the full
  manifest: 130 of 130.
- **Tube 37 settled by data.** Its payload arrived; the auditor gives 13.1 / 13.6 for
  sections 03 / 04 against 0.88 / 0.87 for 01 / 02. `far_label` confirmed.
- **Tests froze cohort facts that legitimately change.** Payload membership, the
  three-section set and the needs-confirmation set are now derived from disk, and the
  unresolved-row tests assert the *mechanism* (excluded and reported) rather than a zero
  count. Red builds caused by new data arriving teach people to ignore red.

**Noted, not acted on.** Atomic writes, `verified_ok` being written but never read,
staleness detection, `section_notes` fail-open modes, and config-not-wired-to-code are
real but do not change a number today. Cage/litter structure, dosing-start age and the
multiplicity question are study-design items the PI has ruled out of scope for the image
pipeline; they are recorded for the analysis stage rather than solved here.

## ADR-0017 — GFAP near band reverted to >0–25 µm; the bleed-through reasoning was wrong

**Date:** 2026-08-07 · **Status:** accepted · **Supersedes the 10 µm floor**

The near band was moved to start at 10 µm on the premise that GFAP peaking *inside* the
deposit is an optical bleed-through signature, because genuine astrogliosis would dip at
the amyloid core and peak in a ring outside. **Both halves were wrong.**

*The premise is contradicted by the literature.* Chromogenic DAB GFAP — where optical
bleed is physically impossible — shows the same filled, monotonic, plaque-centred pattern
(Mandybur & Chuirazzi 1990). In vivo two-photon imaging reports explicitly that there are
no plaque-centred concentric rings (Galea et al. 2015, PNAS). Olabarria et al. 2010 show
hypertrophy at plaques and atrophy far from them. The measured profile *is* the published
biology.

*The direction was implausible anyway.* Cy3 (em ~566 nm) bleeding into a FITC band
(~515–530 nm) requires anti-Stokes emission. The plausible direction is FITC → Cy3, which
contaminates the **primary Aβ endpoint** — and that had never been tested.

It has now been tested: far from any deposit (>50 µm), Cy3 in the brightest 5 % of GFAP
tissue against GFAP-dim tissue gives 1.07, 1.10, 1.11 and 1.32 across four animals. Small,
not separable from tissue-property covariance at this resolution, and an order of magnitude
below the peri-plaque GFAP elevation (up to 2.85×) — so it cannot manufacture that signal.

Band restored to `>0–25 µm`. Absolute GFAP percent area is now reported alongside the
enrichment contrast, because normalising to the >50 µm band makes a *global* astrogliosis
change invisible by construction — and that is exactly what rapamycin is reported to do.

Single-stain controls remain the right answer and would settle both directions with one
number. Declining them on antibody quality was a category error: bleed-through is a
fluorophore and filter property, not an antibody property.

## ADR-0018 — 82E1 detects βCTF, on the pathway rapamycin manipulates

**Date:** 2026-08-07 · **Status:** stated limitation, validation proposed

Retracting the earlier pyroglutamate concern was right, but for a better reason than the
one given: it is a human-AD fact wrongly imported into this mouse. Cryo-EM of App NL-G-F
filaments shows the fibril core running from **D1** — the 82E1 epitope (Yang et al. 2023),
and Aβ3pE-42 is ~0.1 % of total Aβ42 in this line against >40 % in human AD (Iwata et al.
2024).

**But the non-differential justification does not hold.** 82E1 recognises the free Asp1
neoepitope, which **C99/βCTF also carries** — it is used for in situ C99 detection
(Lauritzen et al. 2016). Rapamycin acts on the autophagic-lysosomal pathway governing βCTF
catabolism, so the 82E1 signal contains a **treatment-sensitive component that is not
plaque Aβ**. That is a differential bias on its face.

Compounding it, 82E1 stains a 2.5–4.5× larger area than 6E10 in App NL-F and NL-G-F
specifically, with filamentous immunoreactivity beyond plaque cores (Araki et al. 2026) —
so this endpoint is dominated by non-core material, the compartment most sensitive to
autophagy status.

Not a code change. Record as a limitation, and on a subset pair 82E1 with a fibril-selective
dye and/or 6E10: if the treatment effect holds for both, the worry is closed empirically.

---

## ADR-0019 — Corrections after a follow-up review of the review fixes

**Date:** 2026-08-08 · **Status:** accepted

The reviewer who raised the original findings checked whether the fixes actually landed.
Most did. Three things did not, and one of them was a false claim in this very log.

**ADR-0016 claimed a seed guard that did not exist.** It stated "a guard now rejects
date-shaped and short seeds". The patch that was supposed to add it failed silently and
was never verified — `seed=20260807`, `seed=42` and `seed='password'` were all accepted.
The key in use had been drawn correctly, so nothing was at risk, but the next person to
re-code a cohort had no protection and a decision log telling them they did.

**A documented control with no test is how this happens.** The guard is now implemented
*and* tested — date-shaped, short-numeric and wordlist seeds each have a test asserting
refusal, plus one asserting the CLI never passes the `allow_weak_seed` escape hatch. If
the guard disappears again, the suite goes red.

**Staleness was live, not latent.** ADR-0016 listed it among things that "do not change a
number today". That was wrong at the time of writing: the QuPath project the delineator
opens held 22 animals while 29 had pixels. Someone would have drawn regions on a cohort
short by seven animals with nothing to tell them. `./ihc doctor` now fingerprints the
inputs (index files, which have pixels, and the two records that decide condition) and
fails when the derived artefacts no longer match. Project rebuilt: 113 images, 29 animals.

**`doctor` greenlit a failed transfer.** Tube 58's payload folder exists but contains only
error text — no pixel data. `doctor` globbed the folder, counted it, and reported all
checks passed. That is exactly the silent transfer failure this project names as its first
gotcha, and the health check was endorsing it. It now requires a payload folder to contain
at least one `.ets` and fails by name when it does not.

**Figures in the spec had drifted from the data as the cohort filled in.** Corrected in
`CLAUDE_v1.3.md` with the measurement date, rather than quietly: the between-box gap ratio
is 1.33–4.00× (not 1.7–2.4×; 20 of 31 slides fall outside the old range), and tile sparsity
is 1.4–23.6 % (not 5–11 %; two thirds fall outside).

**The one that matters scientifically:** the positive/negative separation margin has
**narrowed from 8.7× to 1.92×** on the normalised index as the cohort grew from 30 to 122
sections. There is still no overlap and the record and pixels agree 122/122, so the check
is sound — but the documented figure overstated the safety by more than fourfold. An
animal landing between 2 and 4 in future is a case for human review, not an automatic call.

---

## ADR-0020 — The staining cross-check was measured on the padded tile grid; re-measured, and GFAP is the discriminator

**Date:** 2026-08-08 · **Status:** accepted · **Supersedes the figures in ADR-0019**

### Context

Payloads for tubes 31 and 58 arrived, completing the cohort at 31 animals with
pixels. Their recorded staining condition had never been checked against pixels,
because until now they had none. Re-running the check meant re-running the ad-hoc
script that produced the separation figure quoted in `CLAUDE_v1.3.md` — and that
script had never been committed, so it had to be rewritten.

Rewriting it surfaced two defects in the original.

**1. It measured on the padded tile grid.** ETS stores tiles as whole 512 px
blocks, so the last row and column are padded out with zeros. At the coarse pyramid
level this check uses, that pad is enormous: a 17,920 px image occupies 280 px of a
single 512 px tile at level 6, so more than two thirds of the assembled plane is
pad. Those zeros sit below every real pixel and pull the background percentile down,
which inflates every signal measured against it. This is the same padded-dimensions
defect `verify.py` was fixed for in Stage 1 (ADR-0002) — reintroduced in a second
place, which is the recurring shape of this bug.

**2. Its tissue mask could be captured by a single dust speck.** Otsu on raw DAPI is
not robust to extreme outliers. Tube 32 s2 peaks at 37,612 counts against a normal
maximum near 670; Otsu put the threshold at 8,889, so 0.03% of the section was
called "tissue" and every statistic then described the artefact. Nineteen sections
were affected.

### Decision

Cropping the pad and clipping the DAPI histogram at its 99.5th percentile before
Otsu fixes both. The check now lives in `src/ihc/qc/stain_check.py` with tests
pinning each defect, and a `scripts/audit_stain_condition.py` that runs it.

Two consequences for the recorded figures:

**The old separation numbers are not reproducible and are withdrawn.** Both 8.7×
and the 1.92× that replaced it were measured on the padded grid. Re-measured over
all 130 sections: GFAP index 1.25–18.98 in positives against 0.27–1.09 in
negatives, AUC 1.000, medians 9.6× apart, extreme-value margin 1.14×.

**Use GFAP, not Aβ.** Measuring both showed the Aβ index does not separate reliably
(AUC 0.966, 7 of 128 misassigned) while GFAP separates perfectly (AUC 1.000, 0 of
128). The reason is biological, not arithmetic: the negative control is
secondary-only, so both markers drop together, but Aβ is sparse and focal and its
abundance genuinely varies between animals by disease burden — so a low Aβ index is
ambiguous in a way a low GFAP index is not. The earlier check happened to use GFAP;
it is now used deliberately, and `abeta_index` is reported as context only.

**Result:** 128 of 128 measurable sections agree with `slides.csv`, across all 31
animals, including 31 and 58 on their first check. Two sections (tube 49 s1 at 0.6%
tissue, tube 52 s4 at 2.4%) have too little tissue in frame to judge and are
reported rather than scored.

### Consequences

- A figure that has now been wrong twice is finally under test. The pattern in both
  ADR-0016 and here is the same: a number produced once in a terminal, quoted in a
  document, and never re-derived. Anything load-bearing gets a module and a test.
- The QuPath project build had the same problem — no committed entrypoint, rebuilt
  by hand — and `doctor` told the user to "rebuild the QuPath project" without
  saying how. It is now `scripts/build_qupath_project.py`.
- `src/ihc/qc/` is no longer empty, which the follow-up review had flagged.

---

## ADR-0021 — A prior QuPath classifier exists; keep it as a comparator, retrain for measurement

**Date:** 2026-08-10 · **Status:** accepted · **Informs D-5 (resolution), D-2 (artefact method), D-3 (positive class)**

### Context

An earlier manual analysis was carried out in QuPath before this pipeline existed, and
its project was handed over: 130 images, one trained Aβ pixel classifier, a run script,
and a partial measurement export covering 68 images. This finally answers a question open
since the first planning round — whether the lab-standard analysis was ever run. It was,
as a working draft.

Reading the classifier and the script settles several things we had been arguing from
theory.

**Resolution: 2.6 µm/px** (`inputResolution` 2.6000734903780036 µm), i.e. downsample 8
from the native 0.325 µm/px. Our candidate list was `[1, 2, 4]` — the prior work sat well
outside it, four times coarser than the coarsest option we were considering.

**Feature set: three features.** One Gaussian filter, σ = 1.0, applied to each of DAPI,
FITC and Cy3. No multiscale, no edge, no texture features.

**Model: `ANN_MLP` with `layer_sizes [3, 3]`** — three inputs, three outputs, no hidden
layer. Combined with the feature set, this is close to a smoothed three-channel intensity
threshold rather than a learned texture rule.

**No intensity normalisation:** the feature preprocessor has `offsets [0,0,0]` and
`scales [1,1,1]`, so it operates on raw counts.

**Denominator: the whole image.** The run script calls `createFullImageAnnotation(true)`,
so "Aβ %" is measured against the entire frame, glass included.

**Trained on a single animal.**

### Decision

**Retrain rather than reuse, and keep the prior classifier as a comparator.** It is a
legitimate lab-standard result to report the new pipeline against, and the painting effort
behind it carries over. It should not become the frozen classifier, for reasons that are
mechanical rather than a judgement on how carefully it was made:

1. **Raw counts, no normalisation, across a cohort whose exposure varies up to 12.6×**
   (ADR-0004) — and 2.73× between Control IP and Rapamycin IP on the Aβ channel
   specifically. An intensity-keyed rule therefore behaves differently per animal, in a
   way partly aligned with treatment on the channel carrying the primary claim.
2. **The whole-image denominator** is dominated by how much empty slide is in frame, and
   section framing varies. It is not comparable between animals. This is what the
   anatomical ROI work replaces.
3. Median exported Aβ area is **13.1%**, with values reaching 100% on small annotations.
   That is far above any plausible immunoreactive area for this material and indicates
   over-detection, consistent with 1–3 above.

**Resolution: this is the strongest evidence we have for D-5, and it points finer, not
coarser.** At 2.6 µm/px an intracellular Aβ deposit spans roughly one pixel and a 20 µm
plaque about eight. Extracellular plaques survive that sampling; intracellular signal
cannot. The difficulty reported with intracellular Aβ is therefore predicted by the
resolution alone and is not a training deficiency — more painting at 2.6 µm/px will not
fix it. Recommend locking **downsample 2 (0.65 µm/px)**, or downsample 1 if intracellular
Aβ is to be a reported class at all. Candidate list unchanged; the PI decision stands.

**Positive class (D-3):** the same evidence bears on it. If intracellular signal cannot be
resolved at the working resolution *and* 82E1 also detects βCTF on the pathway rapamycin
manipulates (ADR-0018), then defining the positive class as extracellular parenchymal
deposits with a size floor removes a measurement problem and a confound together.
Recommended, still the PI's to lock.

### Consequences

- Three PENDING decisions now have evidence attached rather than argument.
- The prior project's measurement export is a working draft, not a result: it mixes
  training annotations with full-section runs and covers 68 of 130 images. It should not
  be cited as a completed analysis.
- The handed-over project sits under `RawData/`, which is meant to hold read-only source
  images only. It is derived work and belongs elsewhere; noted rather than moved, since
  moving another person's files is theirs to agree to.

---

## ADR-0022 — Decisions from the 2026-08-11 bench meeting, and two corrections to ADR-0021

**Date:** 2026-08-11 · **Status:** accepted · **Corrects ADR-0021**

### Context

The PI and the person who did the staining and imaging went through the open
decisions together, with the prior QuPath project open on screen. Several things
were settled; two of them correct claims made in ADR-0021 on the basis of bench
knowledge this project did not have.

### Corrections to ADR-0021

**1. Resolution does not make intracellular Aβ measurable. ADR-0021 implied it
would, and that was wrong.**

ADR-0021 said the reported difficulty with intracellular Aβ was "predicted by the
resolution alone and is not a training deficiency", and recommended downsample 1 if
intracellular Aβ were to be a reported class. Both parts need correcting.

The difficulty was never visibility. Intracellular signal *is* visible at 2.6 µm/px
— a nucleus with signal around it — and it was deliberately painted out during
training rather than missed. The problem is that it sits close to what would
otherwise be called background.

More importantly, **no downsample fixes this**, because the limit is not lateral
sampling. These are 5 µm sections on a widefield scanner, so every pixel integrates
through the whole section thickness. "Inside this cell" is not establishable from
such an image at any pixel size; it needs optical sectioning — confocal at 60× or
better. Finer sampling makes intracellular signal easier to *see* and no easier to
*attribute*.

Consequence: the restrictive positive class is a **signal-based exclusion**, not a
compartment assignment, and the methods text must say so. A companion neuronal
marker (e.g. NeuN) does not rescue it either — it yields a silhouette, not a
compartment, since the marker channel is subject to the same through-thickness
integration.

**2. Resolution is not a neutral knob.** Finer resolution tightens the boundary
drawn around each deposit and therefore shrinks measured area — the diffuse corona
around a dense core is progressively excluded as resolution increases. Resolution
and `abeta.positive_class` are entangled, not independent, and a sweep must report
them together.

### Decisions

**Regions (D-4), now settled.** Hippocampus = Allen **HIP** (CA1–3 + dentate gyrus),
excluding **RHP**. Cortex = isocortex, full thickness, pia to white matter,
white matter excluded. A **third region** was added: the whole section minus glass.

**Resolution (D-5): measure it, do not argue it.** Keep one set of painted
annotations, retrain at each candidate resolution, compare percent area across a
cohort subset. Retraining on existing annotations is a few minutes' work, so this
is cheaper than the argument. `abeta_locked`/`gfap_locked` move from
PENDING_PI_DECISION to PENDING_SWEEP.

**Positive class (D-3): two classifiers, roles declared in advance.** Restrictive
(extracellular) is primary and reported; generous (all immunoreactivity) is a
pre-declared sensitivity analysis. Declaring the roles now is what keeps this a
robustness check rather than two attempts at the same question.

**Delineation method:** hand-drawing now, ABBA deferred. Atlas registration remains
attractive for subregions (CA1 vs CA3) but is not needed for two large regions.

**Order of work:** Aβ classifier first, GFAP after it is working.

**Δ (D-11) withdrawn as posed.** It was presented as an effect size, which is an
inference question and out of scope for this pipeline. The pipeline-relevant
question is narrower — how much human/machine disagreement is tolerable before the
classifier is trusted — and it can be answered with a conventional agreement
threshold without reference to any expected treatment effect. Proposed default:
ICC ≥ 0.8 at animal level.

**GFAP fallback accepted in advance.** Whole-region intensity, not a trained
classifier, is what most of the field uses for GFAP. If the classifier proves
unworkable, falling back to region intensity is a legitimate standard method and
not a failure. Agreeing this now avoids sunk-cost persistence later.

### Findings recorded

**No photobleaching between the original scans and the rescans ~10 days later.**
The rescans of tubes 51 and 60 look brighter, which was raised as a worry since
bleaching should make them dimmer. Most of it is exposure: those two slides had Cy3
corrected from 146 and 398 ms to 1840 ms, i.e. 12.6× and 4.6×. After normalising for
exposure the rescans still carry **1.24–1.46× more signal per millisecond**, so
there is no evidence of degradation in either direction that matters. Caveat: framing
and DAPI exposure also differ, and the 146 ms original sits near the noise floor
where counts per millisecond stop being linear. No further rescans are needed; a
clean bleaching test would require identical exposure.

**Folds are essentially absent, by design.** Slides were selected at the bench to
avoid folds, and a search through the project found none. Guidance to include a
folded section in classifier training is unachievable and has been changed to use
whatever artefacts exist — a stray hair, a torn edge, tissue displaced into a
ventricle.

**Painting technique.** Closed outlines around deposits teach the classifier shape
rather than signal, and make it seek circular objects. Short strokes only. This was
established independently at the bench and matches the written guidance.

### Consequences

- `config.yaml` now records region definitions, both positive classes, and the sweep.
- Both classifier guides and the delineation guide gain the third region and lose
  the instruction to find a fold.
- The whole-section region is hand-drawn for now. Deriving it from the counterstain
  would save roughly two hours and be more reproducible; recorded as available work,
  not done, because the decision in the room was to draw it.
