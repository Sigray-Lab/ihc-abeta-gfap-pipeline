# IHC analysis pipeline: Aβ and GFAP in App NL-G-F mouse brain

Quantifies amyloid-beta and GFAP immunoreactive percent area, and local GFAP enrichment
around Aβ deposits, from Olympus VS200 whole-slide immunofluorescence.

A human trains pixel classifiers in QuPath. This pipeline handles everything around that:
reading the vendor format, verifying the raw containers, building an auditable manifest,
blinding the cohort, applying frozen classifiers at scale, and measuring.

**Current status of pipeline/project: still in pre-measurement / pre-analysis.** Ingest, verification, manifest and blinding are built and
tested. Classifier training and measurement are not.

---

## Install and run

Requires Python 3.11 with numpy, scipy, scikit-image, imagecodecs, pandas and pyyaml
(`env/environment.yml`), plus QuPath 0.7.0 for the classifier and delineation steps.

```bash
./ihc doctor          # environment, paths, and a check for cloud-dehydrated files
./ihc check-config    # parse the YAML, validate slides.csv, list undecided parameters
./ihc verify          # structural verification of the raw VSI/ETS containers
./ihc meta            # metadata sweep over the .vsi index files
./ihc manifest        # build the section manifest — one row per tissue section
./ihc blind           # split into a private and a blinded manifest
```

Every path resolves through `config/paths.yaml`, so the storage layout changes in one
edit. Raw data is read-only; everything else regenerates from raw + code + config, with
one exception noted under Blinding.

---

## Data

One slide per animal, three or four coronal sections per slide arranged in two PAP-pen
wells. Three channels: DAPI, FITC (GFAP), Cy3 (Aβ), 16-bit, 0.325 µm/px, 280–410
megapixels per section.

A `.vsi` file is an index only — pixels live in a sibling `_Image_NN_` folder. Copying
the `.vsi` alone yields ~1.7 MB against ~1.4 GB of real data, and the failure is silent.
That is why verification is a gate rather than a nicety.

### Format notes

- **The tile grid is sparse.** 4–18 % of tile positions are never acquired, because the
  scanner uses a sample mask. Those regions are *missing support* — not background, and
  not tissue. They must be excluded from any area denominator.
- **Tiles are lossy JPEG2000** (9/7 irreversible wavelet, 2 quality layers). Irreversible,
  and a stated limitation. Compression artefacts are signal-dependent, so validation is
  stratified by exposure and burden rather than assuming the codec is neutral.
- **A missing `.ets` stack fails silently.** Bio-Formats reports every named series and
  exits 0, because names and dimensions come from the index alone. The verifier asserts
  a stack inventory instead of trusting the series list.
- **Bio-Formats returns almost nothing from an index-only `.vsi`** (14 KB of XML against
  745 KB with the payload). The pipeline therefore includes its own dependency-free tag
  reader, so metadata for the whole cohort is available before the bulk data transfers.

---

## Caveat / NB

Each slide holds two PAP-pen wells. **The two sections in one well share a staining
condition** — one well received primary antibody, the other DAPI and secondary only. Which
well is positive varies per slide, and getting it wrong means quantifying negative
controls as data, with every downstream number in a plausible range.

Three things make that harder than it looks, all verified in this cohort:

1. **The section number is acquisition order, not slide position.** Two slides run out of
   order. Well membership is derived from **stage coordinates** — sort by stage X, split
   at the largest gap (the between-well gap is 11–23 mm; within-well spacing is 6–8 mm).
2. **Layouts are not uniform.** 2+2, 2+1, 1+2 and 4+0 all occur. Four slides have both
   wells stained and therefore no negative control at all.
3. **Condition comes only from the wet-lab record** (`config/slides.csv`), never from
   pixels. A pixel-based check runs alongside as an *auditor* and has agreed on every
   checkable section — but it can never write the condition column. Deriving condition
   from images would silently reclassify a failed stain as a control, and would make the
   negative-control check circular, since the controls would have been selected for the
   property being tested.

Where the record contradicts itself, the section is marked `unresolved`, excluded from
the analysis manifest, and reported. It is never guessed.

---

## Location 

The repository sits inside the project folder, beside the data it describes:

```
IHC_analysis_pipeline/
  pipeline/          this repository
  RawData/           read-only source images
  BlindingKey/       the code key — gitignored, never committed
  CLAUDE_v*.md       specification
  EXECUTION_PLAN_v*.md
```

Derived artefacts are written outside the project folder (see `config/paths.yaml`),
because they are high-churn and would otherwise be re-synced continuously.

## Repo layout

```
ihc                       single entry command
config/                   paths, analysis parameters, channel map, wet-lab records
  slides.csv              per-slide positive-well assignment — the authority on condition
  section_notes.csv       per-section facts from the bench (artefacts, scan preference)
src/ihc/ingest/           vsi_meta, verify, manifest, blinding, qupath_export
src/ihc/util/             config loading
qupath/scripts/           Groovy for building the blinded project
env/                      environment and pinned tool versions
docs/decisions.md         ADR log — every non-obvious choice and why
tests/                    619 tests
```

---

## Manifest / Meta-data master file

One row per tissue section, joining scanner metadata to the wet-lab record. Carries
condition, well, stage position, exposure, calibration, dimensions, checksums, bench
notes and provenance.

**Filter on `use_for_measurement`.** It is the single correct filter: condition resolved,
not bench-excluded, and the preferred scan for that physical section. Preference is
resolved *per section*, not per animal, because a second imaging session may be a
replacement for some sections and a supplement for others. An invariant fails the build
if any animal would vanish, lose its positive sections, or be double-counted.

---

## Blinding procedure

Coded IDs from a seeded random permutation. Treatment groups run in contiguous blocks of
tube ID, so any order-preserving scheme — sequential, hash, arithmetic — reproduces the
group structure exactly and blinds nobody.

- The key lives outside the repository and is gitignored. Nothing in the repo maps a code
  to a group.
- The key is **append-only**: an existing key is reused and the run refuses if any issued
  code would change. Without this, adding one late animal reshuffles nearly the whole
  cohort and orphans every measurement made under the old codes.
- Blinded artefacts exclude the slide-label image (which carries the animal ID as printed
  text *and* a barcode), internal series names, file paths, acquisition timestamps and raw
  exposure values — each verified as a live leak vector in this data, not a hypothetical.
- Blinded QuPath projects show **DAPI only**, with the marker channels switched off, so
  anatomical boundaries cannot be drawn toward or away from visible signal.

Blinding here is procedural, not cryptographic. Roughly four animals in thirty remain
identifiable to someone holding the staining record, because a three-section or
double-stained slide is structurally distinctive. That is documented rather than
engineered around; see `docs/decisions.md`.

The key is the one artefact that cannot be regenerated from raw data plus code.

---

## Tests

```bash
python3 -m pytest tests/ -q
```

The load-bearing tests are the well-assignment regressions (including both out-of-order
slides) and a synthetic-corruption battery asserting that verification *fails* on a
removed stack, a cloud conflicted-copy file, a truncated tile table, a wrong-companion
prefix trap, and a label image substituted into a tissue stack.

Tests that depend on the cohort derive their facts from disk rather than freezing them.
Payload folders arrive in batches and open questions get answered, so a frozen expectation
turns legitimate progress into a red build — which teaches people to ignore red.

---

## Design notes

- **Percent-area denominator:** anatomical ROI ∩ acquired support ∩ tissue − artefact.
  Every component is exported separately so an alternative denominator can be recomputed
  without re-running segmentation.
- **The numerator is measured inside that same region, not on the image frame.** This
  looks too obvious to state, which is exactly why it went wrong: the evaluators measured
  the classifier over the whole rectangle and divided by tissue area, and it survived two
  internal review rounds before an external reviewer caught it. On one control section
  96 % of the locally-normalised classifier's "amyloid" was glass. Measurement region is
  supplied by `scripts/make_tissue_rois.py` and enforced by `eval_v4.groovy`, which fails
  rather than returning zero. **ADR-0029.**
- **One frozen classifier per marker**, applied identically everywhere, released only
  after validation stratified by region and burden.
- **Classifier input channels are isolated:** the Aβ classifier may not use the GFAP
  channel and vice versa, or a treatment effect on one marker leaks into the other and
  the peri-plaque measure becomes partly circular.
- **Validation partitions by animal, not by image.** Sections from one animal share
  staining, handling and anatomy.
- **The endpoint is Aβ-immunoreactive area, not plaque load.** A pixel classifier cannot
  separate parenchymal deposits from vascular or intracellular signal without object
  rules, so the honest name is the one that describes the measurement.

Rationale for every non-obvious choice, including several that were reversed after review,
is in `docs/decisions.md`.
