# Measured Aβ area tracks section brightness — investigation, open

**Started 2026-08-14. Status: finding established, fix under test, not yet resolved.**

Read this before trusting any Aβ percent-area number, and before anyone proposes
"just turn on normalisation" — that phrase names a control that does not work here,
for a reason worth understanding once (§4).

---

## 1. The finding

Running the delivered Cy3-only classifier over all 121 sections at 1.30 µm/px, the
measured Aβ percent area correlates with how bright the section is:

| "Brightness" metric | Stained sections | Secondary-only controls |
|---|---|---|
| Cy3 p99.5 in tissue | +0.74 | +0.79 |
| **Cy3 median in tissue** | **+0.84** | **+0.69** |

(log–log Pearson, n = 61 stained / 44 controls, degenerate reads excluded.)

A 1.9× difference in section brightness produces a 3.6× difference in measured area.

## 2. Why this is not simply biology

The obvious objection is that a section with more amyloid *should* be brighter, so the
correlation is the measurement working. Two things rule that out:

**The controls.** Secondary-only sections received no primary antibody and contain no
amyloid-dependent signal at all. Measured "Aβ" in them still tracks brightness at
r ≈ +0.69 to +0.79. There is nothing there for brightness to be a proxy of.

**The metric.** The first analysis used the 99.5th percentile of Cy3 in tissue, which
in a stained section is *dominated by plaques* — so that correlation was partly
tautological, and the objection was fair. Repeating it with the **median** Cy3 in
tissue, a background-level measure that cannot be plaque signal, made the correlation
in stained sections **stronger** (+0.84), not weaker.

**It is not exposure.** 29 of 31 animals were imaged at identical Cy3 exposure
(1840 ms); the two exceptions were re-imaged. This is section-level variation in
background, staining and optics — the residual that remains even under standardised
staining and acquisition, and that is conventionally normalised away.

## 3. A correction to the first analysis

The first write-up reported rescan sections reading **30×** higher than originals.
**That comparison was confounded**: all 10 rescan sections are *positive*, while the
originals are a mixture of positive and negative. Condition-matched, positives only,
it is **2.55% (originals) vs 33.6% (rescans) ≈ 13×** — still large, but the 30×
figure overstated it and should not be quoted.

## 4. Why the obvious fix is not a fix

QuPath's classifier training offers a **feature-preprocessor** normalisation — the
setting that writes `offsets`/`scales` into the classifier JSON. It is the natural
thing to reach for, and it does nothing here.

It fits **one offset and scale per feature** from the pooled training data and applies
them identically to every image, so it cannot respond to one section being brighter
than another. Worse, the trees are trained *after* it, so their split thresholds are
already expressed in normalised units. **A fixed monotonic per-feature rescale cannot
change a decision tree's output at all.**

The consequence matters: turning it on changes the JSON, so a check of the form "are
`offsets`/`scales` non-trivial?" **passes**, while every measured number stays
identical. That is a false-pass test, and one was very nearly adopted.

## 5. What is actually being tested

`ImageOps.Normalize.localNormalization(sigma, sigmaAfter)` — per-image, spatially
local: each pixel is scaled by variation measured in its own neighbourhood, inside its
own image. That is the operation that can remove section-level brightness.

Four classifiers were retrained headlessly from the same annotations, identical in
every other respect:

| Name | Normalisation |
|---|---|
| `abeta_rt_nonorm` | none — control |
| `abeta_rt_norm25` | local, σ = 25 px |
| `abeta_rt_norm50` | local, σ = 50 px |
| `abeta_rt_norm100` | local, σ = 100 px |

**This is not assumed to work.** Local normalisation with σ near plaque size would
normalise plaques away. The falsification criterion is explicit:

> If positive/negative separation collapses under normalisation, normalisation is the
> wrong tool and the unnormalised classifier stands.

Success requires all three: correlation near zero **on the controls**, the
rescan/original gap closing, and separation preserved.

## 6. Two things found while importing the training data

**The annotation-to-image mapping is not one-to-one.** The training files are named by
original tube (`Image_NNr_s01`); this project uses coded names, and the rescans need a
*section* remap as well: **a rescan's series `01` is physical section 03**, so its
annotations belong on the image ending `_s03`, not `_s01`. For at least one animal the
`_s01` image is a *negative control*, so name-matching does not merely misfile the
annotations, it trains on the wrong staining condition. Build the mapping from the
manifest. The mapping itself is part of the blinding key and lives in the custodian
directory, not here.

**The training set is thinner than it appears.** Of the seven annotation files:
`Image_34_s02` and `Image_49_s01` contain **one `Ignore*` object each** — no Abeta, no
Negative — and `Image_51r_s01` carries **8 unclassified objects** that contribute
nothing. Effectively four images carry real training data, not seven. Worth asking
whether annotations were lost on export; it may matter more than the normalisation.

## 7. Reproducing

Training and evaluation are fully scriptable — no GUI needed:

```bash
QuPath script -p <project.qpproj> -i <any-image> \
  --args "<sigma>,<outName>,<umPerPx>,<annotationDir>,<mappingFile>,<outputDir>" \
  scripts/retrain.groovy
```

The mapping file lives in the custodian directory, since it encodes tube → code pairs.

`PixelClassifierTraining` → `OpenCVClassifiers.createStatModel(RTrees)` →
`PixelClassifiers.createClassifier` → `writeClassifier`. Training takes 1–2 minutes per
classifier; evaluation across the cohort is the slow part, and local normalisation
roughly doubles it because of the padding the Gaussian requires.

## 8. Status

Finding established and reproduced with two independent brightness metrics. The fix is
under test. **No Aβ percent-area number should be reported until this closes**, because
the ranking of animals by burden is currently partly a ranking by section brightness.
