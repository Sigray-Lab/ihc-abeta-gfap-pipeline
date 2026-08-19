# Acquisition-brightness normalisation

Whole-slide fluorescence sections vary in brightness for reasons that are not biological:
staining strength, section thickness, autofluorescence, lamp output, acquisition session.
A pixel classifier trained on absolute intensity inherits that variation, and the result
tracks how bright the photograph was rather than how much marker is present.

This document describes the three normalisation options implemented here and how to run
them. It makes no claim about which is correct for a given dataset — that has to be
measured, and §5 gives the tests.

## Diagnosing the problem

The diagnostic that matters is **secondary-only control sections**: tissue processed
without primary antibody, so any classifier-positive area is false by construction. If
measured area in those sections correlates with acquisition brightness, the classifier is
partly measuring the camera.

    r(log background intensity, log measured area)   in controls   should be ~0

Do not use stained sections for this. There, brightness and burden can correlate for real
reasons, and the test cannot separate them.

**Exposure time is not sufficient to correct this.** In the dataset this pipeline was built
for, every endpoint section shares one exposure and still spans 4.1× in tissue background.
Check before assuming exposure explains the variation; if it does not, dividing by it
corrects nothing.

## Option 1 — none

`scripts/retrain.groovy` with `sigma = 0`. The baseline. Retain it as a sensitivity
analysis even when it is not the primary endpoint: it shows how large the acquisition
effect is.

## Option 2 — global, one offset and one gain per image

`scripts/retrain_global.groovy`, then `scripts/make_global_classifiers.py`.

    I' = (I - glass) / (tissue_p50 - glass)

`glass` is the robust background outside the tissue mask, away from the tissue edge:
detector offset plus stray light plus mounting-medium autofluorescence, an additive term
that does not scale with staining. `tissue_p50` is the median inside the tissue mask.
Glass maps to 0 and typical non-marker tissue to 1 in every image.

**Why the tissue anchor is not circular** when marker burden is a small percentage of
tissue area: the median is then marker-free by construction and barely moves if burden
doubles. This is what separates it from whole-histogram matching, which forces the high
tail — the endpoint itself — to agree between images.

**Two constraints.**

*The gain must be large enough to divide by.* A section whose tissue barely exceeds glass
produces a tiny divisor and amplifies noise. `make_global_classifiers.py` refuses below a
configurable floor.

*The anchor must not differ by experimental group.* If it does, normalising to it removes
the effect rather than revealing it. Test it explicitly before using the output.

Training pools per-image-corrected samples before fitting one model, because QuPath bakes a
single operation chain into a classifier and a per-image constant cannot live in a shared
chain. Application writes one small classifier per image with that image's two constants
substituted — the model is identical in every copy. No pixels are rewritten.

Note that a file containing a single annotation class yields no training data, since there
is no contrast to learn; per-image training therefore drops such files where pooled
training would have used them.

## Option 3 — global, per region

`scripts/make_region_classifiers.py`, with the gain estimated inside the region being
measured rather than across the whole section. The offset still comes from the section's
glass, because there is no glass inside an anatomical region.

Per-region normalisation is often warned against, and correctly so **when regions are
compared to each other** — equalising each region's background flattens exactly the
differences between them. That objection does not apply when the endpoint is a comparison
of experimental groups *within* one region. It is a real constraint on interpretation, not
a reason to avoid the method, and it must be stated wherever the output appears.

## Option 4 — local normalisation

`scripts/retrain.groovy` with `sigma > 0` and a variance sigma, giving QuPath's local mean
and standard-deviation normalisation. Strong brightness robustness, at a cost: it asks
whether a pixel exceeds its own neighbourhood, so a deposit broader than that neighbourhood
partly cancels itself. QuPath's own documentation advises avoiding it by default.

**Never quote an absence of broad or diffuse signal from a locally normalised classifier
without the test in §5.** A method that cannot represent low-spatial-frequency signal will
report its absence, and that report is indistinguishable from a true absence.

## 5. Choosing between them

Measure, do not argue. Every test below can be run without new annotation.

| test | what it detects |
|---|---|
| brightness correlation in antibody-free controls | residual acquisition dependence |
| false-positive floor and upper tail in those controls | specificity, and whether the error is constant or scales |
| separation of stained from control sections | whether the two distributions overlap at all |
| repeat acquisitions of the same physical tissue | reproducibility, holding biology exactly constant |
| classified coverage of the measurement region | whether the classifier abstained rather than measured |
| synthetic broad-signal recovery | whether the method can represent diffuse signal at all |

**Classified coverage deserves emphasis.** A classifier that assigns most of the region to
an ignored class reports a small number for a reason that has nothing to do with biology.
Export `(marker + negative) / region area` on every row and gate on it. A near-total
abstention must be recorded as a failure, never as a zero.

**On repeat acquisitions:** if the two differ by an approximately affine intensity map, a
method that is affine-invariant by construction cannot fail that test. Establish what the
acquisitions actually differ by before treating agreement as evidence.

## Anchor files

Both global options take a CSV of per-image anchors:

| column | meaning |
|---|---|
| `image` | image name, matching the project |
| `glass_p50` | robust intensity outside the tissue mask, away from the edge |
| `t_p50` | median intensity inside the tissue mask (option 2) |
| `region`, `r_p50` | region name and its median (option 3) |

Compute these at a coarse pyramid level; they are summary statistics, not morphometry.
