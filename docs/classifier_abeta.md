---
output:
  word_document: default
  html_document: default
  pdf_document: default
---
# Training the Aβ classifier — step by step

You are teaching QuPath one rule: **which pixels are Aβ deposit and which are not.** Once it
is trained you freeze it, and the identical rule runs on all 121 images without further
judgement. That is the whole point — the machine is not smarter than you, it is just
perfectly consistent.

## Where to do this

**You do not need a copy of the main analysis project.** Train on your own raw data, then
send back one small file.

A trained classifier is about **5 kB of JSON**. It contains the rule, not the images — it
refers to channels by name (`DAPI`, `FITC`, `Cy3`) and to resolution in microns per pixel.
So a classifier trained in your project applies unchanged in ours.

Set up, once:

1. QuPath **0.7.0** (the pinned version — please check yours matches).
2. **File → Project… → Create project**, into an empty folder.
3. Drag in **6 `.vsi` files** — see step 3 for how to choose them. When Bio-Formats asks
   which series, pick the **20× tissue** series, not Label, Overview or Macro.
4. Set image type to **Fluorescence** when prompted.

When you are finished, send back **only** the classifier file:
`<your project>/classifiers/pixel_classifiers/Abeta_v1.json`. Not the project, not the
images.

One small thing: choose your six by **how they look**, not by tube number. The measurement
runs on a coded copy at our end, and picking training images by animal identity would put
that choice into the frozen rule.

Budget about 2 hours. It is not a long job, but it is the job everything downstream
inherits, so do it in one sitting rather than across several days — your eye shifts
between days and the classifier records that shift.

---

## 1. Settings

Use exactly these. Each one differs deliberately from the earlier classifier, and the
reasons are in `docs/decisions.md` ADR-0021.

| Setting | Value | Why |
|---|---|---|
| **Classifier** | `Random Trees (RTrees)` | The old one had no hidden layer, so it was effectively a brightness threshold. Random Trees can use texture and shape. |
| **Resolution** | start at **≈0.65 µm/px**, then sweep — see below | 4× finer than the old 2.6 µm/px. But this is the one setting to *measure* rather than assume. |
| **Channels** | **Cy3 only** | The Aβ classifier must not see the GFAP channel, or the two endpoints stop being independent measurements. |
| **Features** | Gaussian · **Laplacian of Gaussian** · Gradient magnitude · Weighted deviation | LoG is a blob detector — it is what actually finds plaques. The old classifier had one Gaussian and nothing else. |
| **Scales** | 1, 2, 4, 8 | ≈0.65–5 µm. A "scale" is how far around each pixel the classifier looks; four of them lets one rule handle both a 5 µm speck and a 50 µm plaque. |
| **Output** | Classification | |
| **Region** | Everywhere | |

### The resolution sweep

Rather than pick a resolution by argument, **measure it.** Once you have painted your
annotations, retrain the same annotations at each resolution — `0.325`, `0.65`, `1.3`,
`2.6` — and compare the percent-area results across a handful of animals. Retraining on
existing annotations takes minutes, so this is cheaper than the debate.

Two things to expect, because they are not artefacts:

- **Finer resolution gives smaller areas.** The diffuse corona around a dense core gets
  progressively excluded as the boundary tightens. So resolution and "what counts as a
  plaque" are the same question asked twice — report them together.
- **Intracellular Aβ will not resolve at any of these settings.** Going finer makes it
  easier to *see* and no easier to *attribute*. These are 5 µm sections on a widefield
  scanner, so every pixel integrates through the whole thickness; establishing that signal
  is *inside a cell* needs optical sectioning (confocal, 60×). Excluding it is a
  signal-based judgement, not a compartment assignment, and the methods will say so.

Whichever resolution is chosen, it is then **frozen** — a later change means retraining
from scratch, not adjusting.

---

## 2. Set up three classes

In the **Annotations** tab, add exactly these, spelled this way:

| Class | Paint it on |
|---|---|
| `Abeta` | the deposits |
| `Negative` | tissue with no deposit — parenchyma, cell layers, neuropil |
| `Ignore*` | folds, tears, edges, debris, out-of-focus patches |

**`Ignore*` is not optional.** Without it a bright fold has nowhere to go, and the
classifier learns it as plaque. Folds are bright in Cy3 too, so they are learnable from
this channel alone.

---

## 3. Choose 6 images to train on

Pick them to span what the classifier will meet:

- two that look **heavily loaded** with deposits
- two **sparse**
- two in the **middle**
- at least one carrying **whatever artefact you can find**

On that last point: the slides were chosen at the bench to avoid folds, and a search
through the project found essentially none — which is a compliment to the sectioning, not
a problem. So do not hunt for a fold. Use what is actually there: a stray hair, a torn
edge, tissue displaced into a ventricle, an out-of-focus patch. `Ignore*` needs *an*
example of "not tissue signal", not specifically a fold.

Write the six names down. They get recorded with the frozen classifier, and they should
not later be used as your check images.

Six is enough. Every image you train on is one you can no longer use as an independent
check.

---

### A note on brightness, since it comes up

Camera exposure on the red channel varies up to **12.6×** across the raw files, which would
matter a great deal: counts scale with exposure, so a fixed rule finds far more "plaque" in
a longer-exposed image than a shorter one, for no biological reason.

**It does not affect this work, and here is why.** 29 of the 31 animals were imaged at the
same Cy3 exposure (1840 ms). The two that were not — tubes 51 and 60 — had their positive
sections re-imaged, and the rescans are what the analysis uses. **Every positive section in
the cohort now sits at the same exposure.**

The only leftovers are the **negative control sections of tubes 51 and 60**, which are still
at the short exposure. So:

- **Do not use tube 51 or 60 sections as training images.**
- **Do not use their negatives for the check in step 6** — they are dimmer than they should
  be, which would make your false-positive rate look better than it is. There are 27 other
  animals with negatives.

Brightness differences you see *within* the standard 1840 ms group are real — section
thickness, staining, tissue quality — and painting across the range is exactly how the
classifier learns to cope with them.

---

## 4. Paint — small strokes, not filled regions

Roughly **10–20 short brush strokes per image**. You are giving examples, not segmenting.

- **`Abeta`** — several deposits across the full range: big and obvious, small, and
  importantly **the dim ones you would still call a plaque**. If you only paint the
  brightest, the classifier will only find the brightest.
- **`Negative`** — parenchyma right next to a deposit (this is what teaches the edge),
  plus cell layers and empty tissue.
- **`Ignore*`** — every fold, tear, bright edge or debris speck you can find.

**Turn Live prediction off while painting** — it recalculates constantly and will crawl.
Turn it on when you want to look.

---

## 5. Look, correct, repeat

Switch Live prediction on and hunt for the two failure modes:

| What you see | What to do |
|---|---|
| Bright non-plaque called `Abeta` | paint `Negative` or `Ignore*` on it |
| Dim real deposits missed | paint `Abeta` on them |

Add a few strokes, look again. **Three or four rounds is normal.** Stop when you are
correcting cosmetics rather than mistakes.

---

## 5b. There is already a baseline to compare against

A plain intensity thresholder has been built and calibrated automatically, and it is
sitting in the project as **`Abeta_threshold_900`**. Load it the same way you would load
any classifier.

It says only "Cy3 above 900 counts is Abeta". No texture, no shape, no size filter. The
threshold was chosen objectively: it is the lowest value that holds the worst
negative-control section below 0.05 % area, so it keeps as much real signal as possible
while still rejecting background.

**Use it as a sanity check on your trained classifier, not as a competitor.** If your
classifier and this one disagree wildly on the same image, one of them is wrong and it is
worth finding out which before freezing.

For reference, on the calibration set it gives a median of **0.89 %** Aβ area inside
tissue, ranging 0.11–5.65 % across animals.

**And here is exactly why the trained classifier is still needed.** At that threshold one
negative control still reports **0.23 %** "Aβ" — a bright speck or edge that a threshold
has no way to reject, and which lands *above* the weakest genuine positive at 0.11 %. A
threshold cannot tell a plaque from any other bright thing. Learning to reject artefacts
is precisely what `Ignore*` and the trained features add, and it is the part you cannot
automate.

---

## 6. The objective check — use the negative controls

This is the part that turns opinion into evidence, and it is better than eyeballing.

The project contains **52 negative-control sections** — sections stained without the
primary antibody. There is no real Aβ signal in them at all. So:

> **A correctly trained classifier should find close to zero Aβ area on a negative
> control.**

Open two or three negative sections you did not train on and run the classifier. If it
reports appreciable Aβ area there, your rule is too permissive — go back and paint more
`Negative` on whatever it is picking up.

This is a direct measurement of your false-positive rate. Do it before you freeze.

---

## 7. Save

Name it with a version: **`Abeta_v1`**. If you retrain, make it `Abeta_v2` — never edit a
classifier in place once anything has been measured with it.

Then open an image you did **not** train on and look at it. If it holds up there, it is
ready.

---

## 8. Measuring with it

Measurement happens **inside the drawn regions**, not the whole image. The earlier
analysis measured against the entire frame including glass, which is why it reported a
median of 13% Aβ.

```groovy
def classifier = loadPixelClassifier('Abeta_v1')
def regions = getAnnotationObjects().findAll {
    it.getPathClass() in [getPathClass('Hippocampus'),
                          getPathClass('Isocortex'),
                          getPathClass('Section')]
}
if (regions.isEmpty()) { println 'NO REGIONS: ' + getProjectEntry().getImageName(); return }
selectObjects(regions)
addPixelClassifierMeasurements(classifier, 'Abeta')
println 'Done: ' + getProjectEntry().getImageName()
```

This needs the region outlines to exist first, so it runs after delineation.

---

## What NOT to do

- **Do not turn on the FITC channel** while training this one. Cy3 only.
- **Do not keep tuning after you have seen group differences.** The rule is frozen before
  measurement, not adjusted until the answer looks right.
- **Do not train on all 121 images.** You need untrained images left over to check against.
- **Do not edit a classifier in place** after measuring with it — new version, new name.

---

## Quick reference

| | |
|---|---|
| Channel | Cy3 only |
| Resolution | ≈0.65 µm/px |
| Model | Random Trees |
| Features | Gaussian, LoG, Gradient magnitude, Weighted deviation |
| Scales | 1, 2, 4, 8 |
| Classes | `Abeta`, `Negative`, `Ignore*` |
| Training images | 6, spanning burden and quality |
| Check | ≈zero Aβ area on negative controls |
| Save as | `Abeta_v1` |
