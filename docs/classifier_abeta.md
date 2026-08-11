# Training the Aβ classifier — step by step

You are teaching QuPath one rule: **which pixels are Aβ deposit and which are not.** Once it
is trained you freeze it, and the identical rule runs on all 121 images without further
judgement. That is the whole point — the machine is not smarter than you, it is just
perfectly consistent.

Work in the blinded project (`~/ihc_work/qupath`). QuPath **0.7.0**.

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
| **Resolution** | the entry reading **≈0.65 µm/px** | The old one ran at 2.6 µm/px, where an intracellular deposit is one pixel. This is 4× finer. |
| **Channels** | **Cy3 only** | The Aβ classifier must not see the GFAP channel, or the two endpoints stop being independent measurements. |
| **Features** | Gaussian · **Laplacian of Gaussian** · Gradient magnitude · Weighted deviation | LoG is a blob detector — it is what actually finds plaques. The old classifier had one Gaussian and nothing else. |
| **Scales** | 1, 2, 4, 8 | ≈0.65–5 µm. A "scale" is how far around each pixel the classifier looks; four of them lets one rule handle both a 5 µm speck and a 50 µm plaque. |
| **Output** | Classification | |
| **Region** | Everywhere | |

> **Do not start until the resolution is confirmed.** If it is later locked at 0.325 µm/px
> instead, the classifier must be retrained from scratch — it cannot be adjusted.

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
- at least one with a **fold, tear or bad edge** — this is what makes `Ignore*` learnable

Write the six names down. They get recorded with the frozen classifier, and they should
not later be used as your check images.

Six is enough. Every image you train on is one you can no longer use as an independent
check.

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
    it.getPathClass() in [getPathClass('Hippocampus'), getPathClass('Isocortex')]
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
