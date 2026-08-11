---
output:
  word_document: default
  html_document: default
  pdf_document: default
---
# Training the GFAP classifier — step by step

Do the Aβ classifier first (`classifier_abeta.md`). The mechanics are the same and this
document assumes you have been through them once — including **Where to do this**: same
arrangement, your own project, and at the end you send back only
`classifiers/pixel_classifiers/GFAP_v1.json`, about 5 kB.

You can reuse the same QuPath project you built for Aβ. The six training images for GFAP
should be chosen on **staining strength**, though, which is a different spread from the
plaque-burden spread you used before — see step 5.

**This one is harder, and it is worth knowing why before you start.**

---

## 1. Why GFAP is not like Aβ

Aβ deposits are **objects**: discrete, bright, with an edge. You can point at one and say
where it stops. Getting the boundary slightly wrong changes the answer slightly.

GFAP is a **network**: astrocyte cell bodies with thin processes branching out and
fading continuously into background. There is no edge. It is a gradient.

The consequence matters:

> Because the signal fades out gradually, **where you put the boundary decides the answer.**
> Slightly more generous, and every thin process gets a little wider and a lot more faint
> ones get included — percent area can move substantially for a decision you made without
> noticing you were making it.

This is the single risk in this job. Nearly everything below exists to control it.

**What you are measuring:** GFAP-immunoreactive **area**. Not astrocyte number, not
"is this astrocyte reactive". Do not try to judge activation state — that is a different
measurement this panel cannot support. Your only question per pixel is *"is there GFAP
staining here?"*

---

## 2. Settings

| Setting | Value | Why |
|---|---|---|
| **Classifier** | `Random Trees (RTrees)` | as for Aβ |
| **Resolution** | the entry reading **≈0.65 µm/px** | Matches the Aβ classifier. They must share a grid, because the peri-plaque GFAP analysis measures distance from Aβ deposits — two different grids cannot be combined. |
| **Channels** | **FITC only** | The GFAP classifier must not see the Aβ channel, or "GFAP near plaques" becomes partly circular. |
| **Features** | Gaussian · Gradient magnitude · **Hessian determinant** · **Hessian eigenvalues** · Structure tensor coherence | These detect *elongated ridge-like* structures. Aβ used Laplacian of Gaussian because it finds blobs; processes are not blobs, they are fibres. This is the main difference between the two classifiers. |
| **Scales** | 1, 2, 4 | Smaller than Aβ's 1–8. Processes are thin, so large scales blur them away rather than help. |
| **Output** | Classification | |
| **Region** | Everywhere | |

---

## 3. Classes

| Class | Paint it on |
|---|---|
| `GFAP` | anything showing GFAP staining — cell bodies **and** processes |
| `Negative` | tissue with no GFAP signal |
| `Ignore*` | folds, tears, edges, debris, out of focus |

---

## 4. Before painting: decide where faint stops counting

**Do this once, deliberately, and write it down.** It is the most important five minutes
in this job.

1. Open one mid-range image.
2. Find a **clearly stained process** — obvious, no argument. That is `GFAP`.
3. Find **clearly empty parenchyma**. That is `Negative`.
4. Now find the **faintest thing you would still call a process.** Look at it properly.
   Decide yes or no. **This is your reference point.**
5. Note the image name and roughly where it is — a screenshot is ideal. Keep it open in
   a second window while you train.

Every faint judgement afterwards gets compared to that one reference, not to whatever you
decided ten minutes ago. Without this, the boundary drifts across a session and the
classifier averages your drift.

**The drift is real and it has a direction.** As you look at faint staining for an hour,
you start seeing more of it. Check yourself against the reference every image or two.

---

## 5. Choose 6 images

Same principle as Aβ — span the range:

- two that look **strongly stained**
- two **weakly stained**
- two **middling**
- at least one carrying any artefact you can find — a stray hair, a torn edge, an
  out-of-focus patch. Folds are essentially absent in this material, so do not go
  looking for one

Two useful facts about this cohort:

- **GFAP brightness is fairly consistent.** The green channel was imaged at only two
  exposure settings across all 31 animals, a 1.65× range — unlike the red channel, which
  varies 12.6×. So you are not fighting large brightness differences here.
- **Only 4 sections in the whole cohort remain at a non-standard exposure** — the two
  negative sections of **tube 51** and the two of **tube 60**, the slides that were
  re-imaged. The rescans replaced their positive sections, so **no measurement section is
  affected**. Two consequences: do not pick those as training images, and **do not use
  their negatives for the check in step 8** — being dimmer, they would flatter your
  false-positive rate. Use any of the other 27 animals with negatives.

---

## 6. Painting

Roughly **15–25 strokes per image** — a little more than Aβ, because the boundary is
harder and needs more examples.

**`GFAP`** — deliberately cover the range:
- astrocyte **cell bodies** (bright, easy)
- **thick proximal processes**
- **thin distal processes**, including ones near your faintness reference

If you only paint the bright cell bodies, the classifier will find only cell bodies and
your percent area will be a fraction of the truth.

**`Negative`** — and this is the part people skip:
- parenchyma **immediately next to a stained process**. This is what teaches the edge, and
  the edge is the whole measurement. Paint plenty of these.
- empty neuropil away from any staining

**`Ignore*`** — folds, tears, bright edges, debris.

### Normal anatomy that is genuinely GFAP-bright

Some regions are strongly GFAP-positive in a completely healthy brain:

- the **surface of the brain** (glia limitans, just under the pia)
- **white matter** tracts
- the **hippocampal fissure**
- tissue **around the ventricles**

**Paint these as `GFAP`, because that is what they are.** Do not exclude them and do not
treat them as artefact. They are real staining, they are present in every animal, and the
region outlines already exclude the white matter tracts. Trying to hand-exclude "normal"
GFAP would mean deciding what counts as pathological, by eye, differently in every
section — exactly the inconsistency we are trying to avoid.

---

## 7. Look, correct, repeat

| What you see | What to do |
|---|---|
| Processes visibly present but not detected | paint `GFAP` on them |
| Background haze being called GFAP | paint `Negative` on it — **this is the common one** |
| Detected area spreading well beyond the visible process | boundary too generous; more `Negative` beside processes |
| Folds or edges detected | paint `Ignore*` |

Expect more rounds than the Aβ classifier — five or six is normal. GFAP takes longer to
settle because you are teaching a boundary rather than an object.

---

## 8. The objective check — negative controls

Same tool as for Aβ, and here it matters more.

The project has **52 negative-control sections across 27 animals**, stained without primary
antibody. **There is no real GFAP signal in any of them.**

> Run `GFAP_v1` on two or three negative sections you did not train on. It should report
> **close to zero** GFAP area.

If it reports meaningful area on a secondary-only section, your boundary is too permissive
— it is detecting background, and it will detect that background in every real section
too, inflating every number. Go back and paint `Negative` on whatever it found.

This is the only measurement in the whole job that tells you objectively whether your
faintness threshold is right. **Do not skip it, and do not freeze before it passes.**

Four animals (35, 38, 45, 53) have no negative control, so use one of the other 27.

---

## 9. Save

Name it **`GFAP_v1`**. Then run it on two or three images you did not train on and look at
them properly before calling it done.

---

## 10. Measuring with it

```groovy
def classifier = loadPixelClassifier('GFAP_v1')
def regions = getAnnotationObjects().findAll {
    it.getPathClass() in [getPathClass('Hippocampus'),
                          getPathClass('Isocortex'),
                          getPathClass('Section')]
}
if (regions.isEmpty()) { println 'NO REGIONS: ' + getProjectEntry().getImageName(); return }
selectObjects(regions)
addPixelClassifierMeasurements(classifier, 'GFAP')
println 'Done: ' + getProjectEntry().getImageName()
```

The peri-plaque GFAP analysis — how much GFAP sits within 25 µm of a deposit — is a
separate downstream step and is **not** your job here. It reuses this same frozen
classifier, which is why the resolution has to match Aβ's.

---

## What NOT to do

- **Do not turn on the Cy3 channel** while training this one. If you can see the plaques
  while deciding what counts as GFAP, "GFAP is elevated near plaques" stops being a
  finding and starts being an assumption you painted in.
- **Do not judge whether an astrocyte is reactive.** Only whether there is staining.
- **Do not let the faintness boundary drift.** Check your reference image regularly.
- **Do not freeze before the negative-control check passes.**
- **Do not train across several days.** Your eye changes; the classifier records it.

---

## Quick reference

| | |
|---|---|
| Channel | FITC only |
| Resolution | ≈0.65 µm/px — same as Aβ, deliberately |
| Model | Random Trees |
| Features | Gaussian, Gradient magnitude, Hessian determinant + eigenvalues, Structure tensor coherence |
| Scales | 1, 2, 4 |
| Classes | `GFAP`, `Negative`, `Ignore*` |
| Before painting | pick a faintness reference and keep it open |
| Paint `Negative` | right beside processes — this teaches the edge |
| Normal bright GFAP | brain surface, white matter, fissure, ventricles — paint as `GFAP` |
| Check | ≈zero GFAP area on negative controls |
| Save as | `GFAP_v1` |
