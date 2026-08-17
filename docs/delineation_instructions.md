---
output:
  word_document: default
  html_document: default
---
# Drawing the brain regions — a step-by-step guide

You do not need to have used QuPath before, and you do not need to know
anything about the code. This is the whole job.

**What you are doing:** on each picture of a brain section, draw three outlines — around
the hippocampus, around the cortex, and around the whole piece of tissue. That is it.
The computer measures what is inside your outlines afterwards.

**How long it takes:** about 5–15 minutes per section once you are used to it. There are
**121 images**. Do them in batches; there is no rush and no deadline inside a single sitting.

---

## 1. Open the project

**The easy way:** double-click **`Open QuPath project.command`** in the project folder.
QuPath opens with the project already loaded. The first launch takes 10–20 seconds and
leaves a small black terminal window behind, which you can close.

**If that does not work:** open QuPath yourself, then **File → Project… → Open Project**,
and pick `project.qpproj` inside `ihc_work/qupath/qupath/`.

Either way, a list of images appears in the panel on the left. Each is named something
like **`K07_s01`** — a code, not a mouse number. Double-click the first to open it.

If QuPath asks you to "set image type", something has gone wrong — close it without
saving and tell the person who set this up. It should already be set.

---

## 2. What you are looking at

The picture will be **blue-grey, showing cell nuclei only**. It will look like a plain
outline of the brain with brighter bands where cells are packed tightly together — the
dentate gyrus and the CA layers of the hippocampus stand out clearly, which is exactly
what you need to find the boundaries.

**The green and red channels are switched off on purpose.** You are not missing
anything and nothing is broken.

The reason matters, so here it is plainly. If you can see the plaques and the GFAP while
you are drawing, your hand drifts. Not deliberately — nobody does it deliberately — but
an edge gets nudged a little to take in an interesting patch, or a little away from an
empty one, and it happens more in some sections than others. Since the number we
finally report is "how much signal, divided by how much area you outlined", moving the
edge changes the answer. Drawing on nuclei only makes that impossible. The outline ends
up where the anatomy is, not where the signal is.

The image names are codes. `K07_s01` means "animal K07, section 01" — the code is not the
tube number, so just take the names as given.

---

## 3. Draw the three regions

You do this three times per image: hippocampus, cortex, and the whole section.

**Set up the classes once, the first time only:**

1. In the **Annotations** tab on the left, find the class list.
2. Click the **`...`** (or the cog) → **Add/Remove…** → **Add class**.
3. Add exactly these three, spelled exactly like this:
   - `Hippocampus`
   - `Isocortex`
   - `Section`

Spelling matters — the software matches on these names. Capitals as shown, no spaces.

**Then, for each image:**

1. Choose the **Polygon** tool from the toolbar (the icon that looks like a shape made
   of straight edges). The **Brush** tool also works if you prefer painting; either is
   fine.
2. Click your way around the hippocampus, one click per corner. Double-click to close
   the shape. You do not need hundreds of points — follow the boundary, and use more
   points where it curves tightly.
3. With the new shape still selected, click **`Hippocampus`** in the class list, then
   the **Set class** button. The outline changes colour.
4. Do the same for the cortex (**`Isocortex`**) and for the whole piece of tissue
   (**`Section`**).

**How accurate does it need to be?** Follow the anatomical boundary as you see it. A
point every so often along a smooth edge is fine; do not agonise over single pixels. Be
**consistent** rather than perfect — the same judgement applied the same way to every
section is worth far more than heroic precision on a few.

**Where the boundaries go:**

- **Hippocampus** — **CA1, CA2, CA3 and the dentate gyrus.** In Allen atlas terms this is
  **HIP**, and it deliberately does **not** include the subiculum or entorhinal cortex
  (Allen calls those **RHP**, the retrohippocampal region).

  The cell layers give you the edge: follow the outside of the dark packed bands. The one
  boundary worth knowing is where **CA1 becomes subiculum** — the pyramidal cell layer
  suddenly widens and loses its sharp edge. Stop there. That transition is the single
  most common place for two people to disagree, so it is the one place worth slowing down.

  Do not include the thalamus below, the cortex above, or the white matter tracts
  (fimbria, alveus) wrapping the outside.

  **One hippocampus per section, not two.** If more than one hippocampal profile is
  visible, outline only the one that is easier to define by eye, and leave the other. Do
  not draw both. Quick and consistent beats complete here — the full extent comes from the
  atlas later, and because the atlas labels every pixel we can trim it back to match what
  you drew whenever the two need comparing.

- **Isocortex** — the full thickness of cortex at this level, from the outer surface of
  the brain (the pia) down to the white matter, but **not including** the white matter
  itself. Stop where the cortical layers stop.
- **Section** — the whole piece of tissue, everything inside the outer edge of the
  brain. This one does not need care: follow the tissue boundary roughly and keep the
  glass out. It is the denominator for a whole-slide burden number, so only the
  tissue/glass border matters, not anatomy. A quick Brush pass is fine.
- **One hemisphere only.** These sections are single hemisphere. If a second, partial
  piece of tissue is present, leave it alone.
- If part of a region runs off the edge of the picture, outline the part that is there.
  Do not guess where the missing part would have been.

**On coronal level.** The sections are all cut at roughly hippocampal level, but not
identically. Do **not** spend time working out an exact position for each one — that is
slow and we do not need it. Just glance at whether the hippocampus looks well developed
with both blades of the dentate gyrus visible. If a section is obviously much further
forward or back than the others, write that in the description box and carry on. We check
level-matching afterwards from the outlined areas, which is quicker and more consistent
than judging it by eye.

**If you are unsure where a boundary is:** draw your best judgement and add a note (see
section 5). Do not skip the section, and do not open the other channels to help you
decide.

---

## 4. Save

Press **Ctrl+S** (Windows) or **Cmd+S** (Mac) before you move to the next image. Or use
**File → Save**.

QuPath does not save automatically. If you close an image without saving, the outlines
are gone.

When you finish a session, that is all — there is nothing to export and nothing to send.
the person who set this up reads the project directly.

---

## 5. If a section looks damaged

Damage is **rare in this material** — the slides were chosen at the bench to avoid folds,
and a search through the project found essentially none. But tears, holes, bubbles, stray
hairs and out-of-focus patches do turn up. **Do not throw the section away and do not try
to fix it.**

Do this instead:

- **Small damage inside a region** (a fold, a tear, a bubble): draw the region outline
  normally, right around the outside as usual. Then draw a second shape *around the
  damaged part* and give it the class **`Artefact`** (add that class the same way you
  added the other two). The software subtracts it.
- **A region is more than roughly half destroyed**: outline what is genuinely there,
  and add a note.
- **The whole section is unusable**: do not draw anything. Add a note.

**To add a note:** in the image list on the left, right-click the image → **Edit
description** (or use the Description box), and write in plain words what is wrong —
"big fold across the hippocampus", "out of focus lower half", "section torn, cortex
missing on the left". Then save.

Those notes are genuinely used. A section you flagged and drew is far more useful than a
section quietly skipped, because we can decide later what to do with it — but only if we
know.

If something looks wrong in a way this section does not cover, write the note and ask.
Asking is always the right call. Nothing here is urgent enough to guess about.

---

## 6. What NOT to do

**Do not turn the other channels on.** If you can see the plaques while drawing, the edge
drifts towards them — not deliberately, but it does, and it cannot be detected afterwards.
Since the final number is signal divided by the area you outlined, that changes the answer.
If you turn them on by accident, close the image **without saving**, reopen, carry on.

**Do not rename the images.** The names are the only link between an outline and the
mouse it came from. A renamed image is an outline that belongs to nobody.

**Do not add other classes**, and do not delete or move images in and out of the project.

**Do not redraw a region days later** unless something is clearly wrong. A fresher eye
quietly makes that section different from the rest, and consistency is what we are after.

---

## 7. Quick reference

| | |
|---|---|
| Open | File → Project… → Open Project → `project.qpproj` |
| Draw | Polygon tool, click round the edge, double-click to finish |
| Label | Select the shape → click the class → **Set class** |
| Classes | `Hippocampus`, `Isocortex`, `Section`, plus `Artefact` for damage |
| Save | Ctrl+S / Cmd+S, **every image, before moving on** |
| Damaged | Draw what is there, mark damage as `Artefact`, write a note |
| Hippocampus | CA1–3 + dentate gyrus. **Not** subiculum or entorhinal |
| Never | turn on the green/red channels · rename images |

Anything unclear, anything that looks odd, anything not covered here: ask the person who set this up before
carrying on. There is no such thing as a silly question about this, and a wrong guess
here is expensive to find later.
