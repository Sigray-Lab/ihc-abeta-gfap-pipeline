# Plaque Coder — manual assessment of the Aβ classification (2026-11-05)

**Start:** double-click `PlaqueCoder.command` on the Desktop. Type your name. The program opens in **Rate** mode on the
20-box rating set.

## Rate mode (the main job today)
Each box shows one candidate delineation at a time, by name, always in the same order:
v4 GLOBAL (as is) → v4 + 1.5× brightness gate → v4 + 2× gate → v4 + 3× gate → v6 (v4 minus DAPI veto) → v5
(retrained with control negatives). Look at the magenta outline against the raw Cy3 image and press

- `g` too overinclusive  (outline takes in tissue that is not deposit)
- `a` about right
- `t` too underinclusive (outline misses deposit or cuts it down to the core)
- `u` unsure

Optional comment in the box before pressing the key. The program advances to the next method, then to the next box.
`o` switches the outline off and on to see the raw image; `←`/`→` move between methods; `1`–`4` channels; `+`/`−` zoom;
contrast slider for faint deposits. Everything is saved immediately to `~/ihc_work/manual_coding/coders/<name>/ratings.csv`,
so you can stop and resume at any time. Rate alone or together; if together, use one name.

**The set:** 10 stained sections × (1 hippocampus box + 1 isocortex box) = 20 boxes, chosen to span brightness,
load and scanning session: T18_s03, J09_s04, C18_s04, Z23_s02, W23_s03, F15_s04, V07_s02, Y06_s01 (July originals,
2 per brightness×load cell) and R22b_s03, A09b_s03 (August rescans). Box sides 250–500 µm depending on how much room
the region outline has. `python3 tabulate_ratings.py` prints the verdict table per method, per coder and per
stained/control, with the comments.

## Practice mode (teaching)
Boxes marked "(practice)" in the list: the other hippocampal boxes of the four first sections plus two antibody-free
controls (Q07_s03, almost nothing flagged; C18_s01, dense parenchymal false positives). "spot ▶" steps through spots the
classifiers flagged, in teaching order (typical deposit, deposit with wide halo, small bright deposit, dim and v4-only,
DAPI-bright autofluorescence suspect, control false positive); decide, then Reveal.

## Code mode (dot counting, optional)
Left click = deposit, right click = uncertain, click a dot to remove, `z` undo, `s` save. Per coder and box; Compare
mode shows the other coder's dots and the agreement counts; `compare_coders.py` writes the table.

## Files
`boxes_index.csv` (coordinates, set = rating / practice), `prep_boxes.py` (rebuilds `data/`), `overviews/` (where each box
sits on its section), `v4_objects/`, `v5_objects/` (classifier exports used for the outlines).
