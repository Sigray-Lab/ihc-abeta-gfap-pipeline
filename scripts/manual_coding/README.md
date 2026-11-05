# Manual coding of Aβ deposits (2026-11-05)

**Start:** double-click `PlaqueCoder.command` on the Desktop. Type your name. Pick a box in the list.

**Code mode.** Left click = deposit (red). Right click (or ctrl-click) = uncertain (yellow). Click a dot to remove it.
`z` undo, `s` save, `n`/`p` next/previous box, `1`–`4` channels (Cy3 = the Aβ antibody; DAPI; GFAP; merge),
`+`/`−` zoom, contrast slider for faint deposits. Tick "finished" when a box is done. Saving is automatic when you
change box or quit; files go to `~/ihc_work/manual_coding/coders/<name>/`.

**Practice mode.** "spot ▶" walks through places the classifiers flagged in the box, in teaching order: typical
deposits, deposits with a wide halo, small bright deposits, dim objects only v4 calls amyloid, DAPI-bright
autofluorescence suspects, and (on the control sections) pure false positives. Decide first, then "Reveal what the machine says". The overlays
(v4 outline, 3× gate, DAPI veto) can be switched on at any time, also while coding, but code with them off.

**Compare mode.** Shows the other coder's dots in cyan and the counts agreed / only you / only them (15 µm match).
`compare_coders.py` writes the full table to `~/ihc_work/manual_coding/comparison.csv`.

**Boxes.** Hippocampus: T18_s03 (high load, mid brightness, 3 × 500 µm), J09_s04 (low load, bright, 3 × 300 µm),
C18_s04 (moderate load, dim, 3 × 400 µm), Z23_s02 (very low load, dim, 3 × 300 µm). Controls for practice: Q07_s03
(hippocampus, 3 × 500 µm, almost nothing flagged — a test of over-clicking) and C18_s01 (the worst standard-exposure
control, 3 × 500 µm placed where v4's false positives are densest — what rims, autofluorescent cells and debris look like).
`boxes_index.csv` has the coordinates; `prep_boxes.py` rebuilds the data folder.

Coding rules to agree on before starting (write the answers into this file): what counts as one deposit when
several touch; whether vessel-wall signal counts; whether intracellular puncta count; minimum size.
