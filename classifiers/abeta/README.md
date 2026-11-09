# Frozen Aβ classifiers (ADR-0040)

The two classifiers behind the presented Aβ results. Frozen: do not retrain or edit; a change is a new name and a new ADR.

| File | Approach | SHA-256 |
|---|---|---|
| `v2_nonorm_s1.json` | Unnormalised | `c4e837287390bf21db8e58f456651e8fab2ec76085ddf25a501cdcae1330ff60` |
| `v4_global.json` | v4 GLOBAL, base forest | `859edf092a1872cbe9f54749ac5a3c77550488d8be0ec5d7934bbcaaa0757986` |

Both: QuPath 0.7.0 pixel classifier, random trees, Cy3 channel only, 1.3 µm/px working resolution, features Gaussian,
Laplacian of Gaussian, weighted standard deviation and gradient magnitude at sigma 1.3, 2.6 and 5.2 µm, classes
Abeta / Ignore* / Negative. Trained on the wet-lab scientist's 12 annotation files (202 pen strokes; 10 files contribute,
4,782 samples for v4). Trainers: `scripts/retrain.groovy` (unnormalised) and `scripts/retrain_global.groovy` (v4; ADR-0032
fixed the label mapping). The training annotations are not in this repository.

**v4 GLOBAL is applied per section.** `v4_global.json` is never loaded directly: `scripts/make_global_classifiers.py` writes one
copy per section, `gx_<image>.json`, with that section's two brightness anchors from `scripts/make_section_anchors.py`.
`v2_nonorm_s1.json` is applied as is. See `docs/RUNBOOK_abeta.md`.
