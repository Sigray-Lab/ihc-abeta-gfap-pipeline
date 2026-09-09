# Archived run scripts

These are the QuPath scripts that actually produced the reported measurements. Until
2026-08-19 they existed **only in `/tmp`** and would not have survived a reboot, so the
measurement stage of the pipeline was not reproducible from this repository. Copied here
verbatim, unmodified, on 2026-08-19.

`eval_v5.groovy` produced the current section-level table. Note that table also carries
`haze` and `core_um2` columns which `report_in_tissue.py` does not write, so at least one
further untracked stage sits between pixels and the reported numbers. Reconcile before
publishing anything from it.

`eval_v4.groovy` produced the full-cohort run. It exports `Abeta/(Abeta+Negative)` under
the column name `abeta_pct_in_tissue` -- the mislabelling ADR-0029 exists to prevent.
`report_in_tissue.py` renames it defensively, which works, but the name is a trap.

The maintained, version-controlled evaluation script is
`pipeline/qupath/scripts/eval_in_tissue.groovy`. Prefer it for all new work.
