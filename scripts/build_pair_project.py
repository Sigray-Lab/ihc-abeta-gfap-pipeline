#!/usr/bin/env python3
"""Build a QuPath project holding BOTH acquisitions of every re-scanned section.

Why this exists
---------------
Six slides were put back on the microscope, and for nine physical sections the
result is two independent acquisitions of **the same piece of tissue**. The wet-lab
note on two of those slides says why in as many words: *"REDO different exposure"*.

That is a controlled exposure experiment that already happened. The external review
asked us to go back to the microscope and acquire one; we do not need to, because
the biology is held exactly constant across a deliberate change in acquisition
brightness, and the images are already on disk.

The main measurement project deliberately keeps only one acquisition per physical
section (`use_for_measurement`, ADR-0019) -- correct there, because counting the
same tissue twice would inflate n. This project deliberately keeps **both**, and is
therefore a diagnostic project that must never be used for measurement.

    python3 scripts/build_pair_project.py [--dry-run]

Nothing here re-derives the pairing rule: `physical_section_label` is the manifest's
own record of which piece of tissue a series came from, and it is the only safe key.
Matching on `section_label` would pair the wrong tissue -- on tube 51 the rescan's
series `01` is physical section `03` -- and matching on dimensions would pair
whatever happened to be the same size.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import pandas as pd  # noqa: E402

from ihc.ingest.qupath_export import build_project_spec, write_project_spec  # noqa: E402
from ihc.util.config import load_paths  # noqa: E402


def select_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """Rows for physical sections acquired twice -- both acquisitions, or neither."""
    sec = df[df["row_kind"] == "section"].copy()
    sec["phys"] = sec["physical_section_label"].astype(str)
    keep = []
    for (tube, phys), g in sec.groupby(["tube_id", "phys"]):
        if set(g["scan"]) == {"original", "rescan"}:
            keep.append(g)
    if not keep:
        return sec.iloc[0:0]
    return pd.concat(keep).sort_values(["tube_id", "phys", "scan"])


def main(argv: list[str] | None = None) -> int:
    paths = load_paths()
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--custodian-dir", type=Path, default=paths["custodian_root"])
    ap.add_argument("--raw-root", type=Path, default=paths["raw_root"])
    ap.add_argument("--out-dir", type=Path,
                    default=Path(paths["work_root"]) / "qupath_pairs")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    prov = args.custodian_dir / "provenance_manifest.csv"
    key = args.custodian_dir / "blinding_key.json"
    for p in (prov, key):
        if not p.exists():
            print(f"  FAIL  missing {p}", file=sys.stderr)
            return 1

    df = pd.read_csv(prov)
    codes = {int(k): v for k, v in json.loads(key.read_text())["codes"].items()}
    sel = select_pairs(df)

    n_pairs = sel.groupby(["tube_id", "phys"]).ngroups
    print(f"  ok    {n_pairs} physical section(s) acquired twice -> {len(sel)} image(s)")
    for (tube, phys), g in sel.groupby(["tube_id", "phys"]):
        cond = g["condition"].iloc[0]
        print(f"          tube {tube} physical {phys} ({cond}): "
              + " + ".join(f"{r.scan}/series{r.section_label}" for r in g.itertuples()))

    # Both members of a pair share a condition by construction -- they are the same
    # tissue. If that ever fails, the manifest is wrong and nothing downstream is safe.
    bad = [k for k, g in sel.groupby(["tube_id", "phys"]) if g["condition"].nunique() != 1]
    if bad:
        print(f"  FAIL  same tissue recorded with two conditions: {bad}", file=sys.stderr)
        return 1

    spec = build_project_spec(sel, codes, args.raw_root, args.out_dir,
                              project_name="rescan_pairs")
    c = spec["counts"]
    print(f"  ok    {c['images']} image(s) across {c['animals']} animal(s)")
    for s in spec["skipped"]:
        print(f"  WARN  skipped {s.get('code')}_{s.get('section_label')}: {s.get('reason')}")
    for w in spec["warnings"]:
        print(f"  WARN  {w}")

    if args.dry_run:
        print("\n  --dry-run: nothing written")
        return 0
    out = write_project_spec(spec, args.out_dir)
    print(f"\n  ok    {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
