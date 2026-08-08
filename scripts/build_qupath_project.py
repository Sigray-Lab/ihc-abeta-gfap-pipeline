#!/usr/bin/env python3
"""Build the blinded QuPath delineation project from the custodian's manifest.

This is a **custodian action**, which is why it lives here and not behind `./ihc`.
It needs the tube_id -> code mapping to join coded sections back to raw files, so
whoever runs it can see the blinding key by construction. The person who will draw
the regions must not run it; they receive the built project.

Run it whenever `./ihc doctor` reports the derived artefacts have gone stale --
which is what happens when payloads arrive for animals that previously had only an
index file. That has now happened twice, both times leaving the project quietly
short of animals, which is the failure this script exists to make one command long.

    ./ihc manifest && ./ihc blind --seed <recorded seed> --force
    python3 scripts/build_qupath_project.py

Rows are selected on `use_for_measurement` -- the single column manifest.py derives
for this purpose. Do not re-derive that rule here: `analysis_include` and
`scan_is_preferred` each look like the right filter and each is wrong on its own
(ADR-0019), and having a second copy of the rule is how they drift apart.
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


def main(argv: list[str] | None = None) -> int:
    paths = load_paths()
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--custodian-dir", type=Path, default=paths["custodian_root"])
    ap.add_argument("--raw-root", type=Path, default=paths["raw_root"])
    ap.add_argument("--out-dir", type=Path, default=Path(paths["work_root"]) / "qupath")
    ap.add_argument("--dry-run", action="store_true",
                    help="build and report the spec, write nothing")
    args = ap.parse_args(argv)

    prov = args.custodian_dir / "provenance_manifest.csv"
    key = args.custodian_dir / "blinding_key.json"
    for p in (prov, key):
        if not p.exists():
            print(f"  FAIL  missing {p}\n        run ./ihc blind first", file=sys.stderr)
            return 1

    print(f"  reads:  {prov}")
    print(f"          {key}")
    print(f"          {args.raw_root}")
    print(f"  writes: {args.out_dir}" if not args.dry_run else "  writes: nothing (--dry-run)")
    print()

    df = pd.read_csv(prov)
    codes = {int(k): v for k, v in json.loads(key.read_text())["codes"].items()}

    selected = df[df["use_for_measurement"].astype(bool)].copy()
    print(f"  ok    {len(selected)} of {len(df)} section(s) selected on use_for_measurement")

    spec = build_project_spec(selected, codes, args.raw_root, args.out_dir)

    c = spec["counts"]
    print(f"  ok    {c['images']} image(s): {c['positive']} positive, "
          f"{c['negative']} negative, across {c['animals']} animal(s)")

    if spec["skipped"]:
        print(f"  WARN  {len(spec['skipped'])} section(s) skipped:")
        for s in spec["skipped"][:10]:
            print(f"          {s.get('code')}_{s.get('section_label')}: {s.get('reason')}")
    if spec["excluded"]:
        print(f"  WARN  {len(spec['excluded'])} section(s) excluded")
    for w in spec["warnings"]:
        print(f"  WARN  {w}")

    if args.dry_run:
        print("\n  --dry-run: nothing written")
        return 0

    out = write_project_spec(spec, args.out_dir)
    print(f"\n  ok    {out}")
    print(f"  ok    coded symlinks in {spec['images_dir']}")
    print("\n  Next: open QuPath and run the import script against that spec.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
