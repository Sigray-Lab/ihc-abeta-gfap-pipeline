#!/usr/bin/env python3
"""Cross-check every recorded staining condition against the pixels.

Reports whether `condition` in the manifest agrees with what the images show, using
the GFAP index (see `ihc.qc.stain_check` for why GFAP and not Aß). It does **not**
change the manifest and must never be used to assign `condition` -- a genuinely
failed stain would be silently relabelled a control.

    python3 scripts/audit_stain_condition.py

Run it after `./ihc manifest`, and whenever payloads arrive for animals that
previously had only an index file -- their condition assignment has never been
checked against pixels until their pixels exist.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ihc.qc.stain_check import audit_manifest  # noqa: E402
from ihc.util.config import load_paths  # noqa: E402

# Sections whose tissue mask is this degenerate are not evidence about staining --
# there is barely a section there to measure. They are reported, not judged.
MIN_TISSUE_FRACTION = 0.05


def main(argv: list[str] | None = None) -> int:
    paths = load_paths()
    work = Path(paths["work_root"])
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", type=Path, default=work / "manifest" / "manifest.csv")
    ap.add_argument("--out", type=Path, default=work / "qc" / "stain_audit.csv")
    args = ap.parse_args(argv)

    if not args.manifest.exists():
        print(f"  FAIL  no manifest at {args.manifest}\n        run ./ihc manifest first",
              file=sys.stderr)
        return 1

    print(f"  reads:  {args.manifest}")
    print(f"  writes: {args.out}\n")

    df = pd.read_csv(args.manifest)
    res = audit_manifest(df)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.out, index=False)

    measured = res[res.gfap_index.notna()]
    print(f"  ok    measured {len(measured)} section(s) across {measured.tube_id.nunique()} animal(s)")

    thin = measured[measured.tissue_fraction <= MIN_TISSUE_FRACTION]
    if len(thin):
        print(f"  WARN  {len(thin)} section(s) have almost no tissue and cannot be judged:")
        for _, r in thin.iterrows():
            print(f"          tube {r.tube_id} {r.scan} s{r.section_label}: "
                  f"{r.tissue_fraction:.1%} tissue")

    good = measured[measured.tissue_fraction > MIN_TISSUE_FRACTION]
    pos = good[good.condition == "positive"]["gfap_index"]
    neg = good[good.condition == "negative"]["gfap_index"]
    if not len(pos) or not len(neg):
        print("  WARN  not enough of both conditions to compare")
        return 0

    disagree = good[
        ((good.condition == "positive") & (good.gfap_index < neg.max()))
        | ((good.condition == "negative") & (good.gfap_index > pos.min()))
    ]

    print(f"  ok    positive (n={len(pos)}): {pos.min():.2f} - {pos.max():.2f}")
    print(f"  ok    negative (n={len(neg)}): {neg.min():.2f} - {neg.max():.2f}")
    print(f"  ok    medians {np.median(pos)/np.median(neg):.1f}x apart")

    if len(disagree) == 0:
        print(f"  ok    NO OVERLAP - margin {pos.min()/neg.max():.2f}x, "
              f"{len(good)}/{len(good)} agree with the record")
        return 0

    print(f"  FAIL  {len(disagree)} section(s) where the pixels and the record disagree.")
    print("        This is for a human to settle at the bench. Do NOT edit condition")
    print("        to match the pixels -- a failed stain looks exactly like this.")
    for _, r in disagree.iterrows():
        print(f"          tube {r.tube_id} {r.scan} s{r.section_label}: "
              f"recorded {r.condition}, gfap_index {r.gfap_index:.2f}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
