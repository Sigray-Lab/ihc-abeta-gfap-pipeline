#!/usr/bin/env python3
"""One classifier per region per section, carrying that region's own offset and gain.

Within-region normalisation, agreed 2026-08-19. Same construction as the whole-section
global arm, except the gain is estimated inside the region being measured:

    I' = (I - glass_section) / (region_p50 - glass_section)

The offset still comes from the section's off-tissue glass, because there is no glass
inside a hippocampus. Only the gain is regional.

Why this is not the per-region normalisation that reviewers warn against: that objection is
that equalising every region destroys comparisons BETWEEN regions. Here regions are never
compared to each other -- the endpoint is arm versus arm inside one region -- so the cost
does not apply. The risk that does carry over is an anchor that moves with treatment, and
that is checked separately before these numbers are used.

    python3 scripts/make_region_classifiers.py --anchors /tmp/region_anchors.csv
"""
from __future__ import annotations

import argparse, json
from pathlib import Path

import pandas as pd

DEFAULT_CLS = Path.home() / "ihc_work/qupath/qupath/classifiers/pixel_classifiers"
MIN_GAIN = 50.0     # dividing by a smaller gain amplifies noise beyond usefulness


def find_ops(node, out):
    if isinstance(node, dict):
        if node.get("type") in ("op.core.subtract", "op.core.divide"):
            out.append(node)
        for v in node.values():
            find_ops(v, out)
    elif isinstance(node, list):
        for v in node:
            find_ops(v, out)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--classifier", default="v3_global")
    ap.add_argument("--anchors", type=Path, required=True,
                    help="CSV with columns image, region, glass_p50, r_p50")
    ap.add_argument("--dir", type=Path, default=DEFAULT_CLS)
    ap.add_argument("--prefix", default="rx")
    a = ap.parse_args()

    base = a.dir / f"{a.classifier}.json"
    if len(find_ops(json.loads(base.read_text()), [])) != 2:
        print("  FAIL reference chain does not hold exactly one subtract and one divide")
        return 1

    d = pd.read_csv(a.anchors)
    d["gain"] = d.r_p50 - d.glass_p50
    skip = d[(d.gain < MIN_GAIN) | ~d.glass_p50.notna()]
    for _, r in skip.iterrows():
        print(f"  SKIP {r.image}/{r.region}: gain {r.gain:.0f} below {MIN_GAIN:.0f}")
    d = d.drop(skip.index)

    for _, r in d.iterrows():
        j = json.loads(base.read_text())
        sub, div = find_ops(j, [])
        sub["values"] = [float(r.glass_p50)]
        div["values"] = [float(r.gain)]
        (a.dir / f"{a.prefix}_{r.region}_{r.image}.json").write_text(json.dumps(j))

    print(f"  wrote {len(d)} region classifiers")
    for reg, g in d.groupby("region"):
        print(f"    {reg:<12} n={len(g):>3}  gain {g.gain.min():.0f}-{g.gain.max():.0f} "
              f"(median {g.gain.median():.0f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
