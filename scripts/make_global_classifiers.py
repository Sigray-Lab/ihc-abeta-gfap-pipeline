#!/usr/bin/env python3
"""Write one classifier per section, carrying that section's global correction.

QuPath bakes a single op chain into a classifier, so a per-image constant cannot live in
a shared chain. Rather than rewrite 42 GB of pixels, this copies the trained reference
classifier once per section and rewrites the two leading constants:

    I' = (I - glass) / (tissue_p50 - glass)

The forest itself is byte-identical in every copy -- only the two numbers change. That is
the whole of the global arm: one offset and one gain per section, nothing regional, and
every spatial frequency left intact.

    python3 scripts/make_global_classifiers.py --classifier v3_global
"""
from __future__ import annotations

import argparse, json, shutil
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DEFAULT_CLS = Path.home() / "ihc_work/qupath/qupath/classifiers/pixel_classifiers"


def find_ops(node, out):
    """Collect every dict that looks like a subtract/divide op, in chain order."""
    if isinstance(node, dict):
        t = node.get("type", "")
        if t in ("op.core.subtract", "op.core.divide"):
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
                    help="CSV with columns image, glass_p50, t_p50 (see docs/normalisation.md)")
    ap.add_argument("--dir", type=Path, default=DEFAULT_CLS)
    ap.add_argument("--prefix", default="gx")
    a = ap.parse_args()

    base = a.dir / f"{a.classifier}.json"
    ref = json.loads(base.read_text())
    probe = find_ops(ref, [])
    if len(probe) != 2:
        print(f"  FAIL expected exactly one subtract and one divide in the reference "
              f"chain, found {len(probe)}")
        return 1

    anchors = pd.read_csv(a.anchors)
    anchors["gain"] = anchors.t_p50 - anchors.glass_p50
    bad = anchors[anchors.gain <= 0]
    if len(bad):
        print(f"  FAIL non-positive gain on {list(bad.image)} -- tissue not above glass")
        return 1

    written = []
    for _, r in anchors.iterrows():
        j = json.loads(base.read_text())
        sub, div = find_ops(j, [])
        if sub["type"] != "op.core.subtract" or div["type"] != "op.core.divide":
            print(f"  FAIL op order unexpected on {r.image}")
            return 1
        sub["values"] = [float(r.glass_p50)]
        div["values"] = [float(r.gain)]
        out = a.dir / f"{a.prefix}_{r.image}.json"
        out.write_text(json.dumps(j))
        written.append((r.image, r.glass_p50, r.gain))

    print(f"  wrote {len(written)} per-section classifiers to {a.dir}")
    g = pd.DataFrame(written, columns=["image", "glass", "gain"])
    print(f"  glass {g.glass.min():.0f}-{g.glass.max():.0f}   "
          f"gain {g.gain.min():.0f}-{g.gain.max():.0f} ({g.gain.max()/g.gain.min():.1f}x)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
