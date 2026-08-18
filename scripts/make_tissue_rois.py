#!/usr/bin/env python3
"""Export the DAPI tissue mask of every section as a GeoJSON annotation for QuPath.

Why this exists
---------------
External review 2026-08-18 identified a stop-level defect: the evaluation scripts called
``createFullImageAnnotation`` and measured the classifier over the whole rectangular frame,
then divided by DAPI tissue area. The endpoint computed was

    100 * (Abeta anywhere in the frame) / (DAPI tissue area)

which is not a percent area of anything -- numerator and denominator referred to different
spatial domains. It was not a small effect for the locally-normalised classifier: its
classified area (Abeta + Negative) came to 1.169x the DAPI tissue area, above 1.0 on every
image measured, because local variance normalisation amplifies noise on flat glass and the
classifier then assigns that glass to a measured class.

Writing the mask out as an annotation makes the numerator and the denominator the same
region by construction, for every classifier.

Coordinates are full-resolution pixels, which is what QuPath's GeoJSON reader expects.
Holes are preserved: a ventricle excluded from the mask stays excluded, rather than being
handed back by the polygon conversion.

    python3 scripts/make_tissue_rois.py [--out DIR] [--project SPEC.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import cv2  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ihc.qc.stain_check import CH_DAPI, read_pyramid_plane  # noqa: E402
from ihc.qc.tissue_mask import tissue_mask  # noqa: E402
from ihc.util.config import load_paths  # noqa: E402

# Reject rings by AREA, never by vertex count. CHAIN_APPROX_SIMPLE reduces a rectangle to
# four corners, so a vertex threshold discards large simple shapes -- verified: a plain
# 160x160 px square exported as zero polygons under the previous 12-vertex rule, outer ring
# included. It did not corrupt the 2026-08-18 run (all 121 ROIs exported 1-4 polygons, and
# the never-acquired tiles in 15 images are all outside the tissue outline, so there were no
# interior holes to lose) but the rule was unsafe. External review round 2, §3.
MIN_RING_AREA_PX = 4.0    # a ring smaller than this cannot represent real support
MIN_RING_VERTICES = 4     # a polygon needs three distinct corners plus closure


def mask_to_geojson(mask: np.ndarray, scale: float, name: str = "Tissue") -> dict:
    """Mask -> one MultiPolygon feature in full-resolution pixel coordinates.

    ``RETR_CCOMP`` gives outer boundaries and holes as separate contours with a
    parent/child hierarchy, so enclosed background stays enclosed.
    """
    m = np.ascontiguousarray(mask.astype(np.uint8))
    contours, hierarchy = cv2.findContours(m, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    polys: list[list] = []
    if hierarchy is not None:
        hierarchy = hierarchy[0]
        def keep(ring) -> bool:
            return (len(ring) >= MIN_RING_VERTICES
                    and cv2.contourArea(ring) >= MIN_RING_AREA_PX)

        for i, c in enumerate(contours):
            if hierarchy[i][3] != -1 or not keep(c):
                continue                                    # child ring, handled below
            rings = [(c[:, 0, :] * scale).tolist()]
            child = hierarchy[i][2]
            while child != -1:
                h = contours[child]
                if keep(h):
                    rings.append((h[:, 0, :] * scale).tolist())
                child = hierarchy[child][0]
            for r in rings:
                r.append(r[0])                              # GeoJSON rings must close
            polys.append(rings)
    return {
        "type": "Feature",
        "geometry": {"type": "MultiPolygon", "coordinates": polys},
        "properties": {"objectType": "annotation",
                       "classification": {"name": name, "colorRGB": -16776961}},
    }


def main(argv: list[str] | None = None) -> int:
    paths = load_paths()
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--project", type=Path,
                    default=Path(paths["work_root"]) / "qupath" / "project_spec.json")
    ap.add_argument("--custodian-dir", type=Path, default=paths["custodian_root"])
    ap.add_argument("--out", type=Path, default=Path(paths["work_root"]) / "tissue_rois")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    man = pd.read_csv(args.custodian_dir / "provenance_manifest.csv")
    man = man[man.row_kind == "section"]
    spec = json.loads(args.project.read_text())
    rows = []
    for e in spec["images"]:
        sel = man[(man.code == e["code"]) & (man.scan == e["scan"]) &
                  (man.section_label.astype(str).str.zfill(2) == str(e["section_label"]).zfill(2))]
        if len(sel) != 1:
            print(f"  FAIL  {e['image_name']}: {len(sel)} manifest rows", file=sys.stderr)
            return 1
        r = sel.iloc[0]
        plane, lv = read_pyramid_plane(r.ets_path, CH_DAPI, pixel_size_um=float(r.pixel_size_um),
                                       true_size=(int(r.width_px), int(r.height_px)))
        mask = tissue_mask(plane)
        scale = float(2 ** lv)
        um2 = mask.sum() * (float(r.pixel_size_um) * scale) ** 2
        gj = {"type": "FeatureCollection", "features": [mask_to_geojson(mask, scale)]}
        (args.out / f"{e['image_name']}.geojson").write_text(json.dumps(gj))
        rows.append(dict(image=e["image_name"], level=lv, area_mm2=um2 / 1e6,
                         n_poly=len(gj["features"][0]["geometry"]["coordinates"]),
                         missing_support_px=int((~np.isfinite(plane)).sum())))
        print(f"  {e['image_name']:<12} {um2/1e6:6.2f} mm2  "
              f"{rows[-1]['n_poly']:>3} polygon(s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(args.out / "tissue_roi_index.csv", index=False)
    print(f"\nwrote {len(df)} ROIs to {args.out}")
    print(f"  area: median {df.area_mm2.median():.2f} mm2, range "
          f"{df.area_mm2.min():.2f}-{df.area_mm2.max():.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
