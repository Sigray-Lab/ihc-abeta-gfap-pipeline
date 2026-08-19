#!/usr/bin/env python3
"""Carry hand-drawn regions from an original acquisition onto its rescan.

Where a slide was re-imaged, regions may have been drawn on the ORIGINAL acquisition while
the measurement uses the rescan, leaving those sections with no regional endpoint. Rather
than lose the animals or ask for a redraw, the outlines are transformed across.

This is well posed because the two acquisitions are the SAME physical section -- verified
independently by tissue-mask shape agreement, Dice 0.989-0.996, against 0.84 for genuinely
different sections from the same animal. The difference between them is a stage
repositioning: measured offsets are 159-608 um with no reflection or scaling.

Method. Registration is done on DAPI, not on the outlines: DAPI carries strong interior
structure (dentate gyrus, CA layers) whereas an outline offers only a smooth boundary and
fits unstably. A small rotation sweep is combined with phase correlation for translation,
at coarse resolution, and the resulting rigid transform is applied to the polygon vertices.

Acceptance test, and it can fail. The rescan tissue mask is transformed into the original's
frame and scored by Dice against the original tissue mask. Registration must BEAT the
unregistered agreement already measured for that pair; if it does not, the fit is rejected
and the section is left without regions rather than silently mis-assigned.

    python3 scripts/register_rescan_regions.py --src <unzipped regions> [--dry-run]
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.registration import phase_cross_correlation
from shapely.geometry import shape, mapping, MultiPolygon, Polygon
from shapely.affinity import affine_transform

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
from ihc.qc.stain_check import read_pyramid_plane, CH_DAPI          # noqa: E402
from ihc.qc.tissue_mask import tissue_mask                          # noqa: E402

CANON = {"tissue": "Tissue", "tisssue": "Tissue", "hippocampus": "Hippocampus",
         "isocortex": "Isocortex", "isocortes": "Isocortex", "cortex": "Isocortex"}
ANGLES = np.arange(-4.0, 4.01, 0.25)
REGION_COLOR = {"Hippocampus": -65281, "Isocortex": -16776961}


def dice(a: np.ndarray, b: np.ndarray) -> float:
    s = a.sum() + b.sum()
    return float(2 * (a & b).sum() / s) if s else 0.0


def fit_rigid(fixed: np.ndarray, moving: np.ndarray):
    """Best (angle_deg, dy, dx) mapping moving onto fixed, both boolean masks."""
    best = (None, -1.0)
    for ang in ANGLES:
        rot = ndi.rotate(moving.astype(float), ang, reshape=False, order=1) > 0.5
        sh, _, _ = phase_cross_correlation(fixed.astype(float), rot.astype(float),
                                           upsample_factor=4)
        moved = ndi.shift(rot.astype(float), sh, order=1) > 0.5
        d = dice(fixed, moved)
        if d > best[1]:
            best = ((ang, sh[0], sh[1]), d)
    return best


def _polygonal(geom):
    """Keep only the Polygon parts of a geometry.

    A transformed outline clipped to tissue can return a GeometryCollection -- polygons
    plus zero-width line or point contacts where boundaries graze. QuPath's GeoJSON reader
    rejects GeometryCollection outright, and the degenerate parts carry no area, so drop
    them. Matches the same fix in make_region_rois.py.
    """
    if geom.geom_type == "Polygon":
        return geom
    parts = []
    if geom.geom_type in ("MultiPolygon", "GeometryCollection"):
        for g in geom.geoms:
            q = _polygonal(g)
            if q.is_empty:
                continue
            parts.extend(q.geoms if q.geom_type == "MultiPolygon" else [q])
    return MultiPolygon(parts) if parts else Polygon()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--src", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path.home() / "ihc_work/region_rois")
    ap.add_argument("--tissue-rois", type=Path, default=Path.home() / "ihc_work/tissue_rois")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    man = pd.read_csv(REPO.parent / "BlindingKey" / "provenance_manifest.csv")
    man = man[man.row_kind == "section"].copy()
    man["stem"] = ("Image_" + man.tube_id.astype(str)
                   + np.where(man.scan == "rescan", "r", "")
                   + "_s" + man.section_label.astype(str).str.zfill(2))
    lab = np.where(man.scan == "rescan", man.physical_section_label,
                   man.section_label.astype(str).str.zfill(2))
    man["image"] = (man.code + np.where(man.scan == "rescan", "b", "")
                    + "_s" + pd.Series(lab).astype(str).values)
    man["phys"] = pd.Series(lab).astype(str).values

    # pairs: same tube + same physical section, present as both original and rescan
    pairs = []
    for (tube, phys), g in man.groupby(["tube_id", "phys"]):
        if set(g.scan) == {"original", "rescan"}:
            pairs.append((g[g.scan == "original"].iloc[0], g[g.scan == "rescan"].iloc[0]))

    results = []
    for orig, resc in pairs:
        src = next((f for f in a.src.rglob(f"{orig.stem}_annotations.geojson")
                    if "MACOSX" not in str(f)), None)
        if src is None:
            continue
        if not (a.tissue_rois / f"{resc.image}.geojson").exists():
            continue

        px = float(orig.pixel_size_um)
        do, lv = read_pyramid_plane(orig.ets_path, CH_DAPI, pixel_size_um=px,
                                    true_size=(int(orig.width_px), int(orig.height_px)))
        dr, _ = read_pyramid_plane(resc.ets_path, CH_DAPI, pixel_size_um=float(resc.pixel_size_um),
                                   true_size=(int(resc.width_px), int(resc.height_px)))
        mo, mr = tissue_mask(do), tissue_mask(dr)
        h = max(mo.shape[0], mr.shape[0]); w = max(mo.shape[1], mr.shape[1])
        pad = lambda m: np.pad(m, ((0, h - m.shape[0]), (0, w - m.shape[1])))
        mo, mr = pad(mo), pad(mr)

        before = dice(mo, mr)
        (ang, dy, dx), after = fit_rigid(mo, mr)          # maps rescan -> original
        scale = 2 ** lv                                   # coarse px -> full-res px
        results.append(dict(tube=int(orig.tube_id), phys=orig.phys, orig=orig.image,
                            rescan=resc.image, angle_deg=ang,
                            dy_px=dy * scale, dx_px=dx * scale,
                            dice_before=before, dice_after=after,
                            accepted=after > before))
        print(f"  tube {orig.tube_id} phys {orig.phys}: {orig.image} -> {resc.image}   "
              f"rot {ang:+.2f}deg  shift ({dy*scale:+.0f},{dx*scale:+.0f}) px   "
              f"Dice {before:.4f} -> {after:.4f}  {'ACCEPT' if after > before else 'REJECT'}",
              flush=True)
        if a.dry_run or after <= before:
            continue

        # invert: we fitted rescan->original, so original outlines need the inverse
        th = np.deg2rad(ang)
        cy, cx = (h - 1) / 2 * scale, (w - 1) / 2 * scale
        cos, sin = np.cos(th), np.sin(th)
        # forward (rescan->orig): rotate about centre by ang, then shift by (dy,dx)
        # inverse (orig->rescan): shift by (-dy,-dx), then rotate about centre by -ang
        def to_rescan(g):
            g = affine_transform(g, [1, 0, 0, 1, -dx * scale, -dy * scale])
            return affine_transform(g, [cos, sin, -sin, cos,
                                        cx - cos * cx + sin * cy,
                                        cy - sin * cx - cos * cy])

        d = json.loads(src.read_text())
        feats = d.get("features", d) if isinstance(d, dict) else d
        tj = json.loads((a.tissue_rois / f"{resc.image}.geojson").read_text())
        tis = shape(tj["features"][0]["geometry"]).buffer(0)
        for ft in feats:
            name = CANON.get((ft.get("properties", {}).get("name") or "").strip().lower())
            if name not in ("Hippocampus", "Isocortex"):
                continue
            g = shape(ft["geometry"]).buffer(0)
            g = to_rescan(g).intersection(tis)
            if g.is_empty or g.area * px ** 2 / 1e6 < 0.5:
                print(f"      {name}: empty or degenerate after transform, skipped")
                continue
            sub = a.out / name
            sub.mkdir(parents=True, exist_ok=True)
            (sub / f"{resc.image}.geojson").write_text(json.dumps({
                "type": "FeatureCollection", "features": [{
                    "type": "Feature", "id": f"{resc.image}_{name}",
                    "geometry": mapping(_polygonal(g)),
                    "properties": {"objectType": "annotation",
                                   "classification": {"name": name,
                                                      "colorRGB": REGION_COLOR[name]}}}]}))
            print(f"      {name}: {g.area*px**2/1e6:.2f} mm2 written for {resc.image}")

    r = pd.DataFrame(results)
    if len(r):
        r.to_csv(a.out / "rescan_registration.csv", index=False)
        print(f"\n  {r.accepted.sum()} of {len(r)} pairs accepted "
              f"(Dice {r.dice_before.median():.4f} -> {r.dice_after.median():.4f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
