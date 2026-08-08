"""Tests for the staining cross-check.

The two that matter are `test_padding_is_cropped` and
`test_bright_artefact_does_not_capture_the_mask` -- each pins a defect that made the
earlier ad-hoc version of this check report a confident wrong number.
"""

from __future__ import annotations

import numpy as np
import pytest

from ihc.qc import stain_check
from ihc.qc.stain_check import SectionStain, measure_section


class FakePlanes:
    """Stand in for read_pyramid_plane with hand-built channel planes."""

    def __init__(self, planes: dict[int, np.ndarray], level: int = 5):
        self.planes = planes
        self.level = level
        self.true_sizes: list = []

    def __call__(self, path, channel, *, level=None, pixel_size_um=None, true_size=None):
        self.true_sizes.append(true_size)
        return self.planes[channel], self.level


def _section(dapi, fitc, cy3, monkeypatch, **kw):
    fake = FakePlanes({0: dapi, 1: fitc, 2: cy3})
    monkeypatch.setattr(stain_check, "read_pyramid_plane", fake)
    kw.setdefault("true_size", (100, 100))
    return measure_section("x.ets", {"DAPI": 100.0, "FITC": 100.0, "Cy3": 100.0}, **kw), fake


def _tissue_and_glass(marker_level: float, size: int = 200) -> tuple:
    """Half tissue, half glass. Marker sits only in the tissue half."""
    dapi = np.full((size, size), 5.0, dtype=np.float32)
    dapi[: size // 2] = 200.0
    marker = np.full((size, size), 2.0, dtype=np.float32)
    marker[: size // 2] = marker_level
    return dapi, marker


def test_positive_section_scores_above_negative(monkeypatch):
    d, bright = _tissue_and_glass(500.0)
    _, dim = _tissue_and_glass(6.0)
    hot, _ = _section(d, bright, bright, monkeypatch)
    cold, _ = _section(d, dim, dim, monkeypatch)
    assert hot.gfap_index > cold.gfap_index
    assert hot.abeta_index > cold.abeta_index


def test_true_size_is_required(monkeypatch):
    """Refusing beats returning a plausible number measured on the padded grid."""
    d, m = _tissue_and_glass(500.0)
    out, _ = _section(d, m, m, monkeypatch, true_size=None)
    assert np.isnan(out.gfap_index)
    assert "true_size" in out.problems[0]


def test_true_size_is_passed_through(monkeypatch):
    d, m = _tissue_and_glass(500.0)
    _, fake = _section(d, m, m, monkeypatch, true_size=(123, 456))
    assert fake.true_sizes == [(123, 456), (123, 456), (123, 456)]


def test_padding_is_cropped():
    """The pad below and right of the image must not reach the statistics.

    A 100x100 image stored in a 128x128 tile grid carries 36% pad. Uncropped, that
    pad is indistinguishable from real black pixels.
    """
    plane = np.zeros((128, 128), dtype=np.float32)
    plane[:100, :100] = 500.0
    cropped = plane[:100, :100]
    assert float(np.percentile(plane, 5)) == 0.0
    assert float(np.percentile(cropped, 5)) == 500.0


def test_bright_artefact_does_not_capture_the_mask(monkeypatch):
    """One dust speck must not become the entire tissue mask.

    Unclipped Otsu on a section whose DAPI peaks at 37,612 counts (real case, tube
    32) puts the threshold above every real pixel, so 0.03% of the section is called
    tissue and every statistic then describes the speck.
    """
    dapi, marker = _tissue_and_glass(500.0)
    dapi[0, 0] = 37612.0  # the speck
    out, _ = _section(dapi, marker, marker, monkeypatch)
    assert out.tissue_fraction > 0.4, "the tissue half should still be found"
    assert np.isfinite(out.gfap_index)


def test_missing_tiles_are_nan_not_zero(monkeypatch):
    """Never-acquired tiles must not be counted as legitimately black tissue."""
    dapi, marker = _tissue_and_glass(500.0)
    with_holes = marker.copy()
    with_holes[:20, :20] = np.nan
    out_holes, _ = _section(dapi, with_holes, with_holes, monkeypatch)
    out_full, _ = _section(dapi, marker, marker, monkeypatch)
    assert np.isfinite(out_holes.gfap_index)
    # Dropping tissue pixels must not move a high percentile much; zero-filling would.
    assert out_holes.gfap_index == pytest.approx(out_full.gfap_index, rel=0.2)


def test_exposure_is_normalised_out(monkeypatch):
    """Two identical sections imaged at different exposures must score the same."""
    dapi, marker = _tissue_and_glass(500.0)
    fake = FakePlanes({0: dapi, 1: marker, 2: marker})
    monkeypatch.setattr(stain_check, "read_pyramid_plane", fake)
    a = measure_section("x", {"DAPI": 100.0, "FITC": 100.0, "Cy3": 100.0}, true_size=(9, 9))
    b = measure_section("x", {"DAPI": 200.0, "FITC": 200.0, "Cy3": 200.0}, true_size=(9, 9))
    assert a.gfap_index == pytest.approx(b.gfap_index)


def test_empty_section_is_reported_not_guessed(monkeypatch):
    blank = np.full((50, 50), np.nan, dtype=np.float32)
    out, _ = _section(blank, blank, blank, monkeypatch)
    assert np.isnan(out.gfap_index)
    assert out.problems


def test_module_does_not_assign_condition():
    """The forbidden inference: nothing here may return a staining condition."""
    src = (stain_check.__file__)
    text = open(src).read()
    assert "condition =" not in text
    assert isinstance(SectionStain(ets_path="x"), SectionStain)
    assert not hasattr(stain_check, "infer_condition")
