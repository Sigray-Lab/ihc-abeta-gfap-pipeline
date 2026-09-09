"""The rescan->original forward model must match scipy, and the original->rescan transform
must be its exact inverse. Guards the sign/centre conventions in register_rescan_regions.py
(external review 2026-09-09, F6: the transfer applied the forward rotation where the
inverse was needed)."""
import importlib.util, sys
from pathlib import Path
import numpy as np
import pytest
from scipy import ndimage as ndi
from shapely.geometry import Point
from shapely.affinity import affine_transform

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "register_rescan_regions.py"


def _mod():
    spec = importlib.util.spec_from_file_location("rrr", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    sys.modules["rrr"] = m
    spec.loader.exec_module(m)
    return m


def _observed(ang, dy, dx, h, w, r0, c0):
    """Where scipy actually puts a Gaussian blob centred at (r0, c0) after fit_rigid's forward model."""
    rr, cc = np.mgrid[0:h, 0:w]
    img = np.exp(-((rr - r0) ** 2 + (cc - c0) ** 2) / (2 * 2.0 ** 2))
    rot = ndi.rotate(img, ang, reshape=False, order=1)
    moved = ndi.shift(rot, (dy, dx), order=1)
    return ndi.center_of_mass(moved)          # (row, col)


@pytest.mark.parametrize("ang,dy,dx", [(3.0, 4.0, -6.0), (-2.5, -3.0, 5.0), (0.0, 2.0, 2.0)])
def test_forward_matrix_matches_scipy(ang, dy, dx):
    m = _mod()
    h, w, r0, c0 = 120, 150, 45.0, 55.0
    obs_r, obs_c = _observed(ang, dy, dx, h, w, r0, c0)
    M = m.forward_matrix_xy(ang, dy, dx, h, w, scale=1.0)
    pred = M @ np.array([c0, r0, 1.0])
    assert abs(pred[0] - obs_c) < 0.5 and abs(pred[1] - obs_r) < 0.5, (pred, (obs_c, obs_r))


@pytest.mark.parametrize("ang,dy,dx,scale", [(3.0, 4.0, -6.0, 8.0), (-2.5, -3.0, 5.0, 1.0)])
def test_original_to_rescan_is_exact_inverse(ang, dy, dx, scale):
    m = _mod()
    h, w = 120, 150
    obs_r, obs_c = _observed(ang, dy, dx, h, w, 45.0, 55.0)
    # observed landmark lives in ORIGINAL coords; the inverse must return it to the rescan source
    landmark = Point(obs_c * scale, obs_r * scale)
    back = affine_transform(landmark, m.original_to_rescan_params(ang, dy, dx, h, w, scale))
    assert abs(back.x - 55.0 * scale) < 0.6 * scale and abs(back.y - 45.0 * scale) < 0.6 * scale, (back.x, back.y)


def test_old_forward_applied_as_inverse_would_fail():
    """The pre-fix code rotated by +ang (forward) instead of -ang: at 3 deg on a 150-px-wide
    frame that misplaces a landmark by several pixels. Keep this so the failure mode is named."""
    m = _mod()
    h, w = 120, 150
    obs_r, obs_c = _observed(3.0, 4.0, -6.0, h, w, 45.0, 55.0)
    M = m.forward_matrix_xy(3.0, 4.0, -6.0, h, w, 1.0)
    wrong = M @ np.array([obs_c, obs_r, 1.0])           # forward applied twice
    assert abs(wrong[0] - 55.0) > 2 or abs(wrong[1] - 45.0) > 2
