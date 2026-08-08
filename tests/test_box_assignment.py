"""Regression tests for PAP-pen box assignment (`ihc.ingest.vsi_meta.assign_boxes`).

WHAT CLASS OF BUG THIS PROTECTS AGAINST
---------------------------------------
Each slide carries three or four sections in two PAP-pen boxes.  The two
sections in *one box* share a staining condition: one box got primary antibody
(positive), the other got DAPI + secondary only (negative control).  Box
membership therefore decides which pixels are data and which are a blank.

A wrong box assignment is the worst failure mode in this pipeline because it is
**silent**: negative-control sections get quantified as data, every downstream
number stays in a plausible range, and no other test in the suite can catch it.

The specific bug guarded here is the one CLAUDE_v1.2.md corrects from v1.1: the
belief that sections `_01`+`_02` share a box and `_03`+`_04` share a box.  The
section number is **acquisition order, not slide position**.  Box membership
must be derived from stage X (distance from the label end of the slide) by
sorting and splitting at the largest gap -- never from the section number.

These tests are entirely synthetic (no file I/O) so they run in milliseconds and
stay meaningful on a machine with no data attached.
"""

from __future__ import annotations

import inspect
import itertools
import re

import pytest

from ihc.ingest.vsi_meta import BoxAssignmentError, SeriesMeta, assign_boxes


# ---------------------------------------------------------------------------
# building synthetic SeriesMeta objects
# ---------------------------------------------------------------------------
# `SeriesMeta`'s constructor signature is not part of the documented interface --
# only its attributes are.  We therefore try the real class first and fall back
# to a duck-typed stand-in, so that these tests exercise `assign_boxes` rather
# than the dataclass's argument order.

_DEFAULTS = {
    "name": "20x_DAPI, FITC, Cy3_01",
    "section_label": "01",
    "stage_x_um": 0.0,
    "stage_y_um": 12000.0,
    "exposure_ms": (128.547, 397.927, 1839.999),
    "pixel_size_um": 0.325,
    "width_px": 16384,
    "height_px": 20480,
}


class _StandInSeries:
    """Minimal duck type carrying the documented SeriesMeta attributes."""

    def __init__(self, **fields):
        for key, value in fields.items():
            setattr(self, key, value)

    def __repr__(self):  # keeps assertion output readable
        return f"<series {self.section_label} x={self.stage_x_um:.1f}um>"


def _construct(fields: dict):
    try:
        return SeriesMeta(**fields)
    except TypeError:
        pass
    try:
        params = inspect.signature(SeriesMeta).parameters
        obj = SeriesMeta(**{k: v for k, v in fields.items() if k in params})
    except (TypeError, ValueError):
        return _StandInSeries(**fields)
    for key, value in fields.items():
        if getattr(obj, key, None) is None:
            try:
                setattr(obj, key, value)
            except Exception:
                return _StandInSeries(**fields)
    if getattr(obj, "stage_x_um", None) is None:
        return _StandInSeries(**fields)
    return obj


def make_series(section_label: str, stage_x_mm: float, **overrides):
    """A synthetic section at `stage_x_mm` millimetres from the label end."""
    fields = dict(_DEFAULTS)
    fields["section_label"] = section_label
    fields["name"] = f"20x_DAPI, FITC, Cy3_{section_label}"
    fields["stage_x_um"] = float(stage_x_mm) * 1000.0
    fields.update(overrides)
    return _construct(fields)


def make_series_um(section_label: str, stage_x_um: float, **overrides):
    """A synthetic section at `stage_x_um` micrometres from the label end."""
    return make_series(section_label, stage_x_um / 1000.0, **overrides)


# ---------------------------------------------------------------------------
# reading the result back without depending on what the lists contain
# ---------------------------------------------------------------------------
_LABEL_RE = re.compile(r"(0[1-4])$")


def _label_of(item) -> str:
    """Normalise whatever `assign_boxes` puts in its lists to '01'..'04'."""
    for attr in ("section_label", "label", "section"):
        value = getattr(item, attr, None)
        if value is not None:
            return _label_of(value)
    if isinstance(item, dict):
        for key in ("section_label", "label", "section", "name"):
            if key in item:
                return _label_of(item[key])
    if isinstance(item, bool):
        pytest.fail(f"box list contained a bool, not a section: {item!r}")
    if isinstance(item, int):
        return f"{item:02d}"
    if isinstance(item, str):
        match = _LABEL_RE.search(item)
        if match:
            return match.group(1)
        if item.isdigit():
            return f"{int(item):02d}"
    pytest.fail(
        "cannot read a section label out of a box entry: "
        f"{item!r} (type {type(item).__name__}). Box lists must contain "
        "SeriesMeta objects, section labels, or dicts carrying one."
    )


def _labels(box) -> set[str]:
    assert isinstance(box, (list, tuple)), (
        f"box must be a list or tuple, got {type(box).__name__}"
    )
    return {_label_of(entry) for entry in box}


def _boxes(result) -> tuple[set[str], set[str]]:
    assert isinstance(result, dict), (
        f"assign_boxes must return a dict, got {type(result).__name__}"
    )
    for key in ("near_label", "far_label", "gap_ratio", "within_gaps_mm", "split_gap_mm"):
        assert key in result, (
            f"assign_boxes result is missing the documented key {key!r}; "
            f"got keys {sorted(result)}"
        )
    return _labels(result["near_label"]), _labels(result["far_label"])


def _expect_refusal(series, scenario: str):
    """assign_boxes must raise BoxAssignmentError rather than guess a split."""
    try:
        result = assign_boxes(series)
    except BoxAssignmentError:
        return
    pytest.fail(
        f"assign_boxes accepted {scenario} instead of raising "
        f"BoxAssignmentError. It returned: {result!r}"
    )


# ---------------------------------------------------------------------------
# the ordinary case
# ---------------------------------------------------------------------------
def test_standard_two_plus_two_layout():
    """26 / 33 mm in one box, 45 / 52 mm in the other: a textbook 2+2 slide."""
    series = [
        make_series("01", 26.0),
        make_series("02", 33.0),
        make_series("03", 45.0),
        make_series("04", 52.0),
    ]
    result = assign_boxes(series)
    near, far = _boxes(result)

    assert near == {"01", "02"}, f"near_label box should be the low-stage-X pair, got {near}"
    assert far == {"03", "04"}, f"far_label box should be the high-stage-X pair, got {far}"


def test_near_label_box_is_the_low_stage_x_box():
    """'near_label' means nearest the slide label, i.e. the SMALLEST stage X.

    Stage X measures distance from the label end.  If this orientation is ever
    flipped, every slide in the cohort silently swaps positives and negatives.
    """
    series = [
        make_series("01", 26.0),
        make_series("02", 33.0),
        make_series("03", 45.0),
        make_series("04", 52.0),
    ]
    result = assign_boxes(series)

    near_x = [s.stage_x_um for s in series if s.section_label in _labels(result["near_label"])]
    far_x = [s.stage_x_um for s in series if s.section_label in _labels(result["far_label"])]
    assert max(near_x) < min(far_x), (
        "every near_label section must sit at a smaller stage X than every "
        f"far_label section; near={sorted(near_x)} far={sorted(far_x)}"
    )


def test_every_input_section_lands_in_exactly_one_box():
    """No section may be dropped or duplicated -- a dropped one silently vanishes
    from the animal's measurement, a duplicated one gets counted twice."""
    series = [
        make_series("01", 26.0),
        make_series("02", 33.0),
        make_series("03", 45.0),
        make_series("04", 52.0),
    ]
    result = assign_boxes(series)
    near, far = _boxes(result)

    assert near & far == set(), f"a section was placed in both boxes: {near & far}"
    assert near | far == {"01", "02", "03", "04"}, (
        f"boxes must partition the input; missing {{'01','02','03','04'}} - {near | far}"
    )
    assert len(result["near_label"]) + len(result["far_label"]) == len(series), (
        "a section was duplicated: "
        f"{len(result['near_label'])} + {len(result['far_label'])} != {len(series)}"
    )


# ---------------------------------------------------------------------------
# THE REGRESSION CASE -- tube 49
# ---------------------------------------------------------------------------
# CLAUDE_v1.2.md, section 2, "[CORRECTED from v1.1]":
#
#     v1.1 stated that "sections 01 and 02 always share a condition and 03 and
#     04 always share a condition".  THIS IS FALSE.  In animal 49 the sections
#     sharing a condition are _01 and _04.
#
# Stage X for tube 49 (micrometres, ground truth from Bio-Formats 8.5.0):
#     _01 = 46994.60   _02 = 25600.00   _03 = 32349.60   _04 = 54409.60
#
# Sorted by stage X the slide reads:  [_02 _03]  <-- 14.645 mm gap -->  [_01 _04]
#     within-box gaps: 32349.6-25600.0 = 6.7496 mm
#                      54409.6-46994.6 = 7.4150 mm
#     between-box gap: 46994.6-32349.6 = 14.6450 mm   (1.98-2.17x the within gaps)
#
# The operator scanned the physically-third section first, so the section number
# is acquisition order and carries no positional information whatsoever.
#
# If this test fails, `config/slides.csv` for tube 49 ("far_label" = the
# positives = _01 and _04) will be applied to the wrong two sections and the
# negative controls will be quantified as data.
# ---------------------------------------------------------------------------

TUBE_49_STAGE_X_UM = {
    "01": 46994.60,
    "02": 25600.00,
    "03": 32349.60,
    "04": 54409.60,
}


def _tube_49_series():
    return [make_series_um(label, x) for label, x in TUBE_49_STAGE_X_UM.items()]


def test_tube_49_boxes_are_02_03_and_01_04():
    """REGRESSION: tube 49 proves section number != slide position."""
    result = assign_boxes(_tube_49_series())
    near, far = _boxes(result)

    assert near == {"02", "03"}, (
        "tube 49 near_label box must be {'02','03'} -- the two sections at "
        f"25.60 and 32.35 mm -- got {near}. The 01+02 / 03+04 pairing rule is "
        "WRONG and was corrected in CLAUDE_v1.2.md section 2."
    )
    assert far == {"01", "04"}, (
        "tube 49 far_label box must be {'01','04'} -- the two sections at "
        f"46.99 and 54.41 mm -- got {far}. slides.csv marks tube 49's positives "
        "as far_label, so getting this wrong quantifies the negative controls."
    )


def test_tube_49_split_gap_is_the_14_65_mm_gap():
    """The split must fall on the 14.645 mm between-box gap, not a within-box one."""
    result = assign_boxes(_tube_49_series())

    assert result["split_gap_mm"] == pytest.approx(14.645, abs=0.01), (
        "tube 49's between-box gap is 14.645 mm; assign_boxes reported "
        f"{result['split_gap_mm']!r} mm (units must be millimetres)"
    )


def test_tube_49_within_box_gaps_are_reported_in_mm():
    """Both within-box gaps (6.75 and 7.415 mm) must be reported, in millimetres."""
    result = assign_boxes(_tube_49_series())
    gaps = sorted(float(g) for g in result["within_gaps_mm"])

    assert len(gaps) == 2, (
        f"tube 49 has two within-box gaps (one per box), got {len(gaps)}: {gaps}"
    )
    assert gaps[0] == pytest.approx(6.7496, abs=0.01), (
        f"smaller within-box gap should be 6.7496 mm, got {gaps[0]}"
    )
    assert gaps[1] == pytest.approx(7.4150, abs=0.01), (
        f"larger within-box gap should be 7.4150 mm, got {gaps[1]}"
    )


def test_tube_49_gap_ratio_is_in_the_documented_range():
    """The between/within gap ratio is 1.7-2.4x in every animal checked."""
    result = assign_boxes(_tube_49_series())
    ratio = float(result["gap_ratio"])

    assert ratio > 1.3, (
        f"gap_ratio {ratio} is below the 1.3 ambiguity floor, yet tube 49's "
        "split is unambiguous (14.645 mm vs 6.75/7.42 mm)"
    )
    assert 1.7 <= ratio <= 2.4, (
        f"gap_ratio {ratio} is outside the documented 1.7-2.4x range for real "
        "slides; check whether the ratio is taken against the mean, max or min "
        "within-box gap"
    )


@pytest.mark.parametrize("permutation", list(itertools.permutations(["01", "02", "03", "04"])))
def test_input_order_never_changes_the_answer(permutation):
    """Shuffling the input must not move a single section.

    Tube 49's stage coordinates are deliberately non-monotonic in section
    number, so any implementation that leans on list order, on the section
    label, or on stack order will disagree with itself across these 24 orderings.
    """
    series = [make_series_um(label, TUBE_49_STAGE_X_UM[label]) for label in permutation]
    result = assign_boxes(series)
    near, far = _boxes(result)

    assert near == {"02", "03"}, (
        f"input order {permutation} changed the near_label box to {near}; "
        "assign_boxes must sort by stage X and ignore input order entirely"
    )
    assert far == {"01", "04"}, (
        f"input order {permutation} changed the far_label box to {far}"
    )
    assert result["split_gap_mm"] == pytest.approx(14.645, abs=0.01), (
        f"input order {permutation} changed split_gap_mm to {result['split_gap_mm']}"
    )


# ---------------------------------------------------------------------------
# three-section slides (tubes 30, 33, 34, 42, 53, 54)
# ---------------------------------------------------------------------------
def test_three_sections_split_two_plus_one():
    """2+1: two sections near the label, one far.  Must not be rejected."""
    series = [
        make_series("01", 26.0),
        make_series("02", 33.0),
        make_series("03", 45.0),
    ]
    result = assign_boxes(series)
    near, far = _boxes(result)

    assert near == {"01", "02"}, f"expected near_label {{'01','02'}}, got {near}"
    assert far == {"03"}, f"expected far_label {{'03'}}, got {far}"
    assert result["split_gap_mm"] == pytest.approx(12.0, abs=0.01)


def test_three_sections_split_one_plus_two():
    """1+2: one section near the label, two far.  Also legal."""
    series = [
        make_series("01", 26.0),
        make_series("02", 38.0),
        make_series("03", 45.0),
    ]
    result = assign_boxes(series)
    near, far = _boxes(result)

    assert near == {"01"}, f"expected near_label {{'01'}}, got {near}"
    assert far == {"02", "03"}, f"expected far_label {{'02','03'}}, got {far}"
    assert result["split_gap_mm"] == pytest.approx(12.0, abs=0.01)


def test_three_section_split_ignores_input_order():
    """The 2+1 answer must survive every input ordering too."""
    labels_x = {"01": 26.0, "02": 33.0, "03": 45.0}
    for permutation in itertools.permutations(labels_x):
        series = [make_series(label, labels_x[label]) for label in permutation]
        near, far = _boxes(assign_boxes(series))
        assert near == {"01", "02"}, f"order {permutation} gave near_label {near}"
        assert far == {"03"}, f"order {permutation} gave far_label {far}"


# ---------------------------------------------------------------------------
# layouts that must be refused
# ---------------------------------------------------------------------------
# `assign_boxes` must raise rather than guess.  A guessed split is exactly the
# silent failure the whole design is trying to avoid: it produces a plausible
# answer with no signal that anything went wrong.
# ---------------------------------------------------------------------------

def test_raises_on_three_plus_one_split():
    """3+1 is not a legal PAP-pen layout -- a box holds at most two sections.

    Gaps here are 3.0 / 3.0 / 18.0 mm, so the split point itself is
    unambiguous (ratio 6.0).  This isolates the cardinality check.
    """
    series = [
        make_series("01", 26.0),
        make_series("02", 29.0),
        make_series("03", 32.0),
        make_series("04", 50.0),
    ]
    _expect_refusal(series, "a 3+1 split (sections at 26.0/29.0/32.0/50.0 mm)")


def test_raises_on_one_plus_three_split():
    """1+3 is the mirror image of 3+1 and equally illegal."""
    series = [
        make_series("01", 26.0),
        make_series("02", 44.0),
        make_series("03", 47.0),
        make_series("04", 50.0),
    ]
    _expect_refusal(series, "a 1+3 split (sections at 26.0/44.0/47.0/50.0 mm)")


def test_raises_on_four_sections_in_a_single_box():
    """4+0: all four sections clustered inside one box, no second box at all.

    Physically this is a slide where only one PAP-pen box was used, so there is
    no between-box gap to split on.  A largest-gap split will still return
    *some* partition of these four points, and that partition is meaningless --
    it must be refused, not reported.  Sections span 26.0-27.5 mm here, far
    tighter than the ~7 mm within-box spacing seen in real slides.
    """
    series = [
        make_series("01", 26.0),
        make_series("02", 26.5),
        make_series("03", 27.0),
        make_series("04", 27.5),
    ]
    _expect_refusal(
        series,
        "a 4+0 layout (all four sections clustered at 26.0-27.5 mm, no second box)",
    )


def test_raises_on_ambiguous_gap():
    """Evenly spaced sections have no box boundary to find (gap_ratio < 1.3).

    Gaps are 7.0 / 7.2 / 6.8 mm.  The largest-gap split would return a tidy
    2+2, so the cardinality check passes and only the gap ratio (1.03-1.06)
    can catch it.  This isolates the ambiguity check.
    """
    series = [
        make_series("01", 26.0),
        make_series("02", 33.0),
        make_series("03", 40.2),
        make_series("04", 47.0),
    ]
    _expect_refusal(
        series,
        "four evenly spaced sections (gaps 7.0/7.2/6.8 mm, gap ratio ~1.04)",
    )


def test_raises_on_ambiguous_gap_with_three_sections():
    """Three evenly spaced sections are ambiguous for the same reason."""
    series = [
        make_series("01", 26.0),
        make_series("02", 33.0),
        make_series("03", 40.1),
    ]
    _expect_refusal(
        series, "three evenly spaced sections (gaps 7.0/7.1 mm, gap ratio ~1.01)"
    )


def test_box_assignment_error_is_an_exception_subclass():
    """Callers catch BoxAssignmentError specifically; it must be catchable."""
    assert issubclass(BoxAssignmentError, Exception), (
        "BoxAssignmentError must derive from Exception"
    )


@pytest.mark.parametrize("n_sections", [0, 1, 2])
def test_raises_on_too_few_sections(n_sections):
    """Fewer than three sections cannot form two boxes of the documented shape.

    Every animal in the cohort has three or four tissue series; anything less
    means the series filter dropped something and must not be papered over.
    """
    labels = ["01", "02"][:n_sections]
    series = [make_series(label, 26.0 + 19.0 * i) for i, label in enumerate(labels)]
    _expect_refusal(series, f"a slide with only {n_sections} tissue section(s)")


# ---------------------------------------------------------------------------
# reported geometry must be internally consistent
# ---------------------------------------------------------------------------
def test_split_gap_exceeds_every_within_box_gap():
    """The reported split must genuinely be the largest gap on the slide."""
    series = [
        make_series("01", 26.0),
        make_series("02", 33.0),
        make_series("03", 45.0),
        make_series("04", 52.0),
    ]
    result = assign_boxes(series)
    split = float(result["split_gap_mm"])
    within = [float(g) for g in result["within_gaps_mm"]]

    assert within, "within_gaps_mm must not be empty for a 2+2 slide"
    assert split > max(within), (
        f"split_gap_mm ({split}) must exceed every within-box gap {within}"
    )


def test_gap_ratio_is_consistent_with_the_reported_gaps():
    """gap_ratio must be derivable from the gaps the function itself reports.

    Any of mean / max / min within-box gap is a defensible denominator; a value
    consistent with none of them means the ratio is being computed from
    something other than the reported geometry.  Tube 49's geometry is used
    because its two within-box gaps differ (6.75 vs 7.42 mm), so the three
    candidate definitions give three different answers and the test can
    actually distinguish them.
    """
    result = assign_boxes(_tube_49_series())
    split = float(result["split_gap_mm"])
    within = [float(g) for g in result["within_gaps_mm"]]
    ratio = float(result["gap_ratio"])

    candidates = {
        "mean": split / (sum(within) / len(within)),
        "max": split / max(within),
        "min": split / min(within),
    }
    assert any(ratio == pytest.approx(v, rel=1e-6) for v in candidates.values()), (
        f"gap_ratio {ratio} matches none of the defensible definitions "
        f"{ {k: round(v, 4) for k, v in candidates.items()} } derived from "
        f"split_gap_mm={split} and within_gaps_mm={within}"
    )


def test_gaps_are_millimetres_not_micrometres():
    """A units slip here turns a 12 mm gap into 12000 and defeats every threshold."""
    series = [
        make_series("01", 26.0),
        make_series("02", 33.0),
        make_series("03", 45.0),
        make_series("04", 52.0),
    ]
    result = assign_boxes(series)

    assert 1.0 < float(result["split_gap_mm"]) < 75.0, (
        f"split_gap_mm={result['split_gap_mm']} is not a plausible millimetre "
        "value on a 75 x 25 mm slide -- check for a micrometre/millimetre mix-up"
    )
    for gap in result["within_gaps_mm"]:
        assert 0.0 < float(gap) < 75.0, (
            f"within-box gap {gap} is not a plausible millimetre value on a "
            "75 x 25 mm slide"
        )
