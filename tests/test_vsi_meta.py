"""Ground-truth tests for `ihc.ingest.vsi_meta.read_vsi_meta` on the real cohort.

WHAT CLASS OF BUG THIS PROTECTS AGAINST
---------------------------------------
`read_vsi_meta` walks the Olympus VSI tag block by hand, because Bio-Formats
returns almost nothing from a `.vsi` whose pixel payload is absent -- and 23 of
31 animals currently have no payload on this machine.  Hand-rolled binary
parsing fails in ways that look like data rather than like bugs:

*   a tag walker that follows `nextField` relative to the wrong base drifts and
    starts reading a neighbouring tag's bytes as a stage coordinate;
*   a series-name filter written as an equality test silently drops all four of
    tube 60's series, because tube 60 alone is named `60_20x_DAPI, FITC, Cy3_01`;
*   an exposure reader that collapses values to a set cannot tell "a channel is
    missing" from "two channels share a value";
*   a micro/milli/nano unit slip in exposure or a micro/milli slip in stage
    position produces numbers that are wrong by exactly 1000x and still plot.

Every expectation below is pinned to an externally established fact: Bio-Formats
8.5.0 output for stage positions, the validated cohort exposure sweep, and the
physical section counts recorded in CLAUDE_v1.2.md section 2.

These tests need only the small `.vsi` index files (1.37-1.91 MB each), which
are always present -- no pixel payload required.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

from ihc.ingest.vsi_meta import read_vsi_meta


def _cohort_constants():
    """Load the cohort constants from conftest, whatever the pytest import mode."""
    try:
        from conftest import ALL_TUBES, THREE_SECTION_TUBES  # noqa: PLC0415

        return ALL_TUBES, THREE_SECTION_TUBES
    except ImportError:
        spec = importlib.util.spec_from_file_location(
            "_ihc_tests_conftest", Path(__file__).with_name("conftest.py")
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module.ALL_TUBES, module.THREE_SECTION_TUBES


ALL_TUBES, THREE_SECTION_TUBES = _cohort_constants()

pytestmark = pytest.mark.requires_data


# ---------------------------------------------------------------------------
# ground truth
# ---------------------------------------------------------------------------
# PositionX in micrometres, read by Bio-Formats 8.5.0 from the payload-bearing
# datasets.  Note tube 49: stage X is NOT monotonic in section number, which is
# the fact that killed the 01+02 / 03+04 pairing rule.
GROUND_TRUTH_STAGE_X_UM = {
    29: {"01": 26996.49, "02": 34637.33, "03": 53462.17, "04": 61301.52},
    41: {"01": 26311.95, "02": 33301.10, "03": 44999.89, "04": 51915.11},
    49: {"01": 46994.60, "02": 25600.00, "03": 32349.60, "04": 54409.60},
    51: {"01": 26448.86, "02": 33171.04, "03": 44452.26, "04": 51299.02},
    60: {"01": 25731.46, "02": 33455.81, "03": 47027.49, "04": 55125.60},
}
STAGE_X_TOLERANCE_UM = 0.5

# Exposure per tissue series, (DAPI, FITC, Cy3) in milliseconds.  Tag 100002
# stores MICROseconds, so a reader that forgets to divide by 1000 fails here.
STANDARD_EXPOSURE_MS = (128.547, 397.927, 1839.999)
DEVIANT_EXPOSURE_MS = {
    51: (60.517, 240.822, 145.677),
    60: (128.547, 397.927, 397.927),
}
EXPOSURE_TOLERANCE_MS = 0.002  # 2 us; the underlying values are integers in us

NOMINAL_PIXEL_SIZE_UM = 0.325
PIXEL_SIZE_TOLERANCE_UM = 1e-4

VALID_SECTION_LABELS = ("01", "02", "03", "04")


def expected_series_count(tube: int) -> int:
    return 3 if tube in THREE_SECTION_TUBES else 4


# ---------------------------------------------------------------------------
# small normalisers, so a defensible representation choice is not a failure
# ---------------------------------------------------------------------------
def _label(series) -> str:
    value = getattr(series, "section_label", None)
    assert value is not None, f"series {series!r} has no .section_label"
    if isinstance(value, int):
        return f"{value:02d}"
    text = str(value).strip()
    return f"{int(text):02d}" if text.isdigit() else text


def _box_label(entry) -> str:
    """Normalise a `assign_boxes` box entry to '01'..'04'.

    The boxes may hold SeriesMeta objects or bare section labels; both are
    defensible, so neither is treated as a failure here.
    """
    if isinstance(entry, str):
        text = entry.strip()
        return f"{int(text):02d}" if text.isdigit() else text
    return _label(entry)


CHANNEL_ORDER = ("DAPI", "FITC", "Cy3")


def _exposure_triplet(series) -> tuple[float, float, float]:
    """`.exposure_ms` as (DAPI, FITC, Cy3), from either a mapping or a sequence.

    A channel-name mapping is the preferred shape -- CLAUDE_v1.2.md section 5
    says to "index by channel name, not position" -- but a positional triplet in
    acquisition order is also well defined, so both are accepted here.
    """
    value = getattr(series, "exposure_ms", None)
    assert value is not None, f"series {series!r} has no .exposure_ms"
    assert not isinstance(value, (str, bytes)), (
        f".exposure_ms must hold three per-channel exposures, got the string {value!r}"
    )

    if isinstance(value, dict):
        lookup = {str(k).strip().upper(): v for k, v in value.items()}
        missing = [c for c in CHANNEL_ORDER if c.upper() not in lookup]
        assert not missing, (
            f".exposure_ms is keyed by channel but is missing {missing}; "
            f"got keys {sorted(value)}. Channels are C0 DAPI, C1 FITC, C2 Cy3."
        )
        assert len(value) == 3, (
            "each tissue series carries exactly three channels; .exposure_ms "
            f"has {len(value)} keys: {sorted(value)}"
        )
        return tuple(float(lookup[c.upper()]) for c in CHANNEL_ORDER)

    try:
        values = tuple(float(v) for v in value)
    except TypeError:
        pytest.fail(
            ".exposure_ms must be a mapping of channel name to exposure, or a "
            f"three-element sequence (DAPI, FITC, Cy3); got the scalar {value!r}"
        )
    assert len(values) == 3, (
        "each tissue series carries exactly three channels (DAPI, FITC, Cy3), "
        f"so .exposure_ms must have three entries; got {len(values)}: {values}"
    )
    return values


# ---------------------------------------------------------------------------
# the index-file inventory
# ---------------------------------------------------------------------------
def test_the_expected_31_index_files_are_present(vsi_paths):
    """Tubes 29-58 and 60.  Mouse 59 was excluded before imaging."""
    assert set(vsi_paths) == set(ALL_TUBES), (
        f"missing tubes {sorted(set(ALL_TUBES) - set(vsi_paths))}, "
        f"unexpected tubes {sorted(set(vsi_paths) - set(ALL_TUBES))}"
    )
    assert len(vsi_paths) == 31, f"expected 31 index files, found {len(vsi_paths)}"


def test_rescan_files_do_not_shadow_the_originals(data_root, vsi_paths):
    """`Rescan/` holds re-acquired 51 and 60 under the SAME filenames.

    A recursive glob would yield `Image_51.vsi` twice with different pixel data
    and different exposure, and whichever came last would silently win.  The
    rescans are a separate acquisition that has to be tracked explicitly, not
    merged into the cohort sweep by accident.
    """
    rescan = data_root / "Rescan"
    if not rescan.is_dir():
        pytest.skip(f"no Rescan folder under {data_root}")

    rescanned = sorted(
        int(p.stem.split("_", 1)[1])
        for p in rescan.glob("Image_*.vsi")
        if p.stem.split("_", 1)[-1].isdigit()
    )
    assert rescanned, f"{rescan} exists but holds no Image_*.vsi"

    for tube, path in vsi_paths.items():
        assert path.parent == data_root, (
            f"tube {tube} resolved to {path}, which is below the data root. "
            "The cohort sweep must not recurse into Rescan/ -- it would shadow "
            f"the original acquisition for tubes {rescanned}."
        )


def test_no_index_file_for_mouse_59(vsi_paths):
    """Mouse 59 was excluded pre-imaging for a mounting fault; it has no scan."""
    assert 59 not in vsi_paths, (
        "Image_59.vsi appeared. Mouse 59 was excluded before imaging (PAP-pen "
        "and antibody leak); if a file now exists the exclusion record is stale."
    )


# ---------------------------------------------------------------------------
# every file parses
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("tube", ALL_TUBES)
def test_index_file_parses(tube, read_meta_cached):
    """All 31 `.vsi` index files must parse without the pixel payload."""
    meta = read_meta_cached(tube)
    assert meta is not None, f"read_vsi_meta returned None for tube {tube}"
    for attr in ("tube_id", "n_tissue_series", "series", "warnings"):
        assert hasattr(meta, attr), (
            f"VsiMeta for tube {tube} is missing the documented attribute .{attr}"
        )


@pytest.mark.parametrize("tube", ALL_TUBES)
def test_warnings_is_a_sequence(tube, read_meta_cached):
    """`.warnings` must always be a list, never None -- callers iterate it."""
    warnings = read_meta_cached(tube).warnings
    assert isinstance(warnings, (list, tuple)), (
        f"tube {tube}: .warnings must be a list or tuple, got "
        f"{type(warnings).__name__} ({warnings!r})"
    )


# ---------------------------------------------------------------------------
# tube identity
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("tube", ALL_TUBES)
def test_in_file_tube_id_matches_the_filename(tube, read_meta_cached, vsi_paths):
    """`Image_NN.vsi` must carry tube NN in its own metadata.

    A mismatch means the file was renamed or copied over another animal's data,
    which silently reassigns a whole slide to the wrong treatment group.
    """
    meta = read_meta_cached(tube)
    assert int(meta.tube_id) == tube, (
        f"{vsi_paths[tube].name} reports tube_id {meta.tube_id!r} but its "
        f"filename says {tube}"
    )


# ---------------------------------------------------------------------------
# series counts and labels
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("tube", ALL_TUBES)
def test_tissue_series_count(tube, read_meta_cached):
    """Tubes 30, 33, 34, 42, 53 and 54 have three sections; the rest have four.

    The count must come from a tolerant regex anchored on `_(0[1-4])$`.  Tube 60
    is named `60_20x_DAPI, FITC, Cy3_01`, so an equality test against
    `20x_DAPI, FITC, Cy3_01` finds zero series there.
    """
    meta = read_meta_cached(tube)
    expected = expected_series_count(tube)
    assert meta.n_tissue_series == expected, (
        f"tube {tube}: expected {expected} tissue series "
        f"(three-section tubes are {sorted(THREE_SECTION_TUBES)}), "
        f"got {meta.n_tissue_series}"
    )


@pytest.mark.parametrize("tube", ALL_TUBES)
def test_n_tissue_series_matches_len_of_series(tube, read_meta_cached):
    """The advertised count and the actual list must not diverge."""
    meta = read_meta_cached(tube)
    assert len(meta.series) == meta.n_tissue_series, (
        f"tube {tube}: .n_tissue_series is {meta.n_tissue_series} but "
        f".series holds {len(meta.series)} entries"
    )


@pytest.mark.parametrize("tube", ALL_TUBES)
def test_section_labels_are_unique_and_contiguous(tube, read_meta_cached):
    """Labels must be `01`..`0N` with no duplicates and no gaps.

    A duplicate means two stacks were mapped to the same section; a gap means one
    was dropped.  Both change the animal's denominator without any other signal.
    """
    meta = read_meta_cached(tube)
    labels = [_label(s) for s in meta.series]
    n = expected_series_count(tube)

    assert len(set(labels)) == len(labels), (
        f"tube {tube}: duplicate section labels {sorted(labels)}"
    )
    assert set(labels) == set(VALID_SECTION_LABELS[:n]), (
        f"tube {tube}: expected labels {list(VALID_SECTION_LABELS[:n])}, "
        f"got {sorted(labels)}"
    )


def test_tube_60_non_uniform_series_naming_is_tolerated(read_meta_cached):
    """Tube 60's series are named `60_20x_DAPI, FITC, Cy3_0N`, unlike every other
    animal.  This is the concrete case that breaks an equality-based name filter.
    """
    meta = read_meta_cached(60)
    assert meta.n_tissue_series == 4, (
        "tube 60 must yield 4 tissue series despite its `60_` name prefix; "
        f"got {meta.n_tissue_series}. Parse names with a tolerant regex "
        "anchored on `_(0[1-4])$`, never an equality test."
    )
    names = [str(getattr(s, "name", "")) for s in meta.series]
    assert any("60_" in name for name in names), (
        f"tube 60's series names should retain their `60_` prefix; got {names}"
    )


# ---------------------------------------------------------------------------
# stage coordinates -- validated against Bio-Formats 8.5.0
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("tube", sorted(GROUND_TRUTH_STAGE_X_UM))
def test_stage_x_matches_bio_formats(tube, read_meta_cached):
    """Per-series stage X must reproduce Bio-Formats to better than 0.5 um.

    This is the single most load-bearing number in the ingest stage: box
    membership -- and therefore positive vs negative-control status -- is derived
    from it.  Tag 2018 is VECTOR_DOUBLE_2 in micrometres.
    """
    meta = read_meta_cached(tube)
    expected = GROUND_TRUTH_STAGE_X_UM[tube]
    actual = {_label(s): float(s.stage_x_um) for s in meta.series}

    assert set(actual) == set(expected), (
        f"tube {tube}: series labels {sorted(actual)} do not match the "
        f"ground-truth labels {sorted(expected)}"
    )
    for label in sorted(expected):
        assert actual[label] == pytest.approx(
            expected[label], abs=STAGE_X_TOLERANCE_UM
        ), (
            f"tube {tube} series _{label}: stage X {actual[label]} um differs "
            f"from the Bio-Formats value {expected[label]} um by "
            f"{abs(actual[label] - expected[label]):.3f} um "
            f"(tolerance {STAGE_X_TOLERANCE_UM} um)"
        )


def test_tube_49_stage_x_is_not_monotonic_in_section_number(read_meta_cached):
    """Guards the fact the whole box-assignment design rests on.

    If a future scanner template made stage X monotonic in section number again,
    this test fires -- not because monotonicity is wrong, but because the
    non-monotonic case is the evidence that section number is acquisition order.
    Its disappearance means the ground truth needs re-checking, not that the
    pairing rule may be reinstated.
    """
    meta = read_meta_cached(49)
    by_label = {_label(s): float(s.stage_x_um) for s in meta.series}
    ordered = [by_label[label] for label in ("01", "02", "03", "04")]

    assert ordered != sorted(ordered), (
        f"tube 49 stage X in section order is {ordered}, which is now "
        "monotonic. Bio-Formats reports 46994.6 / 25600.0 / 32349.6 / 54409.6 "
        "-- re-verify the parser against the ground truth before trusting it."
    )


@pytest.mark.parametrize("tube", ALL_TUBES)
def test_stage_coordinates_are_plausible_slide_positions(tube, read_meta_cached):
    """Stage X and Y in micrometres must land on a 75 x 25 mm slide.

    A micrometre/millimetre slip shows up here as values ~1000x too small.
    """
    meta = read_meta_cached(tube)
    for series in meta.series:
        x = float(series.stage_x_um)
        y = float(series.stage_y_um)
        assert math.isfinite(x) and math.isfinite(y), (
            f"tube {tube} series _{_label(series)}: non-finite stage position "
            f"({x}, {y})"
        )
        # Every ground-truth PositionX falls in 25600-61302 um; the brute-force
        # double scan that first located tag 2018 searched 20000-70000 um and
        # found them all.  These bounds are deliberately wider than that.
        assert 1000.0 <= x <= 90000.0, (
            f"tube {tube} series _{_label(series)}: stage X {x} is not a "
            "plausible micrometre position on a 75 mm slide -- check units "
            "(a millimetre value would land near 26-61 here)"
        )
        # Stage Y is not pinned by ground truth, so this only rules out a
        # units slip or a garbage read, not a specific position.
        assert abs(y) <= 200000.0, (
            f"tube {tube} series _{_label(series)}: stage Y {y} um is far off "
            "any slide -- likely a tag-walk drift or a units error"
        )


# ---------------------------------------------------------------------------
# pixel calibration
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("tube", ALL_TUBES)
def test_pixel_size_is_0_325_um_per_series(tube, read_meta_cached):
    """Read per series, pixel size falls in 0.32500-0.32502 um.

    Asserted with a tolerance rather than bit-equality, as CLAUDE_v1.2.md
    section 5 requires -- but tightly enough that a binning change or a
    magnification mix-up (which would give ~0.65 or ~0.1625) is caught.
    """
    meta = read_meta_cached(tube)
    for series in meta.series:
        value = float(series.pixel_size_um)
        assert value == pytest.approx(
            NOMINAL_PIXEL_SIZE_UM, abs=PIXEL_SIZE_TOLERANCE_UM
        ), (
            f"tube {tube} series _{_label(series)}: pixel size {value} um "
            f"differs from the nominal {NOMINAL_PIXEL_SIZE_UM} um by "
            f"{abs(value - NOMINAL_PIXEL_SIZE_UM):.6f} um "
            f"(tolerance {PIXEL_SIZE_TOLERANCE_UM})"
        )


# ---------------------------------------------------------------------------
# image dimensions
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("tube", ALL_TUBES)
def test_series_dimensions_are_positive_integers(tube, read_meta_cached):
    """True dimensions come from tag 2053, not from the tile grid.

    The tile-grid product overestimates by up to 2.4 %, so a width taken from
    `n_tiles_x * tile_width` inflates every percent-area denominator.
    """
    meta = read_meta_cached(tube)
    for series in meta.series:
        width = series.width_px
        height = series.height_px
        assert isinstance(width, int) and isinstance(height, int), (
            f"tube {tube} series _{_label(series)}: dimensions must be ints, "
            f"got {type(width).__name__} x {type(height).__name__}"
        )
        assert width > 0 and height > 0, (
            f"tube {tube} series _{_label(series)}: non-positive dimensions "
            f"{width} x {height}"
        )


# ---------------------------------------------------------------------------
# exposure -- the assignment-free sweep
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("tube", ALL_TUBES)
def test_all_series_within_an_animal_share_one_exposure_triplet(tube, read_meta_cached):
    """Every tissue series of an animal was acquired at the same settings.

    This is the assignment-free check from EXECUTION_PLAN_v3.md section 2: dump
    every exposure triplet and require them to be identical, so that a single
    deviant series inside an otherwise-standard animal cannot hide.  Collapsing
    to a set -- the flaw in the first sweep -- cannot distinguish a missing
    channel from two channels sharing a value.
    """
    meta = read_meta_cached(tube)
    triplets = [_exposure_triplet(s) for s in meta.series]

    assert len(triplets) == expected_series_count(tube), (
        f"tube {tube}: expected {expected_series_count(tube)} exposure triplets "
        f"(3 records per tissue series), got {len(triplets)}"
    )
    reference = triplets[0]
    for label, triplet in zip((_label(s) for s in meta.series), triplets):
        assert triplet == pytest.approx(reference, abs=EXPOSURE_TOLERANCE_MS), (
            f"tube {tube} series _{label} was acquired at {triplet} ms but "
            f"series _{_label(meta.series[0])} at {reference} ms -- exposure "
            "varies within a single animal"
        )


@pytest.mark.parametrize("tube", ALL_TUBES)
def test_exposure_channels_are_ordered_dapi_fitc_cy3(tube, read_meta_cached):
    """DAPI < FITC < Cy3 holds for every animal including both deviants.

    Catches a reversed or rotated channel order, which would silently attach the
    Cy3 exposure to the GFAP channel when the two slides are exposure-corrected.
    """
    meta = read_meta_cached(tube)
    for series in meta.series:
        dapi, fitc, cy3 = _exposure_triplet(series)
        assert dapi < fitc, (
            f"tube {tube} series _{_label(series)}: DAPI exposure {dapi} ms is "
            f"not below FITC {fitc} ms -- channel order may be wrong"
        )
        assert fitc <= cy3 or tube == 51, (
            f"tube {tube} series _{_label(series)}: FITC exposure {fitc} ms "
            f"exceeds Cy3 {cy3} ms; only tube 51 has that inversion"
        )


def test_exposure_sweep_finds_exactly_two_deviant_slides(read_meta_cached):
    """Across all 31 animals, exactly slides 51 and 60 deviate from standard.

    Standard is DAPI 128.547 / FITC 397.927 / Cy3 1839.999 ms.  Slide 51 was
    acquired at 60.517 / 240.822 / 145.677 and slide 60 at Cy3 397.927 instead
    of 1839.999.  A third deviant appearing means either a parsing regression or
    a genuinely new acquisition problem -- both need a human.
    """
    deviants: dict[int, list[tuple[float, float, float]]] = {}
    for tube in ALL_TUBES:
        meta = read_meta_cached(tube)
        odd = [
            triplet
            for triplet in (_exposure_triplet(s) for s in meta.series)
            if triplet != pytest.approx(STANDARD_EXPOSURE_MS, abs=EXPOSURE_TOLERANCE_MS)
        ]
        if odd:
            deviants[tube] = odd

    assert set(deviants) == set(DEVIANT_EXPOSURE_MS), (
        f"expected exactly tubes {sorted(DEVIANT_EXPOSURE_MS)} to deviate from "
        f"the standard exposure {STANDARD_EXPOSURE_MS} ms; got "
        f"{ {t: v[0] for t, v in deviants.items()} }"
    )


@pytest.mark.parametrize("tube", sorted(DEVIANT_EXPOSURE_MS))
def test_deviant_exposure_values_are_exact(tube, read_meta_cached):
    """The two deviants' actual values must match, not merely differ from standard.

    These numbers drive the numerical exposure correction applied to slides 51
    and 60 as a cross-check against the rescans, so they must be right to the
    microsecond, in milliseconds.
    """
    meta = read_meta_cached(tube)
    expected = DEVIANT_EXPOSURE_MS[tube]
    for series in meta.series:
        triplet = _exposure_triplet(series)
        assert triplet == pytest.approx(expected, abs=EXPOSURE_TOLERANCE_MS), (
            f"tube {tube} series _{_label(series)}: exposure {triplet} ms "
            f"differs from the validated {expected} ms"
        )


@pytest.mark.parametrize("tube", sorted(set(ALL_TUBES) - set(DEVIANT_EXPOSURE_MS)))
def test_standard_exposure_values_are_exact(tube, read_meta_cached):
    """The other 29 animals must all carry the standard triplet exactly."""
    meta = read_meta_cached(tube)
    for series in meta.series:
        triplet = _exposure_triplet(series)
        assert triplet == pytest.approx(
            STANDARD_EXPOSURE_MS, abs=EXPOSURE_TOLERANCE_MS
        ), (
            f"tube {tube} series _{_label(series)}: exposure {triplet} ms is "
            f"not the standard {STANDARD_EXPOSURE_MS} ms"
        )


def test_exposure_is_milliseconds_not_microseconds(read_meta_cached):
    """A missed us->ms conversion makes every exposure 1000x too large.

    Tag 100002 stores microseconds (128547 for DAPI); `.exposure_ms` must be
    128.547.  Both are 'plausible-looking' numbers, so only a magnitude check
    catches the slip.
    """
    dapi, fitc, cy3 = _exposure_triplet(read_meta_cached(29).series[0])
    assert 0.1 < dapi < 10000.0, (
        f"tube 29 DAPI exposure {dapi} is outside 0.1-10000 ms; tag 100002 is "
        "in MICROseconds and must be divided by 1000"
    )
    assert dapi == pytest.approx(128.547, abs=EXPOSURE_TOLERANCE_MS), (
        f"tube 29 DAPI exposure should be 128.547 ms, got {dapi}"
    )
    assert cy3 == pytest.approx(1839.999, abs=EXPOSURE_TOLERANCE_MS), (
        f"tube 29 Cy3 exposure should be 1839.999 ms, got {cy3}"
    )


# ---------------------------------------------------------------------------
# the end-to-end path: metadata -> box assignment, on real coordinates
# ---------------------------------------------------------------------------
def test_every_animal_yields_a_legal_box_split(read_meta_cached):
    """Real stage coordinates must produce a legal 2+2 / 2+1 / 1+2 everywhere.

    `test_box_assignment.py` proves the splitter is correct on synthetic input;
    this proves the cohort's actual geometry is the shape the splitter expects,
    so no animal has to be resolved by hand.
    """
    from ihc.ingest.vsi_meta import BoxAssignmentError, assign_boxes

    failures: dict[int, str] = {}
    layouts: dict[int, tuple[int, int]] = {}
    for tube in ALL_TUBES:
        meta = read_meta_cached(tube)
        try:
            result = assign_boxes(meta.series)
        except BoxAssignmentError as exc:
            failures[tube] = str(exc)
            continue
        layouts[tube] = (len(result["near_label"]), len(result["far_label"]))

    assert not failures, (
        f"assign_boxes refused {len(failures)} real animal(s): {failures}"
    )
    illegal = {t: shape for t, shape in layouts.items() if shape not in {(2, 2), (2, 1), (1, 2)}}
    assert not illegal, (
        f"these animals produced a layout outside 2+2 / 2+1 / 1+2: {illegal}"
    )


def test_tube_49_real_metadata_puts_01_and_04_in_the_same_box(read_meta_cached):
    """The regression case, end to end from the real file rather than synthetics.

    `config/slides.csv` records tube 49's positives as `far_label` with the note
    "CAREFUL 01 and 04 the positives".  If this test fails, that row is being
    applied to sections _02 and _03 -- the negative controls.
    """
    from ihc.ingest.vsi_meta import assign_boxes

    meta = read_meta_cached(49)
    result = assign_boxes(meta.series)
    far = {_box_label(s) for s in result["far_label"]}
    near = {_box_label(s) for s in result["near_label"]}

    assert far == {"01", "04"}, (
        f"tube 49 far_label box must be {{'01','04'}}, got {far}"
    )
    assert near == {"02", "03"}, (
        f"tube 49 near_label box must be {{'02','03'}}, got {near}"
    )
