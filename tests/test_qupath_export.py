"""Tests for the blinded QuPath project export.

The point of this module is blinding, so most of these tests are adversarial: they
assert that identifiers *cannot* reach the artefact rather than that the happy path
works. The happy path is the easy half.

Nothing here needs QuPath, Java or Bio-Formats. The Groovy import script itself has
never been executed — see ADR-0015 and ADR-0017 in ``docs/decisions.md`` for what that
leaves unverified.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from ihc.ingest.qupath_export import (
    BlindingLeakError,
    ProjectSpecError,
    audit_codes,
    build_project_spec,
    scan_for_tube_identifiers,
    write_project_spec,
)

# --------------------------------------------------------------------------- #
# a synthetic raw tree: no real pixels needed, only the file/folder shape
# --------------------------------------------------------------------------- #

CHANNELS_CFG = {
    "expected_order": ["DAPI", "FITC", "Cy3"],
    "channels": [
        {"name": "DAPI", "target": "nuclei"},
        {"name": "FITC", "target": "GFAP"},
        {"name": "Cy3", "target": "Abeta"},
    ],
}


@pytest.fixture
def raw_root(tmp_path, monkeypatch):
    """A raw tree shaped like the real one, with the metadata reader stubbed out.

    ``read_vsi_meta`` and ``read_vsi_index`` are exercised against the real cohort by
    ``test_vsi_meta.py``; re-deriving them here would test the fixture, not the code.
    What this module cares about is what the export does with their output.
    """
    root = tmp_path / "RawData"
    (root / "Rescan").mkdir(parents=True)

    def make(directory: Path, tube: int, sections: list[str], payload: bool = True):
        (directory / f"Image_{tube}.vsi").write_bytes(b"II*\x00" + b"\x00" * 64)
        if payload:
            companion = directory / f"_Image_{tube}_"
            for stack in [1, 10000] + [10002 + 3 * i for i in range(len(sections))]:
                (companion / f"stack{stack}").mkdir(parents=True)

    make(root, 29, ["01", "02", "03", "04"])
    make(root, 30, ["01", "02", "03"])
    make(root, 37, ["01", "02", "03", "04"])
    make(root, 44, ["01", "02", "03", "04"], payload=False)  # index-only animal
    make(root / "Rescan", 51, ["01", "02"])
    make(root, 51, ["01", "02", "03", "04"])

    sizes = {"01": (17882, 19104), "02": (18025, 19677),
             "03": (17882, 19104), "04": (17098, 19344)}

    class _Series:
        def __init__(self, label, name):
            self.section_label = label
            self.name = name
            self.width_px, self.height_px = sizes[label]
            self.pixel_size_um = 0.325005

    class _Meta:
        def __init__(self, labels, prefix):
            self.series = [_Series(s, f"{prefix}20x_DAPI, FITC, Cy3_{s}") for s in labels]
            # A realistic reader warning: it quotes the internal series name, which
            # for tube 60 embeds the tube number. This is the leak _redact() exists for.
            self.warnings = [f"series {self.series[0].name!r}: no pixel size (tag 2019)"]

    def fake_read_vsi_meta(path):
        stem = Path(path).stem
        tube = int(stem.split("_")[1])
        labels = ["01", "02"] if "Rescan" in str(path) else (
            ["01", "02", "03"] if tube == 30 else ["01", "02", "03", "04"])
        prefix = f"{tube}_" if tube == 60 else ""
        return _Meta(labels, prefix)

    def fake_read_vsi_index(path):
        meta = fake_read_vsi_meta(path)
        out = [{"kind": "label", "section": None}, {"kind": "overview", "section": None},
               {"kind": "other", "section": None}]
        for s in meta.series:
            out.append({"kind": "tissue", "section": s.section_label})
            out.append({"kind": "other", "section": None})
        return out

    monkeypatch.setattr("ihc.ingest.qupath_export.read_vsi_meta", fake_read_vsi_meta)
    monkeypatch.setattr("ihc.ingest.qupath_export.read_vsi_index", fake_read_vsi_index)
    return root


CODES = {29: "K07", 30: "K19", 37: "K02", 44: "K31", 51: "K11"}


def rows(*specs):
    """``rows(("K07", "01", "positive"), ...)`` -> manifest rows."""
    out = []
    for spec in specs:
        code, section, condition = spec[:3]
        row = {"coded_id": code, "section_label": section, "condition": condition}
        if len(spec) > 3:
            row.update(spec[3])
        out.append(row)
    return out


def build(tmp_path, raw_root, manifest, **kwargs):
    kwargs.setdefault("channels_cfg", CHANNELS_CFG)
    return build_project_spec(manifest, CODES, raw_root, tmp_path / "out", **kwargs)


# --------------------------------------------------------------------------- #
# the happy path
# --------------------------------------------------------------------------- #


def test_builds_one_entry_per_section(tmp_path, raw_root):
    spec = build(tmp_path, raw_root,
                 rows(("K07", "01", "positive"), ("K07", "03", "negative")))
    assert [e["image_name"] for e in spec["images"]] == ["K07_s01", "K07_s03"]
    assert spec["counts"] == {"images": 2, "positive": 1, "negative": 1, "animals": 1,
                              "skipped": 0, "excluded": 0}


def test_negative_controls_are_included_by_default(tmp_path, raw_root):
    """They are needed: the negative-control QC gate is assessed per region."""
    spec = build(tmp_path, raw_root, rows(("K07", "02", "negative")))
    assert spec["counts"]["negative"] == 1


def test_negative_controls_can_be_left_out_explicitly(tmp_path, raw_root):
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive"), ("K07", "02", "negative")),
                 include_conditions=("positive",))
    assert spec["counts"]["images"] == 1
    assert spec["excluded"][0]["reason"] == "condition_negative_not_requested"


def test_image_type_and_channel_visibility(tmp_path, raw_root):
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive")))
    assert spec["image_type"] == "FLUORESCENCE"
    assert spec["channels"]["visible"] == ["DAPI"]
    assert spec["channels"]["hidden"] == ["FITC", "Cy3"]
    assert spec["images"][0]["channel_visible"] == [True, False, False]


def test_series_index_skips_label_and_overview(tmp_path, raw_root):
    """stack1 (label) is 0 and stack10000 (overview) is 1, so tissue starts at 2."""
    spec = build(tmp_path, raw_root,
                 rows(("K07", "01", "positive"), ("K07", "04", "positive")))
    by_name = {e["image_name"]: e for e in spec["images"]}
    assert by_name["K07_s01"]["series_index"] == 2
    assert by_name["K07_s04"]["series_index"] == 5
    assert by_name["K07_s01"]["series_match_suffix"] == "_01"


def test_dimensions_are_carried_for_the_import_time_cross_check(tmp_path, raw_root):
    spec = build(tmp_path, raw_root, rows(("K07", "02", "positive")))
    assert (spec["images"][0]["width_px"], spec["images"][0]["height_px"]) == (18025, 19677)


def test_entries_are_sorted_by_code_not_by_manifest_order(tmp_path, raw_root):
    """Manifest order is tube order, which is treatment-group order."""
    spec = build(tmp_path, raw_root,
                 rows(("K19", "01", "positive"), ("K07", "01", "positive"),
                      ("K11", "01", "positive")))
    assert [e["code"] for e in spec["images"]] == ["K07", "K11", "K19"]


def test_rescan_gets_a_distinct_name_and_alias(tmp_path, raw_root):
    spec = build(tmp_path, raw_root,
                 rows(("K11", "01", "positive"),
                      ("K11", "01", "positive", {"scan": "rescan"})))
    assert sorted(e["image_name"] for e in spec["images"]) == ["K11_s01", "K11b_s01"]
    assert len({e["vsi_path"] for e in spec["images"]}) == 2


# --------------------------------------------------------------------------- #
# what must be left out
# --------------------------------------------------------------------------- #


def test_needs_confirmation_rows_are_excluded_not_guessed(tmp_path, raw_root):
    """Tube 37 is the live case: box column and bench annotation disagree."""
    spec = build(tmp_path, raw_root,
                 rows(("K02", "01", "positive", {"needs_confirmation": "CONTRADICTION"}),
                      ("K07", "01", "positive")))
    assert [e["image_name"] for e in spec["images"]] == ["K07_s01"]
    # Check the fields that carry meaning, not the whole dict: an exact-equality
    # assertion breaks whenever a new diagnostic field is added, which punishes
    # improving the error report.
    assert len(spec["excluded"]) == 1
    got = spec["excluded"][0]
    assert got["code"] == "K02"
    assert got["section_label"] == "01"
    assert got["reason"] == "needs_confirmation"
    assert "not known" in got["detail"]


def test_missing_payload_is_skipped_not_an_error(tmp_path, raw_root):
    """23 of 31 animals are index-only. Normal state, not a failure."""
    spec = build(tmp_path, raw_root,
                 rows(("K31", "01", "positive"), ("K07", "01", "positive")))
    assert spec["counts"]["images"] == 1
    assert spec["skipped"][0]["reason"] == "payload_absent"


def test_section_absent_from_the_slide_is_skipped(tmp_path, raw_root):
    spec = build(tmp_path, raw_root, rows(("K19", "04", "positive")))
    assert spec["images"] == []
    assert spec["skipped"][0]["reason"] == "section_not_in_file"


def test_unresolved_condition_is_excluded_and_reported_not_guessed(tmp_path, raw_root):
    """manifest.py's third condition state. Spec §2: exclude it and report it."""
    spec = build(tmp_path, raw_root,
                 rows(("K07", "01", "unresolved"), ("K07", "02", "positive")))
    assert [e["image_name"] for e in spec["images"]] == ["K07_s02"]
    assert spec["excluded"][0]["reason"] == "condition_unresolved"


def test_blank_condition_is_treated_as_unresolved(tmp_path, raw_root):
    spec = build(tmp_path, raw_root, rows(("K07", "01", "")))
    assert spec["images"] == []
    assert spec["excluded"][0]["reason"] == "condition_unresolved"


def test_an_uninterpretable_condition_is_refused(tmp_path, raw_root):
    with pytest.raises(ProjectSpecError, match="not one of positive"):
        build(tmp_path, raw_root, rows(("K07", "01", "maybe")))


def test_code_missing_from_the_key_is_refused(tmp_path, raw_root):
    with pytest.raises(ProjectSpecError, match="not in the code mapping"):
        build(tmp_path, raw_root, rows(("K99", "01", "positive")))


def test_colliding_rows_are_excluded_together_never_merged(tmp_path, raw_root):
    """The tubes 51 / 60 case, in miniature.

    Two rows sharing (code, section_label) but disagreeing about condition are not a
    duplicate to be de-duplicated: on the real data they are the original and the
    rescan, which cover *different PAP-pen boxes*, so one is a datum and the other is
    a negative control. Picking either would delineate a section the record describes
    as the other thing. Both go, and the collision is reported.
    """
    spec = build(tmp_path, raw_root,
                 rows(("K07", "01", "positive"), ("K07", "01", "negative"),
                      ("K07", "02", "positive")))
    assert [e["image_name"] for e in spec["images"]] == ["K07_s02"]
    collision = spec["excluded"][0]
    assert collision["reason"] == "ambiguous_duplicate_rows"
    assert collision["conditions_seen"] == ["negative", "positive"]
    assert "DISAGREE" in collision["detail"]
    assert any("K07_s01" in w for w in spec["warnings"])


# --------------------------------------------------------------------------- #
# blinding
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("text", [
    "/RawData/Image_29.vsi", "Image 60", "image-42.vsi", "1007344 - 29",
    "_Image_49_", "60_20x_DAPI, FITC, Cy3_01", "tube 55",
])
def test_identifier_scanner_flags_real_leak_vectors(text):
    assert scan_for_tube_identifiers(text)


@pytest.mark.parametrize("text", [
    "K42_s01", "K42b_s02", "K02.vsi", "_K16b_", "20x_DAPI, FITC, Cy3_01",
    "DAPI", "FITC", "Cy3", "_04", "delineation", "0.325", "K29_s01",
])
def test_identifier_scanner_leaves_coded_names_alone(text):
    assert scan_for_tube_identifiers(text) == []


def test_loose_scan_ignores_ordinary_numbers():
    assert scan_for_tube_identifiers('"images": 34,', strict=False) == []
    assert scan_for_tube_identifiers("2026-08-07T09:34:12+00:00", strict=False) == []


def test_no_string_in_the_spec_carries_an_identifier(tmp_path, raw_root):
    spec = build(tmp_path, raw_root,
                 rows(("K07", "01", "positive"), ("K11", "02", "negative"),
                      ("K31", "01", "positive")))

    def walk(node):
        if isinstance(node, str):
            yield node
        elif isinstance(node, dict):
            for k, v in node.items():
                yield from walk(v)
        elif isinstance(node, list):
            for v in node:
                yield from walk(v)

    # The tmp directory pytest hands out is named `pytest-33`, `pytest-42`, ... so it
    # trips the bare-number rule on its own. That prefix is test scaffolding, not part
    # of the artefact, so it is removed before scanning rather than exempted by key.
    prefix = str(tmp_path)
    for key, value in spec.items():
        if key == "link_plan":
            continue
        for text in walk(value):
            stripped = text.replace(prefix, "<tmp>")
            assert scan_for_tube_identifiers(stripped) == [], f"{key}: {text!r}"


def test_reader_warnings_are_redacted(tmp_path, raw_root):
    """A metadata warning quotes the internal series name, which for tube 60 has the
    tube number in it. The warning must survive; the number must not."""
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive")))
    assert spec["warnings"], "the stub reader emits a warning; it should be carried"
    assert all(scan_for_tube_identifiers(w) == [] for w in spec["warnings"])
    assert any("no pixel size" in w for w in spec["warnings"])


def test_project_name_with_an_identifier_is_refused(tmp_path, raw_root):
    with pytest.raises(BlindingLeakError):
        build(tmp_path, raw_root, rows(("K07", "01", "positive")),
              project_name="delineation_Image_29")


def test_uncoded_paths_need_an_explicit_acknowledgement(tmp_path, raw_root):
    with pytest.raises(ProjectSpecError, match="acknowledge_uncoded_paths"):
        build(tmp_path, raw_root, rows(("K07", "01", "positive")), alias_images=False)


def test_uncoded_paths_are_allowed_once_acknowledged(tmp_path, raw_root):
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive")),
                 alias_images=False, acknowledge_uncoded_paths="PPS, engineering pilot")
    assert spec["blinding"]["paths_are_coded"] is False
    # The path now genuinely contains the tube number; that is the whole point of
    # requiring the acknowledgement, and the Groovy script reports every one of them.
    assert scan_for_tube_identifiers(spec["images"][0]["vsi_path"])


# --------------------------------------------------------------------------- #
# code auditing
# --------------------------------------------------------------------------- #


def test_sequential_codes_are_refused():
    codes = {tube: f"K{i:02d}" for i, tube in enumerate(range(29, 45), start=1)}
    with pytest.raises(ProjectSpecError, match="monotonic"):
        audit_codes(codes)


def test_reversed_codes_are_refused():
    codes = {tube: f"K{i:02d}" for i, tube in enumerate(range(44, 28, -1), start=1)}
    with pytest.raises(ProjectSpecError, match="preserves it|monotonic"):
        audit_codes(codes)


def test_arithmetic_transform_of_the_tube_id_is_refused():
    codes = {tube: f"K{tube + 100}" for tube in range(29, 45)}
    with pytest.raises(ProjectSpecError, match="monotonic"):
        audit_codes(codes)


def test_a_real_permutation_passes():
    import random

    tubes = list(range(29, 59)) + [60]
    labels = [f"K{i:02d}" for i in range(1, len(tubes) + 1)]
    random.Random(20260807).shuffle(labels)
    audit_codes(dict(zip(tubes, labels)))  # must not raise


def test_duplicate_codes_are_refused():
    with pytest.raises(ProjectSpecError, match="duplicate"):
        audit_codes({29: "K01", 30: "K01"})


def test_inverted_mapping_is_diagnosed():
    with pytest.raises(ProjectSpecError, match="other way round"):
        audit_codes({"K01": "29", "K02": "30"})


def test_code_with_an_identifier_in_it_is_refused():
    with pytest.raises(BlindingLeakError):
        audit_codes({29: "Image29x", 30: "K02"})


# --------------------------------------------------------------------------- #
# writing
# --------------------------------------------------------------------------- #


def test_write_creates_coded_symlinks_that_resolve(tmp_path, raw_root):
    out = tmp_path / "out"
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive")))
    path = write_project_spec(spec, out)

    assert path == out / "project_spec.json"
    link = out / "images" / "K07.vsi"
    payload = out / "images" / "_K07_"
    assert link.is_symlink() and link.resolve() == (raw_root / "Image_29.vsi")
    assert payload.is_symlink() and (payload / "stack10002").is_dir()


def test_written_json_does_not_contain_the_link_targets(tmp_path, raw_root):
    out = tmp_path / "out"
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive")))
    text = write_project_spec(spec, out).read_text()
    assert "Image_29" not in text
    assert json.loads(text).get("link_plan") is None
    assert json.loads(text)["blinding"]["link_plan_redacted"] is True
    assert scan_for_tube_identifiers(text, strict=False) == []


def test_write_is_idempotent(tmp_path, raw_root):
    out = tmp_path / "out"
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive")))
    write_project_spec(spec, out)
    write_project_spec(spec, out)  # must not raise


def test_write_refuses_to_repoint_a_reused_code(tmp_path, raw_root):
    out = tmp_path / "out"
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive")))
    write_project_spec(spec, out)
    spec["link_plan"][0]["target_vsi"] = str(raw_root / "Image_30.vsi")
    with pytest.raises(ProjectSpecError, match="already points at"):
        write_project_spec(spec, out)


def test_write_refuses_a_directory_the_spec_does_not_describe(tmp_path, raw_root):
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive")))
    with pytest.raises(ProjectSpecError, match="describes"):
        write_project_spec(spec, tmp_path / "somewhere_else")


def test_written_spec_is_valid_json_with_the_expected_shape(tmp_path, raw_root):
    out = tmp_path / "out"
    spec = build(tmp_path, raw_root, rows(("K07", "01", "positive"), ("K07", "03", "negative")))
    loaded = json.loads(write_project_spec(spec, out).read_text())
    assert loaded["schema"] == "ihc.ingest.qupath_export/1"
    assert loaded["project_dir"] == str(out / "qupath")
    assert Path(loaded["project_dir"]).is_dir()
    for entry in loaded["images"]:
        assert set(entry) >= {"image_name", "vsi_path", "series_index",
                              "series_match_suffix", "channel_names", "channel_visible"}
        assert os.path.isabs(entry["vsi_path"])


# =========================================================================== #
# APPENDIX: adversarial checks against the REAL cohort
# --------------------------------------------------------------------------- #
# Everything above runs on a synthetic raw tree with the metadata readers
# stubbed out, which is the right way to test the export's logic.  This section
# does the opposite: it runs the real chain -- 31 real .vsi index files ->
# `build_manifest` -> `split_manifest` -> `build_project_spec` -- and then walks
# whatever comes out with generic scanners that do not know the spec's shape.
#
# Two failure modes that stubs cannot reach:
#
#   * anything that depends on the real metadata.  Tube 60's internal series
#     names really do read "60_20x_DAPI, FITC, Cy3_01", and the real payload
#     folders really do order their stacks 1, 10000, 10002, 10005, 10008, 10011.
#     A stub asserts the layout it was given.
#   * a NEW key added to the spec later that nobody remembers to test.  The
#     walker below scans every leaf of every nesting level, so a leak in a field
#     invented next month is caught without anyone editing this file.
#
# `link_plan` is deliberately exempt from the annotator-facing scans: it is the
# custodian's materialisation plan and it MUST carry the real targets, which is
# its entire purpose.  What matters is that it never reaches the annotator --
# `test_written_json_does_not_contain_the_link_targets` above is the guard for
# that, and `test_real_spec_confines_tube_ids_to_the_link_plan` below is the
# complement from this direction.
# =========================================================================== #

import re as _re  # noqa: E402

from _stage3_helpers import (  # noqa: E402
    ALL_TUBES,
    GROUP_TOKENS,
    NEEDS_CONFIRMATION_TUBES,
    require_cohort_index_files,
    resolve_column,
    text_forbidden_hits,
    text_tube_id_hits,
)

#: Stack IDs inside an `_Image_NN_` payload folder.  Non-contiguous on purpose:
#: the gaps hold the focus map and the sample mask.
TISSUE_STACKS = (10002, 10005, 10008, 10011)
LABEL_STACK = 1
OVERVIEW_STACK = 10000

#: Subtrees that exist for the custodian and never reach the annotator.
CUSTODIAN_SECTIONS = ("link_plan",)

#: Key names allowed to carry an original path even outside those subtrees.
SOURCE_KEY_TOKENS = ("path", "uri", "url", "file", "source", "location", "dir",
                     "link", "target", "payload")

#: Key names whose string values QuPath shows to the annotator.
NAME_KEY_TOKENS = ("name", "title", "label", "display", "caption", "entry", "code")


def walk_spec(node, path=()):
    """Yield ``(key_path, value)`` for every leaf in a nested dict/list."""
    if isinstance(node, dict):
        for key, value in node.items():
            yield from walk_spec(value, path + (str(key),))
    elif isinstance(node, (list, tuple)):
        for index, value in enumerate(node):
            yield from walk_spec(value, path + (f"[{index}]",))
    else:
        yield path, node


def _key_has(key_path, tokens):
    joined = " ".join(key_path).lower()
    return any(token in joined for token in tokens)


def _is_custodian(key_path):
    return bool(key_path) and key_path[0] in CUSTODIAN_SECTIONS


def annotator_strings(spec):
    """Every string leaf the annotator could conceivably see."""
    return [
        (p, v) for p, v in walk_spec(spec)
        if isinstance(v, str) and not _is_custodian(p)
    ]


def _real_blinded(data_root, *, include_rescans):
    manifest_mod = pytest.importorskip("ihc.ingest.manifest")
    blinding_mod = pytest.importorskip("ihc.ingest.blinding")
    pd = pytest.importorskip("pandas")

    require_cohort_index_files(data_root)

    df = manifest_mod.build_manifest(raw_root=data_root, include_rescans=include_rescans)
    tube_column = resolve_column(df, "tube_id")
    usable = df[~df[tube_column].astype("int64").isin(sorted(NEEDS_CONFIRMATION_TUBES))]
    codes = blinding_mod.generate_codes(
        sorted(set(int(t) for t in usable[tube_column].unique())),
        seed=20260807, allow_weak_seed=True,   # readable fixture seed; see conftest
    )
    _, blinded = blinding_mod.split_manifest(usable, codes)
    assert isinstance(blinded, pd.DataFrame)
    return blinded, codes


@pytest.fixture(scope="module")
def real_chain(data_root):
    """(blinded, codes) for the real cohort, built by the real modules.

    Rescans are excluded here so the leak tests below keep doing their job while
    the defect `test_real_chain_survives_rescans` reports stays open.  That
    defect deserves one clearly-named failure, not nine unrelated blinding
    checks erroring out inside a fixture.
    """
    return _real_blinded(data_root, include_rescans=False)


@pytest.mark.requires_data
def test_real_chain_survives_rescans(data_root, tmp_path):
    """Rescan rows must stay distinguishable all the way to the project.

    Tubes 51 and 60 were re-acquired, so the manifest holds two rows carrying
    the same section label for each of them -- one per scan, with a DIFFERENT
    box and a different condition, because the rescan covers the positive box
    only.  Something has to tell those rows apart downstream.

    If the discriminator is stripped during blinding, two things break at once:
    the blinded manifest shows one section twice with contradictory conditions,
    and the project export cannot name the entries at all.  Rescan status is
    technical metadata, which the blinded manifest is explicitly allowed to
    carry (CLAUDE_v1.2 section 6 step 3) -- on its own it identifies no animal.
    """
    blinded, codes = _real_blinded(data_root, include_rescans=True)
    key = [c for c in ("code", "section_label") if c in blinded.columns]
    assert len(key) == 2, f"blinded manifest lacks code/section_label: {list(blinded.columns)}"

    discriminators = [c for c in ("scan", "scan_kind", "is_rescan") if c in blinded.columns]
    duplicated = blinded.duplicated(subset=key, keep=False)
    assert not duplicated.any() or discriminators, (
        f"{int(duplicated.sum())} blinded rows share a (code, section_label) with "
        f"another row and there is no scan column to tell them apart. Columns "
        f"present: {list(blinded.columns)}\n"
        f"{blinded[duplicated].to_string()}\n"
        "These are the original and the rescan of the same slide -- note they "
        "carry OPPOSITE conditions, so a reader of this manifest cannot tell "
        "which one is the negative control."
    )
    build_project_spec(blinded, codes, data_root, tmp_path / "with_rescans")


@pytest.fixture(scope="module")
def real_spec(real_chain, data_root, tmp_path_factory):
    blinded, codes = real_chain
    out_dir = tmp_path_factory.mktemp("real_qupath")
    spec = build_project_spec(blinded, codes, data_root, out_dir)
    assert isinstance(spec, dict) and spec
    return spec


@pytest.mark.requires_data
def test_real_spec_opens_only_tissue_stacks(real_spec, real_chain, data_root):
    """`series_index` resolved against the ACTUAL payload folder layout.

    The stacks inside `_Image_NN_` sort as 1 (label), 10000 (overview), then
    10002 / 10005 / 10008 / 10011 (tissue).  So a series index of 0 or 1 opens
    the slide label or the whole-slide overview -- the label image carries the
    tube number as printed text *and* as a DataMatrix barcode.

    This walks the index back through the real directory listing rather than
    trusting a number, which is the only way to know the off-by-one is right.
    """
    _, codes = real_chain
    inverse = {str(c): t for t, c in codes.items()}
    entries = real_spec.get("images") or []
    assert entries, "the real spec contains no image entries"

    checked = 0
    for entry in entries:
        tube = inverse[str(entry["code"])]
        payload = Path(data_root) / f"_Image_{tube}_"
        if not payload.is_dir():
            continue
        stacks = sorted(
            int(p.name[5:]) for p in payload.iterdir()
            if p.is_dir() and p.name.startswith("stack") and p.name[5:].isdigit()
        )
        index = entry["series_index"]
        resolved = stacks[index] if 0 <= index < len(stacks) else index
        assert resolved in TISSUE_STACKS, (
            f"{entry['image_name']}: series_index {index} resolves to stack "
            f"{resolved} in {payload.name} (label={LABEL_STACK}, "
            f"overview={OVERVIEW_STACK}, tissue={list(TISSUE_STACKS)})"
        )
        checked += 1
    assert checked, "no entry could be resolved against a real payload folder"


@pytest.mark.requires_data
def test_real_spec_series_suffix_matches_the_section_it_claims(real_spec):
    """The name-suffix fallback must agree with the section it is filed under.

    The import script resolves by name suffix rather than by index, so a suffix
    that disagrees with `section_label` opens a different section of the same
    animal -- a wrong coronal level, silently, with no dimension mismatch to
    give it away.
    """
    for entry in real_spec.get("images") or []:
        assert entry["series_match_suffix"] == f"_{entry['section_label']}", (
            f"{entry['image_name']}: suffix {entry['series_match_suffix']!r} does "
            f"not match section {entry['section_label']!r}"
        )
        assert _re.fullmatch(r"_0[1-4]", entry["series_match_suffix"]), (
            f"{entry['image_name']}: {entry['series_match_suffix']!r} is not a "
            "tissue-section suffix"
        )


@pytest.mark.requires_data
def test_real_spec_never_names_the_label_or_overview_series(real_spec):
    """No image entry may reference the label, overview, mask or focus series."""
    for entry in real_spec.get("images") or []:
        blob = " ".join(str(v) for v in entry.values()).lower()
        for token in ("label", "overview", "macro", "thumbnail", "sample mask",
                      "focusmap", "focuspoints"):
            assert token not in blob, (
                f"{entry.get('image_name')} references the {token!r} series"
            )
        assert not _re.search(r"\bstack1\b", blob)
        assert "stack10000" not in blob


@pytest.mark.requires_data
def test_real_spec_display_names_carry_no_tube_id(real_spec):
    """The tube 60 regression, on the real file.

    Tube 60's internal series names read `"60_20x_DAPI, FITC, Cy3_01"` and
    QuPath displays internal names in preference to filenames.  Renaming the
    file achieves nothing; the entry has to be renamed.
    """
    offenders = []
    for key_path, value in annotator_strings(real_spec):
        if _key_has(key_path, SOURCE_KEY_TOKENS):
            continue
        if not _key_has(key_path, NAME_KEY_TOKENS):
            continue
        hits = text_tube_id_hits(value)
        if hits:
            offenders.append((".".join(key_path), value, hits))
    assert not offenders, f"display names carrying a tube ID: {offenders}"


@pytest.mark.requires_data
def test_real_spec_reuses_no_raw_series_name_or_filename(real_spec):
    offenders = []
    for key_path, value in annotator_strings(real_spec):
        if _key_has(key_path, SOURCE_KEY_TOKENS):
            continue
        lowered = value.lower()
        for token in ("60_", "20x_dapi", "image_"):
            if token in lowered:
                offenders.append((".".join(key_path), value, token))
    assert not offenders, f"spec values reuse a raw series name or filename: {offenders}"


#: An ISO-8601 date-time.  A *build* timestamp is wall-clock noise whose digits
#: land on a tube number roughly one minute in three ("...T08:20:34" contains
#: 34), so scanning it for tube IDs produces a test that fails on the clock.
#: An *acquisition* timestamp is a different animal entirely -- it varies per
#: mouse and is a perfect group proxy -- so it is checked separately, by
#: `test_real_spec_carries_no_per_entry_timestamp`, rather than by digits.
_ISO_DT_RE = _re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}")


@pytest.mark.requires_data
def test_real_spec_confines_tube_ids_to_the_link_plan(real_spec):
    """A tube ID may live only in the custodian's materialisation plan.

    Anywhere else -- a note, a warning, a count, a metadata field -- it is a
    leak.  Reader warnings are the sharp case: they quote the internal series
    name, which for tube 60 embeds the tube number.
    """
    offenders = []
    for key_path, value in walk_spec(real_spec):
        if _is_custodian(key_path) or _key_has(key_path, SOURCE_KEY_TOKENS):
            continue
        if isinstance(value, str):
            if _ISO_DT_RE.search(value):
                continue
            hits = text_tube_id_hits(value)
            if hits:
                offenders.append((".".join(key_path), value[:160], hits))
        elif isinstance(value, int) and not isinstance(value, bool):
            if value in ALL_TUBES and not _key_has(
                key_path, ("series", "stack", "index", "width", "height", "count", "n_")
            ):
                offenders.append((".".join(key_path), value, ["exact"]))
    assert not offenders, f"tube IDs outside the link plan: {offenders}"


@pytest.mark.requires_data
def test_real_spec_carries_no_per_entry_timestamp(real_spec):
    """The complement to the exclusion above, and the leak that matters.

    A single build timestamp is harmless -- it is the same for every animal.  A
    timestamp attached to an image ENTRY is not: the cohort was scanned in
    ascending tube order, which is ascending group order, so per-entry times
    sort the animals straight back into their treatment blocks.
    """
    offenders = []
    for entry in real_spec.get("images") or []:
        for key, value in entry.items():
            if isinstance(value, str) and _ISO_DT_RE.search(value):
                offenders.append((entry.get("image_name"), key, value))
    assert not offenders, (
        f"image entries carry an acquisition-style timestamp: {offenders}"
    )


@pytest.mark.requires_data
def test_real_spec_names_no_treatment_group(real_spec):
    """No arm or group name outside the custodian subtree.

    Note the raw tree lives under a directory named after the drug, so the study
    path contains "Rapamycin" for every animal alike.  That is a constant, not a
    discriminator, and it is confined to `link_plan`; a group name attached to an
    individual entry would be the real thing.
    """
    offenders = []
    for key_path, value in annotator_strings(real_spec):
        if _key_has(key_path, SOURCE_KEY_TOKENS):
            continue
        hits = text_forbidden_hits(value, GROUP_TOKENS)
        if hits:
            offenders.append((".".join(key_path), value[:160], hits))
    assert not offenders, f"the spec names a treatment group: {offenders}"


@pytest.mark.requires_data
def test_real_spec_addresses_every_entry_by_a_coded_id(real_spec, real_chain):
    """The spec must USE the codes, not merely avoid the tube IDs.

    Entries named `entry_0 .. entry_119` leak no tube ID and are still fatal:
    the ordinal is order-preserving, so it hands the tube order -- and therefore
    the treatment blocks -- straight back.
    """
    _, codes = real_chain
    issued = {str(c) for c in codes.values()}
    entries = real_spec.get("images") or []
    assert entries
    for entry in entries:
        assert str(entry["code"]) in issued, (
            f"entry {entry.get('image_name')!r} is not addressed by an issued code"
        )
        assert str(entry["code"]) in str(entry["image_name"]), (
            f"image_name {entry['image_name']!r} does not contain its code"
        )


@pytest.mark.requires_data
def test_real_spec_entry_count_never_exceeds_the_sections_available(real_spec):
    """One entry per section, never one per series in the file.

    Six series exist per `.vsi`; only three or four are tissue.  An entry count
    that tracks the series count is the signature of the label and the overview
    having been imported alongside the sections.
    """
    entries = real_spec.get("images") or []
    seen = set()
    for entry in entries:
        key = (entry["code"], entry["section_label"], entry.get("scan", "original"))
        assert key not in seen, f"duplicate entry for {key}"
        seen.add(key)
    per_code = {}
    for entry in entries:
        per_code[entry["code"]] = per_code.get(entry["code"], 0) + 1
    assert per_code, "no entries"
    assert max(per_code.values()) <= 6, (
        f"an animal has {max(per_code.values())} entries; four sections plus a "
        f"two-section rescan is the maximum: {per_code}"
    )
