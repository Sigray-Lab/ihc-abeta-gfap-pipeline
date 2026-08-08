"""Corruption tests for `ihc.ingest.verify` -- the ingest gate.

WHAT CLASS OF BUG THIS PROTECTS AGAINST
---------------------------------------
A `.vsi` file is an *index*: about 1.7 MB against a real dataset of roughly
1.5 GB.  Series names, dimensions and channel counts all come from the index
alone, so a dataset whose pixels are missing, truncated, duplicated or bound to
the wrong animal still *describes itself perfectly*.  CLAUDE_v1.2.md section 5
records the consequence measured with Bio-Formats: with one `.ets` stack
deleted, all six named series were still reported and conversion exited 0 with
no warning.  Two transfers have already lost payloads silently.

The gate exists to convert those silent losses into loud failures.  So the tests
here are almost all *negative* tests: build a specific, realistic corruption and
assert `verify_dataset` refuses it.  A gate that passes everything is worse than
no gate, because it manufactures confidence.

Corruptions covered:
    * a removed `.ets` stack (directory removed, and file removed with the
      directory left in place)
    * a Dropbox conflicted-copy file sitting beside the real tile file
    * a truncated tile table (last 10 kB chopped; the chunk table sits at the
      very end of every `.ets` in this dataset, so this is a direct hit)
    * the wrong-companion trap: a `.vsi` whose stem is a prefix of another
      animal's payload folder
    * plus a positive control -- an intact dataset must pass, or every negative
      result above is vacuous.

Datasets are assembled in `tmp_path` with the `.vsi` copied and every `.ets`
symlinked, so a corruption costs kilobytes rather than 1.5 GB, and the real
`RawData/` tree is never written to.
"""

from __future__ import annotations

import hashlib
import os
import shutil
from pathlib import Path

import pytest

from ihc.ingest.verify import verify_dataset, verify_directory

pytestmark = [pytest.mark.requires_data, pytest.mark.requires_payload]

TRUNCATE_BYTES = 10 * 1024


# ---------------------------------------------------------------------------
# building synthetic datasets
# ---------------------------------------------------------------------------
def clone_dataset(
    vsi: Path,
    companion: Path,
    dest: Path,
    *,
    vsi_name: str | None = None,
    companion_name: str | None = None,
    omit_stacks: tuple[str, ...] = (),
    omit_files: tuple[str, ...] = (),
    materialise_stacks: tuple[str, ...] = (),
) -> tuple[Path, Path]:
    """Reproduce a real dataset under `dest`, copying the index and symlinking pixels.

    `omit_stacks`     -- stack directories to leave out entirely
    `omit_files`      -- 'stackNNNNN/filename' entries to leave out
    `materialise_stacks` -- stacks whose files are real copies, so they can be
                            mutated (truncated) without touching RawData/
    """
    dest.mkdir(parents=True, exist_ok=True)
    out_vsi = dest / (vsi_name or vsi.name)
    shutil.copy2(vsi, out_vsi)

    out_companion = dest / (companion_name or companion.name)
    out_companion.mkdir(parents=True, exist_ok=True)

    for stack in sorted(p for p in companion.iterdir() if p.is_dir()):
        if stack.name in omit_stacks:
            continue
        stack_dest = out_companion / stack.name
        stack_dest.mkdir(parents=True, exist_ok=True)
        for src in sorted(stack.iterdir()):
            if src.name.startswith(".") or not src.is_file():
                continue
            if f"{stack.name}/{src.name}" in omit_files:
                continue
            target = stack_dest / src.name
            if stack.name in materialise_stacks:
                shutil.copy2(src, target)
            else:
                target.symlink_to(src.resolve())
    return out_vsi, out_companion


def smallest_tile_file(companion: Path) -> Path:
    """The smallest `frame_t*.ets` in a dataset -- cheapest real ETS to corrupt."""
    candidates = [
        p
        for p in companion.rglob("*.ets")
        if p.name.startswith("frame_t") and p.is_file()
    ]
    assert candidates, f"no frame_t*.ets found under {companion}"
    return min(candidates, key=lambda p: p.stat().st_size)


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# tolerant readers for the documented result shape
# ---------------------------------------------------------------------------
DOCUMENTED_KEYS = ("ok", "vsi_sha256", "companion", "stacks", "problems", "warnings")


def check_shape(result) -> dict:
    assert isinstance(result, dict), (
        f"verify_dataset must return a dict, got {type(result).__name__}"
    )
    missing = [key for key in DOCUMENTED_KEYS if key not in result]
    assert not missing, (
        f"verify_dataset result is missing documented key(s) {missing}; "
        f"got {sorted(result)}"
    )
    return result


def assert_failed(result, scenario: str):
    check_shape(result)
    assert result["ok"] is False, (
        f"the gate PASSED {scenario}. This is exactly the silent failure the "
        f"gate exists to stop. Result: {_brief(result)}"
    )
    assert result["problems"], (
        f"the gate reported ok=False for {scenario} but left `problems` empty, "
        "so nothing tells the operator what is wrong"
    )


def _brief(result: dict) -> str:
    return ", ".join(
        f"{k}={result.get(k)!r}" for k in ("ok", "companion", "problems", "warnings")
    )


def _walk_strings(node):
    """Every string anywhere in a nested result structure."""
    if isinstance(node, str):
        yield node
    elif isinstance(node, dict):
        for key, value in node.items():
            yield str(key)
            yield from _walk_strings(value)
    elif isinstance(node, (list, tuple, set)):
        for item in node:
            yield from _walk_strings(item)
    else:
        yield str(node)


def mentions(result, needle: str) -> bool:
    return any(needle in text for text in _walk_strings(result))


def _walk_hex_digests(node):
    """Every 64-character lowercase hex string in a nested result structure."""
    for text in _walk_strings(node):
        candidate = text.strip()
        if len(candidate) == 64:
            try:
                int(candidate, 16)
            except ValueError:
                continue
            yield candidate.lower()


def stack_count(result) -> int:
    stacks = result.get("stacks")
    if stacks is None:
        return 0
    if isinstance(stacks, dict):
        return len(stacks)
    if isinstance(stacks, (list, tuple, set)):
        return len(stacks)
    pytest.fail(
        f"`stacks` must be a list or a mapping keyed by stack name, got "
        f"{type(stacks).__name__}: {stacks!r}"
    )


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def intact(tmp_path, payload_dataset):
    """A complete, uncorrupted clone of a real dataset."""
    vsi, companion = payload_dataset
    return clone_dataset(vsi, companion, tmp_path / "intact")


@pytest.fixture(scope="session")
def real_stack_names(payload_dataset) -> list[str]:
    _, companion = payload_dataset
    return sorted(p.name for p in companion.iterdir() if p.is_dir())


@pytest.fixture(scope="session")
def tissue_stack(payload_dataset) -> str:
    """Name of one fluorescence tissue stack, e.g. 'stack10002'.

    Selected by looking for `frame_t_0.ets` rather than by name: the label
    (stack1) and overview (stack10000) stacks use `frame_t.ets`, and stack10000
    also carries `blob_21_f_Frame#0.ets` and `blob_21_f.meta`.  Losing a tissue
    stack is the documented silent failure, so that is what these tests break.
    """
    _, companion = payload_dataset
    for path in sorted(p for p in companion.iterdir() if p.is_dir()):
        if (path / "frame_t_0.ets").is_file():
            return path.name
    pytest.skip(f"no fluorescence tissue stack (frame_t_0.ets) under {companion}")


# ---------------------------------------------------------------------------
# positive control -- without this, every negative test below is vacuous
# ---------------------------------------------------------------------------
def test_intact_dataset_passes(intact):
    """An untouched dataset must pass cleanly, with no problems reported."""
    vsi_path, companion = intact
    result = check_shape(verify_dataset(vsi_path, hash_ets=False))

    assert result["ok"] is True, (
        "an intact dataset was rejected; the gate is unusable if it cannot "
        f"pass good data. Problems: {result['problems']!r}"
    )
    assert not result["problems"], (
        f"intact dataset reported problems: {result['problems']!r}"
    )


def test_intact_dataset_binds_to_its_own_companion(intact):
    """`companion` must point at the folder actually holding the pixels."""
    vsi_path, companion = intact
    result = verify_dataset(vsi_path, hash_ets=False)

    reported = result["companion"]
    assert reported is not None, "verify_dataset reported companion=None for a good dataset"
    assert Path(str(reported)).resolve() == companion.resolve(), (
        f"companion resolved to {reported!r}, expected {companion}"
    )


def test_intact_dataset_inventories_every_stack(intact, real_stack_names):
    """Every `stackNNNNN` folder present on disk must appear in the report.

    CLAUDE_v1.2.md section 5: "An explicit stack-inventory assertion is
    mandatory."  A stack the gate never looked at is a stack whose loss it
    cannot detect.
    """
    vsi_path, _ = intact
    result = verify_dataset(vsi_path, hash_ets=False)

    unlisted = [name for name in real_stack_names if not mentions(result, name)]
    assert not unlisted, (
        f"these stack folders exist under the companion but are absent from "
        f"the verify report: {unlisted}. Report was: {result['stacks']!r}"
    )
    assert stack_count(result) >= 3, (
        f"the inventory holds {stack_count(result)} entries; a real dataset has "
        f"{len(real_stack_names)} stacks ({real_stack_names})"
    )


# ---------------------------------------------------------------------------
# corruption 1 -- a removed .ets stack
# ---------------------------------------------------------------------------
# CLAUDE_v1.2.md section 5: "a missing `.ets` stack fails silently: with one
# stack removed, Bio-Formats still reported all six named series and conversion
# exited 0 with no warning, because names and dimensions come from the `.vsi`
# index alone. An explicit stack-inventory assertion is mandatory."
# ---------------------------------------------------------------------------

def test_removed_stack_directory_fails(tmp_path, payload_dataset, tissue_stack):
    """Deleting a whole tissue stack folder must fail the gate."""
    vsi, companion = payload_dataset
    vsi_path, _ = clone_dataset(
        vsi, companion, tmp_path / "no_stack_dir", omit_stacks=(tissue_stack,)
    )

    result = verify_dataset(vsi_path, hash_ets=False)
    assert_failed(result, f"a dataset with the whole {tissue_stack}/ folder deleted")
    # Not asserted: that the report names `stack10002` specifically.  The index
    # records how many tissue series exist but not which stackNNNNN folder each
    # one lives in, so the gate can legitimately only report a count mismatch.
    # It must still point the operator at the stack inventory rather than at
    # some unrelated symptom.
    assert any("stack" in str(p).lower() for p in result["problems"]), (
        "the gate failed but no problem mentions the stack inventory, so the "
        "operator is not told that a stack is missing. Problems: "
        f"{result['problems']!r}"
    )


def test_removed_ets_file_with_directory_intact_fails(
    tmp_path, payload_dataset, tissue_stack
):
    """The subtler loss: the stack folder survives but its `.ets` is gone.

    A gate that only counts directories passes this, and Bio-Formats reports the
    series regardless because the name lives in the index.
    """
    vsi, companion = payload_dataset
    tile_names = [
        p.name for p in (companion / tissue_stack).iterdir() if p.suffix == ".ets"
    ]
    assert tile_names, f"{tissue_stack} holds no .ets file to remove"

    vsi_path, out_companion = clone_dataset(
        vsi,
        companion,
        tmp_path / "no_ets_file",
        omit_files=tuple(f"{tissue_stack}/{name}" for name in tile_names),
    )
    assert (out_companion / tissue_stack).is_dir(), "the empty stack folder must remain"

    result = verify_dataset(vsi_path, hash_ets=False)
    assert_failed(
        result, f"a dataset whose {tissue_stack}/ folder is present but empty"
    )


def test_stack_count_drops_when_a_stack_is_removed(
    tmp_path, payload_dataset, intact, tissue_stack
):
    """The inventory must actually shrink, not merely flip an `ok` flag."""
    vsi, companion = payload_dataset
    victim = tissue_stack
    vsi_path, _ = clone_dataset(
        vsi, companion, tmp_path / "count_check", omit_stacks=(victim,)
    )

    good = verify_dataset(intact[0], hash_ets=False)
    bad = verify_dataset(vsi_path, hash_ets=False)
    assert stack_count(bad) == stack_count(good) - 1, (
        f"removing {victim} should drop the inventory from {stack_count(good)} "
        f"to {stack_count(good) - 1}, got {stack_count(bad)}"
    )


# ---------------------------------------------------------------------------
# corruption 2 -- a Dropbox conflicted copy
# ---------------------------------------------------------------------------
CONFLICTED_NAME = "frame_t_0 (user's conflicted copy 2026-07-29).ets"


def test_dropbox_conflicted_copy_is_not_counted_as_a_series(
    tmp_path, payload_dataset, intact, tissue_stack
):
    """A conflicted copy beside the real tile file must not inflate the inventory.

    Dropbox writes these when two machines edit the same path.  The file is a
    valid ETS container, so any inventory built by globbing `*.ets` will count it
    as an extra series or an extra tile source -- and the animal silently gains a
    section that does not exist.
    """
    vsi, companion = payload_dataset
    vsi_path, out_companion = clone_dataset(vsi, companion, tmp_path / "conflicted")
    real_tile = companion / tissue_stack / "frame_t_0.ets"
    (out_companion / tissue_stack / CONFLICTED_NAME).symlink_to(real_tile.resolve())

    result = verify_dataset(vsi_path, hash_ets=False)
    baseline = verify_dataset(intact[0], hash_ets=False)

    assert stack_count(result) == stack_count(baseline), (
        f"the conflicted copy inflated the stack inventory from "
        f"{stack_count(baseline)} to {stack_count(result)}. A conflicted copy "
        "is a duplicate of an existing tile file, never a new series."
    )


def test_dropbox_conflicted_copy_is_reported(tmp_path, payload_dataset, tissue_stack):
    """It must also be surfaced, because it means the transfer was disturbed.

    Whether the gate treats it as a hard problem or a warning is a design call,
    but staying silent is not: the operator has no other way to learn that two
    divergent copies of this dataset existed.
    """
    vsi, companion = payload_dataset
    vsi_path, out_companion = clone_dataset(
        vsi, companion, tmp_path / "conflicted_report"
    )
    real_tile = companion / tissue_stack / "frame_t_0.ets"
    (out_companion / tissue_stack / CONFLICTED_NAME).symlink_to(real_tile.resolve())

    result = check_shape(verify_dataset(vsi_path, hash_ets=False))
    surfaced = mentions(result["problems"], "conflicted") or mentions(
        result["warnings"], "conflicted"
    )
    assert surfaced, (
        "a Dropbox conflicted copy sat beside the real tile file and the gate "
        f"said nothing. problems={result['problems']!r} "
        f"warnings={result['warnings']!r}"
    )


# ---------------------------------------------------------------------------
# corruption 3 -- a truncated tile table
# ---------------------------------------------------------------------------
def test_truncated_ets_fails(tmp_path, payload_dataset):
    """Chopping the last 10 kB off an `.ets` must fail the gate.

    In every `.ets` in this dataset the chunk table sits at the very end of the
    file (`chunk_offset + n_chunks * (20 + 4*n_dim) == file_size` exactly), so
    removing 10 kB removes most of the table: the header still advertises
    `n_chunks` entries and the table now runs past EOF.  A reader that trusts the
    header and seeks blindly produces garbage tiles instead of an error.
    """
    vsi, companion = payload_dataset
    victim_file = smallest_tile_file(companion)
    victim_stack = victim_file.parent.name

    vsi_path, out_companion = clone_dataset(
        vsi,
        companion,
        tmp_path / "truncated",
        materialise_stacks=(victim_stack,),
    )
    target = out_companion / victim_stack / victim_file.name
    original_size = target.stat().st_size
    assert original_size > TRUNCATE_BYTES * 4, (
        f"{target} is only {original_size} bytes; too small for this test"
    )
    assert not target.is_symlink(), (
        "the truncation target must be a real copy, never a symlink into RawData/"
    )
    os.truncate(target, original_size - TRUNCATE_BYTES)

    assert target.stat().st_size == original_size - TRUNCATE_BYTES
    result = verify_dataset(vsi_path, hash_ets=False)
    assert_failed(
        result,
        f"a dataset whose {victim_stack}/{victim_file.name} lost its last "
        f"{TRUNCATE_BYTES} bytes (the tile table)",
    )


def test_truncation_never_touches_the_real_data(payload_dataset, tmp_path):
    """Paranoia check: RawData/ is read-only and must stay byte-identical.

    A corruption test that mutates a symlink target would destroy irreplaceable
    raw data.  This asserts the clone's materialised stack is a real file whose
    parent tree is the tmp_path, not RawData/.
    """
    vsi, companion = payload_dataset
    victim_file = smallest_tile_file(companion)
    before = victim_file.stat().st_size

    _, out_companion = clone_dataset(
        vsi,
        companion,
        tmp_path / "safety",
        materialise_stacks=(victim_file.parent.name,),
    )
    clone = out_companion / victim_file.parent.name / victim_file.name
    os.truncate(clone, before - TRUNCATE_BYTES)

    assert victim_file.stat().st_size == before, (
        f"the real file {victim_file} changed size -- the clone helper is "
        "mutating RawData/, which is read-only"
    )


# ---------------------------------------------------------------------------
# corruption 4 -- the wrong-companion trap
# ---------------------------------------------------------------------------
# The companion folder for `Image_49.vsi` is `_Image_49_`.  A resolver written as
# a prefix match or a glob (`_Image_4*`) binds `Image_4.vsi` to `_Image_49_` and
# silently reports another animal's pixels under this animal's identity -- which
# would put one animal's data into another's treatment group.
# ---------------------------------------------------------------------------

def test_vsi_does_not_bind_to_a_prefix_matching_companion(tmp_path, payload_dataset):
    """`Image_4.vsi` beside `_Image_49_` must NOT claim `_Image_49_`'s pixels."""
    vsi, companion = payload_dataset
    trap = tmp_path / "prefix_trap"
    clone_dataset(vsi, companion, trap)  # the legitimate dataset, e.g. _Image_49_

    stem_prefix = vsi.stem[:-1]  # "Image_49" -> "Image_4"
    decoy = trap / f"{stem_prefix}.vsi"
    shutil.copy2(vsi, decoy)
    assert not (trap / f"_{stem_prefix}_").exists(), (
        "the decoy must have no companion folder of its own"
    )

    result = check_shape(verify_dataset(decoy, hash_ets=False))

    reported = result["companion"]
    if reported is not None:
        assert Path(str(reported)).resolve() != (trap / companion.name).resolve(), (
            f"{decoy.name} bound to {reported!r}, which is {companion.name}'s "
            f"payload. `_{stem_prefix}_` does not exist; the resolver is "
            "prefix-matching instead of requiring an exact companion name."
        )
    assert result["ok"] is False, (
        f"{decoy.name} has no companion folder at all, yet the gate passed it: "
        f"{_brief(result)}"
    )


def test_each_vsi_binds_to_its_own_companion_when_both_exist(tmp_path, payload_dataset):
    """With `_Image_4_` and `_Image_49_` both present, each must bind correctly.

    This is the positive half of the trap: exact-name resolution has to keep
    working when a genuine prefix-sharing neighbour is in the same directory.
    """
    vsi, companion = payload_dataset
    root = tmp_path / "both"
    long_vsi, long_companion = clone_dataset(vsi, companion, root)

    stem_prefix = vsi.stem[:-1]  # "Image_49" -> "Image_4"
    short_vsi, short_companion = clone_dataset(
        vsi,
        companion,
        root,
        vsi_name=f"{stem_prefix}.vsi",
        companion_name=f"_{stem_prefix}_",
    )

    short_result = check_shape(verify_dataset(short_vsi, hash_ets=False))
    long_result = check_shape(verify_dataset(long_vsi, hash_ets=False))

    assert Path(str(short_result["companion"])).resolve() == short_companion.resolve(), (
        f"{short_vsi.name} bound to {short_result['companion']!r}, "
        f"expected {short_companion}"
    )
    assert Path(str(long_result["companion"])).resolve() == long_companion.resolve(), (
        f"{long_vsi.name} bound to {long_result['companion']!r}, "
        f"expected {long_companion}"
    )
    # `ok` is deliberately not asserted here: the decoy carries the real
    # animal's in-file tube ID under a different filename, so a gate that
    # cross-checks the two is entitled to flag it.  Companion binding is the
    # property under test.


def test_missing_companion_folder_fails(tmp_path, payload_dataset):
    """The original silent failure: a `.vsi` copied without its payload folder.

    This has already happened twice in transfer.  1.7 MB arrives, 1.5 GB does
    not, and nothing complains.
    """
    vsi, _ = payload_dataset
    lonely = tmp_path / "lonely"
    lonely.mkdir()
    shutil.copy2(vsi, lonely / vsi.name)

    result = verify_dataset(lonely / vsi.name, hash_ets=False)
    assert_failed(result, "a .vsi copied without its companion payload folder")


# ---------------------------------------------------------------------------
# hashing
# ---------------------------------------------------------------------------
def test_vsi_hash_is_stable_and_matches_hashlib(intact):
    """The index hash must be reproducible and actually be sha256 of the file."""
    vsi_path, _ = intact
    first = verify_dataset(vsi_path, hash_ets=False)["vsi_sha256"]
    second = verify_dataset(vsi_path, hash_ets=False)["vsi_sha256"]

    assert first == second, (
        f"vsi_sha256 changed between two runs of the same file: {first} vs {second}"
    )
    assert first == sha256_of(vsi_path), (
        f"vsi_sha256 is {first} but sha256 of the file is {sha256_of(vsi_path)}"
    )


def test_vsi_hash_differs_between_animals(tmp_path, payload_datasets):
    """Two different animals must not hash the same -- a constant hash is useless."""
    if len(payload_datasets) < 2:
        pytest.skip("needs two payload datasets to compare")
    tubes = sorted(payload_datasets)[:2]
    hashes = {}
    for tube in tubes:
        vsi, companion = payload_datasets[tube]
        vsi_path, _ = clone_dataset(vsi, companion, tmp_path / f"tube{tube}")
        hashes[tube] = verify_dataset(vsi_path, hash_ets=False)["vsi_sha256"]
    assert hashes[tubes[0]] != hashes[tubes[1]], (
        f"tubes {tubes[0]} and {tubes[1]} produced the same vsi_sha256 "
        f"{hashes[tubes[0]]}"
    )


@pytest.mark.slow
def test_ets_hashing_is_stable_across_runs(intact, payload_dataset):
    """With hash_ets=True the same bytes must give the same digests every run.

    Reads the full ~1.5 GB payload twice, hence `slow`.  Deselect with
    `-m 'not slow'`.
    """
    vsi_path, _ = intact
    first = verify_dataset(vsi_path, hash_ets=True)
    second = verify_dataset(vsi_path, hash_ets=True)

    digests_first = sorted(_walk_hex_digests(first))
    digests_second = sorted(_walk_hex_digests(second))

    assert digests_first, (
        "hash_ets=True produced no sha256 digests at all. CLAUDE_v1.2.md "
        "section 6 requires hashing the .vsi AND every companion .ets."
    )
    assert digests_first == digests_second, (
        "hashes changed between two runs of the same unmodified dataset:\n"
        f"  run 1: {digests_first}\n  run 2: {digests_second}"
    )

    _, companion = payload_dataset
    known = sha256_of(smallest_tile_file(companion))
    assert known in digests_first, (
        f"the sha256 of {smallest_tile_file(companion).name} ({known}) does not "
        "appear in the result; companion .ets files are not being hashed"
    )


@pytest.mark.slow
def test_ets_hashing_detects_a_changed_byte(tmp_path, payload_dataset):
    """Flipping one byte in a tile file must change the reported hashes.

    A hash computed over the filename or the header alone would not move.
    """
    vsi, companion = payload_dataset
    victim_file = smallest_tile_file(companion)
    victim_stack = victim_file.parent.name

    vsi_path, out_companion = clone_dataset(
        vsi, companion, tmp_path / "flip", materialise_stacks=(victim_stack,)
    )
    before = sorted(_walk_hex_digests(verify_dataset(vsi_path, hash_ets=True)))

    target = out_companion / victim_stack / victim_file.name
    with open(target, "r+b") as handle:
        handle.seek(target.stat().st_size // 2)
        original = handle.read(1)
        handle.seek(target.stat().st_size // 2)
        handle.write(bytes([original[0] ^ 0xFF]))

    after = sorted(_walk_hex_digests(verify_dataset(vsi_path, hash_ets=True)))
    assert before != after, (
        "flipping one byte inside a tile file did not change any reported "
        "hash -- the .ets files are not being hashed over their full contents"
    )


# ---------------------------------------------------------------------------
# sparsity reporting
# ---------------------------------------------------------------------------
def test_sparsity_fraction_is_reported(intact):
    """5-11 % of tile positions inside the bounding box were never acquired.

    Never-acquired tiles are MISSING SUPPORT, not background, and must be
    excluded from the percent-area denominator.  The gate is where that number
    is established, so it has to appear in the report; downstream code cannot
    recover it once the tiles are mosaicked.
    """
    vsi_path, _ = intact
    result = verify_dataset(vsi_path, hash_ets=False)

    found: list[float] = []

    def scan(node):
        if isinstance(node, dict):
            for key, value in node.items():
                if "spars" in str(key).lower() and isinstance(value, (int, float)):
                    found.append(float(value))
                scan(value)
        elif isinstance(node, (list, tuple)):
            for item in node:
                scan(item)

    scan(result)
    assert found, (
        "no sparsity figure anywhere in the verify_dataset result. The tile "
        "grid is 5-11 % sparse and the verifier must report that fraction; "
        f"result keys were {sorted(result)}"
    )
    for value in found:
        assert 0.0 <= value <= 1.0, (
            f"sparsity {value} is not a fraction in [0, 1] -- report a "
            "fraction, not a percentage, or rename the key"
        )


# ---------------------------------------------------------------------------
# verify_directory
# ---------------------------------------------------------------------------
def test_verify_directory_passes_on_a_clean_directory(tmp_path, payload_dataset):
    """A directory holding one intact dataset must pass."""
    vsi, companion = payload_dataset
    root = tmp_path / "clean"
    clone_dataset(vsi, companion, root)

    result = verify_directory(root, hash_ets=False)
    assert isinstance(result, dict), (
        f"verify_directory must return a dict, got {type(result).__name__}"
    )
    assert result.get("ok") is True, (
        f"a directory of intact datasets was rejected: {result!r}"
    )


def test_verify_directory_fails_when_one_dataset_is_broken(
    tmp_path, payload_dataset, tissue_stack
):
    """One bad dataset among good ones must fail the whole directory.

    A per-directory `ok` that ignores an individual failure is how a broken
    animal reaches the manifest.
    """
    vsi, companion = payload_dataset
    root = tmp_path / "mixed"
    clone_dataset(vsi, companion, root)

    stem_alt = f"{vsi.stem}b"
    clone_dataset(
        vsi,
        companion,
        root,
        vsi_name=f"{stem_alt}.vsi",
        companion_name=f"_{stem_alt}_",
        omit_stacks=(tissue_stack,),
    )

    result = verify_directory(root, hash_ets=False)
    assert result.get("ok") is False, (
        f"a directory containing a dataset missing {tissue_stack} was passed: "
        f"{result!r}"
    )
    assert mentions(result, stem_alt), (
        f"the broken dataset {stem_alt} is not named anywhere in the "
        f"directory report: {result!r}"
    )


def test_verify_directory_reports_every_dataset(tmp_path, payload_datasets):
    """No dataset may be skipped -- a skipped one looks like a passing one."""
    if len(payload_datasets) < 2:
        pytest.skip("needs two payload datasets to compare")
    root = tmp_path / "many"
    tubes = sorted(payload_datasets)[:2]
    for tube in tubes:
        vsi, companion = payload_datasets[tube]
        clone_dataset(vsi, companion, root)

    result = verify_directory(root, hash_ets=False)
    for tube in tubes:
        assert mentions(result, f"Image_{tube}"), (
            f"Image_{tube}.vsi is present under {root} but does not appear in "
            f"the verify_directory report: {result!r}"
        )
