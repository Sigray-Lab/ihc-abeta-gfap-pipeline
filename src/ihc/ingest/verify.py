#!/usr/bin/env python3
"""Ingest gate for Olympus/Evident VSI whole-slide datasets.

Why this module exists
----------------------
A ``.vsi`` file is only an *index*: a couple of megabytes of metadata. The
pixels live in a sibling folder called ``_<name>_`` holding one ``stackNNNNN``
subfolder per image series, each containing an ``.ets`` tile container. Copying
or syncing the ``.vsi`` on its own gives you a file that opens, reports sensible
dimensions, and contains no image data at all. The failure is silent until
somebody tries to measure something.

This module is the gate that has to catch that, and a few other ways the data
can be quietly wrong, *before* any analysis runs. It needs no Java, no
Bio-Formats and no QuPath — only the Python standard library.

What it checks
--------------
1. **Storage is real.** Dropbox "online-only" placeholder files report a
   plausible size from ``stat()`` but stall or fail when read. Detected up
   front from allocated blocks, and failed fast with an actionable message.
2. **The companion folder belongs to *this* slide.** Matching is on the exact
   file stem. (The previous implementation used substring matching, so
   ``Image_5.vsi`` could bind to ``_Image_51_``'s pixels and still report PASS.)
3. **Only genuine tile files count as image series.** Exactly ``frame_t.ets``
   or ``frame_t_0.ets``. A Dropbox conflicted copy such as
   ``frame_t_0 (the PI's conflicted copy 2026-07-29).ets`` is a *failure*,
   not a bonus series.
4. **Content integrity.** SHA-256 of the ``.vsi`` and of every ``.ets``, so
   interior corruption cannot pass silently, plus a check that every tile
   offset lies inside its file (catches truncated transfers).
5. **Stack inventory.** Stack IDs must be a subset of the expected set, and the
   number of tissue stacks must match the expected section count.
6. **True image dimensions.** Reported from VSI tag 2053, which is
   authoritative. The tile-grid product is reported too, clearly labelled, but
   it is *not* an image size: it usually exceeds the true size (up to +5 % by
   area for tissue, +8.4 % for the label) because the last tile row and column
   are partial, and it can also fall *short* of it (-7.0 % for tube 30's first
   section) when whole tile columns at an edge hold no tissue and were never
   acquired.
7. **Sparse-tile accounting.** The scanner acquires only tiles that intersect
   its sample mask, so 4-18 % of tile positions inside the bounding box hold no
   data. That is normal and is *not* an error — but it is reported, because
   those regions are **missing support, not background**, and any denominator
   used downstream must exclude them.

Public API
----------
``verify_dataset(vsi_path, *, hash_ets=True)``
    Verify one ``.vsi`` plus its companion payload. Returns a result dict.
``verify_directory(directory, *, hash_ets=True)``
    Verify every ``.vsi`` directly inside a folder. Returns a summary dict.
``main(argv)``
    Command-line entry point. ``argv`` excludes the program name (pass
    ``sys.argv[1:]``). Returns a process exit code: 0 = all passed.

Command line
------------
::

    python -m ihc.ingest.verify RawData
    python -m ihc.ingest.verify RawData RawData/Rescan --expect-sections 3,4
    python -m ihc.ingest.verify RawData --no-hash --json QC/verify.json

Exit code is 0 only if every dataset passed, so this can be used as a gate.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import glob
import hashlib
import json
import math
import os
import re
import struct
import sys
from typing import Any, Sequence

__all__ = [
    "verify_dataset",
    "verify_directory",
    "format_dataset_report",
    "format_directory_report",
    "main",
    "DehydratedFileError",
    "EtsFormatError",
]

TOOL_VERSION = "1.0.0"
SCHEMA = "ihc.ingest.verify/1"

# --------------------------------------------------------------------------
# Expectations for this cohort. Everything tunable lives here so a reader can
# see in one place what "normal" is assumed to be.
# --------------------------------------------------------------------------

#: Stack folder IDs Olympus writes for these slides. 1 = slide label,
#: 10000 = slide overview, the rest are the tissue sections in acquisition
#: order.
EXPECTED_STACK_IDS: frozenset[int] = frozenset({1, 10000, 10002, 10005, 10008, 10011})
LABEL_STACK_ID = 1
OVERVIEW_STACK_ID = 10000

#: The only file names that count as a tile container. Strict allow-list:
#: anything else in a stack folder is either a known auxiliary file (below) or
#: a problem.
ALLOWED_ETS_NAMES: frozenset[str] = frozenset({"frame_t.ets", "frame_t_0.ets"})

#: Auxiliary files Olympus writes into stack10000 (the scanner's sample mask).
#: Present in every real dataset; not image series.
AUX_FILE_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^blob_\d+_f\.meta$"),
    re.compile(r"^blob_\d+_f_Frame#\d+\.ets$"),
)

#: Filesystem clutter that is never data and never worth mentioning.
IGNORED_FILE_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^\.DS_Store$"),
    re.compile(r"^\._"),
    re.compile(r"^Thumbs\.db$", re.IGNORECASE),
    re.compile(r"^Icon\r?$"),
    re.compile(r"^desktop\.ini$", re.IGNORECASE),
)

#: Marks a Dropbox / OneDrive sync conflict. Ambiguous pixel data: hard fail.
CONFLICT_PATTERN = re.compile(r"conflicted copy|\(\d+\)\.ets$|-\s*Conflict", re.IGNORECASE)

#: Tissue series are named like "20x_DAPI, FITC, Cy3_01". Naming is NOT
#: uniform across the cohort (tube 60 prefixes the tube number), so this is
#: anchored only on the trailing section number — never an equality test.
TISSUE_SERIES_RE = re.compile(r"_(0[1-9])$")

#: Number of tissue sections expected per slide in the main cohort. Rescan
#: slides hold only the re-acquired sections, so they need an explicit
#: ``--expect-sections 2``.
DEFAULT_EXPECTED_SECTION_COUNTS: tuple[int, ...] = (3, 4)

#: Objective pixel size, micrometres. Measured 0.32500-0.32502 across the
#: cohort. Compared with a tolerance, never for equality.
EXPECTED_PIXEL_SIZE_UM = 0.325
PIXEL_SIZE_TOL_UM = 0.002

#: Sparse-tile fraction that is considered normal. Measured 3.6-17.9 % across
#: the eight animals with payloads plus both rescans. Outside this band we warn
#: (it may mean an unusual sample mask), but it is never an error.
SPARSITY_WARN_LO = 0.005
SPARSITY_WARN_HI = 0.30

#: A file whose allocated blocks fall below this fraction of its apparent size
#: is treated as a cloud placeholder rather than real data.
DEHYDRATION_RATIO = 0.10
#: Below this size, block-allocation granularity makes the ratio meaningless.
DEHYDRATION_MIN_BYTES = 1 << 20

_HASH_CHUNK_BYTES = 8 << 20

PIXEL_TYPES: dict[int, tuple[str, int]] = {
    1: ("int8", 1),
    2: ("uint8", 1),
    3: ("int16", 2),
    4: ("uint16", 2),
    5: ("int32", 4),
    6: ("uint32", 4),
    9: ("float32", 4),
    10: ("float64", 8),
}

#: Compression codes. Note that 3 does NOT imply lossless: this cohort is lossy
#: JPEG2000 (9/7 irreversible wavelet, 2 quality layers, quality 98).
COMPRESSION_NAMES: dict[int, str] = {
    0: "uncompressed",
    2: "JPEG",
    3: "JPEG2000",
    5: "raw/lossless",
}

_ETS_CHUNK_FIXED_BYTES = 20  # 4 pad + <Q offset + <i length + 4 pad


class DehydratedFileError(Exception):
    """Raised when a file is a cloud placeholder rather than real data."""


class EtsFormatError(Exception):
    """Raised when an ``.ets`` container cannot be parsed or is truncated."""


# --------------------------------------------------------------------------
# Storage checks
# --------------------------------------------------------------------------


def allocated_bytes(path: str) -> int | None:
    """Return the bytes actually allocated on disk for *path*.

    Returns ``None`` on platforms where ``st_blocks`` is unavailable, in which
    case dehydration cannot be detected and is not claimed.
    """
    st = os.stat(path)
    blocks = getattr(st, "st_blocks", None)
    if blocks is None:
        return None
    return int(blocks) * 512


def check_materialised(path: str) -> None:
    """Fail fast if *path* is a cloud placeholder (Dropbox "online-only").

    A dehydrated file reports a plausible ``st_size`` but has almost no blocks
    allocated. Reading it either stalls for minutes while Dropbox downloads it
    or fails outright, so we refuse to start rather than die halfway through.

    Raises:
        DehydratedFileError: with an instruction the user can act on.
    """
    size = os.path.getsize(path)
    if size < DEHYDRATION_MIN_BYTES:
        return
    alloc = allocated_bytes(path)
    if alloc is None:
        return
    if alloc < DEHYDRATION_RATIO * size:
        raise DehydratedFileError(
            f"{path} looks like a cloud placeholder: it claims {size:,} bytes but only "
            f"{alloc:,} bytes are allocated on disk. Reading it would stall or fail. "
            f"In Finder, right-click the enclosing folder and choose "
            f"'Make Available Offline' (Dropbox) / 'Always keep on this device' "
            f"(OneDrive), wait for the download to finish, then re-run."
        )


def sha256_file(path: str, chunk_bytes: int = _HASH_CHUNK_BYTES) -> str:
    """Return the SHA-256 hex digest of *path*, streamed in chunks.

    Streaming matters: a single tissue series is ~350 MB and a whole animal
    ~1.4 GB, which must not be loaded into memory.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk_bytes)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------
# ETS tile container
# --------------------------------------------------------------------------


def read_ets_index(path: str) -> dict[str, Any]:
    """Parse the header and chunk table of an ``.ets`` tile container.

    Only the index is read — no pixels are decoded — so this is fast even for a
    350 MB file.

    Layout (verified against this cohort):
      * 64-byte header: magic ``SIS\\0`` at 0; ``<II`` (version, n_dim) at 8;
        ``<Q`` extra_offset at 16; ``<Q`` chunk_offset at 32; ``<I`` n_chunks
        at 40.
      * At ``extra_offset``: 40 bytes, magic ``ETS\\0``, then ``<8i`` at +8 =
        (pixel_type, size_c, colorspace, compression, quality,
        tile_x, tile_y, tile_z).
      * Chunk table at ``chunk_offset``, record size ``20 + 4 * n_dim``:
        4 pad bytes, ``n_dim`` int32 coordinates, ``<Q`` file offset,
        ``<i`` byte length, 4 pad bytes.
      * Coordinates: ``[0]`` = tile column, ``[1]`` = tile row,
        ``[2]`` = channel (when ``n_dim >= 4``), ``[-1]`` = pyramid level.

    Raises:
        EtsFormatError: on bad magic, an implausible chunk table, or a chunk
            table that runs past the end of the file (truncated transfer).
    """
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        head = fh.read(64)
        if len(head) < 64 or head[0:4] != b"SIS\x00":
            raise EtsFormatError("not an ETS/SIS tile container (bad magic)")
        version, n_dim = struct.unpack_from("<II", head, 8)
        (extra_offset,) = struct.unpack_from("<Q", head, 16)
        (chunk_offset,) = struct.unpack_from("<Q", head, 32)
        (n_chunks,) = struct.unpack_from("<I", head, 40)

        if not 2 <= n_dim <= 8:
            raise EtsFormatError(f"implausible dimension count ({n_dim})")
        if extra_offset + 40 > size:
            raise EtsFormatError("ETS header block lies past the end of the file")

        fh.seek(extra_offset)
        block = fh.read(40)
        if block[0:4] != b"ETS\x00":
            raise EtsFormatError("ETS header block missing")
        (
            pixel_type,
            size_c,
            colorspace,
            compression,
            quality,
            tile_x,
            tile_y,
            tile_z,
        ) = struct.unpack_from("<8i", block, 8)

        record = _ETS_CHUNK_FIXED_BYTES + 4 * n_dim
        table_end = chunk_offset + n_chunks * record
        if table_end > size:
            raise EtsFormatError(
                f"chunk table runs past the end of the file (needs {table_end:,} bytes, "
                f"file is {size:,}) - the file is truncated"
            )
        fh.seek(chunk_offset)
        table = fh.read(n_chunks * record)

    chunks: list[tuple[tuple[int, ...], int, int]] = []
    for i in range(n_chunks):
        p = i * record
        coords = struct.unpack_from("<" + "i" * n_dim, table, p + 4)
        offset, nbytes = struct.unpack_from("<Qi", table, p + 4 + 4 * n_dim)
        chunks.append((coords, offset, nbytes))

    return {
        "path": path,
        "file_bytes": size,
        "version": version,
        "n_dim": n_dim,
        "n_chunks": n_chunks,
        "pixel_type": pixel_type,
        "samples_per_tile": size_c,
        "colorspace": colorspace,
        "compression": compression,
        "quality": quality,
        "tile": (tile_x, tile_y, tile_z),
        "chunks": chunks,
    }


def summarise_ets(path: str, *, hash_ets: bool = True) -> dict[str, Any]:
    """Describe one ``.ets`` container: geometry, sparsity, integrity, hash.

    Args:
        path: the tile container to inspect.
        hash_ets: whether to compute the SHA-256. Reading ~1.4 GB per animal
            takes a few seconds warm and up to a minute from cold storage.

    Returns:
        A dict with the parsed geometry plus ``problems`` and ``warnings``
        lists. Never raises for data reasons — a parse failure is reported as
        a problem so the caller can carry on with the other stacks.
    """
    problems: list[str] = []
    warnings: list[str] = []
    out: dict[str, Any] = {
        "ets_path": path,
        "ets_bytes": None,
        "sha256": None,
        "problems": problems,
        "warnings": warnings,
    }

    try:
        check_materialised(path)
    except DehydratedFileError as exc:
        problems.append(str(exc))
        return out
    except OSError as exc:
        problems.append(f"cannot stat {path}: {exc}")
        return out

    try:
        index = read_ets_index(path)
    except (EtsFormatError, OSError, struct.error) as exc:
        problems.append(f"{os.path.basename(path)}: {exc}")
        return out

    out["ets_bytes"] = index["file_bytes"]
    dtype, bytes_per_sample = PIXEL_TYPES.get(
        index["pixel_type"], (f"unknown({index['pixel_type']})", 1)
    )
    if index["pixel_type"] not in PIXEL_TYPES:
        problems.append(f"unknown pixel type code {index['pixel_type']}")

    tile_x, tile_y, tile_z = index["tile"]
    if tile_x <= 0 or tile_y <= 0:
        problems.append(f"implausible tile size {tile_x}x{tile_y}")
        return out

    chunks = index["chunks"]
    n_dim = index["n_dim"]

    # Truncation check: every tile must lie inside the file.
    overflowing = sum(
        1 for _c, off, nb in chunks if off < 0 or nb < 0 or off + nb > index["file_bytes"]
    )
    if overflowing:
        problems.append(
            f"{overflowing:,} of {len(chunks):,} tiles point past the end of the file "
            f"- truncated or corrupt"
        )

    levels = sorted({c[0][-1] for c in chunks}) if chunks else []
    if not levels:
        problems.append("no tiles in the chunk table")
        return out
    base_level = min(levels)
    if levels != list(range(base_level, base_level + len(levels))):
        warnings.append(f"pyramid levels are not contiguous: {levels}")

    base = [c for c in chunks if c[0][-1] == base_level]
    channel_axis = 2 if n_dim >= 4 else None
    channels = sorted({c[0][channel_axis] for c in base}) if channel_axis is not None else [0]

    xs = [c[0][0] for c in base]
    ys = [c[0][1] for c in base]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    if min_x < 0 or min_y < 0:
        problems.append(f"negative tile indices (x from {min_x}, y from {min_y})")
    grid_cols = max_x + 1
    grid_rows = max_y + 1

    positions = {(c[0][0], c[0][1]) for c in base}
    occupied = len(positions)
    bbox_positions = grid_cols * grid_rows
    sparsity = 1.0 - occupied / bbox_positions if bbox_positions else None

    # Every channel must cover exactly the same tile positions, otherwise a
    # per-channel measurement would silently use different support.
    if channel_axis is not None and len(channels) > 1:
        per_channel = {
            ch: {(c[0][0], c[0][1]) for c in base if c[0][channel_axis] == ch}
            for ch in channels
        }
        reference = per_channel[channels[0]]
        mismatched = [ch for ch in channels[1:] if per_channel[ch] != reference]
        if mismatched:
            problems.append(
                f"channels {mismatched} do not cover the same tile positions as channel "
                f"{channels[0]} - per-channel support differs"
            )

    if hash_ets:
        try:
            out["sha256"] = sha256_file(path)
        except OSError as exc:
            problems.append(f"could not hash {os.path.basename(path)}: {exc}")

    out.update(
        {
            "dtype": dtype,
            "bytes_per_sample": bytes_per_sample,
            "compression_code": index["compression"],
            "compression": COMPRESSION_NAMES.get(
                index["compression"], f"code {index['compression']}"
            ),
            "quality": index["quality"],
            "tile_width_px": tile_x,
            "tile_height_px": tile_y,
            "tile_depth": tile_z,
            "n_dim": n_dim,
            "n_chunks": index["n_chunks"],
            "pyramid_levels": levels,
            "n_pyramid_levels": len(levels),
            "channels": channels,
            "n_channels": len(channels),
            "samples_per_tile": index["samples_per_tile"],
            "grid_cols": grid_cols,
            "grid_rows": grid_rows,
            # Tile-grid extent. This OVERESTIMATES the true image: the last row
            # and column of tiles are partial. Kept only for diagnostics.
            "grid_width_px": grid_cols * tile_x,
            "grid_height_px": grid_rows * tile_y,
            "tiles_present": occupied,
            "tiles_in_bbox": bbox_positions,
            "sparsity_fraction": sparsity,
        }
    )
    return out


# --------------------------------------------------------------------------
# VSI index metadata
# --------------------------------------------------------------------------

_TAG_INLINE_FLAG = 0x40000000
_TAG_EXTENDED_FLAG = 0x10000000

_TAG_IMAGE = 2001  # one record per image/series, inside the 2000 collection
_TAG_STAGE_POSITION = 2018  # VECTOR_DOUBLE_2, micrometres
_TAG_PIXEL_SIZE = 2019  # VECTOR_DOUBLE_2, micrometres
_TAG_SERIES_NAME = 2030  # UTF-16LE
_TAG_IMAGE_RECT = 2053  # RECT [x, y, width, height] - the TRUE image size


def _tag_block_start(buf: bytes, value_offset: int) -> int | None:
    """Return the offset of the nested tag block that a value holds, if any.

    Some extended fields put the nested block right at the value offset;
    others prefix it with a four-byte count. Both are recognised by looking
    for the 24-byte ``IS`` block header.
    """
    for candidate in (value_offset, value_offset + 4):
        if candidate + 24 > len(buf):
            continue
        if buf[candidate + 2 : candidate + 4] != b"IS":
            continue
        if struct.unpack_from("<H", buf, candidate)[0] != 24:
            continue
        return candidate
    return None


def read_vsi_index(vsi_path: str) -> list[dict[str, Any]]:
    """Extract per-series metadata from a ``.vsi`` index, payload not required.

    This is the whole point of having our own reader: Bio-Formats returns
    almost nothing from a ``.vsi`` whose companion folder is absent, but the
    index alone carries the series names, the true image rectangles, the pixel
    sizes and the stage coordinates.

    Tag-block layout: a 24-byte header (``u16`` headerSize = 24, magic ``IS``,
    ``u32`` version, ``u32`` dataFieldOffset, ``u32`` flags, ``u32`` nTags,
    ``u32`` reserved), then records of (``u32`` fieldType, ``u32`` tag,
    ``u32`` nextField). If ``fieldType & 0x40000000`` the value is four inline
    bytes; otherwise a ``u32`` length precedes the payload. ``nextField`` is a
    byte offset relative to the *block header start*, not the current record.

    Returns:
        One dict per series, in index order, with keys ``name``, ``kind``,
        ``section``, ``width_px``, ``height_px``, ``pixel_size_um``,
        ``stage_position_um``.
    """
    with open(vsi_path, "rb") as fh:
        buf = fh.read()

    series: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    visited_blocks: set[int] = set()

    def walk(block_start: int, depth: int) -> None:
        nonlocal current
        if depth > 24 or block_start + 24 > len(buf) or block_start in visited_blocks:
            return
        visited_blocks.add(block_start)
        _version, data_field_offset, _flags, n_tags, _reserved = struct.unpack_from(
            "<IIIII", buf, block_start + 4
        )
        if n_tags > 200_000:  # nonsense: refuse rather than churn
            return
        pos = block_start + data_field_offset
        seen_records: set[int] = set()
        for _ in range(n_tags):
            if pos + 12 > len(buf) or pos in seen_records:
                return
            seen_records.add(pos)
            field_type, tag, next_field = struct.unpack_from("<III", buf, pos)
            inline = bool(field_type & _TAG_INLINE_FLAG)
            extended = bool(field_type & _TAG_EXTENDED_FLAG)
            if inline:
                value_offset, value_len = pos + 12, 4
            else:
                (value_len,) = struct.unpack_from("<I", buf, pos + 12)
                value_offset = pos + 16
            if value_offset > len(buf):
                return
            # The declared length of a container field can overrun the file by
            # a few bytes, so it is used only to advance, never trusted as a
            # bound: every value read below checks its own extent.
            readable = value_offset + value_len <= len(buf)

            if depth == 1 and tag == _TAG_IMAGE:
                current = {
                    "name": None,
                    "kind": None,
                    "section": None,
                    "width_px": None,
                    "height_px": None,
                    "pixel_size_um": None,
                    "stage_position_um": None,
                }
                series.append(current)
            if current is not None and readable:
                if tag == _TAG_SERIES_NAME and not inline and value_len % 2 == 0:
                    # Each image carries a short numeric 2030 record as well as
                    # its UTF-16LE name, so ignore anything that decodes empty.
                    text = (
                        buf[value_offset : value_offset + value_len]
                        .decode("utf-16le", "replace")
                        .rstrip("\x00")
                        .strip()
                    )
                    if text and not current["name"]:
                        current["name"] = text
                elif tag == _TAG_IMAGE_RECT and value_len == 16:
                    if current["width_px"] is None:
                        _x, _y, width, height = struct.unpack_from("<4i", buf, value_offset)
                        current["width_px"] = width
                        current["height_px"] = height
                elif tag == _TAG_PIXEL_SIZE and value_len == 16:
                    if current["pixel_size_um"] is None:
                        current["pixel_size_um"] = list(
                            struct.unpack_from("<dd", buf, value_offset)
                        )
                elif tag == _TAG_STAGE_POSITION and value_len == 16:
                    if current["stage_position_um"] is None:
                        current["stage_position_um"] = list(
                            struct.unpack_from("<dd", buf, value_offset)
                        )

            if extended:
                nested = _tag_block_start(buf, value_offset)
                if nested is not None:
                    walk(nested, depth + 1)
            pos = block_start + next_field if next_field else value_offset + value_len

    walk(8, 0)  # the root tag block sits immediately after the 8-byte TIFF header

    for entry in series:
        name = entry["name"] or ""
        match = TISSUE_SERIES_RE.search(name)
        if match:
            entry["kind"] = "tissue"
            entry["section"] = match.group(1)
        elif name.lower() == "label":
            entry["kind"] = "label"
        elif name.lower() == "overview":
            entry["kind"] = "overview"
        elif name:
            entry["kind"] = "other"
        else:
            entry["kind"] = "unnamed"
    return series


# --------------------------------------------------------------------------
# Companion folder discovery
# --------------------------------------------------------------------------


def find_companion(vsi_path: str) -> str | None:
    """Return the folder holding this ``.vsi``'s stack folders, or ``None``.

    Matching is on the **exact** file stem. Substring matching is what allowed
    ``Image_5.vsi`` to bind to ``_Image_51_``'s pixel data and still report
    PASS, which is the single most dangerous thing this gate can get wrong.

    Accepted layouts, in order of preference:
      1. ``_<stem>_`` beside the ``.vsi`` (what the scanner writes)
      2. ``_<stem>`` or ``<stem>`` beside the ``.vsi``
      3. the ``.vsi``'s own folder, but *only* if it holds exactly one ``.vsi``
         (someone unzipped the stacks next to it; with several ``.vsi`` present
         this layout is ambiguous and is refused)
    """
    directory = os.path.dirname(os.path.abspath(vsi_path))
    stem = os.path.splitext(os.path.basename(vsi_path))[0]

    for name in (f"_{stem}_", f"_{stem}", stem):
        candidate = os.path.join(directory, name)
        if os.path.isdir(candidate) and _stack_dirs(candidate):
            return candidate

    if _stack_dirs(directory):
        siblings = glob.glob(os.path.join(directory, "*.vsi"))
        if len(siblings) == 1:
            return directory
    return None


def _stack_dirs(directory: str) -> list[str]:
    """Return the ``stackNNNNN`` subfolders of *directory*, sorted by ID."""
    found = []
    try:
        entries = os.listdir(directory)
    except OSError:
        return []
    for entry in entries:
        full = os.path.join(directory, entry)
        if entry.startswith("stack") and entry[5:].isdigit() and os.path.isdir(full):
            found.append(full)
    return sorted(found, key=lambda p: int(os.path.basename(p)[5:]))


def _classify_stack_files(stack_dir: str) -> tuple[list[str], list[str], list[str]]:
    """Split a stack folder's files into (tile containers, auxiliary, unexpected).

    Walks recursively, because unzipping stacks one at a time can produce
    ``stack10002/stack10002/frame_t.ets``. The allow-list is applied to the
    file *name*, so a Dropbox conflicted copy never counts as a tile container.
    """
    tiles: list[str] = []
    aux: list[str] = []
    unexpected: list[str] = []
    for dirpath, _dirnames, filenames in os.walk(stack_dir):
        for name in sorted(filenames):
            full = os.path.join(dirpath, name)
            if any(p.search(name) for p in IGNORED_FILE_PATTERNS):
                continue
            if name in ALLOWED_ETS_NAMES:
                tiles.append(full)
            elif any(p.match(name) for p in AUX_FILE_PATTERNS):
                aux.append(full)
            else:
                unexpected.append(full)
    return sorted(tiles), sorted(aux), sorted(unexpected)


def _find_placeholders(directory: str) -> list[str]:
    """Return one message per cloud-placeholder file found under *directory*.

    Run before any hashing, so a folder that was never downloaded fails in a
    second instead of after several minutes of futile reading.
    """
    messages: list[str] = []
    for dirpath, _dirnames, filenames in os.walk(directory):
        for name in sorted(filenames):
            if any(p.search(name) for p in IGNORED_FILE_PATTERNS):
                continue
            try:
                check_materialised(os.path.join(dirpath, name))
            except DehydratedFileError as exc:
                messages.append(str(exc))
            except OSError:
                continue
    return messages


def _directory_bytes(directory: str) -> int:
    total = 0
    for dirpath, _dirnames, filenames in os.walk(directory):
        for name in filenames:
            try:
                total += os.path.getsize(os.path.join(dirpath, name))
            except OSError:
                pass
    return total


# --------------------------------------------------------------------------
# Dataset verification
# --------------------------------------------------------------------------


def verify_dataset(
    vsi_path: str,
    *,
    hash_ets: bool = True,
    expect_sections: Sequence[int] = DEFAULT_EXPECTED_SECTION_COUNTS,
) -> dict[str, Any]:
    """Verify one ``.vsi`` dataset and its companion payload.

    Args:
        vsi_path: path to the ``.vsi`` index file.
        hash_ets: compute SHA-256 for every ``.ets`` as well as the ``.vsi``.
            Reads ~1.4 GB per animal: a few seconds warm, up to a minute from
            cold storage. Switch off for a quick structural pass. Note that
            hashing does not detect corruption on first sight - it *records*
            the digest, so a later re-verification can prove the bytes have
            not changed.
        expect_sections: acceptable numbers of tissue sections. The main cohort
            has 3 or 4; rescan slides hold only the re-acquired sections and
            need this set explicitly.

    Returns:
        A dict with ``ok``, ``vsi_path``, ``vsi_sha256``, ``companion``,
        ``index``, ``stacks``, ``problems`` and ``warnings``. It never raises
        for data reasons: every failure mode becomes an entry in ``problems``
        and ``ok`` is ``False``.
    """
    problems: list[str] = []
    warnings: list[str] = []
    result: dict[str, Any] = {
        "schema": SCHEMA,
        "tool_version": TOOL_VERSION,
        "verified_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "ok": False,
        "dataset": os.path.splitext(os.path.basename(vsi_path))[0],
        "vsi_path": os.path.abspath(vsi_path),
        "vsi_bytes": None,
        "vsi_sha256": None,
        "companion": None,
        "companion_bytes": None,
        "hash_ets": bool(hash_ets),
        "index": {"series": [], "n_tissue_series": 0, "tissue_sections": []},
        "stacks": [],
        "problems": problems,
        "warnings": warnings,
    }

    # --- the index file itself -------------------------------------------
    if not os.path.isfile(vsi_path):
        problems.append(f"no such file: {vsi_path}")
        return result
    try:
        check_materialised(vsi_path)
    except DehydratedFileError as exc:
        problems.append(str(exc))
        return result

    result["vsi_bytes"] = os.path.getsize(vsi_path)
    try:
        result["vsi_sha256"] = sha256_file(vsi_path)
    except OSError as exc:
        problems.append(f"could not read the index file: {exc}")
        return result

    try:
        series = read_vsi_index(vsi_path)
    except (OSError, struct.error) as exc:
        series = []
        problems.append(f"could not parse the VSI index metadata: {exc}")

    tissue_series = [s for s in series if s["kind"] == "tissue"]
    result["index"] = {
        "series": series,
        "n_series": len(series),
        "n_tissue_series": len(tissue_series),
        "tissue_sections": [s["section"] for s in tissue_series],
    }

    if series and not tissue_series:
        problems.append(
            "the index names no tissue series (nothing matching a trailing _01.._04) "
            "- this may not be a tissue slide"
        )
    if len(tissue_series) > len(set(s["section"] for s in tissue_series)):
        problems.append("duplicate tissue section numbers in the index")

    for entry in tissue_series:
        size = entry.get("pixel_size_um")
        if not size:
            warnings.append(f"series {entry['name']!r}: no pixel size in the index")
            continue
        if any(abs(v - EXPECTED_PIXEL_SIZE_UM) > PIXEL_SIZE_TOL_UM for v in size):
            problems.append(
                f"series {entry['name']!r}: pixel size {size[0]:.5f} x {size[1]:.5f} um is "
                f"outside {EXPECTED_PIXEL_SIZE_UM} +/- {PIXEL_SIZE_TOL_UM} um "
                f"- wrong objective or wrong slide?"
            )

    # --- the companion payload -------------------------------------------
    companion = find_companion(vsi_path)
    if companion is None:
        problems.append(
            "no companion payload folder found for this .vsi (expected a sibling folder "
            f"named _{result['dataset']}_ containing stackNNNNN subfolders). The pixel data "
            "is missing: a .vsi on its own is only an index. Re-transfer at the level of the "
            "parent folder so the .vsi and its companion folder travel together."
        )
        result["ok"] = False
        return result

    result["companion"] = os.path.abspath(companion)
    result["companion_bytes"] = _directory_bytes(companion)

    # Sweep the whole payload for cloud placeholders before doing anything
    # expensive. Hashing a 1.4 GB animal takes a minute; discovering on the
    # last stack that the folder was never downloaded wastes all of it.
    placeholders = _find_placeholders(companion)
    if placeholders:
        problems.extend(placeholders)
        return result

    stack_dirs = _stack_dirs(companion)
    if not stack_dirs:
        problems.append(f"companion folder {os.path.basename(companion)} has no stack subfolders")
        return result

    stack_ids = [int(os.path.basename(d)[5:]) for d in stack_dirs]
    unexpected_ids = sorted(set(stack_ids) - EXPECTED_STACK_IDS)
    if unexpected_ids:
        problems.append(
            f"unexpected stack IDs {unexpected_ids}; expected a subset of "
            f"{sorted(EXPECTED_STACK_IDS)}"
        )
    if LABEL_STACK_ID not in stack_ids:
        warnings.append(f"no stack{LABEL_STACK_ID} (slide label image) in the companion folder")
    if OVERVIEW_STACK_ID not in stack_ids:
        warnings.append(f"no stack{OVERVIEW_STACK_ID} (slide overview) in the companion folder")

    tissue_stack_ids = [i for i in stack_ids if i not in (LABEL_STACK_ID, OVERVIEW_STACK_ID)]
    expected = tuple(int(v) for v in expect_sections)
    if len(tissue_stack_ids) not in expected:
        problems.append(
            f"{len(tissue_stack_ids)} tissue stack(s) {tissue_stack_ids}; expected "
            f"{' or '.join(str(v) for v in expected)}. Rescan slides hold only the "
            f"re-acquired sections - run them with --expect-sections set accordingly."
        )
    if tissue_series and len(tissue_stack_ids) != len(tissue_series):
        problems.append(
            f"the index names {len(tissue_series)} tissue series but the companion folder "
            f"holds {len(tissue_stack_ids)} tissue stack(s) - the companion folder may belong "
            f"to a different slide, or the transfer is incomplete"
        )

    # Map .vsi series onto stack folders: the label and overview by name, the
    # tissue series onto the tissue stacks in ascending stack order (which is
    # acquisition order, the same order the index lists them in).
    by_stack: dict[int, dict[str, Any]] = {}
    label_series = next((s for s in series if s["kind"] == "label"), None)
    overview_series = next((s for s in series if s["kind"] == "overview"), None)
    if label_series:
        by_stack[LABEL_STACK_ID] = label_series
    if overview_series:
        by_stack[OVERVIEW_STACK_ID] = overview_series
    for stack_id, entry in zip(sorted(tissue_stack_ids), tissue_series):
        by_stack[stack_id] = entry

    # --- per-stack inspection --------------------------------------------
    for stack_dir, stack_id in zip(stack_dirs, stack_ids):
        stack_result = _verify_stack(stack_dir, stack_id, by_stack.get(stack_id))
        problems.extend(f"stack{stack_id}: {p}" for p in stack_result.pop("_problems"))
        warnings.extend(f"stack{stack_id}: {w}" for w in stack_result.pop("_warnings"))
        result["stacks"].append(stack_result)

    # Fill in hashes and geometry now that the structure is known.
    for stack_result in result["stacks"]:
        path = stack_result.get("ets_path")
        if not path:
            continue
        summary = summarise_ets(path, hash_ets=hash_ets)
        stack_problems = summary.pop("problems")
        stack_warnings = summary.pop("warnings")
        stack_id = stack_result["id"]
        problems.extend(f"stack{stack_id}: {p}" for p in stack_problems)
        warnings.extend(f"stack{stack_id}: {w}" for w in stack_warnings)
        stack_result.update(summary)
        geometry_problems: list[str] = []
        geometry_warnings: list[str] = []
        _cross_check_geometry(stack_result, geometry_problems, geometry_warnings)
        problems.extend(f"stack{stack_id}: {p}" for p in geometry_problems)
        warnings.extend(f"stack{stack_id}: {w}" for w in geometry_warnings)

    if not any(s.get("role") == "tissue" and s.get("dtype") in ("uint16", "int16")
               for s in result["stacks"]):
        problems.append(
            "no 16-bit tissue series found - the fluorescence stacks are missing "
            "(only the 8-bit label/overview images are present)"
        )

    result["ok"] = not problems
    return result


def _verify_stack(stack_dir: str, stack_id: int, series: dict[str, Any] | None) -> dict[str, Any]:
    """Locate and identify the single tile container inside one stack folder."""
    problems: list[str] = []
    warnings: list[str] = []
    if stack_id == LABEL_STACK_ID:
        role = "label"
    elif stack_id == OVERVIEW_STACK_ID:
        role = "overview"
    else:
        role = "tissue"

    tiles, aux, unexpected = _classify_stack_files(stack_dir)

    for path in unexpected:
        name = os.path.basename(path)
        if CONFLICT_PATTERN.search(name):
            problems.append(
                f"sync conflict file present: {name!r}. The pixel data in this folder is "
                f"ambiguous - resolve the conflict in Dropbox and re-run; do not analyse "
                f"the slide until exactly one tile file remains."
            )
        elif name.lower().endswith(".ets"):
            problems.append(
                f"unrecognised .ets file {name!r} - only {sorted(ALLOWED_ETS_NAMES)} and the "
                f"sample-mask blob files are legitimate tile containers"
            )
        else:
            warnings.append(f"unexpected file in the stack folder: {name!r}")

    ets_path: str | None = None
    if not tiles:
        problems.append(
            f"no tile container found (expected one of {sorted(ALLOWED_ETS_NAMES)}); "
            f"the transfer or unzip did not complete"
        )
    elif len(tiles) > 1:
        problems.append(
            "more than one tile container in this stack folder "
            f"({[os.path.relpath(t, stack_dir) for t in tiles]}) - ambiguous, refusing to guess"
        )
    else:
        ets_path = tiles[0]

    return {
        "id": stack_id,
        "role": role,
        "stack_dir": os.path.abspath(stack_dir),
        "ets_path": ets_path,
        "aux_files": [os.path.basename(p) for p in aux],
        "unexpected_files": [os.path.basename(p) for p in unexpected],
        "series_name": series["name"] if series else None,
        "section": series["section"] if series else None,
        # True dimensions come from VSI tag 2053 and are authoritative.
        "true_width_px": series["width_px"] if series else None,
        "true_height_px": series["height_px"] if series else None,
        "pixel_size_um": series["pixel_size_um"] if series else None,
        "stage_position_um": series["stage_position_um"] if series else None,
        "_problems": problems,
        "_warnings": warnings,
    }


def _cross_check_geometry(
    stack: dict[str, Any], problems: list[str], warnings: list[str]
) -> None:
    """Compare the tile grid against the true image rectangle from the index.

    Two things are established here.

    *Overestimate.* ``grid_cols * tile_width`` exceeds the true width because
    the last tile column is partial. Measured 0.6-8.4 % by area in this cohort.
    Both numbers are kept, clearly labelled, so nobody accidentally uses the
    grid product as an image size.

    *Consistency.* Every tile origin must lie inside the true image. If the
    companion folder belonged to a different slide, or the file were corrupt,
    tiles would fall outside. Note the converse does **not** hold: the grid can
    be *narrower* than the true image when whole tile columns at an edge are
    absent from the sample mask (tube 30's first section does exactly that), so
    "grid >= true" must not be asserted.
    """
    true_w = stack.get("true_width_px")
    true_h = stack.get("true_height_px")
    tile_w = stack.get("tile_width_px")
    tile_h = stack.get("tile_height_px")
    grid_cols = stack.get("grid_cols")
    grid_rows = stack.get("grid_rows")
    if not all(isinstance(v, int) and v > 0 for v in (true_w, true_h, tile_w, tile_h)):
        stack["grid_overestimate_fraction"] = None
        stack["sparsity_vs_true_grid"] = None
        if stack.get("ets_path") and true_w is None:
            warnings.append(
                "no true image rectangle in the index for this stack; only the tile-grid "
                "extent is available, which overestimates the image size"
            )
        return

    grid_w = stack["grid_width_px"]
    grid_h = stack["grid_height_px"]
    stack["grid_overestimate_fraction"] = (grid_w * grid_h) / (true_w * true_h) - 1.0

    if (grid_cols - 1) * tile_w >= true_w or (grid_rows - 1) * tile_h >= true_h:
        problems.append(
            f"tile grid ({grid_cols} x {grid_rows} tiles of {tile_w}x{tile_h}) does not fit "
            f"inside the true image ({true_w} x {true_h} px from VSI tag 2053) - the companion "
            f"folder may belong to a different slide, or the file is corrupt"
        )

    nominal_cols = math.ceil(true_w / tile_w)
    nominal_rows = math.ceil(true_h / tile_h)
    nominal_positions = nominal_cols * nominal_rows
    stack["nominal_grid_cols"] = nominal_cols
    stack["nominal_grid_rows"] = nominal_rows
    present = stack.get("tiles_present")
    stack["sparsity_vs_true_grid"] = (
        1.0 - present / nominal_positions if present and nominal_positions else None
    )

    if stack.get("role") == "tissue":
        sparsity = stack.get("sparsity_fraction")
        if sparsity is not None and not (SPARSITY_WARN_LO <= sparsity <= SPARSITY_WARN_HI):
            warnings.append(
                f"sparse-tile fraction {sparsity:.1%} is outside the usual "
                f"{SPARSITY_WARN_LO:.1%}-{SPARSITY_WARN_HI:.0%} for this cohort. Not an error, "
                f"but check the scanner sample mask before trusting any area denominator."
            )
        _assert_tissue_shape(stack, problems)


# Envelope measured over all 34 tissue stacks in this cohort (8 payload animals plus
# both rescans). A tissue container falling outside this is not tissue.
TISSUE_DTYPE = "uint16"
TISSUE_N_CHANNELS = 3
TISSUE_SAMPLES_PER_TILE = 1
TISSUE_MAX_SPARSITY_VS_TRUE = 0.30      # observed 0.036-0.179
TISSUE_MAX_GRID_OVERESTIMATE = 0.15     # observed -0.070 to +0.050


def _assert_tissue_shape(stack: dict[str, Any], problems: list[str]) -> None:
    """Hard-fail a tissue stack whose container is not shaped like tissue.

    These are *problems*, not warnings, and the distinction is the point of the gate.
    Binding the wrong pixels to a section is worse than a missing payload: the run
    completes, every number lands in a plausible range, and a section's Abeta/GFAP
    percent area is measured off an unrelated image. Dropping one of the three
    channel planes is the same class of error and equally invisible downstream.

    Concretely this catches an 8-bit RGB label or overview container substituted into
    a tissue stack directory, and a tissue container written with only 1 or 2 of its
    3 channel planes.
    """
    dtype = stack.get("dtype")
    if dtype is not None and dtype != TISSUE_DTYPE:
        problems.append(
            f"tissue stack has dtype {dtype!r}, expected {TISSUE_DTYPE!r} - this container is "
            f"not 16-bit fluorescence data (a label or overview image substituted into a "
            f"tissue stack looks exactly like this)"
        )

    n_ch = stack.get("n_channels")
    if isinstance(n_ch, int) and n_ch != TISSUE_N_CHANNELS:
        problems.append(
            f"tissue stack has {n_ch} channel plane(s), expected {TISSUE_N_CHANNELS} "
            f"(DAPI, FITC, Cy3) - a missing plane would silently remove a marker from every "
            f"measurement on this section"
        )

    spt = stack.get("samples_per_tile")
    if isinstance(spt, int) and spt != TISSUE_SAMPLES_PER_TILE:
        problems.append(
            f"tissue stack packs {spt} samples per tile, expected {TISSUE_SAMPLES_PER_TILE} - "
            f"that is an interleaved RGB container, not a 16-bit fluorescence one"
        )

    sparse_true = stack.get("sparsity_vs_true_grid")
    if isinstance(sparse_true, float) and sparse_true > TISSUE_MAX_SPARSITY_VS_TRUE:
        problems.append(
            f"{sparse_true:.1%} of the tile positions implied by the true image rectangle hold "
            f"no data (cohort maximum {TISSUE_MAX_SPARSITY_VS_TRUE:.0%}) - the companion folder "
            f"may belong to a different slide, or the transfer is incomplete"
        )

    over = stack.get("grid_overestimate_fraction")
    if isinstance(over, float) and abs(over) > TISSUE_MAX_GRID_OVERESTIMATE:
        problems.append(
            f"tile grid area differs from the true image rectangle by {over:+.1%} "
            f"(cohort range -7% to +5%) - grid and index disagree about the image size"
        )


# --------------------------------------------------------------------------
# Directory verification
# --------------------------------------------------------------------------


def verify_directory(
    directory: str,
    *,
    hash_ets: bool = True,
    expect_sections: Sequence[int] = DEFAULT_EXPECTED_SECTION_COUNTS,
    progress: bool = False,
) -> dict[str, Any]:
    """Verify every ``.vsi`` sitting directly inside *directory*.

    Subfolders are not searched: ``RawData/Rescan`` is a separate cohort with a
    different expected section count, so it is passed explicitly.

    Args:
        directory: folder to scan for ``*.vsi``.
        hash_ets: see :func:`verify_dataset`.
        expect_sections: see :func:`verify_dataset`.
        progress: write one line per dataset to stderr while running.

    Returns:
        A dict with ``ok``, ``directory``, ``n_datasets``, ``n_passed``,
        ``n_failed``, ``datasets`` (per-dataset results) and top-level
        ``problems`` / ``warnings``.
    """
    problems: list[str] = []
    warnings: list[str] = []
    vsi_files = sorted(glob.glob(os.path.join(directory, "*.vsi")), key=_natural_key)
    if not vsi_files:
        problems.append(f"no .vsi files found in {directory}")

    datasets: list[dict[str, Any]] = []
    for i, path in enumerate(vsi_files, start=1):
        if progress:
            print(f"[{i}/{len(vsi_files)}] {os.path.basename(path)}", file=sys.stderr, flush=True)
        datasets.append(
            verify_dataset(path, hash_ets=hash_ets, expect_sections=expect_sections)
        )

    n_passed = sum(1 for d in datasets if d["ok"])
    return {
        "schema": SCHEMA,
        "tool_version": TOOL_VERSION,
        "verified_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "ok": bool(datasets) and n_passed == len(datasets) and not problems,
        "directory": os.path.abspath(directory),
        "hash_ets": bool(hash_ets),
        "expect_sections": list(int(v) for v in expect_sections),
        "n_datasets": len(datasets),
        "n_passed": n_passed,
        "n_failed": len(datasets) - n_passed,
        "datasets": datasets,
        "problems": problems,
        "warnings": warnings,
    }


def _natural_key(path: str) -> tuple[Any, ...]:
    """Sort ``Image_9.vsi`` before ``Image_10.vsi``."""
    name = os.path.basename(path)
    return tuple(int(p) if p.isdigit() else p for p in re.split(r"(\d+)", name))


# --------------------------------------------------------------------------
# Human-readable report
# --------------------------------------------------------------------------


def _fmt_bytes(n: int | None) -> str:
    if n is None:
        return "?"
    if n >= 1 << 30:
        return f"{n / (1 << 30):.2f} GB"
    if n >= 1 << 20:
        return f"{n / (1 << 20):.1f} MB"
    return f"{n:,} B"


def format_dataset_report(result: dict[str, Any], *, verbose: bool = True) -> str:
    """Render one dataset result as plain text a wet-lab scientist can read."""
    lines: list[str] = []
    lines.append("=" * 88)
    lines.append(f"{os.path.basename(result['vsi_path'])}   [{'PASS' if result['ok'] else 'FAIL'}]")
    sha = result.get("vsi_sha256")
    lines.append(
        f"  index      {result['vsi_path']}  ({_fmt_bytes(result['vsi_bytes'])}"
        + (f", sha256 {sha[:16]}...)" if sha else ")")
    )
    index = result.get("index", {})
    sections = index.get("tissue_sections") or []
    if index.get("n_series"):
        lines.append(
            f"  index says {index['n_tissue_series']} tissue series "
            f"({' '.join('_' + s for s in sections) if sections else 'none'}) "
            f"of {index['n_series']} series total"
        )
    if result.get("companion"):
        lines.append(
            f"  payload    {result['companion']}  "
            f"({len(result['stacks'])} stacks, {_fmt_bytes(result['companion_bytes'])})"
        )
    else:
        lines.append("  payload    -- none found --")

    if verbose:
        for stack in result.get("stacks", []):
            lines.extend(_format_stack(stack))

    for problem in result.get("problems", []):
        lines.append(f"  PROBLEM  {problem}")
    for warning in result.get("warnings", []):
        lines.append(f"  warning  {warning}")
    lines.append(f"  RESULT: {'PASS' if result['ok'] else 'FAIL'}")
    return "\n".join(lines)


def _format_stack(stack: dict[str, Any]) -> list[str]:
    name = stack.get("series_name")
    header = f"  stack{stack['id']:<6} {stack['role']}"
    if stack.get("section"):
        header += f" _{stack['section']}"
    if name:
        header += f'  "{name}"'
    lines = [""]
    lines.append(header)
    if not stack.get("ets_path"):
        lines.append("      (no readable tile container - see problems below)")
        return lines

    true_w, true_h = stack.get("true_width_px"), stack.get("true_height_px")
    if true_w and true_h:
        lines.append(
            f"      true size   {true_w:,} x {true_h:,} px   "
            f"(VSI tag 2053 - AUTHORITATIVE)"
        )
    else:
        lines.append("      true size   unknown (not in the index)")

    over = stack.get("grid_overestimate_fraction")
    if over is None:
        note = "  (NOT the image size)"
    elif over >= 0:
        note = f"  ({over:+.1%} by area - OVERESTIMATE, do not use as image size)"
    else:
        note = (
            f"  ({over:+.1%} by area - SMALLER than the image: whole edge tiles are "
            f"absent from the sample mask. Not the image size either way)"
        )
    lines.append(
        f"      tile grid   {stack['grid_cols']} x {stack['grid_rows']} tiles of "
        f"{stack['tile_width_px']}x{stack['tile_height_px']} -> "
        f"{stack['grid_width_px']:,} x {stack['grid_height_px']:,} px" + note
    )

    sparsity = stack.get("sparsity_fraction")
    if sparsity is not None:
        note = (
            "  <- missing support, NOT background"
            if stack["role"] == "tissue"
            else ""
        )
        line = (
            f"      tiles       {stack['tiles_present']:,} of {stack['tiles_in_bbox']:,} "
            f"bounding-box positions present, {sparsity:.1%} absent{note}"
        )
        lines.append(line)
        nominal = stack.get("sparsity_vs_true_grid")
        if nominal is not None and abs(nominal - sparsity) > 0.005:
            lines.append(
                f"                  {nominal:.1%} absent relative to the full "
                f"{stack['nominal_grid_cols']} x {stack['nominal_grid_rows']} grid implied "
                f"by the true size"
            )
    lines.append(
        f"      pixels      {stack.get('dtype')}, {stack.get('n_channels')} channel plane(s), "
        f"{stack.get('samples_per_tile')} sample(s)/tile, {stack.get('compression')} "
        f"q{stack.get('quality')}, {stack.get('n_pyramid_levels')} pyramid levels"
    )
    px = stack.get("pixel_size_um")
    if px:
        lines.append(f"      pixel size  {px[0]:.5f} x {px[1]:.5f} um")
    pos = stack.get("stage_position_um")
    if pos and stack["role"] == "tissue":
        lines.append(f"      stage XY    {pos[0]:.1f}, {pos[1]:.1f} um")
    lines.append(
        f"      file        {os.path.basename(stack['ets_path'])}  "
        f"({_fmt_bytes(stack.get('ets_bytes'))})"
    )
    if stack.get("sha256"):
        lines.append(f"      sha256      {stack['sha256']}")
    else:
        lines.append("      sha256      (not computed - --no-hash)")
    if stack.get("aux_files"):
        lines.append(f"      aux files   {', '.join(stack['aux_files'])}")
    return lines


def format_directory_report(summary: dict[str, Any], *, verbose: bool = True) -> str:
    """Render a directory summary, one block per dataset plus a tally."""
    blocks = [format_dataset_report(d, verbose=verbose) for d in summary.get("datasets", [])]
    blocks.append("=" * 88)
    for problem in summary.get("problems", []):
        blocks.append(f"PROBLEM  {problem}")
    blocks.append(
        f"{summary['directory']}: {summary['n_passed']} of {summary['n_datasets']} "
        f"dataset(s) passed"
        + ("" if summary["hash_ets"] else "   (hashing disabled)")
    )
    failed = [d["dataset"] for d in summary.get("datasets", []) if not d["ok"]]
    if failed:
        blocks.append("FAILED: " + ", ".join(failed))
    return "\n".join(blocks)


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m ihc.ingest.verify",
        description=(
            "Verify Olympus VSI whole-slide datasets before analysis: companion payload "
            "present and belonging to this slide, tile containers parseable and complete, "
            "content hashed, stack inventory as expected, true image size reported."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Exit code is 0 only if every dataset passed, so this can gate a pipeline.\n"
            "Examples:\n"
            "  python -m ihc.ingest.verify RawData\n"
            "  python -m ihc.ingest.verify RawData/Rescan --expect-sections 2\n"
            "  python -m ihc.ingest.verify RawData --no-hash --json QC/verify.json\n"
        ),
    )
    parser.add_argument(
        "paths",
        nargs="+",
        metavar="PATH",
        help="a .vsi file, or a folder containing .vsi files (not searched recursively)",
    )
    parser.add_argument(
        "--no-hash",
        action="store_true",
        help="skip SHA-256 of the .ets payloads (much faster; loses the integrity check)",
    )
    parser.add_argument(
        "--expect-sections",
        default=",".join(str(v) for v in DEFAULT_EXPECTED_SECTION_COUNTS),
        metavar="N[,N...]",
        help="acceptable tissue-section counts per slide (default: %(default)s; "
        "rescan slides need 2)",
    )
    parser.add_argument("--json", metavar="FILE", help="write the full structured result here")
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="print only the per-dataset verdict lines, not the per-stack detail",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="report progress to stderr while hashing (a full cohort takes minutes)",
    )
    return parser


def main(argv: list[str]) -> int:
    """Command-line entry point.

    Args:
        argv: arguments **excluding** the program name, i.e. ``sys.argv[1:]``.

    Returns:
        0 if every dataset passed, 1 if any failed, 2 on a usage/IO error.
    """
    args = _build_parser().parse_args(argv)
    try:
        expect_sections = tuple(int(v) for v in args.expect_sections.split(",") if v.strip())
    except ValueError:
        print(f"--expect-sections: not a comma-separated list of integers: "
              f"{args.expect_sections!r}", file=sys.stderr)
        return 2
    if not expect_sections:
        print("--expect-sections: at least one count is required", file=sys.stderr)
        return 2

    hash_ets = not args.no_hash
    verbose = not args.quiet
    results: list[dict[str, Any]] = []
    reports: list[str] = []
    exit_code = 0

    for path in args.paths:
        if os.path.isdir(path):
            summary = verify_directory(
                path,
                hash_ets=hash_ets,
                expect_sections=expect_sections,
                progress=args.progress,
            )
            results.append(summary)
            reports.append(format_directory_report(summary, verbose=verbose))
            if not summary["ok"]:
                exit_code = 1
        elif os.path.isfile(path):
            result = verify_dataset(path, hash_ets=hash_ets, expect_sections=expect_sections)
            results.append(result)
            reports.append(format_dataset_report(result, verbose=verbose))
            if not result["ok"]:
                exit_code = 1
        else:
            print(f"no such file or directory: {path}", file=sys.stderr)
            exit_code = 2

    print("\n".join(reports))

    if args.json:
        payload: Any = results[0] if len(results) == 1 else results
        try:
            parent = os.path.dirname(os.path.abspath(args.json))
            if parent:
                os.makedirs(parent, exist_ok=True)
            with open(args.json, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=2, sort_keys=False)
                fh.write("\n")
            print(f"\nStructured result written to {args.json}")
        except OSError as exc:
            print(f"could not write {args.json}: {exc}", file=sys.stderr)
            exit_code = 2

    return exit_code


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main(sys.argv[1:]))
