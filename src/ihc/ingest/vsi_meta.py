"""Read Olympus VS200 ``.vsi`` index metadata without touching the pixel payload.

Why this module exists
----------------------
A ``.vsi`` file is only an *index*. The pixels live in a sibling folder
``_Image_NN_/stackNNNNN/*.ets`` that is ~1.4 GB per animal. Bio-Formats returns
almost nothing from a ``.vsi`` whose payload folder is absent, but everything the
ingest stage needs -- series names, stage coordinates, exposure times, pixel size
and true image dimensions -- is present in the 1.4-1.9 MB index itself.

This module reads that index with nothing but the Python standard library, so the
whole cohort can be inventoried from a laptop before any bulk data is rehydrated
from Dropbox.

What it is used for
-------------------
The scientifically load-bearing output is :func:`assign_boxes`. Each slide carries
3 or 4 brain sections arranged in *two PAP-pen boxes*. The two sections inside one
box share a staining condition: one box received primary antibody, the other
received DAPI + secondary only (the negative control). Which box is which must be
read from ``config/slides.csv`` -- it is a bench fact, never inferred from pixels.

But *box membership itself* has to be derived from the stage coordinates, because
the section number is **acquisition order, not slide position**. Animal 49 proves
it: sorted by stage X its sections run ``_02 _03 | _01 _04``, so the two sections
sharing a condition are ``_01`` and ``_04``. Anything that assumes ``_01``/``_02``
are one box and ``_03``/``_04`` the other is wrong, and gets the biology wrong.

File format notes (verified against all 31 index files in this cohort)
---------------------------------------------------------------------
A ``.vsi`` opens as a little-endian TIFF (``II*\\0``). Immediately after the 8-byte
TIFF header, at offset 8, sits the root Olympus tag block. A tag block is::

    u16  header_size      always 24
    2s   magic            b'IS'
    u32  version
    u64  data_field_offset
    u32  flags            n_tags = flags & 0x0FFFFFFF   (often 0; not trustworthy)
    u32  reserved

Records begin at ``block_offset + header_size`` and are chained by ``next_field``,
which is a byte offset **relative to the start of the block**, not to the current
record. Each record is::

    u32  field_type
    u32  tag
    u32  next_field       0 terminates the chain
    u32  length           (absent when the INLINE flag is set)
    u32  extra            (present only when the EXTRA flag is set)
    ...  value

``field_type`` carries flags in its top bits (``INLINE``, ``ARRAY``, ``EXTEND``,
``EXTRA``, ``VOLUME``) and the value type in ``field_type & 0xFFFFFF``.

Two traps, both hit during development and both handled below:

1. ``data_field_offset`` and ``n_tags`` are legitimately 0 in some blocks that
   nonetheless contain records. Records always start at ``header_size``, and the
   ``next_field`` chain -- not ``n_tags`` -- is what terminates iteration.
2. For ``VOLUME`` records the ``length`` field is *not* a usable byte count (it is
   the block-relative end offset minus 20, and is nonsense for the root record).
   A nested block must therefore be bounded by the *next sibling record*, i.e. by
   ``block_offset + next_field``. Without that bound, an empty nested block runs
   straight off its own end and swallows the following sibling -- which silently
   deletes the tissue-series metadata for every animal.

Tag numbers used here
---------------------
====== ==========================================================
2000   root image collection
2001   one image (tissue series, label, overview, focus map, ...)
2002   dimension block; ``2018``/``2037`` -> ``2053``
2005   image properties block
2007   channel collection; ``2008`` = one channel
2015   acquisition timestamp, Unix epoch seconds
2018   stage origin (X, Y) in micrometres, VECTOR_DOUBLE_2
2019   pixel size (X, Y) in micrometres, VECTOR_DOUBLE_2
2021   channel name  (``2419`` is the same name, used as fallback)
2030   series name
2053   true image RECT ``[0, 0, width, height]`` in pixels
2061   tube / slide ID string (``120635`` duplicates it)
100002 exposure time in MICROseconds
====== ==========================================================

``2053`` is the authority on image size: the tile-grid product from the ``.ets``
payload overestimates by up to 2.4 %.

Public API
----------
:func:`read_vsi_meta`, :func:`assign_boxes`, :class:`SeriesMeta`, :class:`VsiMeta`.
Run the module directly to print a per-animal table for the whole cohort.
"""

from __future__ import annotations

import datetime as _dt
import os
import re
import struct
import sys
from dataclasses import dataclass, field
from typing import Iterable, Sequence

__all__ = [
    "SeriesMeta",
    "VsiMeta",
    "BoxAssignmentError",
    "VsiParseError",
    "read_vsi_meta",
    "assign_boxes",
]


# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

#: Cohort default location of the raw ``.vsi`` index files (used only by __main__).
DEFAULT_RAW_DATA_DIR = (
    "~/Dropbox/longevity_project/projects/Rapamycin_CSF_mice_Per"
    "/IHC_analysis_pipeline/RawData"
)

#: Minimum ratio of the between-box stage gap to the largest within-box gap.
#: Measured range across the 31 cohort slides is 1.33 - 4.00; anything below this
#: means the two PAP-pen boxes are not cleanly separated and must be looked at.
MIN_GAP_RATIO = 1.3

#: Below this the split is real but uncomfortably tight (tube 55 sits at 1.33).
#: Not enforced -- reported by the __main__ table so a human can eyeball it.
MARGINAL_GAP_RATIO = 1.5

#: Expected pixel size for the 20x objective used throughout this cohort, and the
#: tolerance for the warning. Never compare pixel sizes for bit-equality: the
#: scanner writes a slightly different value for every series.
EXPECTED_PIXEL_SIZE_UM = 0.325
PIXEL_SIZE_TOL_UM = 0.001

#: Channels as acquired, in order. Exposure records appear in this order.
EXPECTED_CHANNELS = ("DAPI", "FITC", "Cy3")

# Series names are UTF-16LE and NOT uniform: most read "20x_DAPI, FITC, Cy3_01"
# but tube 60 reads "60_20x_DAPI, FITC, Cy3_01". Anchor on the trailing _0N only.
_SECTION_RE = re.compile(r"_(0[1-4])$")

# Tube ID from the filename, e.g. "Image_49.vsi" -> 49.
_FILENAME_TUBE_RE = re.compile(r"(\d+)")

# --- tag-block binary layout ---------------------------------------------- #

_BLOCK_MAGIC = b"\x18\x00IS"  # header_size=24 followed by b'IS'
_BLOCK_HEADER_SIZE = 24
_RECORD_HEADER_SIZE = 12  # field_type + tag + next_field

# field_type flag bits. EXTEND and ARRAY are recorded for completeness only:
# the parser does not need them, because leaf values are read using the declared
# length and nested blocks are located by their magic bytes.
_FLAG_EXTRA = 0x08000000  # value is preceded by an extra u32
_FLAG_EXTEND = 0x10000000  # noqa: F841 - documents the layout
_FLAG_ARRAY = 0x20000000  # noqa: F841 - documents the layout
_FLAG_INLINE = 0x40000000  # value is 4 bytes stored in place of the length field
_FLAG_VOLUME = 0x80000000  # value is a nested tag block

# Tag numbers (see module docstring).
_TAG_COLLECTION = 2000
_TAG_IMAGE = 2001
_TAG_DIMENSIONS = 2002
_TAG_PROPERTIES = 2005
_TAG_CHANNELS = 2007
_TAG_CHANNEL = 2008
_TAG_TIMESTAMP = 2015
_TAG_STAGE_ORIGIN = 2018
_TAG_PIXEL_SIZE = 2019
_TAG_CHANNEL_NAME = 2021
_TAG_CHANNEL_NAME_ALT = 2419
_TAG_SERIES_NAME = 2030
_TAG_DIM_ALT = 2037
_TAG_IMAGE_RECT = 2053
_TAG_SLIDE_INFO = 2062
_TAG_TUBE_ID = 2061
_TAG_TUBE_ID_ALT = 120635
_TAG_SLIDE_BLOCK = 2004
_TAG_EXPOSURE_US = 100002


# --------------------------------------------------------------------------- #
# Exceptions
# --------------------------------------------------------------------------- #


class VsiParseError(RuntimeError):
    """The file is not a readable ``.vsi`` index.

    Raised only for structural failures -- a truncated file, a missing Olympus
    tag block, an absent image collection, or a Dropbox placeholder that has not
    been downloaded. Files that are merely *odd* (unexpected exposure, an
    unusual number of sections, a tube-ID mismatch) never raise; those problems
    are collected in :attr:`VsiMeta.warnings`.
    """


class BoxAssignmentError(ValueError):
    """The sections could not be split into two plausible PAP-pen boxes."""


# --------------------------------------------------------------------------- #
# Public data model
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SeriesMeta:
    """Metadata for one imaged brain section (one tissue series in the ``.vsi``).

    Attributes:
        name: Full internal series name, e.g. ``"20x_DAPI, FITC, Cy3_01"``.
        section_label: ``"01"`` .. ``"04"``, parsed from the trailing ``_0N``.
            This is **acquisition order, not slide position** -- do not use it to
            infer which PAP-pen box the section sat in.
        stage_x_um: Stage origin X in micrometres. Increases with distance from
            the slide's label end, so it orders the sections along the slide.
        stage_y_um: Stage origin Y in micrometres.
        exposure_ms: Per-channel exposure in milliseconds, keyed by channel name
            (normally ``DAPI``, ``FITC``, ``Cy3``). Converted from the
            microseconds stored in the file.
        pixel_size_um: Pixel size along X in micrometres, or ``None`` if absent.
        width_px: True image width from tag 2053, or ``None``.
        height_px: True image height from tag 2053, or ``None``.
        acquisition_time: Acquisition timestamp as an ISO-8601 UTC string, or
            ``""`` if the file did not record one.
    """

    name: str
    section_label: str
    stage_x_um: float
    stage_y_um: float
    exposure_ms: dict[str, float]
    pixel_size_um: float | None
    width_px: int | None
    height_px: int | None
    acquisition_time: str = ""


@dataclass
class VsiMeta:
    """Everything recoverable from one ``.vsi`` index file.

    Attributes:
        path: Absolute path to the ``.vsi`` that was read.
        tube_id: Animal / tube number parsed from the filename, cross-checked
            against the in-file tag 2061. ``None`` if the filename has no number.
        tube_id_in_file: Raw tube-ID string stored in the file, or ``None``.
        n_tissue_series: Number of usable tissue series -- always
            ``len(series)``. Expected to be 3 or 4 for this cohort.
        series: Tissue series **ordered by ``stage_x_um`` ascending**, i.e. in
            physical order along the slide starting at the label end. Label,
            overview, sample mask and focus-map series are excluded.
        acquisition_times: ISO-8601 UTC timestamps, one per entry of
            ``series`` and in the same order. ``""`` where the file had none.
        warnings: Human-readable notes about anything unexpected. Always check
            this; the reader deliberately does not raise on odd-but-readable
            files.
    """

    path: str
    tube_id: int | None
    tube_id_in_file: str | None
    n_tissue_series: int
    series: list[SeriesMeta]
    acquisition_times: list[str]
    warnings: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# Low-level tag-block parser
# --------------------------------------------------------------------------- #


class _Tag:
    """One record inside an Olympus tag block.

    ``value`` holds the raw little-endian bytes for leaf records and is empty for
    collection (``VOLUME``) records, whose content is in ``children``.
    """

    __slots__ = ("tag", "field_type", "value_type", "offset", "value", "children")

    def __init__(
        self,
        tag: int,
        field_type: int,
        offset: int,
        value: bytes,
        children: list["_Tag"],
    ) -> None:
        self.tag = tag
        self.field_type = field_type
        self.value_type = field_type & 0xFFFFFF
        self.offset = offset
        self.value = value
        self.children = children

    def child(self, tag: int) -> "_Tag | None":
        """Return the first direct child with ``tag``, or ``None``."""
        for c in self.children:
            if c.tag == tag:
                return c
        return None

    def children_with(self, tag: int) -> list["_Tag"]:
        """Return all direct children with ``tag``, in file order."""
        return [c for c in self.children if c.tag == tag]

    def descendants_with(self, tag: int, _acc: list["_Tag"] | None = None) -> list["_Tag"]:
        """Return every descendant with ``tag``, depth-first in file order."""
        acc: list["_Tag"] = [] if _acc is None else _acc
        for c in self.children:
            if c.tag == tag:
                acc.append(c)
            c.descendants_with(tag, acc)
        return acc

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"<_Tag {self.tag} @{self.offset} nchildren={len(self.children)}>"


def _parse_block(
    data: bytes,
    block_offset: int,
    limit: int,
    depth: int = 0,
    _stack: tuple[int, ...] = (),
    max_depth: int = 40,
) -> list[_Tag]:
    """Parse one tag block and return its records.

    Args:
        data: The whole file.
        block_offset: Byte offset of the block header (must start with
            ``b'\\x18\\x00IS'``).
        limit: Exclusive upper bound for this block. Nothing at or beyond this
            offset belongs to the block. This is what stops an empty nested block
            from consuming its parent's next sibling.
        depth: Current recursion depth.
        _stack: Block offsets currently on the recursion stack, to break cycles.
        max_depth: Hard recursion cap.

    Returns:
        The block's records in file order. An empty list if ``block_offset`` does
        not hold a valid block header -- callers treat that as "no children".
    """
    if depth > max_depth or block_offset in _stack:
        return []
    end = min(limit, len(data))
    if block_offset < 0 or block_offset + _BLOCK_HEADER_SIZE > end:
        return []
    if data[block_offset : block_offset + 4] != _BLOCK_MAGIC:
        return []

    _stack = _stack + (block_offset,)
    records: list[_Tag] = []

    # Records always begin at header_size. `data_field_offset` and the n_tags
    # field in the header are 0 in several real blocks that do contain records,
    # so neither is used to drive iteration.
    rel = _BLOCK_HEADER_SIZE
    while True:
        rec_offset = block_offset + rel
        if rec_offset + _RECORD_HEADER_SIZE + 4 > end:
            break
        field_type, tag, next_field = struct.unpack_from("<III", data, rec_offset)

        is_inline = bool(field_type & _FLAG_INLINE)
        is_volume = bool(field_type & _FLAG_VOLUME)
        has_extra = bool(field_type & _FLAG_EXTRA)

        # A record's value may not reach past the next sibling; the last record
        # in a block may not reach past the block's own limit.
        rec_end = block_offset + next_field if next_field > rel else end

        value = b""
        children: list[_Tag] = []
        if is_inline:
            value = data[rec_offset + 12 : rec_offset + 16]
        else:
            value_offset = rec_offset + 16 + (4 if has_extra else 0)
            if is_volume:
                # The declared length is unreliable for collections; find the
                # nested block header just after the record header instead.
                nested = _find_block(data, value_offset, rec_end)
                if nested is not None:
                    children = _parse_block(
                        data, nested, min(rec_end, end), depth + 1, _stack, max_depth
                    )
            else:
                (length,) = struct.unpack_from("<I", data, rec_offset + 12)
                if 0 <= length <= end - value_offset:
                    value = data[value_offset : value_offset + length]

        records.append(_Tag(tag, field_type, rec_offset, value, children))

        if next_field == 0 or next_field <= rel:
            break
        rel = next_field

    return records


def _find_block(data: bytes, start: int, limit: int) -> int | None:
    """Locate the nested tag-block header at or just after ``start``.

    The header normally sits exactly at ``start``, but the ``EXTRA``/``EXTEND``
    flag combinations shift it by a word or two. Probing a short aligned window
    is cheaper and far more robust than modelling every flag permutation.
    """
    for delta in range(0, 25, 4):
        pos = start + delta
        if pos + 4 > min(limit, len(data)):
            return None
        if data[pos : pos + 4] == _BLOCK_MAGIC:
            return pos
    return None


# --- value decoders --------------------------------------------------------- #


def _as_str(node: "_Tag | None") -> str | None:
    """Decode a UTF-16LE string value."""
    if node is None:
        return None
    return node.value.decode("utf-16-le", errors="replace").rstrip("\x00")


def _as_int32(node: "_Tag | None") -> int | None:
    """Decode a signed 32-bit integer value."""
    if node is None or len(node.value) < 4:
        return None
    return struct.unpack_from("<i", node.value)[0]


def _as_int64(node: "_Tag | None") -> int | None:
    """Decode a signed 64-bit integer value."""
    if node is None or len(node.value) < 8:
        return None
    return struct.unpack_from("<q", node.value)[0]


def _as_double_pair(node: "_Tag | None") -> tuple | None:
    """Decode a VECTOR_DOUBLE_2 value as ``(x, y)``."""
    if node is None or len(node.value) < 16:
        return None
    return struct.unpack_from("<2d", node.value)


def _as_rect(node: "_Tag | None") -> tuple | None:
    """Decode a RECT_INT value as ``(x, y, width, height)``."""
    if node is None or len(node.value) < 16:
        return None
    return struct.unpack_from("<4i", node.value)


# --------------------------------------------------------------------------- #
# Storage sanity
# --------------------------------------------------------------------------- #


def _check_materialised(path: str) -> None:
    """Fail fast on a Dropbox online-only placeholder.

    A dehydrated file reports a plausible size from ``stat()`` but reading it
    stalls or errors. On macOS a materialised file has allocated blocks roughly
    matching its size; a placeholder has almost none.
    """
    try:
        st = os.stat(path)
    except OSError as exc:
        raise VsiParseError(f"cannot stat {path}: {exc}") from exc

    if st.st_size == 0:
        raise VsiParseError(f"{path} is empty")

    blocks = getattr(st, "st_blocks", None)
    if blocks is None:  # platform without st_blocks; nothing to check
        return
    allocated = blocks * 512
    if allocated < st.st_size * 0.5:
        raise VsiParseError(
            f"{path} looks like an online-only Dropbox placeholder "
            f"({st.st_size} bytes reported, ~{allocated} bytes on disk). "
            f"Make it available offline before running the pipeline."
        )


# --------------------------------------------------------------------------- #
# Series extraction
# --------------------------------------------------------------------------- #


def _channel_exposures(
    channel_block: _Tag, warnings: list[str], context: str
) -> tuple[str, float | None]:
    """Return ``(channel_name, exposure_ms)`` for one 2008 channel block.

    The exposure (tag 100002, microseconds) is nested several levels below the
    channel inside its device-settings subtree, so it is found by descendant
    search rather than a fixed path. Exactly one record is expected per channel;
    the label and overview images carry their own exposures elsewhere and are
    never reached from here.
    """
    name = _as_str(channel_block.child(_TAG_CHANNEL_NAME))
    if not name:
        name = _as_str(channel_block.child(_TAG_CHANNEL_NAME_ALT))
    if not name:
        name = "?"

    records = channel_block.descendants_with(_TAG_EXPOSURE_US)
    values = [v for v in (_as_int32(r) for r in records) if v is not None]
    if len(values) != 1:
        warnings.append(
            f"{context}: expected exactly 1 exposure record for channel {name!r}, "
            f"found {len(values)}"
        )
    if not values:
        return name, None
    return name, values[0] / 1000.0


def _image_rect(image: _Tag) -> tuple[int, int, int, int] | None:
    """Return the true ``(x, y, width, height)`` from tag 2053, or ``None``.

    Tag 2053 hangs off the 2002 dimension block under either 2018 or 2037; both
    agree on every file in this cohort, and 2018 is preferred.
    """
    dims = image.child(_TAG_DIMENSIONS)
    if dims is None:
        return None
    for sub_tag in (_TAG_STAGE_ORIGIN, _TAG_DIM_ALT):
        sub = dims.child(sub_tag)
        if sub is None:
            continue
        rect = _as_rect(sub.child(_TAG_IMAGE_RECT))
        if rect is not None:
            return rect
    return None


def _build_series(
    image: _Tag, name: str, label: str, warnings: list[str]
) -> SeriesMeta | None:
    """Build a :class:`SeriesMeta` for one tissue series, or ``None`` if unusable."""
    props = image.child(_TAG_PROPERTIES)
    if props is None:
        warnings.append(f"series {name!r}: no properties block (tag 2005); skipped")
        return None

    stage = _as_double_pair(props.child(_TAG_STAGE_ORIGIN))
    if stage is None:
        warnings.append(
            f"series {name!r}: no stage origin (tag 2018); skipped, because box "
            f"assignment is impossible without it"
        )
        return None

    pixel = _as_double_pair(props.child(_TAG_PIXEL_SIZE))
    pixel_size_um = None
    if pixel is None:
        warnings.append(f"series {name!r}: no pixel size (tag 2019)")
    else:
        pixel_size_um = pixel[0]
        if abs(pixel[0] - pixel[1]) > PIXEL_SIZE_TOL_UM:
            warnings.append(
                f"series {name!r}: anisotropic pixel size "
                f"X={pixel[0]:.6f} Y={pixel[1]:.6f} um"
            )
        if abs(pixel_size_um - EXPECTED_PIXEL_SIZE_UM) > PIXEL_SIZE_TOL_UM:
            warnings.append(
                f"series {name!r}: pixel size {pixel_size_um:.6f} um is more than "
                f"{PIXEL_SIZE_TOL_UM} um from the expected "
                f"{EXPECTED_PIXEL_SIZE_UM} um"
            )

    rect = _image_rect(image)
    if rect is None:
        warnings.append(f"series {name!r}: no image rectangle (tag 2053)")
        width_px = height_px = None
    else:
        width_px, height_px = int(rect[2]), int(rect[3])

    exposure_ms: dict[str, float] = {}
    channels = image.child(_TAG_CHANNELS)
    channel_blocks = channels.children_with(_TAG_CHANNEL) if channels else []
    if not channel_blocks:
        warnings.append(f"series {name!r}: no channel blocks (tag 2007/2008)")
    observed: list[str] = []
    for block in channel_blocks:
        channel_name, value = _channel_exposures(block, warnings, f"series {name!r}")
        observed.append(channel_name)
        if value is not None:
            if channel_name in exposure_ms:
                warnings.append(
                    f"series {name!r}: duplicate channel {channel_name!r}; "
                    f"keeping the first exposure"
                )
            else:
                exposure_ms[channel_name] = value
    if observed and tuple(observed) != EXPECTED_CHANNELS:
        warnings.append(
            f"series {name!r}: channels are {tuple(observed)}, expected "
            f"{EXPECTED_CHANNELS}"
        )

    epoch = _as_int64(props.child(_TAG_TIMESTAMP))
    if epoch is None:
        warnings.append(f"series {name!r}: no acquisition timestamp (tag 2015)")
        acquisition_time = ""
    else:
        acquisition_time = _dt.datetime.fromtimestamp(
            epoch, _dt.timezone.utc
        ).isoformat()

    return SeriesMeta(
        name=name,
        section_label=label,
        stage_x_um=float(stage[0]),
        stage_y_um=float(stage[1]),
        exposure_ms=exposure_ms,
        pixel_size_um=pixel_size_um,
        width_px=width_px,
        height_px=height_px,
        acquisition_time=acquisition_time,
    )


def _tube_id_in_file(collection: _Tag) -> str | None:
    """Return the tube / slide ID string stored in the file, or ``None``."""
    slide_block = collection.child(_TAG_SLIDE_BLOCK)
    if slide_block is None:
        return None
    for info in slide_block.children_with(_TAG_SLIDE_INFO):
        for tag in (_TAG_TUBE_ID, _TAG_TUBE_ID_ALT):
            value = _as_str(info.child(tag))
            if value:
                return value
    return None


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #


def read_vsi_meta(path: str) -> VsiMeta:
    """Read metadata from a ``.vsi`` index file, without its pixel payload.

    Args:
        path: Path to a ``.vsi`` index file. The sibling ``_Image_NN_`` payload
            folder does not need to be present.

    Returns:
        A :class:`VsiMeta` whose ``series`` list is sorted by ``stage_x_um``
        ascending, i.e. in physical order along the slide from the label end.
        Anything unexpected but readable lands in ``warnings`` -- always check it.

    Raises:
        VsiParseError: The file is unreadable, is not a VSI container, is a
            Dropbox online-only placeholder, or holds no image collection.

    Example:
        >>> meta = read_vsi_meta("RawData/Image_49.vsi")     # doctest: +SKIP
        >>> [s.section_label for s in meta.series]           # doctest: +SKIP
        ['02', '03', '01', '04']
    """
    path = os.path.abspath(path)
    _check_materialised(path)

    try:
        with open(path, "rb") as handle:
            data = handle.read()
    except OSError as exc:
        raise VsiParseError(f"cannot read {path}: {exc}") from exc

    if len(data) < 64:
        raise VsiParseError(f"{path} is only {len(data)} bytes; not a .vsi index")
    if data[:4] != b"II*\x00":
        raise VsiParseError(
            f"{path} does not start with the little-endian TIFF magic 'II*\\0'"
        )

    warnings: list[str] = []

    roots = _parse_block(data, 8, len(data))
    if not roots:
        raise VsiParseError(
            f"{path}: no Olympus tag block at offset 8; the file is not a VSI index "
            f"or is truncated"
        )

    collections = [node for node in roots if node.tag == _TAG_COLLECTION]
    if not collections:
        raise VsiParseError(f"{path}: no image collection (tag 2000) in the tag tree")
    if len(collections) > 1:
        warnings.append(
            f"{len(collections)} image collections (tag 2000); using the first"
        )
    collection = collections[0]

    images = collection.children_with(_TAG_IMAGE)
    if not images:
        raise VsiParseError(f"{path}: image collection contains no images (tag 2001)")

    series: list[SeriesMeta] = []
    for image in images:
        props = image.child(_TAG_PROPERTIES)
        name = _as_str(props.child(_TAG_SERIES_NAME)) if props else None
        if not name:
            continue
        match = _SECTION_RE.search(name)
        if match is None:
            # Label, Overview, Sample Mask, 20x FocusMap, FocusPoints.
            continue
        built = _build_series(image, name, match.group(1), warnings)
        if built is not None:
            series.append(built)

    series.sort(key=lambda s: s.stage_x_um)

    labels = [s.section_label for s in series]
    if len(set(labels)) != len(labels):
        warnings.append(f"duplicate section labels {labels}")
    if len(series) not in (3, 4):
        warnings.append(
            f"{len(series)} tissue series; this cohort has 3 or 4 per slide "
            f"(rescans legitimately have 2)"
        )

    filename_match = _FILENAME_TUBE_RE.search(os.path.basename(path))
    tube_id = int(filename_match.group(1)) if filename_match else None
    if tube_id is None:
        warnings.append(f"no tube number in filename {os.path.basename(path)!r}")

    in_file = _tube_id_in_file(collection)
    if in_file is None:
        # A truncated .vsi is the danger here, not a merely-odd one. All 31 intact
        # index files carry tag 2061; truncation is what removes it. And a truncated
        # 4-section index degrades to a *plausible* 3-section one -- this cohort
        # legitimately contains six 3-section slides, so the section count alone
        # cannot distinguish them, and the box assignment then comes out confidently
        # wrong (Image_49.vsi truncated to 90% yields far_label=['01'] instead of
        # ['01','04']). Refuse rather than guess.
        raise VsiParseError(
            f"{os.path.basename(path)}: no tube ID (tag 2061) in the index. Every intact "
            f"file in this cohort has one, so this file is most likely TRUNCATED. A "
            f"truncated index silently loses whole sections and produces a confidently "
            f"wrong PAP-pen box assignment. Re-transfer the file.")
    elif tube_id is not None:
        try:
            matches = int(in_file.strip()) == tube_id
        except ValueError:
            matches = False
        if not matches:
            warnings.append(
                f"tube ID mismatch: filename says {tube_id}, file says {in_file!r}"
            )

    return VsiMeta(
        path=path,
        tube_id=tube_id,
        tube_id_in_file=in_file,
        n_tissue_series=len(series),
        series=series,
        acquisition_times=[s.acquisition_time for s in series],
        warnings=warnings,
    )


def assign_boxes(series: Sequence[SeriesMeta]) -> dict:
    """Split sections into the two PAP-pen boxes using the stage-X gap.

    Each slide carries two PAP-pen boxes. The sections within one box sit close
    together; the two boxes are separated by a much larger gap. Sorting by stage
    X and splitting at the largest gap therefore recovers box membership -- which
    the section number does **not**, because it records acquisition order rather
    than slide position.

    ``near_label`` is the low-stage-X box, since stage X measures distance from
    the label end of the slide.

    Args:
        series: :class:`SeriesMeta` objects for one slide, in any order.

    Returns:
        dict with:

        - ``near_label``: section labels in the low-stage-X box, ascending.
        - ``far_label``: section labels in the high-stage-X box, ascending.
        - ``gap_ratio``: between-box gap divided by the **largest** within-box
          gap. Measured range across this cohort is 1.33 - 4.00.
        - ``within_gaps_mm``: within-box gaps in millimetres, in stage-X order.
        - ``split_gap_mm``: the between-box gap in millimetres.

    Raises:
        BoxAssignmentError: Fewer than 3 sections, a split that is not 2+2, 2+1
            or 1+2, or ``gap_ratio`` below :data:`MIN_GAP_RATIO`. A two-section
            rescan legitimately fails here: both of its sections come from a
            single box, so there is no split to find and box membership must be
            carried over from the original scan.

    Example:
        >>> boxes = assign_boxes(read_vsi_meta("Image_49.vsi").series)  # doctest: +SKIP
        >>> boxes["near_label"], boxes["far_label"]                     # doctest: +SKIP
        (['02', '03'], ['01', '04'])
    """
    ordered = sorted(series, key=lambda s: s.stage_x_um)
    if len(ordered) < 3:
        raise BoxAssignmentError(
            f"need at least 3 sections to find a box split, got {len(ordered)}"
        )

    xs = [s.stage_x_um for s in ordered]
    gaps = [xs[i + 1] - xs[i] for i in range(len(xs) - 1)]
    split_index = max(range(len(gaps)), key=lambda i: gaps[i])

    near = [s.section_label for s in ordered[: split_index + 1]]
    far = [s.section_label for s in ordered[split_index + 1 :]]
    shape = (len(near), len(far))
    if shape not in {(2, 2), (2, 1), (1, 2)}:
        raise BoxAssignmentError(
            f"stage-X split gives {shape[0]}+{shape[1]} sections; expected 2+2, "
            f"2+1 or 1+2. Section labels in stage order: "
            f"{[s.section_label for s in ordered]}"
        )

    split_gap = gaps[split_index]
    within = [g for i, g in enumerate(gaps) if i != split_index]
    gap_ratio = split_gap / max(within)
    if gap_ratio < MIN_GAP_RATIO:
        raise BoxAssignmentError(
            f"between-box gap {split_gap / 1000:.3f} mm is only {gap_ratio:.2f}x the "
            f"largest within-box gap (minimum {MIN_GAP_RATIO}); the two PAP-pen "
            f"boxes are not cleanly separated and need a human look"
        )

    return {
        "near_label": sorted(near),
        "far_label": sorted(far),
        "gap_ratio": gap_ratio,
        "within_gaps_mm": [g / 1000.0 for g in within],
        "split_gap_mm": split_gap / 1000.0,
    }


# --------------------------------------------------------------------------- #
# Command-line report
# --------------------------------------------------------------------------- #


def _collect_paths(args: Sequence[str]) -> list[str]:
    """Return the ``.vsi`` paths to report on, sorted by tube number."""
    if args:
        paths: list[str] = []
        for arg in args:
            if os.path.isdir(arg):
                paths.extend(
                    os.path.join(arg, n)
                    for n in os.listdir(arg)
                    if n.lower().endswith(".vsi")
                )
            else:
                paths.append(arg)
    else:
        paths = [
            os.path.join(DEFAULT_RAW_DATA_DIR, n)
            for n in os.listdir(DEFAULT_RAW_DATA_DIR)
            if n.lower().endswith(".vsi")
        ]

    def sort_key(p: str):
        m = _FILENAME_TUBE_RE.search(os.path.basename(p))
        return (int(m.group(1)) if m else 0, p)

    return sorted(paths, key=sort_key)


def _format_exposures(series: Iterable[SeriesMeta]) -> str:
    """Render the per-slide exposure triplet, or note that it varies."""
    triplets = {
        tuple(round(s.exposure_ms.get(c, float("nan")), 2) for c in EXPECTED_CHANNELS)
        for s in series
    }
    if len(triplets) == 1:
        d, f, c = next(iter(triplets))
        return f"{d:8.2f} {f:8.2f} {c:9.2f}"
    return "  VARIES WITHIN SLIDE".ljust(28)


def _main(argv: Sequence[str]) -> int:
    """Print a per-animal metadata and box-assignment table."""
    paths = _collect_paths(argv)
    if not paths:
        print("no .vsi files found", file=sys.stderr)
        return 1

    header = (
        f"{'tube':>4}  {'file':>4}  {'n':>1}  {'stage order':<16}  "
        f"{'near_label':<12} {'far_label':<12} {'ratio':>6}  "
        f"{'DAPI':>8} {'FITC':>8} {'Cy3':>9}  px_um"
    )
    print(header)
    print("-" * len(header))

    problems: list[str] = []
    for path in paths:
        name = os.path.basename(path)
        try:
            meta = read_vsi_meta(path)
        except VsiParseError as exc:
            print(f"{'':>4}  {name:<20} UNREADABLE: {exc}")
            problems.append(f"{name}: {exc}")
            continue

        order = " ".join(s.section_label for s in meta.series)
        try:
            boxes = assign_boxes(meta.series)
            near = ",".join(boxes["near_label"])
            far = ",".join(boxes["far_label"])
            ratio = f"{boxes['gap_ratio']:6.2f}"
            if boxes["gap_ratio"] < MARGINAL_GAP_RATIO:
                ratio += "*"
        except BoxAssignmentError as exc:
            near, far, ratio = "-", "-", "  n/a"
            problems.append(f"{name}: box assignment failed: {exc}")

        pixel_sizes = [s.pixel_size_um for s in meta.series if s.pixel_size_um]
        px = f"{min(pixel_sizes):.5f}-{max(pixel_sizes):.5f}" if pixel_sizes else "-"

        print(
            f"{str(meta.tube_id):>4}  {meta.tube_id_in_file or '-':>4}  "
            f"{meta.n_tissue_series:>1}  {order:<16}  {near:<12} {far:<12} "
            f"{ratio:>7}  {_format_exposures(meta.series)}  {px}"
        )
        for warning in meta.warnings:
            problems.append(f"{name}: {warning}")

    print()
    print("* = between-box gap under " f"{MARGINAL_GAP_RATIO}x the within-box gap")
    print("near_label = low stage X (label end of the slide)")
    print(
        "NOTE: box membership comes from stage coordinates; which box is the "
        "positive one comes from config/slides.csv, never from these numbers."
    )
    if problems:
        print(f"\n{len(problems)} warning(s):")
        for problem in problems:
            print(f"  - {problem}")
    else:
        print("\nno warnings")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
