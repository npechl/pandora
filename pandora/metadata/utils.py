from __future__ import annotations

from pandora.schemas.metadata import MetadataProvenance
from pandora.schemas.structure import Structure

RawRow = dict[str, str | None]

MISSING_VALUES = frozenset({"", ".", "?"})
WATER_COMP_IDS = frozenset({"HOH", "WAT", "DOD"})


def clean(value: str | None) -> str | None:
    """Strip value and normalize mmCIF null tokens ("", ".", "?") to None."""

    if value is None:
        return None
    stripped = value.strip()
    return None if stripped in MISSING_VALUES else stripped


def as_int(value: str | int | None) -> int | None:
    """Parse value as int, tolerating mmCIF null tokens; None if unparseable."""

    if isinstance(value, int):
        return value
    cleaned = clean(value)
    if cleaned is None:
        return None
    try:
        return int(cleaned)
    except ValueError:
        return None


def as_float(value: str | float | None) -> float | None:
    """Parse value as float, tolerating mmCIF null tokens; None if unparseable."""

    if isinstance(value, float):
        return value
    cleaned = clean(value)
    if cleaned is None:
        return None
    try:
        return float(cleaned)
    except ValueError:
        return None


def first_value(*values: str | None) -> str | None:
    """First cleaned, non-null value among values."""

    for value in values:
        cleaned = clean(value)
        if cleaned is not None:
            return cleaned
    return None


def clean_sequence(value: str | None) -> str | None:
    """clean() a sequence value and strip internal whitespace."""

    cleaned = clean(value)
    if cleaned is None:
        return None
    return "".join(cleaned.split())


def normalise_category(category: str) -> str:
    """Ensure category has a leading underscore (mmCIF category name form)."""

    return category if category.startswith("_") else f"_{category}"


def raw_rows(structure: Structure, category: str) -> list[RawRow]:
    """Raw mmCIF rows for category from structure.raw, or [] if absent."""

    return structure.raw.get(normalise_category(category), [])


def first_row(structure: Structure, category: str) -> RawRow:
    """First raw row for category, or {} if absent."""

    rows = raw_rows(structure, category)
    return rows[0] if rows else {}


def provenance(
    category: str,
    record_id: str | None = None,
) -> MetadataProvenance:
    """Build a MetadataProvenance stamp for a raw mmCIF source category/record."""

    return MetadataProvenance(
        source="mmcif",
        source_category=normalise_category(category),
        source_record_id=record_id,
    )
