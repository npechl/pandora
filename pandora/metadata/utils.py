from __future__ import annotations

from pandora.schemas.metadata import MetadataProvenance
from pandora.schemas.structure import Structure


RawRow = dict[str, str | None]

MISSING_VALUES = frozenset({"", ".", "?"})
WATER_COMP_IDS = frozenset({"HOH", "WAT", "DOD"})


def clean(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    return None if stripped in MISSING_VALUES else stripped


def as_int(value: str | int | None) -> int | None:
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
    for value in values:
        cleaned = clean(value)
        if cleaned is not None:
            return cleaned
    return None


def clean_sequence(value: str | None) -> str | None:
    cleaned = clean(value)
    if cleaned is None:
        return None
    return "".join(cleaned.split())


def normalise_category(category: str) -> str:
    return category if category.startswith("_") else f"_{category}"


def raw_rows(structure: Structure, category: str) -> list[RawRow]:
    return structure.raw.get(normalise_category(category), [])


def first_row(structure: Structure, category: str) -> RawRow:
    rows = raw_rows(structure, category)
    return rows[0] if rows else {}


def provenance(
    category: str,
    record_id: str | None = None,
) -> MetadataProvenance:
    return MetadataProvenance(
        source="mmcif",
        source_category=normalise_category(category),
        source_record_id=record_id,
    )
