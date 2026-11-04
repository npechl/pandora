"""Flat per-entry provenance rows, for writing as tables."""

from __future__ import annotations

from pathlib import Path

from pandora.schemas.canonicalisation import canonicalisationProvenance
from pandora.schemas.dataset import CurationProvenance
from pandora.schemas.metadata import MetadataRecord
from pandora.schemas.provenance import (
    CanonicalisationEntryRecord,
    CurationEntryRecord,
)
from pandora.schemas.structure import Structure


def _path(path: str | Path | None) -> str | None:
    """A path as a string, or None."""

    return None if path is None else str(path)


def canonicalisation_entry_record(
    structure: Structure,
    provenance: canonicalisationProvenance,
    structure_path: str | Path | None = None,
) -> CanonicalisationEntryRecord:
    """Flatten one entry's canonicalisation provenance into a table row.

    Args:
        structure: The canonical structure.
        provenance: The provenance `canonicalise_structure` returned
            for it.
        structure_path: Where the canonical structure was written, if
            it was.

    Returns:
        The entry's `CanonicalisationEntryRecord`; `write_records` can
        write a list of them as JSONL or Parquet.
    """

    return CanonicalisationEntryRecord(
        entry_id=structure.entry_id,
        structure_path=_path(structure_path),
        policy_id=provenance.policy_id,
        policy_name=provenance.policy_name,
        policy_version=provenance.policy_version,
        canonicalised_at=provenance.canonicalised_at,
        transforms=provenance.transforms,
        dropped_raw_categories=provenance.dropped_raw_categories,
        warnings=provenance.report.get("warnings"),
        errors=provenance.report.get("errors"),
    )


def curation_entry_record(
    structure: Structure,
    metadata: MetadataRecord | None,
    provenance: CurationProvenance,
    structure_path: str | Path | None = None,
) -> CurationEntryRecord:
    """Flatten one kept entry's curation provenance into a table row.

    Args:
        structure: The curated structure `curate_structure` returned.
        metadata: The entry's metadata, for method and resolution.
        provenance: The provenance `curate_structure` returned.
        structure_path: Where the curated structure was written, if it
            was.

    Returns:
        The entry's `CurationEntryRecord`.
    """

    quality = metadata.quality if metadata else None
    polymer_entities = {e.id for e in structure.entities if e.type == "polymer"}
    return CurationEntryRecord(
        entry_id=structure.entry_id,
        structure_path=_path(structure_path),
        experimental_method=quality.experimental_method if quality else None,
        resolution=quality.resolution if quality else None,
        chain_count=sum(
            1 for a in structure.asym_units if a.entity_id in polymer_entities
        ),
        policy_id=provenance.policy_id,
        policy_name=provenance.policy_name,
        policy_version=provenance.policy_version,
        curated_at=provenance.curated_at,
        dropped_raw_categories=provenance.dropped_raw_categories,
    )
