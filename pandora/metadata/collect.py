from __future__ import annotations

from pandora.metadata.mmcif import (
    extract_entity_metadata,
    extract_entry_metadata,
    extract_ligand_metadata,
    extract_quality,
    extract_taxonomies,
    extract_uniprot_mappings,
)
from pandora.schemas.metadata import MetadataRecord
from pandora.schemas.structure import Structure


def collect_metadata(structure: Structure) -> MetadataRecord:
    """Collect the basic source-backed metadata supported by Pandora today."""

    return MetadataRecord(
        entry_id=structure.entry_id,
        entry=extract_entry_metadata(structure),
        quality=extract_quality(structure),
        taxonomies=extract_taxonomies(structure),
        entities=extract_entity_metadata(structure),
        ligands=extract_ligand_metadata(structure),
        uniprot_mappings=extract_uniprot_mappings(structure),
        raw_categories=sorted(structure.raw.keys()),
    )
