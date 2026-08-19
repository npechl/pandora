from __future__ import annotations

from collections.abc import Sequence
from datetime import datetime, timezone

from pandora import __version__
from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.canonicalisation import canonicalisationProvenance
from pandora.schemas.ingestion import IngestionProvenance
from pandora.schemas.metadata import MetadataProvenance, MetadataRecord
from pandora.schemas.provenance import AnnotationProvenanceRecord, ProvenanceBundle
from pandora.schemas.structure import Structure


def collect_metadata_provenance(
    metadata: MetadataRecord,
) -> list[MetadataProvenance]:
    """Flatten the per-record provenance stamps out of a `MetadataRecord`.

    Args:
        metadata: The collected metadata record.

    Returns:
        One `MetadataProvenance` per populated sub-record (entry, quality,
        each taxonomy/entity/ligand/UniProt mapping), in record order.
    """

    records = [metadata.entry, metadata.quality]
    records += metadata.taxonomies
    records += metadata.entities
    records += metadata.ligands
    records += metadata.uniprot_mappings
    return [record.provenance for record in records if record is not None]


def build_provenance_bundle(
    structure: Structure,
    *,
    ingestion: IngestionProvenance | None = None,
    canonicalisation: canonicalisationProvenance | None = None,
    metadata: MetadataRecord | None = None,
    annotations: Sequence[AnnotationLayer] = (),
) -> ProvenanceBundle:
    """Assemble a `ProvenanceBundle` for one (canonicalised) structure.

    Aggregates whatever provenance the caller already produced by running
    the earlier pipeline stages — this function does not re-derive or
    fetch anything itself. Every argument is optional so the bundle can be
    built at any point in the pipeline (e.g. right after parsing, with
    only `ingestion` set).

    Args:
        structure: The structure to attribute the bundle to.
        ingestion: Provenance from `pandora.ingestion.fetch_mmcif`, if run.
        canonicalisation: Provenance from `canonicalise_structure`, if run.
        metadata: The record from `collect_metadata`, if run.
        annotations: Any `AnnotationLayer`s produced by `annotate_*`
            functions.

    Returns:
        A `ProvenanceBundle` for `structure`.
    """

    return ProvenanceBundle(
        entry_id=structure.entry_id,
        pandora_version=__version__,
        generated_at=datetime.now(timezone.utc).isoformat(),
        ingestion=ingestion,
        canonicalisation=canonicalisation,
        metadata_sources=(
            collect_metadata_provenance(metadata) if metadata else []
        ),
        annotations=[
            AnnotationProvenanceRecord(
                layer_name=layer.layer_name,
                layer_type=layer.layer_type,
                method=layer.method,
                provenance=layer.provenance,
            )
            for layer in annotations
        ],
    )
