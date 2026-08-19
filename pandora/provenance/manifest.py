from __future__ import annotations

from collections.abc import Sequence
from datetime import datetime, timezone

from pandora import __version__
from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.canonicalisation import canonicalisationProvenance
from pandora.schemas.dataset import (
    DatasetCurationPolicy,
    DeduplicationProvenance,
    ExclusionRecord,
)
from pandora.schemas.ingestion import IngestionProvenance
from pandora.schemas.metadata import MetadataProvenance, MetadataRecord
from pandora.schemas.provenance import (
    AnnotationProvenanceRecord,
    DatasetManifest,
    ProvenanceBundle,
)
from pandora.schemas.similarity import ClusteringProvenance, PartitionProvenance
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


def build_dataset_manifest(
    *,
    dataset_id: str,
    dataset_name: str,
    dataset_version: str,
    curation_policy: DatasetCurationPolicy | None = None,
    excluded: Sequence[ExclusionRecord] = (),
    deduplication: DeduplicationProvenance | None = None,
    clustering: ClusteringProvenance | None = None,
    partition: PartitionProvenance | None = None,
    splits: dict[str, list[str]] | None = None,
    structures: Sequence[ProvenanceBundle] = (),
) -> DatasetManifest:
    """Assemble a single-file report of how a dataset was constructed.

    Aggregates whatever provenance the caller already produced by running
    the dataset-construction stages — this function does not re-derive,
    recompute, or validate anything itself. Every argument beyond the
    dataset identity is optional so the manifest can be built with
    whichever stages actually ran (e.g. curation without clustering).

    Args:
        dataset_id: Stable identifier for this dataset build.
        dataset_name: Human-readable name.
        dataset_version: Version of this dataset build.
        curation_policy: The policy passed to `curate_structure`, if run.
        excluded: `ExclusionRecord`s returned by `curate_structure`/
            `deduplicate_structures` for every structure left out.
        deduplication: Provenance from `deduplicate_structures`, if run.
        clustering: Provenance from `cluster_similar_items`, if run.
        partition: Provenance from `partition_dataset`, if run.
        splits: The split assignment returned by `partition_dataset`,
            if run.
        structures: A `ProvenanceBundle` per structure retained in the
            dataset (from `build_provenance_bundle`).

    Returns:
        A `DatasetManifest` capturing the policy, exclusions, dedup/
        clustering/partition parameters, split assignment, and every
        retained structure's own provenance in one object.
    """

    return DatasetManifest(
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        dataset_version=dataset_version,
        pandora_version=__version__,
        generated_at=datetime.now(timezone.utc).isoformat(),
        curation_policy=curation_policy,
        excluded=list(excluded),
        deduplication=deduplication,
        clustering=clustering,
        partition=partition,
        splits=splits or {},
        structures=list(structures),
    )
