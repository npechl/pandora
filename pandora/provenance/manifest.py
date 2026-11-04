from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import gemmi
from pydantic import BaseModel

from pandora import __version__
from pandora._util import now_iso
from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.canonicalisation import (
    canonicalisationPolicy,
    canonicalisationProvenance,
)
from pandora.schemas.dataset import (
    CurationProvenance,
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
    StepSummary,
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
    curation: CurationProvenance | None = None,
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
        curation: Provenance from `curate_structure`, if run.
        metadata: The record from `collect_metadata`, if run.
        annotations: Any `AnnotationLayer`s produced by `annotate_*`
            functions.

    Returns:
        A `ProvenanceBundle` for `structure`.
    """

    return ProvenanceBundle(
        entry_id=structure.entry_id,
        pandora_version=__version__,
        generated_at=now_iso(),
        ingestion=ingestion,
        canonicalisation=canonicalisation,
        curation=curation,
        metadata_sources=(
            collect_metadata_provenance(metadata) if metadata else []
        ),
        annotations=[
            AnnotationProvenanceRecord(
                layer_name=layer.layer_name,
                layer_type=layer.layer_type,
                method=layer.method,
                target_ids=layer.target_ids,
                parameters=layer.parameters,
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
    canonicalisation_policy: canonicalisationPolicy | None = None,
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
        canonicalisation_policy: The policy passed to
            `canonicalise_structure`, if run.
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
        generated_at=now_iso(),
        curation_policy=curation_policy,
        canonicalisation_policy=canonicalisation_policy,
        excluded=list(excluded),
        deduplication=deduplication,
        clustering=clustering,
        partition=partition,
        splits=splits or {},
        structures=list(structures),
    )


def build_step_summary(
    step: str,
    policy: BaseModel,
    *,
    input: str | Path,
    output: str | Path,
    started_at: str,
    counts: Mapping[str, int] | None = None,
    previous_summary: str | Path | None = None,
) -> StepSummary:
    """Build the summary of one step run, stamped with tool versions.

    Call it when the step has finished; `finished_at` is set to now.
    Write it with `pandora.export.write_json`.

    Args:
        step: The step's name, e.g. "canonicalisation".
        policy: The policy the step ran with; embedded in full.
        input: What the step read: a path list file or a folder.
        output: The folder the step wrote to.
        started_at: When the step started, as an ISO 8601 timestamp
            (e.g. from `pandora._util.now_iso()` at the start).
        counts: Named counts to record, e.g. `{"entries_in": 100}`.
        previous_summary: The summary file of the step whose output
            this step read, if any.

    Returns:
        The `StepSummary`.
    """

    return StepSummary(
        step=step,
        pandora_version=__version__,
        gemmi_version=gemmi.__version__,
        policy=policy.model_dump(mode="json"),
        input=str(input),
        previous_summary=(
            None if previous_summary is None else str(previous_summary)
        ),
        output=str(output),
        started_at=started_at,
        finished_at=now_iso(),
        counts=dict(counts or {}),
    )
