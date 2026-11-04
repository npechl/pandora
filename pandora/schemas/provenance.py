from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

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
from pandora.schemas.metadata import MetadataProvenance
from pandora.schemas.similarity import ClusteringProvenance, PartitionProvenance


class AnnotationProvenanceRecord(BaseModel):
    """Record of one annotation layer's method and parameters, for
    reproducing it later.

    Attributes:
        layer_name: The annotation layer's human-readable name.
        layer_type: The annotation layer's machine-readable type
            (e.g. "ligand_contacts").
        method: The identifier of the method/algorithm used to
            compute the layer.
        target_ids: The entry id(s) the layer applies to.
        parameters: The parameters the layer was computed with.
        provenance: Additional provenance details about how the layer
            was computed.
    """

    layer_name: str
    layer_type: str
    method: str
    target_ids: list[str] = Field(default_factory=list)
    parameters: dict[str, Any] = Field(default_factory=dict)
    provenance: dict[str, Any] = Field(default_factory=dict)


class ProvenanceBundle(BaseModel):
    """Every provenance collected for one structure: ingestion,
    canonicalisation, metadata sources, and annotations.

    Attributes:
        entry_id: The structure's entry id.
        pandora_version: The Pandora version that produced this
            bundle.
        generated_at: When this bundle was assembled, as an ISO 8601
            timestamp.
        ingestion: How the structure's raw file was fetched, if known.
        canonicalisation: How the structure was canonicalised, if
            known.
        curation: How the structure was curated, if known.
        metadata_sources: Which sources the structure's metadata came
            from.
        annotations: Every annotation layer computed for this
            structure.
    """

    entry_id: str
    pandora_version: str
    generated_at: str
    ingestion: IngestionProvenance | None = None
    canonicalisation: canonicalisationProvenance | None = None
    curation: CurationProvenance | None = None
    metadata_sources: list[MetadataProvenance] = Field(default_factory=list)
    annotations: list[AnnotationProvenanceRecord] = Field(default_factory=list)


class DatasetManifest(BaseModel):
    """A single-file report of how a dataset was built: policies,
    exclusions, dedup/clustering/partition provenance, splits, and every
    retained structure's `ProvenanceBundle`.

    Attributes:
        dataset_id: A unique identifier for the dataset.
        dataset_name: A human-readable name for the dataset.
        dataset_version: The dataset's version string.
        pandora_version: The Pandora version that produced this
            manifest.
        generated_at: When this manifest was assembled, as an ISO
            8601 timestamp.
        curation_policy: The curation policy applied, if any.
        canonicalisation_policy: The canonicalisation policy applied,
            if any.
        excluded: Every structure excluded during curation or
            deduplication.
        deduplication: How deduplication was applied, if any.
        clustering: How similarity clustering was applied, if any.
        partition: How the train/val/test split was computed, if any.
        splits: The resulting split assignment, keyed by split name.
        structures: The `ProvenanceBundle` for every retained
            structure.
    """

    dataset_id: str
    dataset_name: str
    dataset_version: str
    pandora_version: str
    generated_at: str
    curation_policy: DatasetCurationPolicy | None = None
    canonicalisation_policy: canonicalisationPolicy | None = None
    excluded: list[ExclusionRecord] = Field(default_factory=list)
    deduplication: DeduplicationProvenance | None = None
    clustering: ClusteringProvenance | None = None
    partition: PartitionProvenance | None = None
    splits: dict[str, list[str]] = Field(default_factory=dict)
    structures: list[ProvenanceBundle] = Field(default_factory=list)


class StepSummary(BaseModel):
    """The recipe of one dataset-building step run (e.g. canonicalisation
    or curation over a list of entries): what it ran with, what it read,
    what it wrote, and how many entries went where.

    Attributes:
        step: The step's name, e.g. "canonicalisation".
        pandora_version: The Pandora version that ran the step.
        gemmi_version: The gemmi version that ran the step.
        policy: The full policy the step ran with, as plain data.
        input: What the step read: a path list file or a folder.
        previous_summary: The summary file of the step whose output this
            step read, if any.
        output: The folder the step wrote to.
        started_at: When the step started, as an ISO 8601 timestamp.
        finished_at: When the step finished, as an ISO 8601 timestamp.
        counts: Named counts, e.g. entries in, kept, failed, and
            excluded per reason code.
    """

    step: str
    pandora_version: str
    gemmi_version: str
    policy: dict[str, Any]
    input: str
    previous_summary: str | None = None
    output: str
    started_at: str
    finished_at: str
    counts: dict[str, int] = Field(default_factory=dict)


class CanonicalisationEntryRecord(BaseModel):
    """One canonicalised entry as a flat table row: where its structure
    was written and how it was canonicalised.

    Attributes:
        entry_id: The structure's entry id.
        structure_path: Where the canonical structure was written, if it
            was.
        policy_id: The id of the canonicalisation policy that ran.
        policy_name: The name of that policy.
        policy_version: The version of that policy.
        canonicalised_at: When canonicalisation ran, as an ISO 8601
            timestamp.
        transforms: The transform labels applied.
        dropped_raw_categories: The verbatim mmCIF categories removed.
        warnings: The number of warnings, if a report was recorded.
        errors: The number of errors, if a report was recorded.
    """

    entry_id: str
    structure_path: str | None = None
    policy_id: str
    policy_name: str
    policy_version: str
    canonicalised_at: str
    transforms: list[str] = Field(default_factory=list)
    dropped_raw_categories: list[str] = Field(default_factory=list)
    warnings: int | None = None
    errors: int | None = None


class CurationEntryRecord(BaseModel):
    """One curated (kept) entry as a flat table row: where its structure
    was written, its quality, and how it was curated.

    Attributes:
        entry_id: The structure's entry id.
        structure_path: Where the curated structure was written, if it
            was.
        experimental_method: The entry's experimental method, if known.
        resolution: The entry's resolution in angstroms, if known.
        chain_count: The number of polymer chains kept.
        policy_id: The id of the curation policy that ran.
        policy_name: The name of that policy.
        policy_version: The version of that policy.
        curated_at: When curation ran, as an ISO 8601 timestamp.
        dropped_raw_categories: The verbatim mmCIF categories removed.
    """

    entry_id: str
    structure_path: str | None = None
    experimental_method: str | None = None
    resolution: float | None = None
    chain_count: int
    policy_id: str
    policy_name: str
    policy_version: str
    curated_at: str
    dropped_raw_categories: list[str] = Field(default_factory=list)
