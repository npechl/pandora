from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

from pandora.schemas.canonicalisation import (
    canonicalisationPolicy,
    canonicalisationProvenance,
)
from pandora.schemas.dataset import (
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
