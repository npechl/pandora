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
    layer_name: str
    layer_type: str
    method: str
    target_ids: list[str] = Field(default_factory=list)
    parameters: dict[str, Any] = Field(default_factory=dict)
    provenance: dict[str, Any] = Field(default_factory=dict)


class ProvenanceBundle(BaseModel):
    entry_id: str
    pandora_version: str
    generated_at: str
    ingestion: IngestionProvenance | None = None
    canonicalisation: canonicalisationProvenance | None = None
    metadata_sources: list[MetadataProvenance] = Field(default_factory=list)
    annotations: list[AnnotationProvenanceRecord] = Field(default_factory=list)


class DatasetManifest(BaseModel):
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
