from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

from pandora.schemas.canonicalisation import canonicalisationProvenance
from pandora.schemas.ingestion import IngestionProvenance
from pandora.schemas.metadata import MetadataProvenance


class AnnotationProvenanceRecord(BaseModel):
    layer_name: str
    layer_type: str
    method: str
    provenance: dict[str, Any] = Field(default_factory=dict)


class ProvenanceBundle(BaseModel):
    entry_id: str
    pandora_version: str
    generated_at: str
    ingestion: IngestionProvenance | None = None
    canonicalisation: canonicalisationProvenance | None = None
    metadata_sources: list[MetadataProvenance] = Field(default_factory=list)
    annotations: list[AnnotationProvenanceRecord] = Field(default_factory=list)
