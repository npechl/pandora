from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field

SimilarityType = Literal[
    "sequence_similarity", "structure_similarity", "custom"
]


class SimilarityMethod(BaseModel):
    engine: str
    version: str | None = None
    parameters: dict[str, Any] | None = None


class SimilarityRelationshipProvenance(BaseModel):
    computed_at: str | None = None
    source_dataset_id: str | None = None


class SimilarityRelationship(BaseModel):
    source_id: str
    target_id: str
    # source_id < target_id lexicographically, one record per unordered pair.
    similarity_type: SimilarityType
    score: float
    coverage: float | None = None
    identity: float | None = None
    method: SimilarityMethod
    provenance: SimilarityRelationshipProvenance = Field(
        default_factory=SimilarityRelationshipProvenance
    )


class SimilarityCluster(BaseModel):
    components: list[str]
    n_components: int


class ClusteringProvenance(BaseModel):
    clustered_at: str
    threshold: float
    n_relationships: int
    n_clusters: int


class PartitionProvenance(BaseModel):
    partitioned_at: str
    pct_train: float
    pct_val: float
    pct_test: float
    keep_similar_items: bool
    split_sizes: dict[str, int] = Field(default_factory=dict)
