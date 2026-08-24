from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field

SimilarityType = Literal[
    "sequence_similarity", "structure_similarity", "custom"
]


class SimilarityMethod(BaseModel):
    """Which similarity engine (and version/parameters) produced a
    similarity network.

    Attributes:
        engine: The similarity engine used (e.g. "MMseqs2",
            "Foldseek").
        version: The engine's version string, if determined.
        parameters: The parameters the engine was run with.
    """

    engine: str
    version: str | None = None
    parameters: dict[str, Any] | None = None


class SimilarityRelationshipProvenance(BaseModel):
    """When (and from which source dataset) one similarity relationship
    was computed.

    Attributes:
        computed_at: When this relationship was computed, as an ISO
            8601 timestamp.
        source_dataset_id: The dataset this relationship was computed
            for, if known.
    """

    computed_at: str | None = None
    source_dataset_id: str | None = None


class SimilarityRelationship(BaseModel):
    """One pairwise similarity score between two items
    (source_id < target_id).

    Attributes:
        source_id: The lexicographically smaller of the two item ids.
        target_id: The lexicographically larger of the two item ids.
        similarity_type: What kind of similarity this measures.
        score: The relationship's primary similarity score.
        coverage: The alignment coverage between the two items, if
            reported.
        identity: The sequence/structural identity between the two
            items, if reported.
        method: Which engine/parameters produced this relationship.
        provenance: When and for which dataset this relationship was
            computed.
    """

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
    """One connected-component cluster of similar items.

    Attributes:
        components: Every item id in this cluster.
        n_components: The number of items in this cluster.
    """

    components: list[str]
    n_components: int


class ClusteringProvenance(BaseModel):
    """Record of one clustering run: threshold, relationship/cluster
    counts, and the similarity method used.

    Attributes:
        clustered_at: When this clustering run completed, as an ISO
            8601 timestamp.
        threshold: The similarity score threshold used to connect
            items.
        n_relationships: How many relationships were considered.
        n_clusters: How many clusters resulted.
        similarity_method: Which engine/parameters produced the
            clustered relationships.
    """

    clustered_at: str
    threshold: float
    n_relationships: int
    n_clusters: int
    similarity_method: SimilarityMethod | None = None


class PartitionProvenance(BaseModel):
    """Record of one train/val/test partition run: target fractions,
    mode, and resulting split sizes.

    Attributes:
        partitioned_at: When this partition run completed, as an ISO
            8601 timestamp.
        pct_train: The target fraction of items assigned to train.
        pct_val: The target fraction of items assigned to val.
        pct_test: The target fraction of items assigned to test.
        keep_similar_items: Whether each cluster was kept whole in
            one split (leakage-safe) or divided proportionally.
        split_sizes: The resulting number of items per split.
    """

    partitioned_at: str
    pct_train: float
    pct_val: float
    pct_test: float
    keep_similar_items: bool
    split_sizes: dict[str, int] = Field(default_factory=dict)
