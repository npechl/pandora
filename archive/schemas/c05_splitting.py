from __future__ import annotations

from typing import Annotated, Any, Literal, Union

from pydantic import BaseModel, Discriminator, Field, Tag

from .common import AppliedPolicyRef, Diagnostic, DiagnosticBundle, ResultStatus
from .c04_curation import (
    ChainDataset,
    Dataset,
    InterfaceDataset,
    ResidueDataset,
)


# ── PandoraDataset discriminated union ────────────────────────────────────────

PandoraDataset = Annotated[
    Union[
        Annotated[Dataset, Tag("structure")],
        Annotated[ChainDataset, Tag("chain")],
        Annotated[InterfaceDataset, Tag("interface")],
        Annotated[ResidueDataset, Tag("residue")],
    ],
    Discriminator("granularity"),
]


# ── Similarity relationship ────────────────────────────────────────────────────

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
    similarity_type: SimilarityType
    score: float
    coverage: float | None = None
    identity: float | None = None
    method: SimilarityMethod
    provenance: SimilarityRelationshipProvenance = Field(
        default_factory=SimilarityRelationshipProvenance
    )


# ── Similarity network ─────────────────────────────────────────────────────────


class SimilarityEdge(BaseModel):
    source_id: str
    target_id: str
    score: float
    similarity_type: SimilarityType


class GraphStatistics(BaseModel):
    node_count: int
    edge_count: int
    connected_components: int | None = None
    largest_component_size: int | None = None


class SimilarityNetworkProvenance(BaseModel):
    built_at: str | None = None


class SimilarityNetwork(BaseModel):
    network_id: str
    dataset_id: str
    relationships: list[SimilarityRelationship] = Field(default_factory=list)
    nodes: list[str] = Field(default_factory=list)
    edges: list[SimilarityEdge] = Field(default_factory=list)
    graph_statistics: GraphStatistics
    applied_policy: AppliedPolicyRef
    provenance: SimilarityNetworkProvenance = Field(
        default_factory=SimilarityNetworkProvenance
    )


# ── Clusters ───────────────────────────────────────────────────────────────────


class Cluster(BaseModel):
    cluster_id: str
    members: list[str]
    cluster_size: int
    representative_id: str | None = None


class ClusteringSummary(BaseModel):
    total_clusters: int
    singleton_clusters: int
    largest_cluster_size: int
    mean_cluster_size: float


class SimilarityClustersProvenance(BaseModel):
    clustered_at: str | None = None


class SimilarityClusters(BaseModel):
    clustering_id: str
    dataset_id: str
    clusters: list[Cluster] = Field(default_factory=list)
    clustering_summary: ClusteringSummary
    applied_policy: AppliedPolicyRef
    provenance: SimilarityClustersProvenance = Field(
        default_factory=SimilarityClustersProvenance
    )


# ── Leakage-safe dataset ───────────────────────────────────────────────────────


class Partitions(BaseModel):
    train: list[str] = Field(default_factory=list)
    validation: list[str] = Field(default_factory=list)
    test: list[str] = Field(default_factory=list)


class PartitionSummary(BaseModel):
    train_count: int
    validation_count: int
    test_count: int
    train_fraction_achieved: float
    validation_fraction_achieved: float
    test_fraction_achieved: float


class LeakageSummary(BaseModel):
    max_cross_split_similarity: float | None = None
    leakage_detected: bool = False
    leakage_diagnostics: list[Diagnostic] = Field(default_factory=list)


class LeakageSafeDatasetProvenance(BaseModel):
    split_at: str | None = None
    similarity_engines: list[str] = Field(default_factory=list)


class LeakageSafeDataset(BaseModel):
    dataset_id: str
    dataset_name: str
    dataset_version: str
    granularity: Literal["structure", "chain", "interface", "residue"]
    source_dataset: PandoraDataset
    similarity_network: SimilarityNetwork
    similarity_clusters: SimilarityClusters
    partitions: Partitions = Field(default_factory=Partitions)
    partition_summary: PartitionSummary
    leakage_summary: LeakageSummary = Field(default_factory=LeakageSummary)
    applied_policy: AppliedPolicyRef
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: LeakageSafeDatasetProvenance = Field(
        default_factory=LeakageSafeDatasetProvenance
    )


class LeakageAnalysisBatchSummary(BaseModel):
    total: int
    success: int
    warning: int
    failed: int


class LeakageAnalysisBatchResultItem(BaseModel):
    dataset_id: str
    status: ResultStatus
    leakage_safe_dataset: LeakageSafeDataset | None = None
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)


class LeakageAnalysisBatchResult(BaseModel):
    mode: Literal["sequential", "parallel"]
    summary: LeakageAnalysisBatchSummary
    results: list[LeakageAnalysisBatchResultItem]


# ── Policy schema ──────────────────────────────────────────────────────────────


class SequenceSimilarityRules(BaseModel):
    enabled: bool = True
    engine: str = "MMseqs2"
    threshold: float | None = 0.3
    coverage_threshold: float | None = None


class StructureSimilarityRules(BaseModel):
    enabled: bool = False
    engine: str = "Foldseek"
    threshold: float | None = 0.5
    coverage_threshold: float | None = None


class SimilarityRules(BaseModel):
    sequence_similarity: SequenceSimilarityRules = Field(
        default_factory=SequenceSimilarityRules
    )
    structure_similarity: StructureSimilarityRules = Field(
        default_factory=StructureSimilarityRules
    )


class ClusteringRules(BaseModel):
    enabled: bool = True
    strategy: Literal["connected_components", "single_linkage", "custom"] = (
        "connected_components"
    )
    min_cluster_size: int | None = None


class PartitionRules(BaseModel):
    train_fraction: float = 0.8
    validation_fraction: float = 0.1
    test_fraction: float = 0.1
    keep_similar_items_together: bool = True
    stratify_by_cluster: bool = False


class LeakageRules(BaseModel):
    forbid_cross_split_similarity: bool = False
    max_allowed_cross_split_similarity: float | None = None
    enforce_cluster_isolation: bool = True


class LeakageProvenanceRules(BaseModel):
    record_similarity_method: bool = True
    record_engine_versions: bool = True
    record_thresholds: bool = True
    record_partition_history: bool = True


class LeakagePolicy(BaseModel):
    policy_id: str
    policy_name: str
    policy_version: str
    description: str = ""
    similarity_rules: SimilarityRules = Field(default_factory=SimilarityRules)
    clustering_rules: ClusteringRules = Field(default_factory=ClusteringRules)
    partition_rules: PartitionRules = Field(default_factory=PartitionRules)
    leakage_rules: LeakageRules = Field(default_factory=LeakageRules)
    provenance_rules: LeakageProvenanceRules = Field(
        default_factory=LeakageProvenanceRules
    )
