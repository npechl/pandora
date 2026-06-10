"""Component 05 — Similarity & Leakage-Safe Splitting: public functions."""
from __future__ import annotations

from datetime import datetime, timezone

from pandora.schemas.c04_curation import (
    ChainDataset,
    Dataset,
    InterfaceDataset,
    ResidueDataset,
)
from pandora.schemas.c05_splitting import (
    Cluster,
    ClusteringSummary,
    GraphStatistics,
    LeakagePolicy,
    LeakageSafeDataset,
    LeakageSafeDatasetProvenance,
    LeakageSummary,
    PandoraDataset,
    PartitionSummary,
    Partitions,
    SimilarityClusters,
    SimilarityClustersProvenance,
    SimilarityNetwork,
    SimilarityNetworkProvenance,
    SimilarityRelationship,
)
from pandora.schemas.common import AppliedPolicyRef


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _policy_ref(policy: LeakagePolicy) -> AppliedPolicyRef:
    return AppliedPolicyRef(
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
    )


def _get_item_ids(dataset: PandoraDataset) -> list[str]:
    """Return a flat list of item IDs for any dataset granularity."""
    if isinstance(dataset, Dataset):
        return [
            s.canonical_structure_result.entry_id for s in dataset.structures
        ]
    if isinstance(dataset, ChainDataset):
        return [f"{c.entry_id}:{c.chain_id}" for c in dataset.chains]
    if isinstance(dataset, InterfaceDataset):
        return [i.interface_id for i in dataset.interfaces]
    return [r.residue_id for r in dataset.residues]  # type: ignore[union-attr]


# ── Public API ────────────────────────────────────────────────────────────────

def compute_similarity_relationships(
    dataset: PandoraDataset,
    policy: LeakagePolicy,
) -> list[SimilarityRelationship]:
    # TODO: implement — MMseqs2 subprocess wrapper for sequence similarity,
    #   Foldseek subprocess wrapper for structure similarity,
    #   per policy.similarity_rules
    return []


def build_similarity_network(
    dataset: PandoraDataset,
    relationships: list[SimilarityRelationship],
    policy: LeakagePolicy,
) -> SimilarityNetwork:
    """Construct a similarity graph from pairwise relationship records."""
    # TODO: implement graph construction; compute connected-components statistics
    nodes = _get_item_ids(dataset)
    return SimilarityNetwork(
        network_id=f"{dataset.dataset_id}:network",
        dataset_id=dataset.dataset_id,
        relationships=relationships,
        nodes=nodes,
        edges=[],
        graph_statistics=GraphStatistics(
            node_count=len(nodes),
            edge_count=0,
        ),
        applied_policy=_policy_ref(policy),
        provenance=SimilarityNetworkProvenance(built_at=_now_iso()),
    )


def cluster_by_similarity(
    network: SimilarityNetwork,
    policy: LeakagePolicy,
) -> SimilarityClusters:
    """Cluster nodes by similarity per policy.clustering_rules."""
    # TODO: implement connected_components, single_linkage, custom strategies
    clusters = [
        Cluster(cluster_id=f"c{i}", members=[node], cluster_size=1)
        for i, node in enumerate(network.nodes)
    ]
    n = len(clusters)
    return SimilarityClusters(
        clustering_id=f"{network.network_id}:clusters",
        dataset_id=network.dataset_id,
        clusters=clusters,
        clustering_summary=ClusteringSummary(
            total_clusters=n,
            singleton_clusters=n,
            largest_cluster_size=1 if n else 0,
            mean_cluster_size=1.0 if n else 0.0,
        ),
        applied_policy=_policy_ref(policy),
        provenance=SimilarityClustersProvenance(clustered_at=_now_iso()),
    )


def create_leakage_safe_dataset(
    dataset: PandoraDataset,
    policy: LeakagePolicy,
) -> LeakageSafeDataset:
    """Compute similarity, cluster, and partition a dataset into leakage-safe splits."""
    # TODO: implement cluster-aware assignment to guarantee no cross-split
    #   leakage per policy.leakage_rules
    relationships = compute_similarity_relationships(dataset, policy)
    network = build_similarity_network(dataset, relationships, policy)
    clusters = cluster_by_similarity(network, policy)

    item_ids = _get_item_ids(dataset)
    n = len(item_ids)
    train_end = int(n * policy.partition_rules.train_fraction)
    val_end = train_end + int(n * policy.partition_rules.validation_fraction)

    train = item_ids[:train_end]
    val = item_ids[train_end:val_end]
    test = item_ids[val_end:]

    def _frac(count: int) -> float:
        return count / n if n else 0.0

    return LeakageSafeDataset(
        dataset_id=f"{dataset.dataset_id}:split",
        dataset_name=getattr(dataset, "dataset_name", dataset.dataset_id),
        dataset_version=dataset.dataset_version,
        granularity=dataset.granularity,
        source_dataset=dataset,
        similarity_network=network,
        similarity_clusters=clusters,
        partitions=Partitions(train=train, validation=val, test=test),
        partition_summary=PartitionSummary(
            train_count=len(train),
            validation_count=len(val),
            test_count=len(test),
            train_fraction_achieved=_frac(len(train)),
            validation_fraction_achieved=_frac(len(val)),
            test_fraction_achieved=_frac(len(test)),
        ),
        leakage_summary=LeakageSummary(),
        applied_policy=_policy_ref(policy),
        provenance=LeakageSafeDatasetProvenance(split_at=_now_iso()),
    )
