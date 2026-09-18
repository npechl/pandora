from __future__ import annotations

from pandora._util import now_iso
from pandora.schemas.similarity import (
    ClusteringProvenance,
    SimilarityCluster,
    SimilarityRelationship,
)


def cluster_similar_items(
    item_ids: list[str],
    relationships: list[SimilarityRelationship],
    threshold: float,
) -> tuple[list[SimilarityCluster], ClusteringProvenance]:
    """Group item ids into clusters via connected components.

    Two items land in the same cluster iff connected through a chain of
    relationships each scoring >= threshold. Items with no such edges
    (isolates) form their own singleton cluster, so every id in
    `item_ids` ends up in exactly one cluster.

    Args:
        item_ids: Every item id to place into a cluster.
        relationships: Pairwise similarity relationships between items,
            as returned by `compute_sequence_similarity`/
            `compute_structure_similarity`. Must all share the same
            `.method.engine`.
        threshold: Minimum `SimilarityRelationship.score` for an edge
            to count as a connection between two items.

    Returns:
        `(clusters, provenance)` — the clusters, and a record of the
        threshold applied, how many relationships/clusters resulted, and
        the `SimilarityMethod` of the first relationship — enough to
        reproduce this clustering step given the same input structures.

    Raises:
        ValueError: `relationships` mixes more than one similarity
            engine; clustering assumes a single engine/parameters
            produced the whole network.
    """

    if len({rel.method.engine for rel in relationships}) > 1:
        raise ValueError(
            "cluster_similar_items: relationships use more than one "
            "similarity engine; clustering assumes a single "
            "engine/parameters for the whole network"
        )

    parent = {item_id: item_id for item_id in item_ids}

    def find(item_id: str) -> str:
        """Union-find: root id of item_id's cluster, with path compression."""

        while parent[item_id] != item_id:
            parent[item_id] = parent[parent[item_id]]
            item_id = parent[item_id]
        return item_id

    for rel in relationships:
        if (
            rel.score >= threshold
            and rel.source_id in parent
            and rel.target_id in parent
        ):
            root_source, root_target = find(rel.source_id), find(rel.target_id)
            if root_source != root_target:
                parent[root_target] = root_source

    groups: dict[str, list[str]] = {}
    for item_id in item_ids:
        groups.setdefault(find(item_id), []).append(item_id)

    clusters = [
        SimilarityCluster(components=sorted(members), n_components=len(members))
        for members in sorted(groups.values(), key=lambda members: members[0])
    ]
    provenance = ClusteringProvenance(
        clustered_at=now_iso(),
        threshold=threshold,
        n_relationships=len(relationships),
        n_clusters=len(clusters),
        similarity_method=relationships[0].method if relationships else None,
    )
    return clusters, provenance


def pair_cluster_keys(
    pairs: list[tuple[str, str]],
    clusters: list[SimilarityCluster],
) -> dict[tuple[str, str], tuple[str, str]]:
    """Derive a per-side cluster key for each PPI pair from item-level
    clusters (e.g. Pinder-style `{cluster_id_R, cluster_id_L}`).

    Each cluster's own lexicographically-smallest component id is used
    as that cluster's key, so the same item always maps to the same key
    regardless of clustering order.

    Args:
        pairs: `(item_id_1, item_id_2)` tuples — e.g. the two chain ids
            of a PPI pair. Both ids must appear in `clusters`.
        clusters: Item-level clusters, as returned by
            `cluster_similar_items()`.

    Returns:
        `{(item_id_1, item_id_2): (cluster_id_1, cluster_id_2)}`, one
        entry per input pair.

    Raises:
        KeyError: A pair references an item id not present in any
            cluster.
    """

    cluster_of: dict[str, str] = {
        item_id: cluster.components[0]
        for cluster in clusters
        for item_id in cluster.components
    }
    return {
        (item_1, item_2): (cluster_of[item_1], cluster_of[item_2])
        for item_1, item_2 in pairs
    }
