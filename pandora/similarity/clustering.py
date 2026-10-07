from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

from pandora._util import now_iso
from pandora.schemas.similarity import (
    ClusteringProvenance,
    HitFilter,
    SimilarityCluster,
    SimilaritySearch,
)
from pandora.similarity.hits import iter_edges


def cluster_edges(
    item_ids: list[str],
    edges: Iterable[Sequence[Any]],
) -> tuple[list[SimilarityCluster], ClusteringProvenance]:
    """Group item ids into connected-component clusters.

    Two items land in the same cluster iff connected through a chain
    of edges. Items with no edges form their own singleton cluster, so
    every id in `item_ids` ends up in exactly one cluster. `edges` is
    consumed once and never stored, so it can be a stream.

    Args:
        item_ids: Every item id to place into a cluster.
        edges: Pairs of item ids, as `Edge`s from `iter_edges()` or any
            tuples whose first two items are the ids. Duplicates are
            harmless.

    Returns:
        `(clusters, provenance)` — the clusters (sorted by their
        smallest member), and the edge and cluster counts.
        `provenance.hit_filter` and `.search` are None.
    """

    parent = {item_id: item_id for item_id in item_ids}

    def find(item_id: str) -> str:
        """Union-find: root id of item_id's cluster, with path compression."""

        while parent[item_id] != item_id:
            parent[item_id] = parent[parent[item_id]]
            item_id = parent[item_id]
        return item_id

    n_edges = 0
    n_unknown = 0
    for edge in edges:
        n_edges += 1
        source_id, target_id = edge[0], edge[1]
        if source_id not in parent or target_id not in parent:
            n_unknown += 1
            continue
        root_source, root_target = find(source_id), find(target_id)
        if root_source != root_target:
            parent[root_target] = root_source

    groups: dict[str, list[str]] = {}
    for item_id in item_ids:
        groups.setdefault(find(item_id), []).append(item_id)

    clusters = [
        SimilarityCluster(components=sorted(members), n_components=len(members))
        for members in sorted(groups.values(), key=lambda members: min(members))
    ]
    provenance = ClusteringProvenance(
        clustered_at=now_iso(),
        n_edges=n_edges,
        n_edges_unknown_ids=n_unknown,
        n_clusters=len(clusters),
    )
    return clusters, provenance


def cluster_similar_items(
    item_ids: list[str],
    search: SimilaritySearch,
    hit_filter: HitFilter,
    *,
    interface_residues: dict[str, set[int]] | None = None,
) -> tuple[list[SimilarityCluster], ClusteringProvenance]:
    """Cluster item ids from a search's hit file, filtered by a policy.

    Streams `search`'s hit file through `hit_filter` (`iter_edges()`)
    into `cluster_edges()`; memory is bounded by the number of items,
    not hits.

    Args:
        item_ids: Every item id to place into a cluster.
        search: The search whose hit file to read.
        hit_filter: Which hit rows become edges.
        interface_residues: Item id -> interface residue positions,
            needed when `hit_filter.min_interface_coverage` is set.

    Returns:
        `(clusters, provenance)`, with `provenance.hit_filter` and
        `.search` recorded so the clustering can be reproduced.

    Raises:
        ValueError: The filter can't apply to the search, or the hit
            file is malformed (see `iter_edges()`).
    """

    clusters, provenance = cluster_edges(
        item_ids,
        iter_edges(search, hit_filter, interface_residues=interface_residues),
    )
    return clusters, provenance.model_copy(
        update={"hit_filter": hit_filter, "search": search}
    )


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
