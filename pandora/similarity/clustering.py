from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from typing import Any

from pandora._util import now_iso
from pandora.schemas.similarity import ClusteringProvenance, SimilarityCluster
from pandora.similarity.layer import parse_item_id


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
        edges: Pairs of item ids, as `Edge`s from `similarity_edges()` or any
            tuples whose first two items are the ids. Duplicates are
            harmless.

    Returns:
        `(clusters, provenance)` — the clusters (sorted by their
        smallest member), and the edge and cluster counts.
        `provenance.layers` is empty.
    """

    # ponytail: dict union-find on string ids; integer arrays if memory
    # bites at archive scale (~1M items is a few hundred MB here).
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


def entry_edges(edges: Iterable[Sequence[Any]]) -> Iterator[tuple[str, str]]:
    """Roll chain-item edges up to entry-item edges.

    Each edge between two chains becomes an edge between their entries,
    so clustering the result puts two entries together when any of their
    chains are similar. Edges within one entry are dropped.

    Args:
        edges: Pairs of chain item ids (`{entry}_{chain}`), as `Edge`s
            or any tuples whose first two items are the ids.

    Returns:
        An iterator of `(entry_id_1, entry_id_2)` pairs.

    Raises:
        ValueError: An id is not a chain item id (while iterating).
    """

    for edge in edges:
        entry_1 = parse_item_id("chain", edge[0])[0]
        entry_2 = parse_item_id("chain", edge[1])[0]
        if entry_1 != entry_2:
            yield entry_1, entry_2


def _chain_key(entry_id: str, chain_id: str) -> str:
    """Case-insensitive lookup key of a chain."""

    return f"{entry_id.upper()}_{chain_id}"


def interface_edges(
    edges: Iterable[Sequence[Any]], interface_ids: Iterable[str]
) -> Iterator[tuple[str, str]]:
    """Roll chain-item edges up to interface-item edges.

    Clustering the result puts two interfaces together when any chain
    of one is linked to any chain of the other. Interfaces that share a
    chain are linked too, since a chain matches itself. Edges naming a
    chain that is in no interface are dropped: such a chain must not
    join the interfaces of two chains it is similar to. Each chain edge
    gives one interface edge, so the output grows with the input, not
    with the number of interfaces per chain.

    Args:
        edges: Pairs of chain item ids (`{entry}_{chain}`), as `Edge`s
            or any tuples whose first two items are the ids.
        interface_ids: The interface item ids (`{entry}_{chain}_{chain}`)
            to link.

    Returns:
        An iterator of `(interface_id_1, interface_id_2)` pairs: first
        the links between interfaces sharing a chain, then one per
        chain edge whose chains are in different interfaces.

    Raises:
        ValueError: An id is not an interface or chain item id (while
            iterating).
    """

    first: dict[str, str] = {}
    for interface in interface_ids:
        entry, *chains = parse_item_id("interface", interface)
        for chain in chains:
            key = _chain_key(entry, chain)
            if key in first:
                yield first[key], interface
            else:
                first[key] = interface
    for edge in edges:
        side_1 = first.get(_chain_key(*parse_item_id("chain", edge[0])))
        side_2 = first.get(_chain_key(*parse_item_id("chain", edge[1])))
        if side_1 is not None and side_2 is not None and side_1 != side_2:
            yield side_1, side_2


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
            `cluster_edges()`.

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
