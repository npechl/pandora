"""Pick one representative item per cluster, by a ranking policy."""

from __future__ import annotations

from collections.abc import Collection, Iterable, Mapping, Sequence
from typing import Any

from pandora._util import rank_sort_key
from pandora.schemas.similarity import (
    RepresentativeCluster,
    RepresentativePolicy,
    SimilarityCluster,
)


def member_representatives(
    clusters: Iterable[RepresentativeCluster],
) -> dict[str, str]:
    """Map every member of every cluster to its cluster's representative.

    Args:
        clusters: Clusters with their representatives, as written by
            `pandora representatives`.

    Returns:
        Item id -> representative id, for every member (the
        representative maps to itself).
    """

    return {
        member: cluster.representative
        for cluster in clusters
        for member in cluster.members
    }


def restrict_clusters(
    clusters: Iterable[SimilarityCluster], item_ids: Collection[str]
) -> list[SimilarityCluster]:
    """Keep only the cluster members in `item_ids`.

    Args:
        clusters: Clusters, e.g. read from `pandora cluster`'s output.
        item_ids: The items to keep, e.g. the chains of a step folder.

    Returns:
        New clusters holding only kept members, in input order; clusters
        left empty are dropped.
    """

    kept = []
    for cluster in clusters:
        members = [m for m in cluster.components if m in item_ids]
        if members:
            kept.append(
                SimilarityCluster(components=members, n_components=len(members))
            )
    return kept


def cluster_hit_counts(
    clusters: Iterable[SimilarityCluster], edges: Iterable[Sequence[Any]]
) -> dict[str, int]:
    """Count, per item, the distinct other items of its own cluster it
    shares an edge with ("most connected").

    Args:
        clusters: The clusters; only their members are counted.
        edges: Pairs of item ids, as `Edge`s from `similarity_edges()` or
            any tuples whose first two items are the ids. Self-edges and
            edges to ids outside `clusters` are ignored.

    Returns:
        Item id -> number of partners, for every clustered item (0 when
        it has none).
    """

    cluster_of = {
        member: index
        for index, cluster in enumerate(clusters)
        for member in cluster.components
    }
    partners: dict[str, set[str]] = {item: set() for item in cluster_of}
    for edge in edges:
        a, b = edge[0], edge[1]
        if a == b or a not in cluster_of or b not in cluster_of:
            continue
        if cluster_of[a] == cluster_of[b]:
            partners[a].add(b)
            partners[b].add(a)
    return {item: len(found) for item, found in partners.items()}


def pick_representatives(
    clusters: Iterable[SimilarityCluster],
    values: Mapping[str, Mapping[str, float | int | None]],
    policy: RepresentativePolicy,
    eligible: Collection[str] | None = None,
) -> list[RepresentativeCluster]:
    """Pick the best member of each cluster by `policy.rank`.

    Args:
        clusters: The clusters to pick from.
        values: Item id -> {rank key: value}, for every member that can
            be picked. Null values rank last.
        policy: The ranking keys, best first; the lowest item id breaks
            the last tie.
        eligible: The items that may be picked (e.g. no symmetry copies);
            None means all. A cluster with no eligible member picks from
            all its members.

    Returns:
        One `RepresentativeCluster` per cluster, sorted by representative.

    Raises:
        ValueError: A rank key is missing from a ranked item's values.
    """

    picked = []
    for cluster in clusters:
        members = cluster.components
        pool = [m for m in members if eligible is None or m in eligible]
        try:
            best = min(
                pool or members,
                key=lambda m: rank_sort_key(m, values[m], policy.rank),
            )
        except KeyError as exc:
            raise ValueError(
                f"rank key or item {exc.args[0]!r} has no value"
            ) from exc
        picked.append(
            RepresentativeCluster(
                representative=best,
                members=list(members),
                n_members=len(members),
            )
        )
    return sorted(picked, key=lambda r: r.representative)
