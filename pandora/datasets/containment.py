"""Find entries (biounits) contained in another: X is contained in Y when,
for every cluster, Y has at least as many chains in it as X."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from pandora._util import rank_sort_key
from pandora.schemas.common import RankKey
from pandora.schemas.dataset import ContainmentRecord


def _contains(outer: Mapping[str, int], inner: Mapping[str, int]) -> bool:
    """True if `outer` has at least `inner`'s count in every cluster."""

    return all(outer.get(c, 0) >= n for c, n in inner.items())


def _containers(
    cluster_counts: Mapping[str, Mapping[str, int]],
) -> dict[str, list[str]]:
    """Entry id -> sorted ids of the other entries containing it; entries
    with no chains are left out."""

    index: dict[str, set[str]] = {}
    for entry, counts in cluster_counts.items():
        for cluster, n in counts.items():
            if n > 0:
                index.setdefault(cluster, set()).add(entry)

    containers: dict[str, list[str]] = {}
    for entry, counts in cluster_counts.items():
        clusters = [c for c, n in counts.items() if n > 0]
        if not clusters:
            continue
        rarest = min(clusters, key=lambda c: len(index[c]))
        containers[entry] = sorted(
            other
            for other in index[rarest]
            if other != entry and _contains(cluster_counts[other], counts)
        )
    return containers


def find_contained(
    cluster_counts: Mapping[str, Mapping[str, int]],
    rank_values: Mapping[str, Mapping[str, float | int | None]],
    rank: Sequence[RankKey],
) -> list[ContainmentRecord]:
    """Find the entries to remove because another entry contains them.

    X is contained in Y when, for every cluster, Y has at least as many
    chains in it as X (so a monomer A is contained in a homodimer AA,
    not the other way round). Strictly contained entries are removed.
    Among entries with equal counts, the best by `rank` is kept and the
    others are removed. Entries with no chains are always kept.

    Args:
        cluster_counts: Entry id -> {cluster id: how many of the entry's
            chains are in that cluster}.
        rank_values: Entry id -> {rank key: value}, for every entry.
        rank: Ranking keys deciding ties, best first; the lowest entry
            id breaks the last tie. Nulls rank last.

    Returns:
        One record per removed entry, sorted by entry id. `contained_in`
        is the lowest-id kept entry that contains it.

    Raises:
        KeyError: An entry in a tie lacks a value for a rank key.
    """

    containers = _containers(cluster_counts)
    removed: dict[str, bool] = {}
    for entry, found in containers.items():
        counts = cluster_counts[entry]
        if any(not _contains(counts, cluster_counts[o]) for o in found):
            removed[entry] = False
    for entry, found in containers.items():
        if entry in removed or not found:
            continue
        # Every container left has equal counts: a tie group.
        best = min(
            [entry, *found],
            key=lambda e: rank_sort_key(e, rank_values[e], rank),
        )
        if best != entry:
            removed[entry] = True

    records = []
    for entry in sorted(removed):
        kept = [o for o in containers[entry] if o not in removed]
        assert kept, f"no kept container for {entry}"
        records.append(
            ContainmentRecord(
                entry_id=entry, contained_in=kept[0], equal=removed[entry]
            )
        )
    return records
