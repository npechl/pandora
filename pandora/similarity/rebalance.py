"""Homomer/heteromer rebalancing of dimer groups (PPI v2.1 Step 7)."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any

from pandora.schemas.similarity import (
    DroppedDimer,
    PositiveDimer,
    SimilarityCluster,
)
from pandora.similarity.layer import item_id, parse_item_id


def _link(
    adjacency: dict[str, set[str]],
    side_1: Iterable[str],
    side_2: Iterable[str],
) -> None:
    """Make every dimer of `side_1` adjacent to every one of `side_2`."""

    side_2 = list(side_2)
    for a in side_1:
        for b in side_2:
            if a != b:
                adjacency[a].add(b)
                adjacency[b].add(a)


def dimer_adjacency(
    interface_ids: Iterable[str],
    chain_edges: Iterable[Sequence[Any]],
    dimer_edges: Iterable[Sequence[Any]],
    chain_map: Mapping[str, str] | None = None,
) -> dict[str, set[str]]:
    """Build the dimer-to-dimer similarity graph.

    Two dimers are adjacent when they share a chain, when a chain-level
    edge links a chain of one to a chain of the other, or when an
    interface-level edge links them. With `chain_map`, every chain is
    first replaced by the id it maps to (e.g. its 90% representative),
    so dimers whose chains share a representative are adjacent too.

    Unlike `interface_edges`, which keeps only what connectivity needs,
    this gives every adjacency, so degrees can be counted.

    Args:
        interface_ids: The dimers (interface item ids); each becomes a
            key of the result.
        chain_edges: Pairs of chain item ids from chain-level layers.
        dimer_edges: Pairs of interface item ids from interface-level
            layers; pairs naming an unknown dimer are ignored.
        chain_map: Chain item id -> the id to key it by; chains missing
            from it keep their own id.

    Returns:
        Dimer -> its neighbours (never itself).

    Raises:
        ValueError: An id is not an interface or chain item id.
    """

    mapped = chain_map or {}

    def key(chain: str) -> str:
        return mapped.get(chain, chain)

    adjacency: dict[str, set[str]] = {}
    holders: dict[str, set[str]] = defaultdict(set)
    for dimer in interface_ids:
        adjacency[dimer] = set()
        entry, *chains = parse_item_id("interface", dimer)
        for chain in chains:
            holders[key(item_id("chain", entry, chain))].add(dimer)

    # ponytail: all-pairs per chain; fine for dimer groups, quadratic
    # if one chain sits in thousands of dimers.
    for dimers in holders.values():
        _link(adjacency, dimers, dimers)
    for edge in chain_edges:
        _link(
            adjacency,
            holders.get(key(edge[0]), ()),
            holders.get(key(edge[1]), ()),
        )
    for edge in dimer_edges:
        if edge[0] in adjacency and edge[1] in adjacency:
            _link(adjacency, [edge[0]], [edge[1]])
    return adjacency


def _positive(dimer: str, group: str, kind: str, degree: int) -> PositiveDimer:
    """A PositiveDimer record for `dimer`."""

    entry, chain_1, chain_2 = parse_item_id("interface", dimer)
    return PositiveDimer(
        dimer=dimer,
        entry_id=entry,
        chain_id_1=chain_1,
        chain_id_2=chain_2,
        group=group,
        kind=kind,
        degree=degree,
    )


def rebalance(
    groups: Sequence[SimilarityCluster],
    adjacency: Mapping[str, set[str]],
    homomers: set[str],
    sort_key: Callable[[str], Any],
) -> tuple[list[PositiveDimer], list[DroppedDimer]]:
    """Pick the positive dimers of each group.

    A group of homomers only keeps its most connected dimer. A group
    with heteromers drops its homomers, then repeatedly keeps the
    heteromer with the fewest remaining heteromer neighbours and drops
    those neighbours; degrees are recounted every round.

    Args:
        groups: The dimer groups (connected components of `adjacency`).
        adjacency: Dimer -> neighbours, e.g. from `dimer_adjacency`.
        homomers: The dimers whose two chains are the same protein.
        sort_key: Breaks ties between equally connected dimers: the
            smallest key wins (e.g. a `RepresentativePolicy` ranking).

    Returns:
        The positives, group by group in pick order, and the dropped
        dimers with their reasons.

    Raises:
        ValueError: An adjacency links two different groups, so the
            groups were not built from the same graph.
    """

    group_of = {d: g.components[0] for g in groups for d in g.components}
    for dimer, group in group_of.items():
        for other in adjacency.get(dimer, ()):
            if group_of.get(other, group) != group:
                raise ValueError(
                    f"{dimer} and {other} are similar but in different "
                    "groups; were the groups built from the same layers?"
                )

    positives: list[PositiveDimer] = []
    dropped: list[DroppedDimer] = []
    for g in groups:
        group = g.components[0]
        heteromers = {d for d in g.components if d not in homomers}
        if not heteromers:
            degree = {d: len(adjacency.get(d, ())) for d in g.components}
            kept = min(g.components, key=lambda d: (-degree[d], sort_key(d)))
            positives.append(_positive(kept, group, "homomer", degree[kept]))
            dropped += [
                DroppedDimer(dimer=d, group=group, reason="not_most_connected")
                for d in g.components
                if d != kept
            ]
            continue
        dropped += [
            DroppedDimer(dimer=d, group=group, reason="homomer_in_mixed_group")
            for d in g.components
            if d in homomers
        ]
        remaining = set(heteromers)
        # ponytail: O(n^2) per group; a heap if groups reach thousands.
        while remaining:
            degree = {
                d: len(adjacency.get(d, set()) & remaining) for d in remaining
            }
            kept = min(remaining, key=lambda d: (degree[d], sort_key(d)))
            positives.append(_positive(kept, group, "heteromer", degree[kept]))
            neighbours = adjacency.get(kept, set()) & remaining
            dropped += [
                DroppedDimer(
                    dimer=d, group=group, reason="similar_to_kept", kept=kept
                )
                for d in sorted(neighbours)
            ]
            remaining -= neighbours | {kept}
    return positives, dropped
