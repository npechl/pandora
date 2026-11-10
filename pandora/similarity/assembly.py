"""The final PPI pairs table: positives and negatives with their splits
and C-chain leak flags."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from pandora.schemas.similarity import DatasetPair, NegativePair, PositiveDimer
from pandora.similarity.layer import item_id


def _split_of(
    positives: Sequence[PositiveDimer], splits: Mapping[str, Sequence[str]]
) -> dict[str, str]:
    """Positive dimer -> its split; ValueError on a positive in no split
    or two, or a split item that is not a positive."""

    known = {p.dimer for p in positives}
    found: dict[str, str] = {}
    for name, members in splits.items():
        for dimer in members:
            if dimer not in known:
                raise ValueError(f"split item {dimer} is not a positive")
            if dimer in found:
                raise ValueError(f"positive {dimer} is in two splits")
            found[dimer] = name
    missing = sorted(known - found.keys())
    if missing:
        raise ValueError(f"positive {missing[0]} is in no split")
    return found


def _rep(representative_of: Mapping[str, str], chain: str) -> str:
    """`chain`'s representative; ValueError if it has none."""

    try:
        return representative_of[chain]
    except KeyError:
        raise ValueError(f"chain {chain} has no representative") from None


def assemble_pairs(
    positives: Sequence[PositiveDimer],
    negatives: Sequence[NegativePair],
    splits: Mapping[str, Sequence[str]],
    representative_of: Mapping[str, str],
    tm_edges: Iterable[Sequence[Any]],
) -> list[DatasetPair]:
    """Build the final pairs table with splits and C-chain leak flags.

    Each negative takes its positive's split, kind and group. A negative
    (A, C) is flagged `c_leak_positive` when C's representative is
    structurally similar to the representative of a positive's chain in
    another split, and `c_leak_negative` when it is similar to another
    split's negative C.

    Args:
        positives: The positive dimers (from `rebalance`).
        negatives: The negative pairs (from `sample_negatives`).
        splits: Split name -> positive dimer ids (from `partition`).
        representative_of: Chain item id -> representative id.
        tm_edges: Pairs of chain item ids that are structurally similar
            (already filtered). Ends missing from `representative_of`
            are kept as they are.

    Returns:
        Positives (sorted by dimer), then negatives (sorted by positive,
        then C).

    Raises:
        ValueError: A positive in no split or in two, a split item that
            is not a positive, a negative whose positive is unknown, or
            a pair chain with no representative.
    """

    split_of = _split_of(positives, splits)
    by_dimer = {p.dimer: p for p in positives}
    for negative in negatives:
        if negative.positive not in by_dimer:
            raise ValueError(
                f"negative {negative.chain_1}–{negative.chain_2} names "
                f"unknown positive {negative.positive}"
            )

    tm: dict[str, set[str]] = defaultdict(set)
    for edge in tm_edges:
        a = representative_of.get(edge[0], edge[0])
        b = representative_of.get(edge[1], edge[1])
        tm[a].add(b)
        tm[b].add(a)
    positive_reps: dict[str, set[str]] = defaultdict(set)
    rows: list[DatasetPair] = []
    for p in sorted(positives, key=lambda p: p.dimer):
        chains = [
            item_id("chain", p.entry_id, c)
            for c in (p.chain_id_1, p.chain_id_2)
        ]
        for chain in chains:
            positive_reps[split_of[p.dimer]].add(_rep(representative_of, chain))
        rows.append(
            DatasetPair(
                chain_a=chains[0],
                chain_b=chains[1],
                label=1,
                kind=p.kind,
                cluster=p.group,
                split=split_of[p.dimer],
                positive=p.dimer,
            )
        )
    ordered = sorted(negatives, key=lambda n: (n.positive, n.chain_2))
    c_reps: dict[str, set[str]] = defaultdict(set)
    for n in ordered:
        _rep(representative_of, n.chain_1)
        c_reps[split_of[n.positive]].add(_rep(representative_of, n.chain_2))
    names = set(splits)

    def others(sets: Mapping[str, set[str]], split: str) -> set[str]:
        return set().union(*(sets[s] for s in names - {split}))

    for n in ordered:
        p = by_dimer[n.positive]
        split = split_of[p.dimer]
        similar = tm[representative_of[n.chain_2]]
        rows.append(
            DatasetPair(
                chain_a=n.chain_1,
                chain_b=n.chain_2,
                label=0,
                kind=p.kind,
                cluster=p.group,
                split=split,
                positive=p.dimer,
                c_leak_positive=bool(similar & others(positive_reps, split)),
                c_leak_negative=bool(similar & others(c_reps, split)),
            )
        )
    return rows
