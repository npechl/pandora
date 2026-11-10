"""Seeded negative chain-pair sampling (PPI v2.1 Step 8)."""

from __future__ import annotations

import random
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from pandora.schemas.similarity import NegativePair, PositiveDimer
from pandora.similarity.layer import item_id


def _strict(representative_of: Mapping[str, str], chain: str) -> str:
    """`chain`'s representative; ValueError if it has none."""

    try:
        return representative_of[chain]
    except KeyError:
        raise ValueError(f"chain {chain} has no representative") from None


def _neighbours(
    pairs: Iterable[tuple[str, str]], loops: bool
) -> dict[str, set[str]]:
    """Undirected adjacency of `pairs`; self-pairs kept only if `loops`."""

    graph: dict[str, set[str]] = defaultdict(set)
    for a, b in pairs:
        if a != b or loops:
            graph[a].add(b)
            graph[b].add(a)
    return graph


def sample_negatives(
    positives: Sequence[PositiveDimer],
    representative_of: Mapping[str, str],
    tm_edges: Iterable[Sequence[Any]],
    binding_pairs: Iterable[tuple[str, str]],
    seed: int,
) -> tuple[list[NegativePair], list[str]]:
    """Draw one negative partner per positive dimer.

    Every chain is looked up through its representative. The pool is
    every representative minus those of the positives' chains. For each
    positive, in dimer-id order: pick one of its two chains, A, at
    random; exclude for this pick set 1 (A and the chains structurally
    similar to it), set 2 (the chains binding set 1) and set 3 (the
    chains similar to set 2); draw C from the rest and remove it from
    the pool for good.

    Args:
        positives: The positive dimers, e.g. from `rebalance`.
        representative_of: Chain item id -> its representative's id
            (from `member_representatives`); its values are the pool.
        tm_edges: Pairs of chain item ids that are structurally similar
            (already filtered, e.g. TM-score >= 0.5). Ends missing from
            `representative_of` are kept as they are.
        binding_pairs: Pairs of chain item ids in contact (any contact).
        seed: Seed of the one random generator; same inputs and seed
            give the same output.

    Returns:
        The negative pairs, and the dimer ids of the positives the pool
        ran out for.

    Raises:
        ValueError: A positive's chain or a binding pair's chain has no
            representative.
    """

    def loose(chain: str) -> str:
        return representative_of.get(chain, chain)

    tm = _neighbours(
        ((loose(e[0]), loose(e[1])) for e in tm_edges), loops=False
    )
    binders = _neighbours(
        (
            (_strict(representative_of, a), _strict(representative_of, b))
            for a, b in binding_pairs
        ),
        loops=True,
    )
    ordered = sorted(positives, key=lambda p: p.dimer)
    sides = {
        p.dimer: [
            item_id("chain", p.entry_id, p.chain_id_1),
            item_id("chain", p.entry_id, p.chain_id_2),
        ]
        for p in ordered
    }
    pool = set(representative_of.values()) - {
        _strict(representative_of, c) for pair in sides.values() for c in pair
    }

    rng = random.Random(seed)
    negatives: list[NegativePair] = []
    unpaired: list[str] = []
    for positive in ordered:
        chain_a = rng.choice(sides[positive.dimer])
        a = representative_of[chain_a]
        set_1 = {a} | tm[a]
        set_2 = set().union(*(binders[r] for r in set_1))
        set_3 = set().union(*(tm[r] for r in set_2))
        candidates = sorted(pool - set_1 - set_2 - set_3)
        if not candidates:
            unpaired.append(positive.dimer)
            continue
        chain_c = rng.choice(candidates)
        pool.discard(chain_c)
        negatives.append(
            NegativePair(
                chain_1=chain_a, chain_2=chain_c, positive=positive.dimer
            )
        )
    return negatives, unpaired
