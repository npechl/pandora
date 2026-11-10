"""Interface-vs-interface similarity: per-side Jaccard of interface
residues between interfaces of the same two proteins."""

from __future__ import annotations

import itertools
import re
from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping

import gemmi

from pandora.schemas.annotation import LayerRow, ResidueContactRow
from pandora.schemas.similarity import SimilarityCluster
from pandora.similarity.layer import item_id, parse_item_id

_COUNT_KEYS = (
    "rows_ignored",
    "interfaces",
    "interfaces_without_cluster",
    "interfaces_without_sequence",
    "interfaces_sequence_mismatch",
    "pairs_compared",
    "pairs_below_identity",
    "pairs_swapped",
    "pairs_zero_score",
    "alignments",
)

_Key = tuple[str, str, str]  # (upper-cased entry, smaller chain, other chain)
_Alignment = tuple[dict[int, int] | None, float]  # (None = same positions)


def interface_similarity(
    contacts: Iterable[ResidueContactRow],
    chain_clusters: Iterable[SimilarityCluster],
    sequences: Mapping[tuple[str, str], str],
    min_identity: float = 0.9,
) -> tuple[list[LayerRow], dict[str, int]]:
    """Score pairs of interfaces of the same two proteins by how much
    their interface residues overlap, side by side.

    An interface is an entry's unordered chain pair with at least one
    contact row where both residues are polymer residues (`seq_id`
    set) on different chains, as in `interface_items`. Each side is
    the set of `seq_id`s of its chain in contact with the other chain.

    Two interfaces are compared when their chains fall in matching
    `chain_clusters`, in either chain order ("direct": smaller chain
    ids paired; "swapped": crossed). An order counts only if each
    paired chain is at least `min_identity` identical to its
    counterpart on a pairwise gemmi alignment (identical sequences
    skip the alignment). Residues map across the alignment; a residue
    it cannot place counts as not shared. Per side,
    `jaccard = shared / (|side 1| + |side 2| - shared)`, and `score`
    is the smaller of the two sides. With both orders possible
    (homomers) the higher score wins, "direct" on a tie.

    `seq_id` must be the position in the chain's sequence (structures
    that keep `label_seq_id`). An interface whose contact rows name a
    residue that disagrees with the sequence (its `comp_id`'s one-letter
    code, or a `seq_id` past the end) is skipped and counted.

    Args:
        contacts: The rows of one `residue_contacts` layer.
        chain_clusters: Chain clusters (components `ENTRY_CHAIN`), e.g.
            sequence clusters at 90% identity.
        sequences: One-letter chain sequence by
            `(entry_id.upper(), chain_id)`.
        min_identity: The lowest pairwise identity, per chain pair,
            for two interfaces to be the same proteins.

    Returns:
        The rows, one per compared pair, sorted: `item_id_1` <
        `item_id_2` (interface ids), `score`, `jaccard_1` and
        `identity_1` (the smaller chain of `item_id_1` against its
        counterpart), `jaccard_2`, `chain_order` ("direct" or
        "swapped"), `identity_2`; and counts: `rows_ignored`,
        `interfaces`, `interfaces_without_cluster`,
        `interfaces_without_sequence`, `interfaces_sequence_mismatch`,
        `pairs_compared`, `pairs_below_identity`, `pairs_swapped`,
        `pairs_zero_score` and `alignments` (alignments computed).
    """

    counts = dict.fromkeys(_COUNT_KEYS, 0)
    sides, names, mismatched = _collect_sides(contacts, sequences, counts)
    cluster_of = _cluster_index(chain_clusters)
    groups = _group(sides, cluster_of, sequences, mismatched, counts)

    cache: dict[tuple[str, str], _Alignment] = {}

    def align(chain_1: tuple[str, str], chain_2: tuple[str, str]) -> _Alignment:
        """Cached alignment of two chains' sequences."""

        query, target = sequences[chain_1], sequences[chain_2]
        if query == target:
            return None, 1.0
        if (query, target) not in cache:
            cache[(query, target)] = _align(query, target)
            counts["alignments"] += 1
        return cache[(query, target)]

    rows = []
    for group in groups.values():
        # ponytail: all-vs-all per cluster pair is quadratic; stream rows
        # in batches if a group gets huge.
        for key_1, key_2 in itertools.combinations(
            sorted(group, key=names.__getitem__), 2
        ):
            best = _best_order(
                key_1,
                key_2,
                _orders(key_1, key_2, cluster_of),
                align,
                sides,
                names,
                min_identity,
            )
            if best is None:
                counts["pairs_below_identity"] += 1
                continue
            counts["pairs_swapped"] += best.chain_order == "swapped"
            counts["pairs_zero_score"] += best.score == 0
            rows.append(best)
    rows.sort(key=lambda r: (r.item_id_1, r.item_id_2))
    counts["pairs_compared"] = len(rows)
    return rows, counts


def _best_order(
    key_1: _Key,
    key_2: _Key,
    orders: list[tuple[str, list[tuple[tuple[str, str], tuple[str, str]]]]],
    align: Callable[[tuple[str, str], tuple[str, str]], _Alignment],
    sides: Mapping[_Key, Mapping[str, set[int]]],
    names: Mapping[_Key, str],
    min_identity: float,
) -> LayerRow | None:
    """The row of the best-scoring chain order passing the identity
    gate, or None if no order passes."""

    best = None
    for order, pairs in orders:
        alignments = [align(x, y) for x, y in pairs]
        if any(ident < min_identity for _, ident in alignments):
            continue
        jaccards = [
            _jaccard(sides[key_1][x[1]], mapping, sides[key_2][y[1]])
            for (x, y), (mapping, _) in zip(pairs, alignments, strict=True)
        ]
        if best is None or min(jaccards) > best.score:
            best = LayerRow(
                item_id_1=names[key_1],
                item_id_2=names[key_2],
                score=min(jaccards),
                jaccard_1=jaccards[0],
                jaccard_2=jaccards[1],
                chain_order=order,
                identity_1=alignments[0][1],
                identity_2=alignments[1][1],
            )
    return best


def _group(
    sides: Mapping[_Key, object],
    cluster_of: Mapping[tuple[str, str], int],
    sequences: Mapping[tuple[str, str], str],
    mismatched: set[_Key],
    counts: dict[str, int],
) -> dict[frozenset[int], list[_Key]]:
    """Comparable interfaces by their unordered cluster pair; counts
    the ones left out."""

    groups: dict[frozenset[int], list[_Key]] = defaultdict(list)
    for key in sides:
        chains = [(key[0], key[1]), (key[0], key[2])]
        if not all(c in cluster_of for c in chains):
            counts["interfaces_without_cluster"] += 1
        elif not all(c in sequences for c in chains):
            counts["interfaces_without_sequence"] += 1
        elif key in mismatched:
            counts["interfaces_sequence_mismatch"] += 1
        else:
            groups[frozenset(cluster_of[c] for c in chains)].append(key)
    return groups


def _collect_sides(
    contacts: Iterable[ResidueContactRow],
    sequences: Mapping[tuple[str, str], str],
    counts: dict[str, int],
) -> tuple[dict[_Key, dict[str, set[int]]], dict[_Key, str], set[_Key]]:
    """Each interface's residue sets per chain, its item id, and the
    interfaces whose rows disagree with a chain's sequence."""

    sides: dict[_Key, dict[str, set[int]]] = {}
    names: dict[_Key, str] = {}
    mismatched: set[_Key] = set()
    for row in contacts:
        if (
            row.seq_id_1 is None
            or row.seq_id_2 is None
            or row.chain_id_1 == row.chain_id_2
        ):
            counts["rows_ignored"] += 1
            continue
        entry = row.entry_id.upper()
        key = (entry, *sorted((row.chain_id_1, row.chain_id_2)))
        if key not in sides:
            sides[key] = {key[1]: set(), key[2]: set()}
            names[key] = item_id("interface", row.entry_id, *key[1:])
        for chain, seq_id, comp_id in (
            (row.chain_id_1, row.seq_id_1, row.comp_id_1),
            (row.chain_id_2, row.seq_id_2, row.comp_id_2),
        ):
            sides[key][chain].add(seq_id)
            sequence = sequences.get((entry, chain))
            if sequence is not None and not _residue_fits(
                sequence, seq_id, comp_id
            ):
                mismatched.add(key)
    counts["interfaces"] = len(sides)
    return sides, names, mismatched


def _cluster_index(
    chain_clusters: Iterable[SimilarityCluster],
) -> dict[tuple[str, str], int]:
    """Cluster number by `(entry_id.upper(), chain_id)`."""

    cluster_of: dict[tuple[str, str], int] = {}
    for index, cluster in enumerate(chain_clusters):
        for component in cluster.components:
            entry, chain = parse_item_id("chain", component)
            cluster_of[(entry.upper(), chain)] = index
    return cluster_of


def _orders(
    key_1: _Key, key_2: _Key, cluster_of: Mapping[tuple[str, str], int]
) -> list[tuple[str, list[tuple[tuple[str, str], tuple[str, str]]]]]:
    """The chain orders whose clusters line up, each as its two
    `(chain of key_1, chain of key_2)` pairings."""

    (e1, a, b), (e2, c, d) = key_1, key_2
    candidates = (("direct", (c, d)), ("swapped", (d, c)))
    orders = []
    for name, (x, y) in candidates:
        pairs = [((e1, a), (e2, x)), ((e1, b), (e2, y))]
        if all(cluster_of[p] == cluster_of[q] for p, q in pairs):
            orders.append((name, pairs))
    return orders


def _residue_fits(sequence: str, seq_id: int, comp_id: str | None) -> bool:
    """Whether `seq_id` is inside `sequence` and agrees with `comp_id`."""

    if not 1 <= seq_id <= len(sequence):
        return False
    if comp_id is None:
        return True
    info = gemmi.find_tabulated_residue(comp_id)
    code = info.one_letter_code.strip() if info is not None else ""
    return (code.upper() or "X") == sequence[seq_id - 1].upper()


def _align(query: str, target: str) -> tuple[dict[int, int], float]:
    """Position map (1-based, query -> target) and identity over
    aligned columns, from a gemmi global alignment."""

    result = gemmi.align_string_sequences(list(query), list(target), [])
    mapping: dict[int, int] = {}
    identical = 0
    i = j = 0
    for length, op in re.findall(r"(\d+)([MID])", result.cigar_str()):
        for _ in range(int(length)):
            if op == "M":
                mapping[i + 1] = j + 1
                identical += query[i] == target[j]
            i += op in "MI"
            j += op in "MD"
    return mapping, identical / len(mapping) if mapping else 0.0


def _jaccard(
    side_1: set[int], mapping: dict[int, int] | None, side_2: set[int]
) -> float:
    """Jaccard of two sides; residues `mapping` can't place still count
    in the union."""

    mapped = (
        side_1
        if mapping is None
        else {mapping[i] for i in side_1 if i in mapping}
    )
    shared = len(mapped & side_2)
    union = len(side_1) + len(side_2) - shared
    return shared / union if union else 0.0
