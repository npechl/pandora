from __future__ import annotations

import copy

from pandora.schemas.annotation import ResidueContactRow
from pandora.schemas.similarity import SimilarityCluster
from pandora.similarity import interface_similarity

_SEQ = ("ACDEFGHIKLMNPQRSTVWY" * 3)[:50]

_COUNT_KEYS = {
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
}


def _row(entry, a, sa, b, sb, ca=None, cb=None):
    return ResidueContactRow(
        entry_id=entry,
        chain_id_1=a,
        seq_id_1=sa,
        comp_id_1=ca,
        chain_id_2=b,
        seq_id_2=sb,
        comp_id_2=cb,
    )


def _iface(entry, a, side_a, b, side_b):
    """Contact rows giving chain `a` residues `side_a` against `b`."""

    return [_row(entry, a, x, b, y) for x in side_a for y in side_b]


def _cl(*items):
    return SimilarityCluster(components=list(items), n_components=len(items))


def _seqs(*chains, seq=_SEQ):
    return {(entry, chain): seq for entry, chain in chains}


_ALL = _seqs(
    ("1AAA", "A"),
    ("1AAA", "B"),
    ("1AAA", "C"),
    ("1AAA", "D"),
    ("2BBB", "A"),
    ("2BBB", "B"),
    ("3CCC", "A"),
    ("3CCC", "B"),
)
_HETERO = [_cl("1AAA_A", "2BBB_A"), _cl("1AAA_B", "2BBB_B")]
_HOMO = [_cl("1AAA_A", "1AAA_B", "2BBB_A", "2BBB_B")]


def test_identical_heterodimers_score_one():
    contacts = _iface("1AAA", "A", {1, 2, 3}, "B", {4, 5}) + _iface(
        "2BBB", "A", {1, 2, 3}, "B", {4, 5}
    )
    rows, counts = interface_similarity(contacts, _HETERO, _ALL)
    assert [
        (r.item_id_1, r.item_id_2, r.score, r.chain_order) for r in rows
    ] == [("1AAA_A_B", "2BBB_A_B", 1.0, "direct")]
    assert rows[0].identity_1 == rows[0].identity_2 == 1.0
    assert counts["pairs_compared"] == 1
    assert counts["alignments"] == 0
    assert set(counts) == _COUNT_KEYS


def test_row_columns_in_order():
    contacts = _iface("1AAA", "A", {1}, "B", {4}) + _iface(
        "2BBB", "A", {1}, "B", {4}
    )
    rows, _ = interface_similarity(contacts, _HETERO, _ALL)
    assert list(rows[0].model_dump()) == [
        "item_id_1",
        "item_id_2",
        "score",
        "jaccard_1",
        "jaccard_2",
        "chain_order",
        "identity_1",
        "identity_2",
    ]


def test_score_is_min_of_sides():
    contacts = _iface("1AAA", "A", {1, 2, 3}, "B", {1, 2}) + _iface(
        "2BBB", "A", {1, 2, 3}, "B", {10, 11}
    )
    rows, counts = interface_similarity(contacts, _HETERO, _ALL)
    assert rows[0].jaccard_1 == 1.0
    assert rows[0].jaccard_2 == 0.0
    assert rows[0].score == 0.0
    assert counts["pairs_zero_score"] == 1


def test_homodimer_keeps_swapped_order():
    contacts = _iface("1AAA", "A", {1, 2}, "B", {20, 21}) + _iface(
        "2BBB", "A", {20, 21}, "B", {1, 2}
    )
    rows, counts = interface_similarity(contacts, _HOMO, _ALL)
    assert rows[0].chain_order == "swapped"
    assert rows[0].score == 1.0
    assert counts["pairs_swapped"] == 1


def test_homodimer_tie_is_direct():
    contacts = _iface("1AAA", "A", {1, 2}, "B", {1, 2}) + _iface(
        "2BBB", "A", {1, 2}, "B", {1, 2}
    )
    rows, counts = interface_similarity(contacts, _HOMO, _ALL)
    assert rows[0].chain_order == "direct"
    assert counts["pairs_swapped"] == 0


def test_same_entry_interfaces_compared():
    contacts = _iface("1AAA", "A", {1, 2}, "B", {3}) + _iface(
        "1AAA", "C", {1, 2}, "D", {3}
    )
    clusters = [_cl("1AAA_A", "1AAA_C"), _cl("1AAA_B", "1AAA_D")]
    rows, _ = interface_similarity(contacts, clusters, _ALL)
    assert [(r.item_id_1, r.item_id_2) for r in rows] == [
        ("1AAA_A_B", "1AAA_C_D")
    ]


def test_reversed_rows_name_one_interface():
    # 1AAA's rows are given as (B, A); its sides must still land on A, B
    contacts = _iface("1AAA", "B", {4, 5}, "A", {1, 2, 3}) + _iface(
        "2BBB", "A", {1, 2, 3}, "B", {4, 5}
    )
    rows, counts = interface_similarity(contacts, _HETERO, _ALL)
    assert rows[0].item_id_1 == "1AAA_A_B"
    assert rows[0].jaccard_1 == 1.0
    assert rows[0].jaccard_2 == 1.0
    assert counts["interfaces"] == 2


def test_insertion_maps_through_alignment():
    # 2BBB chain A has one residue inserted after position 3
    inserted = _SEQ[:3] + "W" + _SEQ[3:]
    seqs = _ALL | {("2BBB", "A"): inserted}
    contacts = _iface("1AAA", "A", {5, 6}, "B", {1, 2}) + _iface(
        "2BBB", "A", {6, 7}, "B", {1, 2}
    )
    rows, counts = interface_similarity(contacts, _HETERO, seqs)
    assert rows[0].jaccard_1 == 1.0
    assert rows[0].identity_1 == 1.0
    assert counts["alignments"] >= 1


def test_unplaced_residue_counts_against_overlap():
    # 2BBB chain A lacks position 10 of 1AAA chain A
    deleted = _SEQ[:9] + _SEQ[10:]
    seqs = _ALL | {("2BBB", "A"): deleted}
    contacts = _iface("1AAA", "A", {9, 10}, "B", {1, 2}) + _iface(
        "2BBB", "A", {9}, "B", {1, 2}
    )
    rows, _ = interface_similarity(contacts, _HETERO, seqs)
    assert rows[0].jaccard_1 == 0.5


def _mutated(seq, every=5):
    """`seq` with every `every`-th residue replaced: 80% identity."""

    return "".join(
        ("W" if c != "W" else "Y") if i % every == 0 else c
        for i, c in enumerate(seq)
    )


def _below_identity():
    seqs = _ALL | {("2BBB", "A"): _mutated(_SEQ)}
    contacts = _iface("1AAA", "A", {2, 3}, "B", {1, 2}) + _iface(
        "2BBB", "A", {2, 3}, "B", {1, 2}
    )
    return contacts, seqs


def test_below_identity_writes_no_row():
    contacts, seqs = _below_identity()
    rows, counts = interface_similarity(contacts, _HETERO, seqs)
    assert rows == []
    assert counts["pairs_below_identity"] == 1
    assert counts["pairs_compared"] == 0


def test_min_identity_zero_lets_it_through():
    contacts, seqs = _below_identity()
    rows, _ = interface_similarity(contacts, _HETERO, seqs, min_identity=0.0)
    assert len(rows) == 1
    assert rows[0].identity_1 == 0.8


def test_different_clusters_not_compared_and_no_cluster_counted():
    contacts = (
        _iface("1AAA", "A", {1}, "B", {2})
        + _iface("2BBB", "A", {1}, "B", {2})
        + _iface("3CCC", "A", {1}, "B", {2})
    )
    clusters = [_cl("1AAA_A", "2BBB_A"), _cl("1AAA_B"), _cl("2BBB_B")]
    rows, counts = interface_similarity(contacts, clusters, _ALL)
    assert rows == []
    assert counts["interfaces"] == 3
    assert counts["interfaces_without_cluster"] == 1


def _guard_counts(row, seqs=_ALL):
    clusters = [_cl("1AAA_A"), _cl("1AAA_B")]
    _, counts = interface_similarity([row], clusters, seqs)
    return counts


def test_sequence_guard_rejects_wrong_residue():
    # position 1 of _SEQ is A, not G
    counts = _guard_counts(_row("1AAA", "A", 1, "B", 1, "GLY", "ALA"))
    assert counts["interfaces_sequence_mismatch"] == 1


def test_sequence_guard_accepts_selenomethionine():
    # position 11 of _SEQ is M; gemmi's one-letter code of MSE is "m"
    counts = _guard_counts(_row("1AAA", "A", 11, "B", 1, "MSE", "ALA"))
    assert counts["interfaces_sequence_mismatch"] == 0


def test_sequence_guard_rejects_out_of_range():
    counts = _guard_counts(_row("1AAA", "A", 51, "B", 1))
    assert counts["interfaces_sequence_mismatch"] == 1


def test_missing_sequence_counted():
    seqs = {("1AAA", "A"): _SEQ}
    counts = _guard_counts(_row("1AAA", "A", 1, "B", 1), seqs)
    assert counts["interfaces_without_sequence"] == 1
    assert counts["interfaces_sequence_mismatch"] == 0


def test_ignored_rows():
    contacts = [
        _row("1AAA", "A", None, "B", 1),
        _row("1AAA", "A", 1, "A", 2),
    ]
    _, counts = interface_similarity(contacts, _HOMO, _ALL)
    assert counts["rows_ignored"] == 2
    assert counts["interfaces"] == 0


def test_inputs_not_mutated():
    contacts = _iface("1AAA", "A", {1, 2}, "B", {20, 21}) + _iface(
        "2BBB", "A", {20, 21}, "B", {1, 2}
    )
    clusters = copy.deepcopy(_HOMO)
    seqs = dict(_ALL)
    before = (copy.deepcopy(contacts), copy.deepcopy(clusters), dict(seqs))
    interface_similarity(contacts, clusters, seqs)
    assert (contacts, clusters, seqs) == before
