from __future__ import annotations

import pytest

from pandora.schemas.similarity import SimilarityCluster
from pandora.similarity import dimer_adjacency, rebalance


def _group(*ids):
    return SimilarityCluster(components=sorted(ids), n_components=len(ids))


def _graph(*pairs, nodes=()):
    adjacency = {n: set() for pair in pairs for n in pair} | {
        n: set() for n in nodes
    }
    for a, b in pairs:
        adjacency[a].add(b)
        adjacency[b].add(a)
    return adjacency


def test_adjacency_shared_chain():
    adjacency = dimer_adjacency(["1A_A_B", "1A_B_C", "1A_D_E"], [], [])
    assert adjacency == {
        "1A_A_B": {"1A_B_C"},
        "1A_B_C": {"1A_A_B"},
        "1A_D_E": set(),
    }


def test_adjacency_expands_chain_edge():
    dimers = ["1A_A_B", "1A_A_C", "2B_C_D", "3C_E_F"]
    adjacency = dimer_adjacency(dimers, [("1A_A", "2B_C")], [])
    assert adjacency["2B_C_D"] == {"1A_A_B", "1A_A_C"}
    assert adjacency["3C_E_F"] == set()


def test_adjacency_interface_edge():
    adjacency = dimer_adjacency(
        ["1A_A_B", "2B_C_D"], [], [("1A_A_B", "2B_C_D")]
    )
    assert adjacency["1A_A_B"] == {"2B_C_D"}


def test_adjacency_chain_map():
    # 2B_C shares 1A_A's representative; the copy 1A_A-2 maps to 1A_A.
    chain_map = {"1A_A": "1A_A", "1A_A-2": "1A_A", "2B_C": "1A_A"}
    adjacency = dimer_adjacency(
        ["1A_A_A-2", "2B_C_D", "3C_E_F"], [], [], chain_map
    )
    assert adjacency["1A_A_A-2"] == {"2B_C_D"}
    assert adjacency["3C_E_F"] == set()


def test_homomer_only_group_keeps_most_connected():
    dimers = ["1A_A_B", "2B_A_B", "3C_A_B"]
    adjacency = _graph(("1A_A_B", "2B_A_B"), ("2B_A_B", "3C_A_B"))
    positives, dropped = rebalance(
        [_group(*dimers)], adjacency, set(dimers), lambda d: d
    )
    assert [(p.dimer, p.kind, p.degree) for p in positives] == [
        ("2B_A_B", "homomer", 2)
    ]
    assert {(d.dimer, d.reason) for d in dropped} == {
        ("1A_A_B", "not_most_connected"),
        ("3C_A_B", "not_most_connected"),
    }


def test_mixed_group_drops_homomers():
    adjacency = _graph(("1A_A_B", "2B_C_D"))
    positives, dropped = rebalance(
        [_group("1A_A_B", "2B_C_D")], adjacency, {"1A_A_B"}, lambda d: d
    )
    assert [p.dimer for p in positives] == ["2B_C_D"]
    assert [(d.dimer, d.reason) for d in dropped] == [
        ("1A_A_B", "homomer_in_mixed_group")
    ]


def test_heteromer_path_keeps_both_ends():
    # X - Y - Z: degrees 1, 2, 1. Keep X, drop Y; Z is then isolated.
    adjacency = _graph(("1X_A_B", "2Y_A_B"), ("2Y_A_B", "3Z_A_B"))
    positives, dropped = rebalance(
        [_group("1X_A_B", "2Y_A_B", "3Z_A_B")], adjacency, set(), lambda d: d
    )
    assert [(p.dimer, p.degree) for p in positives] == [
        ("1X_A_B", 1),
        ("3Z_A_B", 0),
    ]
    assert [(d.dimer, d.reason, d.kept) for d in dropped] == [
        ("2Y_A_B", "similar_to_kept", "1X_A_B")
    ]
    assert {p.group for p in positives} == {"1X_A_B"}
    assert positives[0].entry_id == "1X"
    assert (positives[0].chain_id_1, positives[0].chain_id_2) == ("A", "B")


def test_tie_uses_sort_key():
    adjacency = _graph(nodes=("1A_A_B", "2B_A_B"))
    order = {"2B_A_B": 0, "1A_A_B": 1}
    positives, _ = rebalance(
        [_group("1A_A_B", "2B_A_B")], adjacency, set(), order.__getitem__
    )
    assert [p.dimer for p in positives] == ["2B_A_B", "1A_A_B"]


def test_single_dimer_group():
    positives, dropped = rebalance(
        [_group("1A_A_B")], _graph(nodes=("1A_A_B",)), set(), lambda d: d
    )
    assert [p.dimer for p in positives] == ["1A_A_B"]
    assert dropped == []


def test_cross_group_edge_raises():
    adjacency = _graph(("1A_A_B", "2B_C_D"))
    with pytest.raises(ValueError, match="1A_A_B"):
        rebalance(
            [_group("1A_A_B"), _group("2B_C_D")], adjacency, set(), lambda d: d
        )
