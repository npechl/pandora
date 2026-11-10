import pytest

from pandora.schemas.common import RankKey
from pandora.schemas.similarity import RepresentativePolicy, SimilarityCluster
from pandora.similarity import (
    cluster_hit_counts,
    pick_representatives,
    restrict_clusters,
)


def _c(*ids):
    return SimilarityCluster(components=list(ids), n_components=len(ids))


POLICY = RepresentativePolicy()  # missing_residues min, n_hits max


def _v(missing, hits):
    return {"missing_residues": missing, "n_hits": hits}


def test_fewest_missing_residues_wins():
    reps = pick_representatives(
        [_c("A_A", "B_A")], {"A_A": _v(5, 9), "B_A": _v(1, 0)}, POLICY
    )
    assert reps[0].representative == "B_A"
    assert reps[0].members == ["A_A", "B_A"] and reps[0].n_members == 2


def test_n_hits_then_lowest_id_break_ties():
    values = {"A_A": _v(1, 1), "B_A": _v(1, 2), "C_A": _v(1, 2)}
    reps = pick_representatives([_c("A_A", "B_A", "C_A")], values, POLICY)
    assert reps[0].representative == "B_A"


def test_null_missing_residues_ranks_last():
    values = {"A_A": _v(None, 5), "B_A": _v(9, 0)}
    reps = pick_representatives([_c("A_A", "B_A")], values, POLICY)
    assert reps[0].representative == "B_A"


def test_copy_is_never_picked():
    values = {"A_A": _v(3, 0), "A_A-2": _v(0, 0)}
    reps = pick_representatives(
        [_c("A_A", "A_A-2")], values, POLICY, eligible={"A_A"}
    )
    assert reps[0].representative == "A_A"
    assert reps[0].members == ["A_A", "A_A-2"]


def test_copies_only_cluster_falls_back_to_all_members():
    reps = pick_representatives(
        [_c("A_A-2")], {"A_A-2": _v(0, 0)}, POLICY, eligible=set()
    )
    assert reps[0].representative == "A_A-2"


def test_results_sorted_by_representative():
    values = {i: _v(0, 0) for i in ("A_A", "B_A", "C_A")}
    reps = pick_representatives([_c("C_A"), _c("A_A", "B_A")], values, POLICY)
    assert [r.representative for r in reps] == ["A_A", "C_A"]


def test_unknown_key_is_an_error():
    policy = RepresentativePolicy(rank=[RankKey(key="nope")])
    with pytest.raises(ValueError, match="nope"):
        pick_representatives([_c("A_A")], {"A_A": _v(0, 0)}, policy)


def test_restrict_drops_members_and_empty_clusters():
    got = restrict_clusters([_c("A_A", "B_A"), _c("C_A")], {"A_A"})
    assert [c.components for c in got] == [["A_A"]]
    assert got[0].n_components == 1


def test_hit_counts_stay_inside_clusters():
    clusters = [_c("A_A", "B_A", "C_A"), _c("D_A")]
    edges = [
        ("A_A", "B_A"),
        ("B_A", "A_A"),
        ("A_A", "C_A"),
        ("A_A", "D_A"),
        ("A_A", "A_A"),
        ("A_A", "Z_A"),
    ]
    assert cluster_hit_counts(clusters, edges) == {
        "A_A": 2,
        "B_A": 1,
        "C_A": 1,
        "D_A": 0,
    }


def test_inputs_not_mutated():
    clusters = [_c("A_A", "B_A")]
    restrict_clusters(clusters, {"A_A"})
    assert clusters[0].components == ["A_A", "B_A"]
