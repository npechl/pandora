from __future__ import annotations

import pytest

from pandora.schemas.similarity import PositiveDimer
from pandora.similarity import sample_negatives

# Eight one-chain proteins R1..R8, each its own representative.
REPS = {f"{i}R_A": f"{i}R_A" for i in range(1, 9)}


def _dimer(entry, c1, c2, reps):
    """A positive in `entry` with chains c1, c2, mapped in `reps`."""

    return PositiveDimer(
        dimer=f"{entry}_{c1}_{c2}",
        entry_id=entry,
        chain_id_1=c1,
        chain_id_2=c2,
        group=f"{entry}_{c1}_{c2}",
        kind="heteromer",
        degree=0,
    )


def _reps_with(extra):
    return REPS | extra


def test_partner_never_negative():
    # Positive 9P_A–9P_B: A's rep 1R_A, B's rep 2R_A. Even if A is picked
    # as "A", B (or A) is never C; both reps leave the pool.
    reps = _reps_with({"9P_A": "1R_A", "9P_B": "2R_A"})
    positive = _dimer("9P", "A", "B", reps)
    for seed in range(20):
        negatives, unpaired = sample_negatives(
            [positive], reps, [], [("9P_A", "9P_B")], seed
        )
        assert unpaired == []
        assert negatives[0].chain_2 not in {"1R_A", "2R_A"}
        assert negatives[0].chain_1 in {"9P_A", "9P_B"}
        assert negatives[0].positive == "9P_A_B"


def test_sets_one_to_three_excluded():
    # A -> 1R; TM 1R–3R (set 1), 3R binds 4R (set 2), TM 4R–5R (set 3).
    reps = _reps_with({"9P_A": "1R_A", "9P_B": "1R_A"})
    positive = _dimer("9P", "A", "B", reps)
    tm = [("1R_A", "3R_A"), ("4R_A", "5R_A")]
    binding = [("3R_A", "4R_A")]
    for seed in range(20):
        negatives, _ = sample_negatives([positive], reps, tm, binding, seed)
        assert negatives[0].chain_2 in {"2R_A", "6R_A", "7R_A", "8R_A"}


def test_c_used_once():
    reps = {
        "1R_A": "1R_A",
        "2R_A": "2R_A",
        "3R_A": "3R_A",
        "4R_A": "4R_A",
        "9P_A": "1R_A",
        "9P_B": "1R_A",
        "9Q_A": "2R_A",
        "9Q_B": "2R_A",
    }
    positives = [_dimer("9P", "A", "B", reps), _dimer("9Q", "A", "B", reps)]
    negatives, unpaired = sample_negatives(positives, reps, [], [], 0)
    assert unpaired == []
    assert sorted(n.chain_2 for n in negatives) == ["3R_A", "4R_A"]


def test_pool_exhaustion_records_unpaired():
    reps = {
        "1R_A": "1R_A",
        "3R_A": "3R_A",
        "9P_A": "1R_A",
        "9P_B": "1R_A",
        "9Q_A": "1R_A",
        "9Q_B": "1R_A",
    }
    positives = [_dimer("9P", "A", "B", reps), _dimer("9Q", "A", "B", reps)]
    negatives, unpaired = sample_negatives(positives, reps, [], [], 0)
    assert [n.chain_2 for n in negatives] == ["3R_A"]
    assert unpaired == ["9Q_A_B"]


def test_same_seed_same_output():
    reps = _reps_with({"9P_A": "1R_A", "9P_B": "2R_A"})
    positives = [_dimer("9P", "A", "B", reps)]
    first = sample_negatives(positives, reps, [], [], 7)
    assert sample_negatives(positives, reps, [], [], 7) == first


def test_reversed_and_duplicate_binding_rows_count_once():
    reps = _reps_with({"9P_A": "1R_A", "9P_B": "1R_A"})
    positives = [_dimer("9P", "A", "B", reps)]
    once = sample_negatives(positives, reps, [], [("1R_A", "3R_A")], 3)
    many = sample_negatives(
        positives,
        reps,
        [],
        [("1R_A", "3R_A"), ("3R_A", "1R_A"), ("1R_A", "3R_A")],
        3,
    )
    assert many == once


def test_tm_edge_outside_reps_mapped_through_members():
    # The TM layer names 8X_A, a member of 3R_A, not a representative.
    reps = _reps_with({"9P_A": "1R_A", "9P_B": "1R_A", "8X_A": "3R_A"})
    positives = [_dimer("9P", "A", "B", reps)]
    for seed in range(20):
        negatives, _ = sample_negatives(
            positives, reps, [("9P_A", "8X_A")], [], seed
        )
        assert negatives[0].chain_2 != "3R_A"


def test_unmapped_positive_chain_raises():
    positives = [_dimer("9P", "A", "B", REPS)]
    with pytest.raises(ValueError, match="9P_A"):
        sample_negatives(positives, REPS, [], [], 0)


def test_empty_positives():
    assert sample_negatives([], REPS, [], [], 0) == ([], [])
