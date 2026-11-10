from __future__ import annotations

import pytest

from pandora.schemas.similarity import NegativePair, PositiveDimer
from pandora.similarity import assemble_pairs


def _pos(entry, kind="heteromer", group=None):
    dimer = f"{entry}_A_B"
    return PositiveDimer(
        dimer=dimer,
        entry_id=entry,
        chain_id_1="A",
        chain_id_2="B",
        group=group or dimer,
        kind=kind,
        degree=0,
    )


POSITIVES = [_pos("1P", kind="homomer", group="1P_A_B"), _pos("2Q")]
NEGATIVES = [
    NegativePair(chain_1="1P_A", chain_2="7C_A", positive="1P_A_B"),
    NegativePair(chain_1="2Q_B", chain_2="8D_A", positive="2Q_A_B"),
]
SPLITS = {"train": ["1P_A_B"], "val": [], "test": ["2Q_A_B"]}
REPS = {c: c for c in ["1P_A", "1P_B", "2Q_A", "2Q_B", "7C_A", "8D_A"]}


def _rows(tm=(), negatives=NEGATIVES, splits=SPLITS, reps=REPS):
    return assemble_pairs(POSITIVES, negatives, splits, reps, tm)


def test_negatives_inherit_split_kind_cluster():
    rows = _rows()
    assert [(r.chain_a, r.chain_b, r.label, r.split) for r in rows] == [
        ("1P_A", "1P_B", 1, "train"),
        ("2Q_A", "2Q_B", 1, "test"),
        ("1P_A", "7C_A", 0, "train"),
        ("2Q_B", "8D_A", 0, "test"),
    ]
    negative = rows[2]
    assert (negative.kind, negative.cluster, negative.positive) == (
        "homomer",
        "1P_A_B",
        "1P_A_B",
    )
    assert rows[0].c_leak_positive is None
    assert rows[0].c_leak_negative is None


def test_c_leak_positive():
    rows = _rows(tm=[("7C_A", "2Q_A")])
    assert rows[2].c_leak_positive is True
    assert rows[2].c_leak_negative is False
    assert rows[3].c_leak_positive is False


def test_c_leak_negative():
    rows = _rows(tm=[("7C_A", "8D_A")])
    assert rows[2].c_leak_negative is True
    assert rows[3].c_leak_negative is True
    assert rows[2].c_leak_positive is False


def test_same_split_link_is_not_a_leak():
    rows = _rows(tm=[("7C_A", "1P_B")])
    assert rows[2].c_leak_positive is False


def test_no_tm_edges_no_leaks():
    rows = _rows()
    assert [r.c_leak_positive for r in rows[2:]] == [False, False]
    assert [r.c_leak_negative for r in rows[2:]] == [False, False]


def test_tm_edge_through_members():
    reps = REPS | {"9M_A": "8D_A"}
    rows = _rows(tm=[("7C_A", "9M_A")], reps=reps)
    assert rows[2].c_leak_negative is True


def test_positive_in_no_split_raises():
    with pytest.raises(ValueError, match="2Q_A_B"):
        _rows(splits={"train": ["1P_A_B"], "val": [], "test": []})


def test_positive_in_two_splits_raises():
    with pytest.raises(ValueError, match="1P_A_B"):
        _rows(
            splits={"train": ["1P_A_B"], "val": ["1P_A_B"], "test": ["2Q_A_B"]}
        )


def test_split_item_not_a_positive_raises():
    splits = SPLITS | {"val": ["9Z_A_B"]}
    with pytest.raises(ValueError, match="9Z_A_B"):
        _rows(splits=splits)


def test_negative_with_unknown_positive_raises():
    negatives = [NegativePair(chain_1="X_A", chain_2="7C_A", positive="5X_A_B")]
    with pytest.raises(ValueError, match="5X_A_B"):
        _rows(negatives=negatives)


def test_unmapped_chain_raises():
    reps = {k: v for k, v in REPS.items() if k != "8D_A"}
    with pytest.raises(ValueError, match="8D_A"):
        _rows(reps=reps)


def test_no_negatives_keeps_positives():
    rows = _rows(negatives=[])
    assert [(r.label, r.split) for r in rows] == [(1, "train"), (1, "test")]
