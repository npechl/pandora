from collections import Counter

from pandora.datasets import find_contained
from pandora.schemas.common import RankKey

RES = [RankKey(key="resolution")]


def _counts(**entries):
    """Entry id -> chain counts per cluster: "AA" = two chains in a."""

    return {e: Counter(s.lower()) for e, s in entries.items()}


def _run(counts, res=None):
    res = res or {}
    values = {e: {"resolution": res.get(e)} for e in counts}
    return {
        r.entry_id: (r.contained_in, r.equal)
        for r in find_contained(counts, values, RES)
    }


def test_monomer_contained_in_homodimer():
    assert _run(_counts(E1="A", E2="AA")) == {"E1": ("E2", False)}


def test_homodimer_not_contained_in_monomer():
    got = _run(_counts(E1="AA", E2="A"), {"E2": 1.0, "E1": 3.0})
    assert got == {"E2": ("E1", False)}


def test_heterodimer_contained_in_larger_complex():
    assert _run(_counts(E1="AB", E2="AAB")) == {"E1": ("E2", False)}


def test_unrelated_entries_both_kept():
    assert _run(_counts(E1="AB", E2="AC")) == {}


def test_equal_counts_keep_best_resolution():
    got = _run(_counts(E1="AB", E2="AB"), {"E1": 2.5, "E2": 1.5})
    assert got == {"E1": ("E2", True)}


def test_equal_counts_and_resolution_keep_lowest_id():
    got = _run(_counts(E2="AB", E1="AB"), {"E1": 2.0, "E2": 2.0})
    assert got == {"E2": ("E1", True)}


def test_missing_resolution_loses_a_tie():
    got = _run(_counts(E1="AB", E2="AB"), {"E2": 3.0})
    assert got == {"E1": ("E2", True)}


def test_three_equal_entries_keep_one():
    got = _run(_counts(E1="A", E2="A", E3="A"), {"E1": 3, "E2": 1, "E3": 2})
    assert got == {"E1": ("E2", True), "E3": ("E2", True)}


def test_contained_in_names_a_kept_entry():
    got = _run(_counts(E1="A", E2="AB", E3="ABC"))
    assert got == {"E1": ("E3", False), "E2": ("E3", False)}


def test_strict_container_in_a_tie_group():
    # E1 is contained in both tied entries; it names the kept one.
    got = _run(_counts(E1="A", E2="AB", E3="AB"), {"E2": 2.0, "E3": 1.0})
    assert got == {"E1": ("E3", False), "E2": ("E3", True)}


def test_entry_without_chains_is_kept():
    assert _run(_counts(E1="", E2="AA")) == {}


def test_inputs_not_mutated():
    counts = _counts(E1="A", E2="AA")
    before = {k: dict(v) for k, v in counts.items()}
    _run(counts)
    assert {k: dict(v) for k, v in counts.items()} == before
