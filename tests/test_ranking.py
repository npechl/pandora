from pandora._util import rank_sort_key
from pandora.schemas.common import RankKey

MIN = [RankKey(key="m")]
MAX = [RankKey(key="m", order="max")]


def _best(values, rank):
    return min(values, key=lambda i: rank_sort_key(i, values[i], rank))


def test_min_order_prefers_smaller():
    assert _best({"B": {"m": 1}, "A": {"m": 2}}, MIN) == "B"


def test_max_order_prefers_larger():
    assert _best({"B": {"m": 1}, "A": {"m": 2}}, MAX) == "A"


def test_null_ranks_last_for_both_orders():
    values = {"A": {"m": None}, "B": {"m": 5}}
    assert _best(values, MIN) == "B"
    assert _best(values, MAX) == "B"


def test_ties_fall_back_to_lowest_id():
    assert _best({"B": {"m": 1}, "A": {"m": 1}}, MIN) == "A"


def test_later_keys_break_earlier_ties():
    rank = [RankKey(key="m"), RankKey(key="n", order="max")]
    values = {"A": {"m": 1, "n": 1}, "B": {"m": 1, "n": 3}}
    assert _best(values, rank) == "B"
