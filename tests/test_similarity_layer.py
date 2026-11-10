import pyarrow as pa
import pytest

from pandora.schemas.annotation import LayerRow
from pandora.similarity.layer import (
    chain_items,
    check_similarity_columns,
    entry_items,
    interface_items,
    item_id,
    item_labels,
    match_similarity_batches,
    parse_item_id,
)


def _batch(**columns):
    return pa.RecordBatch.from_pydict(columns)


def _match(batch, items):
    return next(match_similarity_batches([batch], items))


def test_item_id_round_trip_and_interface_sorted():
    assert item_id("interface", "1ABC", "B", "A") == "1ABC_A_B"
    assert item_id("chain", "1ABC", "A") == "1ABC_A"
    assert item_id("entry", "1ABC") == "1ABC"
    assert parse_item_id("chain", "pdb_00001abc_A") == ("pdb_00001abc", "A")
    assert parse_item_id("interface", "1ABC_A_B") == ("1ABC", "A", "B")
    with pytest.raises(ValueError):
        parse_item_id("interface", "1ABC_A")
    with pytest.raises(ValueError):
        parse_item_id("chain", "1ABC_")
    with pytest.raises(ValueError):
        item_id("chain", "1ABC")


def test_item_id_rejects_chain_ids_with_underscores():
    # `1AYI_A_2` would parse back as entry `1AYI_A`, chain `2`.
    with pytest.raises(ValueError, match="_"):
        item_id("chain", "1AYI", "A_2")
    with pytest.raises(ValueError, match="_"):
        item_id("interface", "1AYI", "A", "A_2")


def test_chain_match_is_case_insensitive_on_entry_only():
    items = chain_items([("104M", "A"), ("112M", "A")])
    batch = _batch(
        item_id_1=["104m_A", "104M_a", "999X_A", "112M"],
        item_id_2=["112M_A"] * 4,
        score=[1, 1, 1, 1],
    )
    matched, unmatched = _match(batch, items)
    assert matched.column("item_id_1").to_pylist() == ["104M_A"]
    assert matched.schema.field("score").type == pa.float64()
    assert unmatched.column("reason").to_pylist() == [
        "chain_not_found",
        "entry_not_found",
        "bad_item_id",
    ]


def test_side_2_reason_when_side_1_matches():
    items = chain_items([("104M", "A")])
    matched, unmatched = _match(
        _batch(item_id_1=["104M_A"], item_id_2=["104M_B"], score=[0.5]),
        items,
    )
    assert matched.num_rows == 0
    assert unmatched.column("reason").to_pylist() == ["chain_not_found"]
    assert unmatched.column("score").to_pylist() == ["0.5"]


def test_interface_items_come_from_contacts_only():
    rows = [
        dict(
            entry_id="1A2C",
            chain_id_1="B",
            chain_id_2="A",
            seq_id_1=5,
            seq_id_2=7,
        ),
        dict(
            entry_id="1A2C",
            chain_id_1="B",
            chain_id_2="E",
            seq_id_1=5,
            seq_id_2=None,
        ),
        dict(
            entry_id="1A2C",
            chain_id_1="C",
            chain_id_2="C",
            seq_id_1=1,
            seq_id_2=9,
        ),
    ]
    items = interface_items(rows)
    batch = _batch(
        item_id_1=["1A2C_B_A", "1A2C_B_C", "1A2C_B_E", "9ZZZ_A_B"],
        item_id_2=["1a2c_A_B"] * 4,
        score=[0.9] * 4,
    )
    matched, unmatched = _match(batch, items)
    assert matched.column("item_id_1").to_pylist() == ["1A2C_A_B"]
    assert matched.column("item_id_2").to_pylist() == ["1A2C_A_B"]
    assert unmatched.column("reason").to_pylist() == [
        "interface_not_found",
        "interface_not_found",
        "entry_not_found",
    ]


def test_entry_items_and_bad_value():
    items = entry_items(["104M", "112M"])
    batch = _batch(
        item_id_1=["104M", "104M", "104M"],
        item_id_2=["112M", "112M", "112M"],
        score=["0.5", "x", None],
        keep=["a", "b", "c"],
    )
    matched, unmatched = _match(batch, items)
    assert matched.column("score").to_pylist() == [0.5]
    assert matched.column("keep").to_pylist() == ["a"]
    assert unmatched.column("reason").to_pylist() == ["bad_value", "bad_value"]


def test_recognised_columns_are_cast():
    items = entry_items(["104M", "112M"])
    batch = _batch(
        item_id_1=["104M"],
        item_id_2=["112M"],
        score=[1],
        coverage_1=["0.75"],
        start_1=["3"],
    )
    matched, _ = _match(batch, items)
    assert matched.schema.field("coverage_1").type == pa.float64()
    assert matched.schema.field("start_1").type == pa.int64()
    assert matched.to_pylist()[0]["start_1"] == 3


def test_self_hits_and_repeats_are_kept():
    items = chain_items([("104M", "A")])
    batch = _batch(
        item_id_1=["104M_A"] * 2, item_id_2=["104M_A"] * 2, score=[1.0, 1.0]
    )
    matched, _ = _match(batch, items)
    assert matched.num_rows == 2


def test_check_similarity_columns_names_missing():
    with pytest.raises(ValueError, match="score"):
        check_similarity_columns({"item_id_1", "item_id_2"})
    check_similarity_columns({"item_id_1", "item_id_2", "score", "x"})


def _chain_label(entry, chain, label):
    return LayerRow(entry_id=entry, chain_id=chain, label=label)


def _interface_label(entry, chain_1, chain_2, label):
    return LayerRow(
        entry_id=entry, chain_id_1=chain_1, chain_id_2=chain_2, label=label
    )


def test_item_labels_chains():
    rows = [_chain_label("1ABC", "A", "x"), _chain_label("1ABC", "B", "y")]
    assert item_labels(rows, "chain_labels") == {"1ABC_A": "x", "1ABC_B": "y"}


def test_item_labels_interface_order():
    rows = [_interface_label("1ABC", "B", "A", "het")]
    assert item_labels(rows, "interface_labels") == {"1ABC_A_B": "het"}


def test_item_labels_extended_entry():
    rows = [{"entry_id": "pdb_00001abc", "chain_id": "A", "label": "x"}]
    assert item_labels(rows, "chain_labels") == {"pdb_00001abc_A": "x"}


def test_item_labels_conflict_raises():
    rows = [
        _interface_label("1ABC", "B", "A", "x"),
        _interface_label("1ABC", "A", "B", "y"),
        _interface_label("1ABC", "A", "C", "y"),
        _interface_label("1ABC", "A", "C", "y"),
    ]
    with pytest.raises(ValueError, match="1ABC_A_B"):
        item_labels(rows, "interface_labels")


def test_item_labels_wrong_template_raises():
    with pytest.raises(ValueError, match="residue_labels"):
        item_labels([], "residue_labels")
