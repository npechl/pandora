from __future__ import annotations

import pyarrow as pa
import pytest

from pandora.schemas.similarity import Edge, HitFilter
from pandora.similarity.hits import (
    read_hit_tsv,
    read_similarity_table,
    similarity_edges,
)


def _foldseek_row(
    query,
    target,
    *,
    fident=0.5,
    qcov=0.9,
    tcov=0.9,
    alntm=0.6,
    qtm=0.6,
    ttm=0.6,
    qstart=1,
    qend=100,
    tstart=1,
    tend=100,
):
    values = (
        query,
        target,
        fident,
        100,
        qcov,
        tcov,
        alntm,
        qtm,
        ttm,
        qstart,
        qend,
        tstart,
        tend,
    )
    return "\t".join(str(value) for value in values)


def _batches(tmp_path, rows, engine="Foldseek"):
    path = tmp_path / "hits.tsv"
    path.write_text("".join(row + "\n" for row in rows))
    return list(read_hit_tsv(path, engine))


def _edges(tmp_path, rows, hit_filter=None, **kwargs):
    batches = _batches(tmp_path, rows)
    return list(similarity_edges(batches, hit_filter or HitFilter(), **kwargs))


def _table(**columns):
    return pa.RecordBatch.from_pydict(columns)


def test_similarity_edges_sorts_pair_and_skips_self_hits(tmp_path):
    edges = _edges(tmp_path, [_foldseek_row("b", "a"), _foldseek_row("a", "a")])

    assert edges == [Edge("a", "b", 0.6, 0.5, 0.9, None)]


def test_similarity_edges_keeps_every_passing_row(tmp_path):
    edges = _edges(tmp_path, [_foldseek_row("a", "b"), _foldseek_row("b", "a")])

    assert [(e.source_id, e.target_id) for e in edges] == [
        ("a", "b"),
        ("a", "b"),
    ]


@pytest.mark.parametrize(
    ("field", "row_kwargs"),
    [
        ("min_score", {"alntm": 0.6}),
        ("min_identity", {"fident": 0.6}),
        ("min_coverage", {"qcov": 0.6, "tcov": 0.6}),
    ],
)
def test_similarity_edges_thresholds_are_inclusive(tmp_path, field, row_kwargs):
    rows = [_foldseek_row("a", "b", **row_kwargs)]

    kept = _edges(tmp_path, rows, HitFilter(**{field: 0.6}))
    dropped = _edges(tmp_path, rows, HitFilter(**{field: 0.61}))

    assert len(kept) == 1
    assert dropped == []


def test_similarity_edges_any_passing_row_makes_an_edge(tmp_path):
    # Best-scoring alignment fails coverage; a second alignment passes.
    rows = [
        _foldseek_row("a", "b", alntm=0.9, qcov=0.5, tcov=0.5),
        _foldseek_row("a", "b", alntm=0.55, qcov=0.85, tcov=0.85),
    ]

    edges = _edges(tmp_path, rows, HitFilter(min_score=0.5, min_coverage=0.8))

    assert [e.score for e in edges] == [0.55]


@pytest.mark.parametrize(
    ("coverage_of", "expected"),
    [("both", 0.4), ("side_1", 0.9), ("side_2", 0.4), ("either", 0.9)],
)
def test_similarity_edges_coverage_of(tmp_path, coverage_of, expected):
    rows = [_foldseek_row("a", "b", qcov=0.9, tcov=0.4)]

    (edge,) = _edges(tmp_path, rows, HitFilter(coverage_of=coverage_of))

    assert edge.coverage == expected


def test_similarity_edges_score_column(tmp_path):
    rows = [_foldseek_row("a", "b", alntm=0.1, qtm=0.7)]

    (edge,) = _edges(
        tmp_path, rows, HitFilter(min_score=0.6, score_column="qtmscore")
    )

    assert edge.score == 0.7


def test_similarity_edges_interface_coverage(tmp_path):
    rows = [
        _foldseek_row("a", "b", qstart=1, qend=15, tstart=1, tend=10),
        _foldseek_row("a", "c"),
    ]
    interfaces = {"a": {10, 20}, "b": {5}}

    every = _edges(tmp_path, rows, interface_residues=interfaces)
    kept = _edges(
        tmp_path,
        rows,
        HitFilter(min_interface_coverage=0.5),
        interface_residues=interfaces,
    )
    dropped = _edges(
        tmp_path,
        rows,
        HitFilter(min_interface_coverage=0.6),
        interface_residues=interfaces,
    )

    # a covers 1 of {10, 20} -> 0.5; b covers {5} -> 1.0; min -> 0.5.
    assert [e.interface_coverage for e in every] == [0.5, None]
    assert [(e.source_id, e.target_id) for e in kept] == [("a", "b")]
    assert dropped == []


def test_similarity_edges_minimal_layer():
    batch = _table(item_id_1=["b"], item_id_2=["a"], score=[0.9])

    edges = list(similarity_edges([batch], HitFilter(min_score=0.5)))

    assert edges == [Edge("a", "b", 0.9, None, None, None)]


def test_similarity_edges_null_coverage_fails_threshold():
    batch = _table(
        item_id_1=["a", "a"],
        item_id_2=["b", "c"],
        score=[0.9, 0.9],
        coverage_1=pa.array([None, 0.5], pa.float64()),
        coverage_2=[0.5, 0.5],
    )

    edges = list(similarity_edges([batch], HitFilter(min_coverage=0.1)))

    assert [e.target_id for e in edges] == ["c"]


@pytest.mark.parametrize(
    ("hit_filter", "missing"),
    [
        (HitFilter(min_identity=0.1), "identity"),
        (HitFilter(min_coverage=0.1), "coverage_1"),
        (HitFilter(min_coverage=0.1, coverage_of="side_2"), "coverage_2"),
        (HitFilter(min_score=0.1, score_column="qtmscore"), "qtmscore"),
    ],
)
def test_similarity_edges_threshold_on_absent_column_raises(
    hit_filter, missing
):
    batch = _table(item_id_1=["a"], item_id_2=["b"], score=[0.9])

    with pytest.raises(ValueError, match=missing):
        next(similarity_edges([batch], hit_filter))


def test_similarity_edges_interface_filter_without_residues_raises(tmp_path):
    batches = _batches(tmp_path, [_foldseek_row("a", "b")])

    with pytest.raises(ValueError, match="interface_residues"):
        next(similarity_edges(batches, HitFilter(min_interface_coverage=0.5)))


def test_similarity_edges_interface_residues_need_ranges():
    batch = _table(item_id_1=["a"], item_id_2=["b"], score=[0.9])

    with pytest.raises(ValueError, match="start_1"):
        next(
            similarity_edges(
                [batch], HitFilter(), interface_residues={"a": {1}}
            )
        )


@pytest.mark.parametrize(
    ("field", "column"),
    [
        ("min_score", "score"),
        ("min_identity", "identity"),
        ("min_coverage", "coverage_1"),
    ],
)
def test_similarity_edges_nan_fails_threshold(field, column):
    columns = dict(
        item_id_1=["a"],
        item_id_2=["b"],
        score=[0.9],
        identity=[0.9],
        coverage_1=[0.9],
        coverage_2=[0.9],
    )
    columns[column] = [float("nan")]

    edges = list(
        similarity_edges([_table(**columns)], HitFilter(**{field: 0.1}))
    )

    assert edges == []


def test_similarity_edges_batching_does_not_change_edges(tmp_path):
    rows = [_foldseek_row(f"q{i}", f"t{i}", alntm=i / 10) for i in range(9)]
    hit_filter = HitFilter(min_score=0.3)
    one = list(similarity_edges(_batches(tmp_path, rows), hit_filter))
    table = pa.Table.from_batches(_batches(tmp_path, rows))

    three = list(
        similarity_edges(table.to_batches(max_chunksize=3), hit_filter)
    )

    assert one == three
    assert len(one) == 6


def test_read_hit_tsv_maps_foldseek_columns(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text(_foldseek_row("a_A", "b_A", alntm=0.7, qtm=0.4) + "\n")

    batch = next(read_hit_tsv(path, "Foldseek"))

    row = batch.to_pylist()[0]
    assert row["item_id_1"] == "a_A"
    assert row["item_id_2"] == "b_A"
    assert row["score"] == row["alntmscore"] == 0.7
    assert row["qtmscore"] == 0.4
    assert row["identity"] == 0.5
    assert (row["coverage_1"], row["start_2"], row["end_2"]) == (0.9, 1, 100)
    assert "query" not in row


def test_read_hit_tsv_mmseqs_score_is_identity(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("a_A\tb_A\t0.3\t10\t0.9\t0.8\n")

    row = next(read_hit_tsv(path, "MMseqs2")).to_pylist()[0]

    assert row["score"] == row["identity"] == 0.3
    assert (row["coverage_1"], row["coverage_2"]) == (0.9, 0.8)


def test_read_hit_tsv_unknown_engine(tmp_path):
    with pytest.raises(ValueError, match="engine"):
        list(read_hit_tsv(tmp_path / "x.tsv", "BLAST"))


def test_read_similarity_table_tsv_reads_values_as_text(tmp_path):
    path = tmp_path / "edges.tsv"
    path.write_text("item_id_1\titem_id_2\tscore\tnote\n104M\t112M\t0.5\tk\n")

    row = next(read_similarity_table(path)).to_pylist()[0]

    assert row == {
        "item_id_1": "104M",
        "item_id_2": "112M",
        "score": "0.5",
        "note": "k",
    }


def test_read_similarity_table_bad_extra_column_names_file(tmp_path):
    path = tmp_path / "edges.tsv"
    lines = ["item_id_1\titem_id_2\tscore\tnote_num"]
    lines += [f"104M_A\t112M_A\t0.5\t{i}" for i in range(20000)]
    lines += ["104M_A\t112M_A\t0.5\tunknown"]  # pyarrow reads n/a as null
    path.write_text("\n".join(lines) + "\n")

    with pytest.raises(ValueError, match="edges.tsv"):
        list(read_similarity_table(path, block_size=1 << 12))


def test_read_similarity_table_rejects_other_suffix(tmp_path):
    with pytest.raises(ValueError, match="parquet"):
        list(read_similarity_table(tmp_path / "edges.csv"))
