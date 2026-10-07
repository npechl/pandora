from __future__ import annotations

import pytest

from pandora.schemas.similarity import Edge, HitFilter, SimilaritySearch
from pandora.similarity.hits import (
    HIT_COLUMNS,
    iter_edges,
    load_similarity_search,
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


def _search(tmp_path, rows, engine="Foldseek", newline="\n"):
    path = tmp_path / "hits.tsv"
    path.write_text("".join(row + newline for row in rows), newline="")
    return SimilaritySearch(
        engine=engine, hits_path=str(path), columns=list(HIT_COLUMNS[engine])
    )


def test_iter_edges_sorts_pair_and_skips_self_hits(tmp_path):
    search = _search(
        tmp_path, [_foldseek_row("b", "a"), _foldseek_row("a", "a")]
    )

    edges = list(iter_edges(search, HitFilter()))

    assert edges == [Edge("a", "b", 0.6, 0.5, 0.9, None)]


def test_iter_edges_keeps_every_passing_row(tmp_path):
    search = _search(
        tmp_path, [_foldseek_row("a", "b"), _foldseek_row("b", "a")]
    )

    edges = list(iter_edges(search, HitFilter()))

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
def test_iter_edges_thresholds_are_inclusive(tmp_path, field, row_kwargs):
    search = _search(tmp_path, [_foldseek_row("a", "b", **row_kwargs)])

    kept = list(iter_edges(search, HitFilter(**{field: 0.6})))
    dropped = list(iter_edges(search, HitFilter(**{field: 0.61})))

    assert len(kept) == 1
    assert dropped == []


def test_iter_edges_any_passing_row_makes_an_edge(tmp_path):
    # Best-scoring alignment fails coverage; a second alignment passes.
    search = _search(
        tmp_path,
        [
            _foldseek_row("a", "b", alntm=0.9, qcov=0.5, tcov=0.5),
            _foldseek_row("a", "b", alntm=0.55, qcov=0.85, tcov=0.85),
        ],
    )

    edges = list(iter_edges(search, HitFilter(min_score=0.5, min_coverage=0.8)))

    assert [e.score for e in edges] == [0.55]


@pytest.mark.parametrize(
    ("coverage_of", "expected"),
    [("both", 0.5), ("query", 0.9), ("target", 0.5), ("either", 0.9)],
)
def test_iter_edges_coverage_of(tmp_path, coverage_of, expected):
    search = _search(tmp_path, [_foldseek_row("a", "b", qcov=0.9, tcov=0.5)])

    (edge,) = iter_edges(search, HitFilter(coverage_of=coverage_of))

    assert edge.coverage == expected


@pytest.mark.parametrize(
    ("normalisation", "expected"),
    [("alignment", 0.5), ("query", 0.7), ("target", 0.6), ("max", 0.7)],
)
def test_iter_edges_tm_normalisation(tmp_path, normalisation, expected):
    search = _search(
        tmp_path, [_foldseek_row("a", "b", alntm=0.5, qtm=0.7, ttm=0.6)]
    )

    (edge,) = iter_edges(search, HitFilter(tm_normalisation=normalisation))

    assert edge.score == expected


def test_iter_edges_interface_coverage_at_read_time(tmp_path):
    search = _search(
        tmp_path,
        [
            _foldseek_row("a", "b", qstart=1, qend=15, tstart=1, tend=10),
            _foldseek_row("a", "c"),
        ],
    )
    interfaces = {"a": {10, 20}, "b": {5}}

    every = list(iter_edges(search, HitFilter(), interface_residues=interfaces))
    kept = list(
        iter_edges(
            search,
            HitFilter(min_interface_coverage=0.5),
            interface_residues=interfaces,
        )
    )
    dropped = list(
        iter_edges(
            search,
            HitFilter(min_interface_coverage=0.6),
            interface_residues=interfaces,
        )
    )

    # a covers 1 of {10, 20} -> 0.5; b covers {5} -> 1.0; min -> 0.5.
    assert [e.interface_coverage for e in every] == [0.5, None]
    assert [(e.source_id, e.target_id) for e in kept] == [("a", "b")]
    assert dropped == []


def test_iter_edges_mmseqs_score_is_identity(tmp_path):
    search = _search(tmp_path, ["a\tb\t0.95\t100\t0.9\t0.8"], engine="MMseqs2")

    edges = list(iter_edges(search, HitFilter()))

    assert edges == [Edge("a", "b", 0.95, 0.95, 0.8, None)]


def test_iter_edges_reads_crlf_and_blank_lines(tmp_path):
    search = _search(
        tmp_path,
        [_foldseek_row("a", "b"), "", _foldseek_row("c", "d")],
        newline="\r\n",
    )

    edges = list(iter_edges(search, HitFilter()))

    assert [(e.source_id, e.target_id) for e in edges] == [
        ("a", "b"),
        ("c", "d"),
    ]


def test_iter_edges_reports_wrong_column_count_with_line(tmp_path):
    search = _search(tmp_path, [_foldseek_row("a", "b"), "a\tb\t0.5"])

    with pytest.raises(ValueError, match=r"hits\.tsv:2: expected 13 columns"):
        list(iter_edges(search, HitFilter()))


def test_iter_edges_reports_bad_number_with_line(tmp_path):
    search = _search(tmp_path, [_foldseek_row("a", "b", qcov="n/a")])

    with pytest.raises(ValueError, match=r"hits\.tsv:1:"):
        list(iter_edges(search, HitFilter()))


def test_iter_edges_rejects_tm_normalisation_on_mmseqs(tmp_path):
    search = _search(tmp_path, [], engine="MMseqs2")

    with pytest.raises(ValueError, match="tm_normalisation"):
        iter_edges(search, HitFilter(tm_normalisation="query"))


def test_iter_edges_rejects_interface_filter_without_residues(tmp_path):
    search = _search(tmp_path, [])

    with pytest.raises(ValueError, match="interface_residues"):
        iter_edges(search, HitFilter(min_interface_coverage=0.5))


def test_iter_edges_rejects_interface_residues_on_mmseqs(tmp_path):
    search = _search(tmp_path, [], engine="MMseqs2")

    with pytest.raises(ValueError, match="Foldseek"):
        iter_edges(search, HitFilter(), interface_residues={"a": {1}})


def test_load_similarity_search_wraps_existing_file(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("a\tb\t0.95\t100\t0.9\t0.8\n")
    parameters = {"max_seqs": 300}

    search = load_similarity_search(
        path, "MMseqs2", version="18.8cc5c", parameters=parameters
    )
    parameters["max_seqs"] = 1  # caller's dict must not be aliased

    assert search.origin == "precomputed"
    assert search.hits_path == str(path)
    assert search.columns == HIT_COLUMNS["MMseqs2"]
    assert search.version == "18.8cc5c"
    assert search.parameters == {"max_seqs": 300}


def test_load_similarity_search_accepts_empty_file(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("")

    assert (
        load_similarity_search(path, "Foldseek").columns
        == (HIT_COLUMNS["Foldseek"])
    )


def test_load_similarity_search_rejects_wrong_columns(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("a\tb\t0.95\n")

    with pytest.raises(ValueError, match=r"hits\.tsv:1: expected 6 columns"):
        load_similarity_search(path, "MMseqs2")


def test_load_similarity_search_rejects_unknown_engine(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("")

    with pytest.raises(ValueError, match="engine"):
        load_similarity_search(path, "BLAST")
