# Similarity Hit Files Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the in-memory, list-returning similarity API with
MMseqs2/Foldseek hit files kept on disk, streamed through a `HitFilter`
policy into clustering, so a whole-PDB all-vs-all search can be
clustered in bounded memory.

**Architecture:**

- **Searches** write the tool's raw TSV (fixed `--format-output`
  columns) to a caller-chosen path and return a small `SimilaritySearch`
  record.
- **`iter_edges`** streams that file line by line, applies a `HitFilter`
  to each row, and yields `Edge` tuples.
- **`cluster_edges`** is the existing union-find over any iterable of
  pairs. `cluster_similar_items` wires the two together and records the
  filter and search in `ClusteringProvenance`.
- **Dependencies:** stdlib only.

**Tech Stack:** Python 3.11+, pydantic v2, pytest, Ruff, uv. MMseqs2 and
Foldseek are external binaries; tests use fake shell scripts.

**Spec:** `design/specs/2026-10-07-similarity-tables-design.md`

## Global Constraints

- **Dependencies:** no new ones; stdlib streaming only, no pyarrow or
  numpy in core.
- **Stored file:** the tool's raw TSV with fixed columns.
  - MMseqs2: `query,target,fident,alnlen,qcov,tcov`
  - Foldseek:
    `query,target,fident,alnlen,qcov,tcov,alntmscore,qtmscore,ttmscore,qstart,qend,tstart,tend`
- **Edges:** a pair is an edge if any single hit row passes every
  threshold. Edges are yielded row by row with `source_id < target_id`;
  duplicates are allowed and self-hits are skipped.
- **Thresholds:** `None` means not applied; every other threshold is
  inclusive (`>=`).
- **Score:**
  - MMseqs2: `fident`.
  - Foldseek: by `tm_normalisation` — `alignment` → `alntmscore`,
    `query` → `qtmscore`, `target` → `ttmscore`, `max` → the larger of
    `qtmscore` and `ttmscore`.
- **Coverage:** by `coverage_of` — `both` → min(qcov, tcov), `query` →
  qcov, `target` → tcov, `either` → max(qcov, tcov).
- **Default `max_seqs`:** MMseqs2 300, Foldseek 1000, always recorded in
  `SimilaritySearch.parameters`.
- **`SimilarityRelationship`** stays as a schema only.
  `pair_cluster_keys`, `partition_dataset`,
  `interface_residues_from_annotation`, `residue_positions` and
  `chain_item_id` are unchanged.
- **Code style** follows CLAUDE.md: start modules with
  `from __future__ import annotations`, add type hints, use Google
  docstrings on public functions and a one-line docstring on `_helpers`,
  keep lines at 80 characters or fewer, never mutate inputs, and catch
  only named exceptions.
- **Formatting:** Ruff also formats Python blocks inside Markdown. Run
  `uv run ruff format <file>` on every edited `.md` before committing.
- Tests use local fixtures only and fake binaries.
- Commit messages end with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Known broken window:** `reproduce_dataset`'s clustering branch is
  inconsistent from Task 2 until Task 4 fixes it. No existing test
  covers that branch.

## Review Focus

1. **A manifest written before this change** (an old
   `ClusteringProvenance` JSON with `threshold` and `n_relationships`)
   must still load, and `reproduce_dataset` must raise the friendly
   "did not record" `ValueError`, not a pydantic `ValidationError`.
   Tested in Task 2 (load) and Task 4 (reproduce).
2. **Hit files with Windows line endings (`\r\n`) or blank lines** must
   parse. Tested in Task 1.
3. **A non-numeric value in a numeric column** must raise `ValueError`
   naming `<path>:<line>`, not a bare `float()` error. Tested in Task 1.
4. **`interface_residues` with an MMseqs2 search** (which has no
   alignment ranges) must raise a clear `ValueError`, not a `KeyError`.
   Tested in Task 1.
5. **Reproducing a clustering that used `min_interface_coverage`** must
   raise a clear `ValueError` (interface residues aren't recorded in
   manifests), not fail deep inside `iter_edges`. Tested in Task 4.

---

## File Structure

| File | Change |
|---|---|
| `pandora/schemas/similarity.py` | Add `SimilarityEngine`, `SimilaritySearch`, `HitFilter`, `Edge`; reshape `ClusteringProvenance` |
| `pandora/similarity/hits.py` (new) | `HIT_COLUMNS`, `load_similarity_search`, `iter_edges`, row helpers, `_interface_coverage` (moved from `structure.py`) |
| `pandora/similarity/clustering.py` | `cluster_edges` (union-find), new `cluster_similar_items` |
| `pandora/similarity/sequence.py`, `structure.py` | Searches write hit files and return `SimilaritySearch` |
| `pandora/similarity/__init__.py` | Export `cluster_edges`, `iter_edges`, `load_similarity_search` |
| `pandora/cli/app.py`, `pandora/cli/README.md` | `similarity` and `cluster` subcommands |
| `pandora/provenance/reproduce.py` | Re-run or re-use the search, cluster with the recorded filter |
| `tests/test_hits.py` (new), `tests/test_clustering.py`, `tests/test_sequence_similarity.py`, `tests/test_structure_similarity.py`, `tests/test_cli.py`, `tests/test_provenance.py` | Tests |
| `docs/usage/similarity.md`, `docs/usage/cli.md`, `docs/usage/provenance.md`, `docs/reference/policies.md`, `docs/policies/leakage.yaml`, `docs/recipes/ppi-01.md`, `examples/ppi_dataset_pipeline.py`, `examples/dataset_pipeline.py` | Docs and examples |

---

### Task 1: Schemas and the hit-file reader

**Files:**

- Modify: `pandora/schemas/similarity.py`
- Create: `pandora/similarity/hits.py`
- Modify: `pandora/similarity/structure.py` (import `_interface_coverage`
  from `hits.py` and delete the local copy)
- Modify: `tests/test_structure_similarity.py` (import `_interface_coverage`
  from `pandora.similarity.hits`)
- Create: `tests/test_hits.py`

**Interfaces:**

- Produces:
  - `SimilarityEngine = Literal["MMseqs2", "Foldseek"]`
  - `SimilaritySearch(engine, version=None, hits_path: str, columns: list[str], parameters: dict[str, Any] = {}, origin: Literal["computed", "precomputed"] = "computed", searched_at: str | None = None)`
  - `HitFilter(min_score=None, min_identity=None, min_coverage=None, coverage_of="both", min_interface_coverage=None, tm_normalisation="alignment")`
  - `Edge(NamedTuple): source_id: str, target_id: str, score: float, identity: float, coverage: float, interface_coverage: float | None`
  - `HIT_COLUMNS: dict[str, list[str]]` in `pandora.similarity.hits`
  - `load_similarity_search(hits_path: str | Path, engine: str, *, version: str | None = None, parameters: dict[str, Any] | None = None) -> SimilaritySearch`
  - `iter_edges(search: SimilaritySearch, hit_filter: HitFilter, *, interface_residues: dict[str, set[int]] | None = None) -> Iterator[Edge]`
    (validates eagerly, then returns a generator)
  - `_interface_coverage(residues: set[int], start: int, end: int) -> float`,
    now in `hits.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_hits.py`:

```python
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
```

In `tests/test_structure_similarity.py`, change the import block at the
top to:

```python
from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.structure import AtomSiteRecord, EntryRecord, Structure
from pandora.similarity.hits import _interface_coverage
from pandora.similarity.structure import (
    chain_item_id,
    interface_residues_from_annotation,
    residue_positions,
)
```

Run `uv run ruff format tests/test_hits.py tests/test_structure_similarity.py`.

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_hits.py tests/test_structure_similarity.py -q`
Expected: collection errors with
`ImportError: cannot import name 'Edge'` (or
`No module named 'pandora.similarity.hits'`).

- [ ] **Step 3: Add the schemas**

In `pandora/schemas/similarity.py`, change the `typing` import to
`from typing import Any, Literal, NamedTuple`. Add these after
`SimilarityCluster`, before `ClusteringProvenance`:

```python
SimilarityEngine = Literal["MMseqs2", "Foldseek"]


class SimilaritySearch(BaseModel):
    """One all-vs-all search: where its hit file is and how it was run.

    The hit file is the tool's raw TSV with Pandora's fixed
    `--format-output` columns; it stays on disk and is streamed by
    `iter_edges()`, never loaded whole.

    Attributes:
        engine: Which tool produced the hits.
        version: The tool's version string, if determined.
        hits_path: Path to the hit TSV.
        columns: The TSV's columns, in order.
        parameters: The search settings (binary, sensitivity,
            max_seqs, extra options, ...), as passed to the search
            function, so the search can be re-run.
        origin: "computed" if Pandora ran the search, "precomputed" if
            an existing file was wrapped with `load_similarity_search()`.
        searched_at: When the search ran, as an ISO 8601 timestamp.
    """

    engine: SimilarityEngine
    version: str | None = None
    hits_path: str
    columns: list[str]
    parameters: dict[str, Any] = Field(default_factory=dict)
    origin: Literal["computed", "precomputed"] = "computed"
    searched_at: str | None = None


class HitFilter(BaseModel):
    """Policy deciding which hit rows become similarity edges.

    A pair of items is an edge if any single hit row passes every
    threshold. `None` thresholds are not applied; the others are
    inclusive (`>=`).

    Attributes:
        min_score: Minimum score: TM-score for Foldseek (see
            `tm_normalisation`), identity for MMseqs2.
        min_identity: Minimum fraction of identical aligned residues.
        min_coverage: Minimum alignment coverage (see `coverage_of`).
        coverage_of: Which coverage is tested: "both" (the smaller of
            query and target), "query", "target", or "either" (the
            larger).
        min_interface_coverage: Minimum fraction of each side's
            interface residues inside the alignment (Foldseek only;
            needs `interface_residues`).
        tm_normalisation: Foldseek TM-score used as the score:
            "alignment" (alntmscore), "query" (qtmscore), "target"
            (ttmscore), or "max" (the larger of the two).
    """

    min_score: float | None = None
    min_identity: float | None = None
    min_coverage: float | None = None
    coverage_of: Literal["both", "query", "target", "either"] = "both"
    min_interface_coverage: float | None = None
    tm_normalisation: Literal["alignment", "query", "target", "max"] = (
        "alignment"
    )


class Edge(NamedTuple):
    """One hit row that passed a `HitFilter` (source_id < target_id).

    A plain tuple, not a pydantic model, so streaming millions of edges
    stays cheap.

    Attributes:
        source_id: The lexicographically smaller item id.
        target_id: The lexicographically larger item id.
        score: The row's score, as selected by the filter.
        identity: The row's fraction of identical aligned residues.
        coverage: The row's coverage, as selected by `coverage_of`.
        interface_coverage: Interface-restricted coverage, or None if
            not computed for this pair.
    """

    source_id: str
    target_id: str
    score: float
    identity: float
    coverage: float
    interface_coverage: float | None
```

(`ClusteringProvenance` is reshaped in Task 2.)

- [ ] **Step 4: Create `pandora/similarity/hits.py`**

```python
from __future__ import annotations

from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

from pandora.schemas.similarity import Edge, HitFilter, SimilaritySearch

HIT_COLUMNS: dict[str, list[str]] = {
    "MMseqs2": ["query", "target", "fident", "alnlen", "qcov", "tcov"],
    "Foldseek": [
        "query",
        "target",
        "fident",
        "alnlen",
        "qcov",
        "tcov",
        "alntmscore",
        "qtmscore",
        "ttmscore",
        "qstart",
        "qend",
        "tstart",
        "tend",
    ],
}

_COVERAGE: dict[str, Callable[[float, float], float]] = {
    "both": min,
    "query": lambda query, _target: query,
    "target": lambda _query, target: target,
    "either": max,
}
_TM_COLUMN = {
    "alignment": "alntmscore",
    "query": "qtmscore",
    "target": "ttmscore",
}


def load_similarity_search(
    hits_path: str | Path,
    engine: str,
    *,
    version: str | None = None,
    parameters: dict[str, Any] | None = None,
) -> SimilaritySearch:
    """Wrap a hit file produced outside Pandora as a `SimilaritySearch`.

    The file must be the tool's TSV with Pandora's fixed
    `--format-output` columns (`HIT_COLUMNS[engine]`). Only the first
    non-empty row is checked here; every row is checked again while
    streaming.

    Args:
        hits_path: Path to the existing hit TSV.
        engine: "MMseqs2" or "Foldseek".
        version: The tool version that produced the file, if known.
        parameters: The settings the search was run with, if known.
            Copied, never aliased.

    Returns:
        A `SimilaritySearch` with `origin="precomputed"`.

    Raises:
        ValueError: `engine` is unknown, or the first row has the wrong
            number of columns.
        FileNotFoundError: `hits_path` does not exist.
    """

    if engine not in HIT_COLUMNS:
        raise ValueError(
            f"engine must be one of {sorted(HIT_COLUMNS)}, got {engine!r}"
        )
    path = Path(hits_path)
    columns = HIT_COLUMNS[engine]
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.strip():
                _split_row(line, path, line_number, columns)
                break
    return SimilaritySearch(
        engine=engine,
        version=version,
        hits_path=str(path),
        columns=list(columns),
        parameters=dict(parameters or {}),
        origin="precomputed",
    )


def iter_edges(
    search: SimilaritySearch,
    hit_filter: HitFilter,
    *,
    interface_residues: dict[str, set[int]] | None = None,
) -> Iterator[Edge]:
    """Stream the hit rows of `search` that pass `hit_filter`.

    Reads the hit file line by line; memory use does not grow with the
    file. Self-hits are skipped, each yielded edge has
    `source_id < target_id`, and a pair is yielded once per passing row
    (so possibly more than once).

    Args:
        search: The search whose hit file to read.
        hit_filter: Which rows become edges.
        interface_residues: Item id -> 1-indexed interface residue
            positions in that item's structure file (Foldseek only).
            When given, each row's interface coverage is computed.

    Returns:
        An iterator of `Edge`s.

    Raises:
        ValueError: The filter can't apply to this search
            (`tm_normalisation` other than "alignment", or
            `interface_residues`, on an MMseqs2 search;
            `min_interface_coverage` without `interface_residues`).
            Raised immediately, before any row is read. While iterating:
            a row has the wrong column count or a non-numeric value
            (the message names `<path>:<line>`).
    """

    if search.engine == "MMseqs2" and hit_filter.tm_normalisation != (
        "alignment"
    ):
        raise ValueError(
            "tm_normalisation must be 'alignment' for an MMseqs2 search "
            f"(it has no TM-scores), got {hit_filter.tm_normalisation!r}"
        )
    if search.engine == "MMseqs2" and interface_residues is not None:
        raise ValueError(
            "interface_residues needs Foldseek alignment ranges "
            "(qstart/qend/tstart/tend); an MMseqs2 search has none"
        )
    if (
        hit_filter.min_interface_coverage is not None
        and interface_residues is None
    ):
        raise ValueError(
            "min_interface_coverage is set but no interface_residues were "
            "given, so every edge would be dropped"
        )
    return _iter_edges(search, hit_filter, interface_residues or {})


def _iter_edges(
    search: SimilaritySearch,
    hit_filter: HitFilter,
    interface_residues: dict[str, set[int]],
) -> Iterator[Edge]:
    """Generator behind `iter_edges` (arguments already validated)."""

    path = Path(search.hits_path)
    index = {name: position for position, name in enumerate(search.columns)}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = _split_row(line, path, line_number, search.columns)
            if fields[index["query"]] == fields[index["target"]]:
                continue
            try:
                edge = _edge(
                    fields, index, search.engine, hit_filter, interface_residues
                )
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc
            if edge is not None:
                yield edge


def _split_row(
    line: str, path: Path, line_number: int, columns: list[str]
) -> list[str]:
    """Tab-split one row, raising ValueError on a wrong column count."""

    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != len(columns):
        raise ValueError(
            f"{path}:{line_number}: expected {len(columns)} columns "
            f"({','.join(columns)}), got {len(fields)}"
        )
    return fields


def _edge(
    fields: list[str],
    index: dict[str, int],
    engine: str,
    hit_filter: HitFilter,
    interface_residues: dict[str, set[int]],
) -> Edge | None:
    """The row's Edge if it passes hit_filter, else None."""

    query = fields[index["query"]]
    target = fields[index["target"]]
    identity = float(fields[index["fident"]])
    coverage = _COVERAGE[hit_filter.coverage_of](
        float(fields[index["qcov"]]), float(fields[index["tcov"]])
    )
    if engine == "MMseqs2":
        score = identity
    elif hit_filter.tm_normalisation == "max":
        score = max(
            float(fields[index["qtmscore"]]), float(fields[index["ttmscore"]])
        )
    else:
        score = float(fields[index[_TM_COLUMN[hit_filter.tm_normalisation]]])

    interface_coverage = None
    if query in interface_residues and target in interface_residues:
        interface_coverage = min(
            _interface_coverage(
                interface_residues[query],
                int(fields[index["qstart"]]),
                int(fields[index["qend"]]),
            ),
            _interface_coverage(
                interface_residues[target],
                int(fields[index["tstart"]]),
                int(fields[index["tend"]]),
            ),
        )

    if not _passes(hit_filter, score, identity, coverage, interface_coverage):
        return None
    source_id, target_id = sorted((query, target))
    return Edge(
        source_id, target_id, score, identity, coverage, interface_coverage
    )


def _passes(
    hit_filter: HitFilter,
    score: float,
    identity: float,
    coverage: float,
    interface_coverage: float | None,
) -> bool:
    """Whether one row's values meet every threshold set in hit_filter."""

    if hit_filter.min_score is not None and score < hit_filter.min_score:
        return False
    if hit_filter.min_identity is not None and (
        identity < hit_filter.min_identity
    ):
        return False
    if hit_filter.min_coverage is not None and (
        coverage < hit_filter.min_coverage
    ):
        return False
    if hit_filter.min_interface_coverage is not None and (
        interface_coverage is None
        or interface_coverage < hit_filter.min_interface_coverage
    ):
        return False
    return True


def _interface_coverage(residues: set[int], start: int, end: int) -> float:
    """Fraction of residues (1-indexed positions) that fall within the
    inclusive [start, end] alignment range. 0.0 if residues is empty."""

    if not residues:
        return 0.0
    return sum(1 for r in residues if start <= r <= end) / len(residues)
```

In `pandora/similarity/structure.py`, delete the `_interface_coverage`
function and add `from pandora.similarity.hits import _interface_coverage`
to the imports. It is still used by the old search code until Task 3.

Run `uv run ruff format pandora/similarity pandora/schemas/similarity.py`.

- [ ] **Step 5: Run the tests and confirm they pass**

Run: `uv run pytest -q`
Expected: all PASS (the new `test_hits.py` plus every existing test).

- [ ] **Step 6: Run lint and commit**

```bash
uv run ruff format --check . && uv run ruff check .
git add pandora/schemas/similarity.py pandora/similarity/hits.py pandora/similarity/structure.py tests/test_hits.py tests/test_structure_similarity.py
git commit -m "Add SimilaritySearch, HitFilter and a streaming hit-file reader

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Clustering from streamed edges, and the CLI `cluster` command

**Files:**

- Modify: `pandora/schemas/similarity.py` (`ClusteringProvenance`)
- Modify: `pandora/similarity/clustering.py`, `pandora/similarity/__init__.py`
- Modify: `pandora/cli/app.py` (`_cmd_cluster`, the `cluster` parser,
  imports)
- Test: `tests/test_clustering.py` (rewrite), `tests/test_cli.py`
  (`test_cluster_with_pairs` plus new tests)

**Interfaces:**

- Consumes: Task 1's `SimilaritySearch`, `HitFilter`, `iter_edges` and
  `HIT_COLUMNS`.
- Produces:
  - `ClusteringProvenance(clustered_at: str, hit_filter: HitFilter | None = None, search: SimilaritySearch | None = None, n_edges: int = 0, n_edges_unknown_ids: int = 0, n_clusters: int)`
  - `cluster_edges(item_ids: list[str], edges: Iterable[Sequence[Any]]) -> tuple[list[SimilarityCluster], ClusteringProvenance]`
  - `cluster_similar_items(item_ids: list[str], search: SimilaritySearch, hit_filter: HitFilter, *, interface_residues: dict[str, set[int]] | None = None) -> tuple[list[SimilarityCluster], ClusteringProvenance]`
  - CLI `cluster`: `--search FILE` (required), `--hit-filter YAML`,
    `--min-score`, `--min-coverage`, `--interface-residues`, `--pairs`,
    `--output`.
- Note: `reproduce_dataset` still calls the old signature until Task 4.
  No test covers that path.

- [ ] **Step 1: Write the failing tests**

Replace `tests/test_clustering.py` entirely:

```python
from __future__ import annotations

from pandora.schemas.similarity import (
    ClusteringProvenance,
    HitFilter,
    SimilaritySearch,
)
from pandora.similarity.clustering import (
    cluster_edges,
    cluster_similar_items,
    pair_cluster_keys,
)
from pandora.similarity.hits import HIT_COLUMNS

ITEM_IDS = ["a", "b", "c", "d", "e"]
EDGES = [("a", "b"), ("b", "c")]  # a-b-c joined transitively, no a-c edge


def test_transitive_merge_and_isolate() -> None:
    clusters, provenance = cluster_edges(ITEM_IDS, EDGES)

    assert [c.components for c in clusters] == [["a", "b", "c"], ["d"], ["e"]]
    assert provenance.n_edges == 2
    assert provenance.n_edges_unknown_ids == 0
    assert provenance.n_clusters == 3
    assert provenance.hit_filter is None
    assert provenance.search is None


def test_every_item_appears_exactly_once() -> None:
    clusters, _ = cluster_edges(ITEM_IDS, EDGES)

    seen = [item for cluster in clusters for item in cluster.components]
    assert sorted(seen) == sorted(ITEM_IDS)


def test_unknown_ids_are_ignored_and_counted() -> None:
    clusters, provenance = cluster_edges(["a", "b"], [("a", "b"), ("a", "z")])

    assert [c.components for c in clusters] == [["a", "b"]]
    assert provenance.n_edges == 2
    assert provenance.n_edges_unknown_ids == 1


def test_duplicate_edges_do_not_change_clusters() -> None:
    once, _ = cluster_edges(ITEM_IDS, EDGES)
    repeated, _ = cluster_edges(ITEM_IDS, EDGES + EDGES + [("b", "a")])

    assert repeated == once


def test_cluster_similar_items_filters_hits_and_records_search(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text(
        "a\tb\t0.95\t100\t0.9\t0.9\n"
        "b\tc\t0.95\t100\t0.9\t0.9\n"
        "d\te\t0.10\t100\t0.9\t0.9\n"
    )
    search = SimilaritySearch(
        engine="MMseqs2",
        hits_path=str(path),
        columns=list(HIT_COLUMNS["MMseqs2"]),
    )
    hit_filter = HitFilter(min_score=0.5)

    clusters, provenance = cluster_similar_items(ITEM_IDS, search, hit_filter)

    assert [c.components for c in clusters] == [["a", "b", "c"], ["d"], ["e"]]
    assert provenance.hit_filter == hit_filter
    assert provenance.search == search
    assert provenance.n_edges == 2


def test_old_clustering_provenance_still_loads() -> None:
    old = {
        "clustered_at": "2026-01-01T00:00:00+00:00",
        "threshold": 0.9,
        "n_relationships": 3,
        "n_clusters": 2,
    }

    provenance = ClusteringProvenance.model_validate(old)

    assert provenance.search is None
    assert provenance.hit_filter is None
    assert provenance.n_clusters == 2


def test_pair_cluster_keys_uses_each_clusters_min_component_id() -> None:
    clusters, _ = cluster_edges(ITEM_IDS, EDGES)

    keys = pair_cluster_keys([("a", "d"), ("c", "e")], clusters)

    # a, c share the {a,b,c} cluster (keyed "a"); d, e are singletons.
    assert keys == {("a", "d"): ("a", "d"), ("c", "e"): ("a", "e")}
```

In `tests/test_cli.py`, add these imports at the top:

```python
from pandora.schemas.similarity import SimilaritySearch
from pandora.similarity.hits import HIT_COLUMNS
```

Replace `test_cluster_with_pairs` with these three tests and a helper:

```python
def _cluster_inputs(tmp_path):
    input_dir = tmp_path / "clustered"
    input_dir.mkdir()
    for name in ("104m", "112m", "118l"):
        (input_dir / f"{name}.cif").touch()
    hits = tmp_path / "hits.tsv"
    hits.write_text("104M\t112M\t0.95\t100\t0.95\t0.95\n")
    search = tmp_path / "hits.tsv.search.json"
    search.write_text(
        SimilaritySearch(
            engine="MMseqs2",
            hits_path=str(hits),
            columns=list(HIT_COLUMNS["MMseqs2"]),
        ).model_dump_json()
    )
    return input_dir, search


def test_cluster_with_pairs(tmp_path):
    input_dir, search = _cluster_inputs(tmp_path)
    pairs = tmp_path / "pairs.json"
    pairs.write_text(json.dumps([["104M", "112M"], ["104M", "118L"]]))
    output = tmp_path / "clusters.json"

    main(
        [
            "cluster",
            "--input-dir",
            str(input_dir),
            "--search",
            str(search),
            "--min-score",
            "0.9",
            "--pairs",
            str(pairs),
            "--output",
            str(output),
        ]
    )

    cluster_pairs = json.loads(
        (output.parent / "cluster_pairs.json").read_text()
    )
    by_pair = {(p["item_id_1"], p["item_id_2"]): p for p in cluster_pairs}
    same = by_pair[("104M", "112M")]
    assert same["cluster_id_1"] == same["cluster_id_2"]
    different = by_pair[("104M", "118L")]
    assert different["cluster_id_1"] != different["cluster_id_2"]


def test_cluster_reads_hit_filter_yaml_and_records_it(tmp_path):
    input_dir, search = _cluster_inputs(tmp_path)
    hit_filter = tmp_path / "filter.yaml"
    hit_filter.write_text("min_score: 0.99\nmin_coverage: 0.5\n")
    output = tmp_path / "clusters.json"

    main(
        [
            "cluster",
            "--input-dir",
            str(input_dir),
            "--search",
            str(search),
            "--hit-filter",
            str(hit_filter),
            "--output",
            str(output),
        ]
    )

    clusters = json.loads(output.read_text())
    assert all(c["n_components"] == 1 for c in clusters)  # 0.95 < 0.99
    provenance = json.loads(
        (output.parent / "cluster_provenance.json").read_text()
    )
    assert provenance["hit_filter"]["min_score"] == 0.99
    assert provenance["hit_filter"]["min_coverage"] == 0.5
    assert provenance["search"]["engine"] == "MMseqs2"


def test_cluster_flag_overrides_hit_filter_yaml(tmp_path):
    input_dir, search = _cluster_inputs(tmp_path)
    hit_filter = tmp_path / "filter.yaml"
    hit_filter.write_text("min_score: 0.99\n")
    output = tmp_path / "clusters.json"

    main(
        [
            "cluster",
            "--input-dir",
            str(input_dir),
            "--search",
            str(search),
            "--hit-filter",
            str(hit_filter),
            "--min-score",
            "0.9",
            "--output",
            str(output),
        ]
    )

    clusters = json.loads(output.read_text())
    assert ["104M", "112M"] in [c["components"] for c in clusters]
```

Run `uv run ruff format tests/test_clustering.py tests/test_cli.py`.

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_clustering.py tests/test_cli.py -q`
Expected:
- `test_clustering.py` fails with
  `ImportError: cannot import name 'cluster_edges'`.
- The `cluster` CLI tests fail with argparse
  `error: unrecognized arguments: --search` (`SystemExit: 2`).

- [ ] **Step 3: Reshape `ClusteringProvenance`**

In `pandora/schemas/similarity.py`, replace `ClusteringProvenance` (its
docstring too) with the following. It must stay after `HitFilter` and
`SimilaritySearch`; move the class below `Edge` if needed:

```python
class ClusteringProvenance(BaseModel):
    """Record of one clustering run: which search and filter produced
    the edges, and how many edges and clusters resulted.

    Attributes:
        clustered_at: When this clustering run completed, as an ISO
            8601 timestamp.
        hit_filter: The filter that turned hits into edges, or None
            for hand-built edges (`cluster_edges`).
        search: The search whose hit file was clustered, or None for
            hand-built edges.
        n_edges: How many edges were read (duplicates included).
        n_edges_unknown_ids: How many of those named an id outside
            `item_ids` and were ignored.
        n_clusters: How many clusters resulted.
    """

    clustered_at: str
    hit_filter: HitFilter | None = None
    search: SimilaritySearch | None = None
    n_edges: int = 0
    n_edges_unknown_ids: int = 0
    n_clusters: int
```

`n_edges` defaults to 0 so manifests written before this change still
load (Review Focus 1). Pydantic ignores their extra `threshold` and
`n_relationships` keys.

- [ ] **Step 4: Rewrite the clustering functions**

In `pandora/similarity/clustering.py`, replace the imports and
`cluster_similar_items` (keep `pair_cluster_keys` unchanged) with:

```python
from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

from pandora._util import now_iso
from pandora.schemas.similarity import (
    ClusteringProvenance,
    HitFilter,
    SimilarityCluster,
    SimilaritySearch,
)
from pandora.similarity.hits import iter_edges


def cluster_edges(
    item_ids: list[str],
    edges: Iterable[Sequence[Any]],
) -> tuple[list[SimilarityCluster], ClusteringProvenance]:
    """Group item ids into connected-component clusters.

    Two items land in the same cluster iff connected through a chain
    of edges. Items with no edges form their own singleton cluster, so
    every id in `item_ids` ends up in exactly one cluster. `edges` is
    consumed once and never stored, so it can be a stream.

    Args:
        item_ids: Every item id to place into a cluster.
        edges: Pairs of item ids, as `Edge`s from `iter_edges()` or any
            tuples whose first two items are the ids. Duplicates are
            harmless.

    Returns:
        `(clusters, provenance)` — the clusters (sorted by their
        smallest member), and the edge and cluster counts.
        `provenance.hit_filter` and `.search` are None.
    """

    parent = {item_id: item_id for item_id in item_ids}

    def find(item_id: str) -> str:
        """Union-find: root id of item_id's cluster, with path compression."""

        while parent[item_id] != item_id:
            parent[item_id] = parent[parent[item_id]]
            item_id = parent[item_id]
        return item_id

    n_edges = 0
    n_unknown = 0
    for edge in edges:
        n_edges += 1
        source_id, target_id = edge[0], edge[1]
        if source_id not in parent or target_id not in parent:
            n_unknown += 1
            continue
        root_source, root_target = find(source_id), find(target_id)
        if root_source != root_target:
            parent[root_target] = root_source

    groups: dict[str, list[str]] = {}
    for item_id in item_ids:
        groups.setdefault(find(item_id), []).append(item_id)

    clusters = [
        SimilarityCluster(components=sorted(members), n_components=len(members))
        for members in sorted(groups.values(), key=lambda members: min(members))
    ]
    provenance = ClusteringProvenance(
        clustered_at=now_iso(),
        n_edges=n_edges,
        n_edges_unknown_ids=n_unknown,
        n_clusters=len(clusters),
    )
    return clusters, provenance


def cluster_similar_items(
    item_ids: list[str],
    search: SimilaritySearch,
    hit_filter: HitFilter,
    *,
    interface_residues: dict[str, set[int]] | None = None,
) -> tuple[list[SimilarityCluster], ClusteringProvenance]:
    """Cluster item ids from a search's hit file, filtered by a policy.

    Streams `search`'s hit file through `hit_filter` (`iter_edges()`)
    into `cluster_edges()`; memory is bounded by the number of items,
    not hits.

    Args:
        item_ids: Every item id to place into a cluster.
        search: The search whose hit file to read.
        hit_filter: Which hit rows become edges.
        interface_residues: Item id -> interface residue positions,
            needed when `hit_filter.min_interface_coverage` is set.

    Returns:
        `(clusters, provenance)`, with `provenance.hit_filter` and
        `.search` recorded so the clustering can be reproduced.

    Raises:
        ValueError: The filter can't apply to the search, or the hit
            file is malformed (see `iter_edges()`).
    """

    clusters, provenance = cluster_edges(
        item_ids,
        iter_edges(search, hit_filter, interface_residues=interface_residues),
    )
    return clusters, provenance.model_copy(
        update={"hit_filter": hit_filter, "search": search}
    )
```

The previous code sorted clusters by `members[0]`. Members are added in
`item_ids` order, so for unsorted `item_ids` that wasn't always the
smallest id. `min(members)` makes the order independent of `item_ids`'
order.

In `pandora/similarity/__init__.py`, add `cluster_edges` to the
clustering import, add
`from pandora.similarity.hits import iter_edges, load_similarity_search`,
and add `"cluster_edges"`, `"iter_edges"` and `"load_similarity_search"`
to `__all__` (keep it sorted).

- [ ] **Step 5: Update the CLI `cluster` command**

In `pandora/cli/app.py`:

1. **Imports.** In the `pandora.schemas.similarity` import, replace
   `SimilarityRelationship` with `HitFilter` and `SimilaritySearch`.
   `SimilarityRelationship` is no longer used in this file; Ruff F401
   will confirm.
2. **Command.** Replace `_cmd_cluster`'s body up to and including the
   `print(...)` line with:

```python
    input_dir = Path(args.input_dir)
    search = _load_json_model(SimilaritySearch, Path(args.search))
    hit_filter = (
        _load_yaml_model(HitFilter, Path(args.hit_filter))
        if args.hit_filter
        else HitFilter()
    )
    overrides = {
        field: value
        for field, value in (
            ("min_score", args.min_score),
            ("min_coverage", args.min_coverage),
        )
        if value is not None
    }
    if overrides:
        hit_filter = hit_filter.model_copy(update=overrides)
    interface_residues = (
        _load_interface_residues(Path(args.interface_residues))
        if args.interface_residues
        else None
    )
    item_ids = sorted(path.stem.upper() for path in input_dir.glob("*.cif"))

    clusters, prov = cluster_similar_items(
        item_ids, search, hit_filter, interface_residues=interface_residues
    )

    output = Path(args.output)
    write_records(clusters, output)
    write_json(prov, output.parent / "cluster_provenance.json")
    print(f"{len(clusters)} clusters from {prov.n_edges} edges -> {output}")
```

   Keep the rest of the function (the `if args.pairs:` block) as it is.
3. **Parser.** In the `cluster` subparser, replace the `--relationships`
   and `--threshold` arguments with:

```python
    p.add_argument(
        "--search",
        required=True,
        help="<hits>.search.json written by `pandora similarity`",
    )
    p.add_argument("--hit-filter", help="HitFilter YAML (thresholds)")
    p.add_argument(
        "--min-score",
        type=float,
        default=None,
        help="overrides --hit-filter's min_score",
    )
    p.add_argument(
        "--min-coverage",
        type=float,
        default=None,
        help="overrides --hit-filter's min_coverage",
    )
    p.add_argument(
        "--interface-residues",
        help=(
            "JSON file of {item_id: [residue position, ...]} for "
            "min_interface_coverage (Foldseek searches only)"
        ),
    )
```

Run `uv run ruff format pandora`.

- [ ] **Step 6: Run the tests and confirm they pass**

Run: `uv run pytest -q`
Expected: all PASS.

- [ ] **Step 7: Run lint and commit**

```bash
uv run ruff format --check . && uv run ruff check .
git add pandora/schemas/similarity.py pandora/similarity/clustering.py pandora/similarity/__init__.py pandora/cli/app.py tests/test_clustering.py tests/test_cli.py
git commit -m "Cluster from streamed, filtered hit edges

cluster_edges is the union-find over any edge stream;
cluster_similar_items reads a SimilaritySearch through a HitFilter and
records both in ClusteringProvenance. The CLI cluster command takes
--search with --hit-filter/--min-score/--min-coverage.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Searches write hit files, and the CLI `similarity` command

**Files:**

- Modify: `pandora/similarity/sequence.py`, `pandora/similarity/structure.py`
- Modify: `pandora/cli/app.py` (`_cmd_similarity`, the `similarity`
  parser, imports)
- Test: `tests/test_sequence_similarity.py` (rewrite),
  `tests/test_structure_similarity.py` (add),
  `tests/test_cli.py` (add)

**Interfaces:**

- Consumes: Task 1's `HIT_COLUMNS`, `SimilaritySearch` and
  `load_similarity_search`.
- Produces:
  - `compute_sequence_similarity(sequences, hits_path: str | Path, *, mmseqs_bin="mmseqs", sensitivity=5.7, max_seqs=300, tmp_dir=None, mmseqs_options: list[str] | None = None) -> SimilaritySearch`,
    where `parameters == {"mmseqs_bin", "sensitivity", "max_seqs", "mmseqs_options"}`.
  - `compute_structure_similarity(structures, hits_path: str | Path, *, foldseek_bin="foldseek", sensitivity=9.5, alignment_type=2, max_seqs=1000, exhaustive_search=False, tmp_dir=None, foldseek_options: list[str] | None = None) -> SimilaritySearch`,
    where `parameters == {"foldseek_bin", "sensitivity", "alignment_type", "max_seqs", "exhaustive_search", "foldseek_options"}`.
  - These `parameters` keys are exactly the functions' keyword
    arguments, so `fn(inputs, hits_path, **search.parameters)` re-runs
    the search (Task 4 relies on this).
  - CLI `similarity`: `--output` is the hit TSV, and
    `<output>.search.json` is the record. New flags: `--max-seqs`,
    `--exhaustive-search`, `--precomputed-hits`. `--input-dir` is
    required unless `--precomputed-hits` is given.

- [ ] **Step 1: Write the failing tests**

Replace `tests/test_sequence_similarity.py` entirely:

```python
import stat

from pandora.schemas.dataset import ChainRecord
from pandora.similarity.hits import HIT_COLUMNS
from pandora.similarity.sequence import compute_sequence_similarity

# Stands in for `mmseqs`: records its arguments and the query FASTA, and
# writes a canned easy-search result (a self-hit plus an a->b hit).
FAKE_MMSEQS = """#!/bin/sh
if [ "$1" = version ]; then echo 15.6f452; exit 0; fi
echo "$@" > "$(dirname "$0")/args.txt"
cp "$2" "$(dirname "$0")/query.fasta"
printf 'a_A\\ta_A\\t1.0\\t10\\t1.0\\t1.0\\n' > "$4"
printf 'a_A\\tb_A\\t0.5\\t10\\t0.9\\t0.8\\n' >> "$4"
"""


def _fake_mmseqs(directory):
    mmseqs = directory / "mmseqs"
    mmseqs.write_text(FAKE_MMSEQS)
    mmseqs.chmod(mmseqs.stat().st_mode | stat.S_IEXEC)
    return mmseqs


def test_search_keeps_hit_file_and_records_parameters(tmp_path):
    mmseqs = _fake_mmseqs(tmp_path)
    records = [
        ChainRecord(
            entry_id=entry_id,
            chain_id="A",
            entity_id="1",
            sequence=sequence,
            residue_count=3,
            atom_count=3,
        )
        for entry_id, sequence in [("a", "MKV"), ("b", "MKL"), ("c", None)]
    ]
    hits_path = tmp_path / "out" / "hits.tsv"  # parent doesn't exist yet

    search = compute_sequence_similarity(
        records, hits_path, mmseqs_bin=str(mmseqs), max_seqs=500
    )

    assert (tmp_path / "query.fasta").read_text() == ">a_A\nMKV\n>b_A\nMKL\n"
    assert hits_path.read_text().count("\n") == 2
    assert search.engine == "MMseqs2"
    assert search.version == "15.6f452"
    assert search.hits_path == str(hits_path)
    assert search.columns == HIT_COLUMNS["MMseqs2"]
    assert search.origin == "computed"
    assert search.searched_at is not None
    assert search.parameters == {
        "mmseqs_bin": str(mmseqs),
        "sensitivity": 5.7,
        "max_seqs": 500,
        "mmseqs_options": [],
    }
    args = (tmp_path / "args.txt").read_text().split()
    assert args[args.index("--max-seqs") + 1] == "500"
    assert args[args.index("--format-output") + 1] == ",".join(
        HIT_COLUMNS["MMseqs2"]
    )
```

Append to `tests/test_structure_similarity.py` (add `import stat` and
`from pandora.similarity.hits import HIT_COLUMNS` to its imports, and
add `compute_structure_similarity` to the
`pandora.similarity.structure` import):

```python
FAKE_FOLDSEEK = """#!/bin/sh
if [ "$1" = version ]; then echo 10.941cd33; exit 0; fi
echo "$@" > "$(dirname "$0")/args.txt"
printf 'A\\tB\\t0.5\\t100\\t0.9\\t0.9\\t0.6\\t0.6\\t0.6\\t1\\t100\\t1\\t100\\n' > "$4"
"""


def _fake_foldseek(directory):
    foldseek = directory / "foldseek"
    foldseek.write_text(FAKE_FOLDSEEK)
    foldseek.chmod(foldseek.stat().st_mode | stat.S_IEXEC)
    return foldseek


def test_structure_search_passes_caps_and_keeps_hits(tmp_path):
    foldseek = _fake_foldseek(tmp_path)
    struct_dir = tmp_path / "structures"
    struct_dir.mkdir()
    hits_path = tmp_path / "hits.tsv"

    search = compute_structure_similarity(
        struct_dir,
        hits_path,
        foldseek_bin=str(foldseek),
        alignment_type=1,
        max_seqs=2000,
        exhaustive_search=True,
    )

    args = (tmp_path / "args.txt").read_text().split()
    assert args[args.index("--max-seqs") + 1] == "2000"
    assert args[args.index("--exhaustive-search") + 1] == "1"
    assert args[args.index("--alignment-type") + 1] == "1"
    assert args[args.index("--format-output") + 1] == ",".join(
        HIT_COLUMNS["Foldseek"]
    )
    assert hits_path.read_text().count("\n") == 1
    assert search.engine == "Foldseek"
    assert search.version == "10.941cd33"
    assert search.columns == HIT_COLUMNS["Foldseek"]
    assert search.parameters == {
        "foldseek_bin": str(foldseek),
        "sensitivity": 9.5,
        "alignment_type": 1,
        "max_seqs": 2000,
        "exhaustive_search": True,
        "foldseek_options": [],
    }


def test_structure_search_omits_exhaustive_flag_by_default(tmp_path):
    foldseek = _fake_foldseek(tmp_path)
    struct_dir = tmp_path / "structures"
    struct_dir.mkdir()

    compute_structure_similarity(
        struct_dir, tmp_path / "hits.tsv", foldseek_bin=str(foldseek)
    )

    args = (tmp_path / "args.txt").read_text().split()
    assert "--exhaustive-search" not in args
    assert args[args.index("--max-seqs") + 1] == "1000"
```

Append to `tests/test_cli.py` (add `import shutil`, `import stat` and
`import pytest` at the top if they're missing):

```python
FAKE_MMSEQS_CLI = """#!/bin/sh
if [ "$1" = version ]; then echo 15.6f452; exit 0; fi
echo "$@" > "$(dirname "$0")/args.txt"
printf '104M\\t112M\\t0.99\\t100\\t0.95\\t0.95\\n' > "$4"
"""


def test_similarity_writes_hits_and_search_record(tmp_path):
    mmseqs = tmp_path / "mmseqs"
    mmseqs.write_text(FAKE_MMSEQS_CLI)
    mmseqs.chmod(mmseqs.stat().st_mode | stat.S_IEXEC)
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    for name in ("104m", "112m"):
        shutil.copy(FIXTURES_DIR / f"{name}.cif", input_dir / f"{name}.cif")
    output = tmp_path / "hits.tsv"

    main(
        [
            "similarity",
            "--input-dir",
            str(input_dir),
            "--engine",
            "mmseqs2",
            "--mmseqs-bin",
            str(mmseqs),
            "--max-seqs",
            "500",
            "--output",
            str(output),
        ]
    )

    record = json.loads((tmp_path / "hits.tsv.search.json").read_text())
    assert output.read_text().startswith("104M\t112M")
    assert record["hits_path"] == str(output)
    assert record["parameters"]["max_seqs"] == 500
    assert record["origin"] == "computed"


def test_similarity_wraps_precomputed_hits(tmp_path):
    hits = tmp_path / "elsewhere.tsv"
    hits.write_text("104M\t112M\t0.99\t100\t0.95\t0.95\n")
    output = tmp_path / "hits.tsv"

    main(
        [
            "similarity",
            "--engine",
            "mmseqs2",
            "--precomputed-hits",
            str(hits),
            "--output",
            str(output),
        ]
    )

    record = json.loads((tmp_path / "hits.tsv.search.json").read_text())
    assert record["origin"] == "precomputed"
    assert record["hits_path"] == str(hits)


def test_similarity_needs_input_dir_without_precomputed_hits(tmp_path):
    with pytest.raises(SystemExit, match="--input-dir"):
        main(
            [
                "similarity",
                "--engine",
                "mmseqs2",
                "--output",
                str(tmp_path / "hits.tsv"),
            ]
        )
```

Run `uv run ruff format tests`.

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_sequence_similarity.py tests/test_structure_similarity.py tests/test_cli.py -q`
Expected:
- The search tests fail with `TypeError` (an unexpected `max_seqs` or
  `hits_path` argument, or a list returned).
- The CLI tests fail with argparse
  `unrecognized arguments: --max-seqs` / `--precomputed-hits`, or with
  `the following arguments are required: --input-dir`.

- [ ] **Step 3: Rewrite `compute_sequence_similarity`**

In `pandora/similarity/sequence.py`, replace the
`pandora.schemas.similarity` import and `_OUTPUT_COLUMNS` with:

```python
from pandora._util import now_iso
from pandora.schemas.dataset import ChainRecord
from pandora.schemas.similarity import SimilaritySearch
from pandora.similarity.hits import HIT_COLUMNS
```

(Keep `_FASTA_GLOBS`.) Replace the whole `compute_sequence_similarity`
function with:

```python
def compute_sequence_similarity(
    sequences: dict[str, str] | list[ChainRecord] | str | Path,
    hits_path: str | Path,
    *,
    mmseqs_bin: str = "mmseqs",
    sensitivity: float = 5.7,
    max_seqs: int = 300,
    tmp_dir: str | Path | None = None,
    mmseqs_options: list[str] | None = None,
) -> SimilaritySearch:
    """All-vs-all sequence similarity via MMseqs2 `easy-search`.

    Writes MMseqs2's hits to `hits_path` as a TSV with the columns in
    `HIT_COLUMNS["MMseqs2"]` and keeps it there; nothing is loaded into
    memory. Read it with `iter_edges()` / `cluster_similar_items()`.

    Args:
        sequences: Mapping of item id -> sequence, a list of `ChainRecord`
            (keyed as "{entry_id}_{chain_id}", records with no sequence are
            skipped), or a path to a directory of existing FASTA files to
            run similarity over directly.
        hits_path: Where to write the hit TSV. Parent directories are
            created.
        mmseqs_bin: Path or name of the MMseqs2 binary.
        sensitivity: MMseqs2 `-s` sensitivity value.
        max_seqs: MMseqs2 `--max-seqs`: hits kept per query after the
            prefilter. The default (MMseqs2's own) can miss similar
            pairs in large families; raise it for leakage control.
        tmp_dir: Working directory for FASTA/temporary files. None =
            system temp.
        mmseqs_options: Additional options passed to MMseqs2.

    Returns:
        A `SimilaritySearch` pointing at `hits_path`. Its `parameters`
        are this function's keyword arguments, so the search can be
        re-run with `compute_sequence_similarity(seqs, path,
        **search.parameters)`.

    Raises:
        RuntimeError: `mmseqs_bin` is not on PATH.
        ValueError: `sequences` is a directory with no FASTA files.
        subprocess.CalledProcessError: MMseqs2 failed.
    """

    if shutil.which(mmseqs_bin) is None:
        raise RuntimeError(
            f"mmseqs binary {mmseqs_bin!r} not found (required for "
            "sequence similarity)"
        )

    if isinstance(sequences, list):
        sequences = {
            f"{record.entry_id}_{record.chain_id}": record.sequence
            for record in sequences
            if record.sequence is not None
        }

    options = list(mmseqs_options or [])
    hits = Path(hits_path)
    hits.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=tmp_dir) as work_dir:
        work = Path(work_dir)
        fasta_path = work / "sequences.fasta"
        if isinstance(sequences, dict):
            _write_fasta(sequences, fasta_path)
        else:
            _concat_fasta_dir(Path(sequences), fasta_path)

        subprocess.run(
            [
                mmseqs_bin,
                "easy-search",
                str(fasta_path),
                str(fasta_path),
                str(hits),
                str(work / "tmp"),
                "-s",
                str(sensitivity),
                "--max-seqs",
                str(max_seqs),
                "--format-output",
                ",".join(HIT_COLUMNS["MMseqs2"]),
                "-v",
                "1",
            ]
            + options,
            check=True,
            capture_output=True,
            text=True,
        )

    return SimilaritySearch(
        engine="MMseqs2",
        version=_mmseqs_version(mmseqs_bin),
        hits_path=str(hits),
        columns=list(HIT_COLUMNS["MMseqs2"]),
        parameters={
            "mmseqs_bin": mmseqs_bin,
            "sensitivity": sensitivity,
            "max_seqs": max_seqs,
            "mmseqs_options": options,
        },
        searched_at=now_iso(),
    )
```

- [ ] **Step 4: Rewrite `compute_structure_similarity`**

In `pandora/similarity/structure.py`:

1. **Imports.** Replace the `pandora.schemas.similarity` import and the
   `_interface_coverage` import (from Task 1), and delete
   `_OUTPUT_COLUMNS`. The import block becomes:

```python
from pandora._util import now_iso
from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.similarity import SimilaritySearch
from pandora.schemas.structure import Structure
from pandora.similarity.hits import HIT_COLUMNS
```

2. **Function.** Replace the whole `compute_structure_similarity`
   function with:

```python
def compute_structure_similarity(
    structures: dict[str, str | Path] | str | Path,
    hits_path: str | Path,
    *,
    foldseek_bin: str = "foldseek",
    sensitivity: float = 9.5,
    alignment_type: int = 2,
    max_seqs: int = 1000,
    exhaustive_search: bool = False,
    tmp_dir: str | Path | None = None,
    foldseek_options: list[str] | None = None,
) -> SimilaritySearch:
    """All-vs-all structural similarity via Foldseek `easy-search`.

    Writes Foldseek's hits to `hits_path` as a TSV with the columns in
    `HIT_COLUMNS["Foldseek"]` (including `qtmscore`/`ttmscore` and the
    alignment ranges) and keeps it there; nothing is loaded into
    memory. Interface-restricted coverage is computed when reading,
    with `iter_edges(..., interface_residues=...)`.

    Args:
        structures: Mapping of item id -> structure file path (PDB/mmCIF,
            optionally gzipped), or a path to a directory of existing
            structure files (ids are then Foldseek's own file-derived
            names).
        hits_path: Where to write the hit TSV. Parent directories are
            created.
        foldseek_bin: Path or name of the Foldseek binary.
        sensitivity: Foldseek `-s` sensitivity value.
        alignment_type: Foldseek `--alignment-type` (0: 3Di alignment,
            1: TM-align, 2: 3Di+AA — Foldseek's own default).
        max_seqs: Foldseek `--max-seqs`: hits kept per query after the
            prefilter. The default (Foldseek's own) can miss similar
            pairs in large families; raise it for leakage control.
        exhaustive_search: Pass `--exhaustive-search 1` (skip the
            prefilter; slower, finds every pair).
        tmp_dir: Working directory for structure/temporary files.
            None = system temp.
        foldseek_options: Additional options passed to Foldseek.

    Returns:
        A `SimilaritySearch` pointing at `hits_path`. Its `parameters`
        are this function's keyword arguments, so the search can be
        re-run with `compute_structure_similarity(structures, path,
        **search.parameters)`.

    Raises:
        RuntimeError: `foldseek_bin` is not on PATH.
        subprocess.CalledProcessError: Foldseek failed.
    """

    if shutil.which(foldseek_bin) is None:
        raise RuntimeError(
            f"foldseek binary {foldseek_bin!r} not found (required for "
            "structure similarity)"
        )

    options = list(foldseek_options or [])
    hits = Path(hits_path)
    hits.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=tmp_dir) as work_dir:
        work = Path(work_dir)
        if isinstance(structures, dict):
            struct_dir = work / "structures"
            struct_dir.mkdir()
            _link_structures(structures, struct_dir)
        else:
            struct_dir = Path(structures)

        command = [
            foldseek_bin,
            "easy-search",
            str(struct_dir),
            str(struct_dir),
            str(hits),
            str(work / "tmp"),
            "-s",
            str(sensitivity),
            "--alignment-type",
            str(alignment_type),
            "--max-seqs",
            str(max_seqs),
            "--format-output",
            ",".join(HIT_COLUMNS["Foldseek"]),
            "-v",
            "1",
        ]
        if exhaustive_search:
            command += ["--exhaustive-search", "1"]
        subprocess.run(
            command + options, check=True, capture_output=True, text=True
        )

    return SimilaritySearch(
        engine="Foldseek",
        version=_foldseek_version(foldseek_bin),
        hits_path=str(hits),
        columns=list(HIT_COLUMNS["Foldseek"]),
        parameters={
            "foldseek_bin": foldseek_bin,
            "sensitivity": sensitivity,
            "alignment_type": alignment_type,
            "max_seqs": max_seqs,
            "exhaustive_search": exhaustive_search,
            "foldseek_options": options,
        },
        searched_at=now_iso(),
    )
```

3. **`interface_residues_from_annotation` docstring.** Replace
   "Pass this straight to
   `compute_structure_similarity(structures=..., interface_residues=...)`"
   with "Pass this to `iter_edges()` / `cluster_similar_items()` as
   `interface_residues`, for a search run on files written by
   `export_chain_mmcif()` with the same ids."

- [ ] **Step 5: Update the CLI `similarity` command**

In `pandora/cli/app.py`:

1. **Import.** Add `load_similarity_search` to the `pandora.similarity`
   import.
2. **Command.** Replace `_cmd_similarity` entirely with:

```python
_CLI_ENGINES = {"mmseqs2": "MMseqs2", "foldseek": "Foldseek"}


def _cmd_similarity(args: argparse.Namespace) -> None:
    """Handle the `similarity` subcommand: run (or wrap a precomputed)
    all-vs-all sequence or structure search, keeping its hit file."""

    output = Path(args.output)

    if args.precomputed_hits:
        search = load_similarity_search(
            args.precomputed_hits, _CLI_ENGINES[args.engine]
        )
    elif args.input_dir is None:
        raise SystemExit(
            "similarity: --input-dir is required unless --precomputed-hits "
            "is given"
        )
    elif args.engine == "mmseqs2":
        structures = _read_structures_dir(Path(args.input_dir))
        sequences = entry_sequences(structures)
        kwargs = {"mmseqs_bin": args.mmseqs_bin}
        if args.sensitivity is not None:
            kwargs["sensitivity"] = args.sensitivity
        if args.max_seqs is not None:
            kwargs["max_seqs"] = args.max_seqs
        search = compute_sequence_similarity(sequences, output, **kwargs)
    else:
        kwargs = {
            "foldseek_bin": args.foldseek_bin,
            "alignment_type": args.alignment_type,
            "exhaustive_search": args.exhaustive_search,
        }
        if args.sensitivity is not None:
            kwargs["sensitivity"] = args.sensitivity
        if args.max_seqs is not None:
            kwargs["max_seqs"] = args.max_seqs
        search = compute_structure_similarity(
            Path(args.input_dir), output, **kwargs
        )

    record = output.with_name(output.name + ".search.json")
    write_json(search, record)
    print(f"hits -> {search.hits_path}; search record -> {record}")
```

   `_load_interface_residues` is still used by `_cmd_cluster`.
3. **Parser.** In the `similarity` subparser, replace the
   `--input-dir`, `--interface-residues` and `--output` arguments with
   the following. Keep `--engine`, `--mmseqs-bin`, `--foldseek-bin`,
   `--sensitivity` and `--alignment-type`:

```python
    p.add_argument(
        "--input-dir", help="required unless --precomputed-hits is given"
    )
    p.add_argument(
        "--max-seqs",
        type=int,
        default=None,
        help="hits kept per query (default: the tool's own, MMseqs2 300 "
        "/ Foldseek 1000; raise it for leakage control)",
    )
    p.add_argument(
        "--exhaustive-search",
        action="store_true",
        help="foldseek only: skip the prefilter",
    )
    p.add_argument(
        "--precomputed-hits",
        help="wrap an existing hit TSV (Pandora's --format-output "
        "columns) instead of running a search",
    )
    p.add_argument(
        "--output",
        required=True,
        help="hit TSV path; the search record goes to <output>.search.json",
    )
```

Run `uv run ruff format pandora`.

- [ ] **Step 6: Run the tests and confirm they pass**

Run: `uv run pytest -q`
Expected: all PASS.

- [ ] **Step 7: Run lint and commit**

```bash
uv run ruff format --check . && uv run ruff check .
git add pandora/similarity/sequence.py pandora/similarity/structure.py pandora/cli/app.py tests/test_sequence_similarity.py tests/test_structure_similarity.py tests/test_cli.py
git commit -m "Write MMseqs2/Foldseek hits to disk and return a SimilaritySearch

Searches expose and record --max-seqs (and Foldseek's
--exhaustive-search); Foldseek output gains qtmscore/ttmscore. The CLI
similarity command writes <output>.search.json and can wrap
precomputed hits.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: `reproduce_dataset`

**Files:**

- Modify: `pandora/provenance/reproduce.py`
- Test: `tests/test_provenance.py`

**Interfaces:**

- Consumes: Task 2's `cluster_similar_items(item_ids, search,
  hit_filter)` and `ClusteringProvenance.search/.hit_filter`; Task 3's
  `compute_*_similarity(inputs, hits_path, **search.parameters)`.
- Produces: `_reproduce_similarity(search: SimilaritySearch, structures:
  dict[str, Structure], output_dir: Path) -> SimilaritySearch`. A
  computed search is re-run into `output_dir / "hits.tsv"`.

- [ ] **Step 1: Write the failing tests**

In `tests/test_provenance.py`, add these imports:

```python
import stat

from pandora.schemas.similarity import (
    ClusteringProvenance,
    HitFilter,
    SimilaritySearch,
)
from pandora.similarity.hits import HIT_COLUMNS
```

Append:

```python
FAKE_MMSEQS = """#!/bin/sh
if [ "$1" = version ]; then echo 15.6f452; exit 0; fi
printf '104M\\t104M\\t1.0\\t10\\t1.0\\t1.0\\n' > "$4"
"""


def _patch_fetch(monkeypatch):
    def fake_fetch_mmcif(
        entry_id, provider, source_uri, output_dir, fetch_options=None
    ):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(MMCIF_PATH, output_dir / f"{entry_id.lower()}.cif")
        return IngestionProvenance(provider=provider, source_uri=source_uri)

    monkeypatch.setattr(
        "pandora.provenance.reproduce.fetch_mmcif", fake_fetch_mmcif
    )


def _manifest_with_clustering(clustering):
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    return build_dataset_manifest(
        dataset_id="d1",
        dataset_name="Dataset One",
        dataset_version="1.0.0",
        clustering=clustering,
        structures=[
            ProvenanceBundle(
                entry_id=structure.entry_id,
                pandora_version="0.0.0",
                generated_at="2026-01-01T00:00:00+00:00",
                ingestion=IngestionProvenance(provider="pdbe"),
            )
        ],
    )


def _clustering(search, hit_filter=None):
    return ClusteringProvenance(
        clustered_at="2026-01-01T00:00:00+00:00",
        hit_filter=hit_filter or HitFilter(min_score=0.9),
        search=search,
        n_clusters=1,
    )


def test_reproduce_reruns_computed_search(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    mmseqs = bin_dir / "mmseqs"
    mmseqs.write_text(FAKE_MMSEQS)
    mmseqs.chmod(mmseqs.stat().st_mode | stat.S_IEXEC)
    search = SimilaritySearch(
        engine="MMseqs2",
        hits_path="/somewhere/else/hits.tsv",
        columns=list(HIT_COLUMNS["MMseqs2"]),
        parameters={
            "mmseqs_bin": str(mmseqs),
            "sensitivity": 5.7,
            "max_seqs": 300,
            "mmseqs_options": [],
        },
    )
    out = tmp_path / "out"

    _, new_manifest = reproduce_dataset(
        _manifest_with_clustering(_clustering(search)), out
    )

    assert (out / "hits.tsv").exists()
    assert new_manifest.clustering.search.hits_path == str(out / "hits.tsv")
    assert new_manifest.clustering.search.parameters == search.parameters
    assert new_manifest.clustering.hit_filter == HitFilter(min_score=0.9)


def test_reproduce_reuses_existing_precomputed_hits(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    hits = tmp_path / "pre.tsv"
    hits.write_text("104M\t104M\t1.0\t10\t1.0\t1.0\n")
    search = SimilaritySearch(
        engine="MMseqs2",
        hits_path=str(hits),
        columns=list(HIT_COLUMNS["MMseqs2"]),
        origin="precomputed",
    )

    _, new_manifest = reproduce_dataset(
        _manifest_with_clustering(_clustering(search)), tmp_path / "out"
    )

    assert new_manifest.clustering.search == search


def test_reproduce_precomputed_hits_missing_raises(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    search = SimilaritySearch(
        engine="MMseqs2",
        hits_path=str(tmp_path / "gone.tsv"),
        columns=list(HIT_COLUMNS["MMseqs2"]),
        origin="precomputed",
    )

    with pytest.raises(ValueError, match="precomputed"):
        reproduce_dataset(
            _manifest_with_clustering(_clustering(search)), tmp_path / "out"
        )


def test_reproduce_old_manifest_clustering_raises(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    old = ClusteringProvenance.model_validate(
        {
            "clustered_at": "2026-01-01T00:00:00+00:00",
            "threshold": 0.9,
            "n_relationships": 1,
            "n_clusters": 1,
        }
    )

    with pytest.raises(ValueError, match="did not record"):
        reproduce_dataset(_manifest_with_clustering(old), tmp_path / "out")


def test_reproduce_interface_filter_raises(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    search = SimilaritySearch(
        engine="Foldseek",
        hits_path=str(tmp_path / "hits.tsv"),
        columns=list(HIT_COLUMNS["Foldseek"]),
    )
    clustering = _clustering(search, HitFilter(min_interface_coverage=0.5))

    with pytest.raises(ValueError, match="interface"):
        reproduce_dataset(
            _manifest_with_clustering(clustering), tmp_path / "out"
        )
```

Run `uv run ruff format tests/test_provenance.py`.

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run pytest tests/test_provenance.py -q`
Expected: the five new tests fail. Reproduce still reads
`similarity_method` and `threshold`, so expect `AttributeError` or
`ValueError` with the old message. The existing tests still PASS.

- [ ] **Step 3: Implement**

In `pandora/provenance/reproduce.py`:

1. **Imports.** Replace
   `from pandora.schemas.similarity import SimilarityRelationship` with
   `from pandora.schemas.similarity import SimilaritySearch`.
2. **`_reproduce_similarity`.** Replace the whole function with:

```python
def _reproduce_similarity(
    search: SimilaritySearch,
    structures: dict[str, Structure],
    output_dir: Path,
) -> SimilaritySearch:
    """Re-run search with its recorded parameters into output_dir, or
    re-use a precomputed hit file that still exists."""

    if search.origin == "precomputed":
        if Path(search.hits_path).exists():
            return search
        raise ValueError(
            "cannot reproduce clustering: its precomputed hit file "
            f"{search.hits_path!r} no longer exists; re-run that search "
            "and load it with load_similarity_search()"
        )
    hits_path = output_dir / "hits.tsv"
    if search.engine == "MMseqs2":
        return compute_sequence_similarity(
            entry_sequences(structures), hits_path, **search.parameters
        )
    paths: dict[str, Path] = {}
    for entry_id, structure in structures.items():
        paths[entry_id] = structure_to_mmcif(
            structure, output_dir / f"{entry_id}.reproduced.cif"
        )
    return compute_structure_similarity(paths, hits_path, **search.parameters)
```

3. **`reproduce_dataset`.** Replace the clustering block, from
   `if manifest.clustering is not None:` up to and including the
   `cluster_similar_items(...)` call, with:

```python
if manifest.clustering is not None:
    recorded = manifest.clustering
    if recorded.search is None or recorded.hit_filter is None:
        raise ValueError(
            "cannot reproduce clustering: the original manifest did "
            "not record which similarity search and hit filter "
            "produced the clusters"
        )
    if recorded.hit_filter.min_interface_coverage is not None:
        raise ValueError(
            "cannot reproduce clustering: its hit filter uses "
            "min_interface_coverage, and the interface residues it "
            "was computed with are not recorded in the manifest"
        )
    search = _reproduce_similarity(recorded.search, structures, output_dir)
    clusters, cluster_prov = cluster_similar_items(
        list(structures), search, recorded.hit_filter
    )
```

   The `if manifest.partition is not None:` block that follows is
   unchanged.

Also update `reproduce_dataset`'s docstring wherever it mentions
`similarity_method`: clustering now needs `ClusteringProvenance.search`
and `.hit_filter`.

Run `uv run ruff format pandora/provenance/reproduce.py`.

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `uv run pytest -q`
Expected: all PASS.

- [ ] **Step 5: Run lint and commit**

```bash
uv run ruff format --check . && uv run ruff check .
git add pandora/provenance/reproduce.py tests/test_provenance.py
git commit -m "Reproduce clustering from the recorded search and hit filter

Computed searches are re-run with their recorded parameters;
precomputed hit files are re-used if still present. Old manifests and
interface-coverage filters raise a clear ValueError.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Examples and docs

**Files:**

- Modify: `examples/ppi_dataset_pipeline.py`, `examples/dataset_pipeline.py`
- Modify: `docs/usage/similarity.md`, `docs/usage/cli.md`,
  `docs/usage/provenance.md`, `docs/reference/policies.md`,
  `docs/policies/leakage.yaml`, `docs/recipes/ppi-01.md`,
  `pandora/cli/README.md`
- Regenerate: `docs/assets/diagrams/*`

**Interfaces:**

- Consumes: everything from Tasks 1–4.

- [ ] **Step 1: Update `examples/ppi_dataset_pipeline.py`**

Add `from pandora.schemas.similarity import HitFilter` to the imports.
Replace:

```python
relationships = compute_sequence_similarity(sequences)
item_ids = sorted(sequences)  # same ids compute_sequence_similarity used
clusters, cluster_prov = cluster_similar_items(
    item_ids, relationships, IDENTITY_THRESHOLD
)
```

with:

```python
search = compute_sequence_similarity(sequences, OUTPUT_DIR / "hits.tsv")
item_ids = sorted(sequences)  # same ids compute_sequence_similarity used
clusters, cluster_prov = cluster_similar_items(
    item_ids, search, HitFilter(min_score=IDENTITY_THRESHOLD)
)
```

Move `OUTPUT_DIR.mkdir(parents=True, exist_ok=True)` (currently line
137) above this block. Then run it for real (MMseqs2 is installed):

Run: `uv run python examples/ppi_dataset_pipeline.py`
Expected: it completes and prints the cluster count line.

- [ ] **Step 2: Update `examples/dataset_pipeline.py`**

Change the import `from pandora.similarity import cluster_similar_items,
partition_dataset` to
`from pandora.similarity import cluster_edges, partition_dataset`.
Replace:

```python
clusters, cluster_prov = cluster_similar_items(
    ENTRY_IDS, relationships, IDENTITY_THRESHOLD
)
```

with:

```python
clusters, cluster_prov = cluster_edges(
    ENTRY_IDS,
    (
        (r.source_id, r.target_id)
        for r in relationships
        if r.score >= IDENTITY_THRESHOLD
    ),
)
```

Update the module docstring's sentence about `cluster_similar_items()`
to say these hand-built `SimilarityRelationship`s are filtered by score
and fed to `cluster_edges()`. Then:

Run: `uv run python examples/dataset_pipeline.py`
Expected: it completes, with the same cluster lines as before the
change (run `git stash` / `git stash pop` around a pre-change run to
compare, or compare with `git show main:examples/dataset_pipeline.py`
output).

- [ ] **Step 3: Rewrite `docs/usage/similarity.md`**

Keep the page's setup block and the "Leakage-safe partitioning" section
as they are. Rewrite the sections from the intro paragraph through
"Paired cluster keys" so they cover the points below. Each code block
must be run against the real MMseqs2/Foldseek on the fixtures, and its
real output pasted into the `#` comment lines. Write each prose
paragraph as one line (docs-build convention).

1. **Intro.** Searches write the tool's hit file to disk and return a
   `SimilaritySearch`. A `HitFilter` decides which hits become edges.
   `cluster_similar_items()` streams the file, so memory doesn't grow
   with the number of hits. `cluster_edges()` clusters any pairs you
   build yourself.
2. **Sequence similarity (MMseqs2).** The library tab:

   ```python
   from pandora.datasets import entry_sequences
   from pandora.similarity import compute_sequence_similarity

   sequences = entry_sequences(structures)
   search = compute_sequence_similarity(
       sequences, "datasets/output/hits/mmseqs.tsv", max_seqs=300
   )
   print(search.engine, search.hits_path, search.parameters["max_seqs"])
   ```

   The CLI tab: `pandora similarity --input-dir deduped/ --engine mmseqs2 --max-seqs 300 --output mmseqs.tsv`,
   writing `mmseqs.tsv` and `mmseqs.tsv.search.json`. Include a
   `!!! warning "max_seqs"` admonition: the tools' defaults (MMseqs2
   300, Foldseek 1000) cap hits per query, so pairs in large families
   (antibodies, kinases) can be silently missing, which becomes
   leakage. Raise `max_seqs`, or use Foldseek's `exhaustive_search`,
   for split-building.
3. **Structural similarity (Foldseek).** The same as today's section,
   with `search = compute_structure_similarity(paths,
   "datasets/output/hits/foldseek.tsv")` and a short table of the
   Foldseek columns.
4. **New "Filtering hits: `HitFilter`" section.** Show
   `iter_edges(search, HitFilter(min_score=0.5, min_coverage=0.8))`
   printing the first edges. Explain:
   - any-row-passes semantics, and that duplicates are possible (use
     `{(e.source_id, e.target_id) for e in edges}` for unique pairs);
   - the `coverage_of` and `tm_normalisation` choices, with a table;
   - that `None` means not applied and the rest are inclusive.
5. **Interface-restricted coverage (PPI pairs).** The same warning and
   the "Deriving `interface_residues` correctly" recipe, but
   `interface_residues` is now passed to
   `iter_edges(..., interface_residues=...)` /
   `cluster_similar_items(..., interface_residues=...)` together with
   `HitFilter(min_interface_coverage=0.5)`. The search no longer takes
   it. The CLI's `--interface-residues` moves to `pandora cluster`.
6. **New "Precomputed searches" section.** Show
   `load_similarity_search("hits.tsv", "Foldseek", version="10.941cd33",
   parameters={"max_seqs": 5000, "exhaustive_search": True})`, the CLI
   `--precomputed-hits`, and the exact `--format-output` string each
   tool must be run with (`",".join(HIT_COLUMNS[engine])`). Note that
   `reproduce_dataset()` re-uses the file if it still exists, and
   otherwise asks you to redo the search.
7. **Clustering.** The library tab is
   `clusters, provenance = cluster_similar_items(list(structures),
   search, HitFilter(min_score=0.9))`, printing the clusters (expect
   the same `['104M', '112M']`, `['118L', '138L']`, `['1AYI']` as
   today; paste the real output). Also show `cluster_edges()` for
   hand-built pairs. The CLI tab:
   `pandora cluster --input-dir deduped/ --search mmseqs.tsv.search.json --min-score 0.9 --output clusters.json`,
   and the `--hit-filter filter.yaml` form.
8. **Paired cluster keys.** Unchanged, except its CLI line uses
   `--search ... --min-score 0.9`.
9. **New "Migrating from the list API (≤ 0.5.6)" section**, as a short
   table:

| Before | After |
|---|---|
| `rels = compute_sequence_similarity(seqs)` | `search = compute_sequence_similarity(seqs, "hits.tsv")` |
| `compute_structure_similarity(paths, interface_residues=ir)` | `compute_structure_similarity(paths, "hits.tsv")`, then `interface_residues=ir` on `iter_edges` / `cluster_similar_items` |
| `cluster_similar_items(ids, rels, threshold=0.9)` | `cluster_similar_items(ids, search, HitFilter(min_score=0.9))` |
| hand-built `SimilarityRelationship` list + threshold | `cluster_edges(ids, ((r.source_id, r.target_id) for r in rels if r.score >= 0.9))` |
| `ClusteringProvenance.threshold` / `.similarity_method` | `.hit_filter` / `.search` |
| CLI `similarity --output relationships.json` | `--output hits.tsv` (+ `hits.tsv.search.json`) |
| CLI `cluster --relationships r.json --threshold 0.9` | `cluster --search hits.tsv.search.json --min-score 0.9` |

- [ ] **Step 4: Update the other docs**

1. **`docs/usage/cli.md`.** In the `## cluster` section, replace the
   example with
   `pandora cluster --input-dir deduped/ --search mmseqs.tsv.search.json --min-score 0.9 --output clusters.json`
   and its real printed line. In the `## similarity` section, update
   the example to `--output mmseqs.tsv` and mention the
   `.search.json` record. In `## reproduce`, replace
   "`ClusteringProvenance.similarity_method` to be set" with
   "`ClusteringProvenance.search` and `.hit_filter` to be set (a
   precomputed hit file is re-used only if it still exists; interface
   coverage filters can't be reproduced automatically)".
2. **`docs/usage/provenance.md`.** Make the same replacement in the
   "two hard requirements" paragraph.
3. **`pandora/cli/README.md`.** Replace the paragraph starting
   "`similarity --engine foldseek` takes an optional
   `--interface-residues`" with:

   > `similarity` writes the tool's hit TSV to `--output` and its `SimilaritySearch` record to `<output>.search.json`; `--max-seqs`/`--exhaustive-search` control how many hits are kept, and `--precomputed-hits <path>` wraps a TSV you ran yourself instead of searching. `cluster` reads that record with `--search`, filters hits with `--hit-filter <yaml>` and/or `--min-score`/`--min-coverage`, takes `--interface-residues <path>` (a JSON `{item_id: [residue position, ...]}` file) for interface-coverage filters, and `--pairs <path>` (a JSON array of `[item_id_1, item_id_2]` pairs) to derive a paired cluster key per pair, written to `cluster_pairs.json` next to `--output`. See `docs/usage/similarity.md`.

4. **`docs/reference/policies.md`.** Add a `## Similarity hit filter`
   section before `## Provenance`, with a field table (`min_score`,
   `min_identity`, `min_coverage`, `coverage_of`,
   `min_interface_coverage`, `tm_normalisation`), their defaults and
   meanings (copy the meanings from the `HitFilter` docstring), and
   this YAML example:

   ```yaml
   min_score: 0.5
   min_coverage: 0.8
   coverage_of: both
   tm_normalisation: max
   ```

5. **`docs/policies/leakage.yaml`.** Under `similarity_rules`, replace
   each engine's `threshold:` and `coverage_threshold:` entries (with
   their comments) by:

   ```yaml
      hit_filter: HitFilter
      # Implemented: pandora.schemas.similarity.HitFilter — see
      # docs/reference/policies.md#similarity-hit-filter.
   ```

6. **`docs/recipes/ppi-01.md`.** Mirror Step 1's example change in its
   code blocks (lines about 108–115 and 245–248). Re-run the example and
   replace the printed manifest excerpt (the `"n_relationships": 19`
   block near line 405) with the real `clustering` block now printed:
   `hit_filter`, `search`, `n_edges`.

Run `uv run ruff format docs/usage/similarity.md docs/usage/cli.md docs/usage/provenance.md docs/reference/policies.md docs/recipes/ppi-01.md examples`.

- [ ] **Step 5: Regenerate the diagrams and build the docs**

Run: `uv run --extra docs python docs/scripts/generate_erd.py`

Then `git diff --stat docs/assets/diagrams`. Commit only diagrams whose
diff includes the changed models (`SimilaritySearch`, `HitFilter` or
`ClusteringProvenance`; check with
`grep -l "HitFilter\|SimilaritySearch" docs/assets/diagrams/*.svg`), and
revert the others with `git checkout --`. On this machine they change
from tool-version layout noise alone.

Then invoke the `docs-build` skill:
`uv run zensical build --clean --strict`, the duplicate-anchor check,
and the nav check.
Expected: a clean build, no duplicate anchors, and the nav check prints
`[]`.

- [ ] **Step 6: Run the tests and commit**

```bash
uv run pytest -q
uv run ruff format --check . && uv run ruff check .
git add examples docs pandora/cli/README.md
git commit -m "Document hit files, HitFilter and the migration from the list API

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Scale verification and final checklist

**Files:**

- Create (throwaway, outside the repo):
  `/private/tmp/claude-502/-Users-pechlivanis-Desktop-work-pandora/17b5cbb9-43b9-4032-ae14-d81eb653c4c2/scratchpad/scale_hits.py`

**Interfaces:** none.

- [ ] **Step 1: Write the scale script**

```python
"""Throwaway: peak RSS of clustering a synthetic Foldseek hit file.

usage: python -I scale_hits.py <dir> <n_rows> <n_items>
"""

import random
import resource
import sys
import time
from pathlib import Path

from pandora.schemas.similarity import HitFilter, SimilaritySearch
from pandora.similarity import cluster_similar_items
from pandora.similarity.hits import HIT_COLUMNS

out_dir, n_rows, n_items = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
path = out_dir / f"synthetic_{n_rows}.tsv"
if not path.exists():
    rng = random.Random(0)
    with path.open("w") as handle:
        for _ in range(n_rows):
            q, t = rng.randrange(n_items), rng.randrange(n_items)
            tm = 0.99 if rng.random() < 0.01 else 0.3  # ~1% pass
            handle.write(
                f"I{q}\tI{t}\t0.5\t100\t0.9\t0.9\t{tm}\t{tm}\t{tm}"
                "\t1\t100\t1\t100\n"
            )
search = SimilaritySearch(
    engine="Foldseek",
    hits_path=str(path),
    columns=list(HIT_COLUMNS["Foldseek"]),
)
item_ids = [f"I{i}" for i in range(n_items)]
start = time.perf_counter()
clusters, prov = cluster_similar_items(
    item_ids, search, HitFilter(min_score=0.5)
)
print(
    f"{n_rows} rows, {n_items} items: {prov.n_edges} edges, "
    f"{len(clusters)} clusters, {time.perf_counter() - start:.1f}s, "
    f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.0f} MB"
)
```

- [ ] **Step 2: Run it at 2M and then 20M rows**

```bash
SP=/private/tmp/claude-502/-Users-pechlivanis-Desktop-work-pandora/17b5cbb9-43b9-4032-ae14-d81eb653c4c2/scratchpad
.venv/bin/python -I $SP/scale_hits.py $SP 2000000 1000000
.venv/bin/python -I $SP/scale_hits.py $SP 20000000 1000000
```

Expected: peak RSS is about the same for both runs (bounded by the 1M
items, not the rows) and **well under 1 GB**. The 20M file is about
1.2 GB on disk. Delete the synthetic files afterwards.

For comparison, report the spec's measured 2.4 KB per pair: the old API
would need about 0.5 GB for the 200k passing pairs at 20M rows, and
about 48 GB if every row had to be held, as the old code did before
filtering. Record both runs' numbers for the PR description.

- [ ] **Step 3: Run the CLAUDE.md checklist**

```bash
uv run pytest -q
uv run ruff format --check .
uv run ruff check .
uv run zensical build --clean --strict
```

Expected: all pass. Then confirm by reading the diff from `main`:

- no input is mutated (`parameters` and `options` lists are copied);
- every new public function (`iter_edges`, `load_similarity_search`,
  `cluster_edges`) has a full docstring and an `__all__` entry;
- no new dependency in `pyproject.toml`.
