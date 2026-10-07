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
    path = Path(hits_path).resolve()
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
