from __future__ import annotations

from collections.abc import Callable, Iterator
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pandora.schemas.similarity import Edge, HitFilter, SimilaritySearch
from pandora.similarity.layer import SIMILARITY_FLOATS, SIMILARITY_INTS

if TYPE_CHECKING:
    import pyarrow as pa

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

# Tool TSV column -> similarity template column.
_TEMPLATE_NAMES = {
    "query": "item_id_1",
    "target": "item_id_2",
    "fident": "identity",
    "qcov": "coverage_1",
    "tcov": "coverage_2",
    "qstart": "start_1",
    "qend": "end_1",
    "tstart": "start_2",
    "tend": "end_2",
}
# Tool TSV column that becomes the template's `score`.
_SCORE = {"MMseqs2": "fident", "Foldseek": "alntmscore"}

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
            `search.columns` lacks a column the engine's hit files have.
            Raised immediately, before any row is read. While iterating:
            a row has the wrong column count or a non-numeric value
            (the message names `<path>:<line>`).
        FileNotFoundError: The hit file doesn't exist (while
            iterating).
    """

    missing = [c for c in HIT_COLUMNS[search.engine] if c not in search.columns]
    if missing:
        raise ValueError(
            f"{search.engine} search is missing columns {missing}; "
            f"expected {HIT_COLUMNS[search.engine]}"
        )
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
    """Whether one row's values meet every threshold set in hit_filter
    (a NaN value fails)."""

    # `not x >= limit` rather than `x < limit`, so NaN fails.
    if hit_filter.min_score is not None and not score >= hit_filter.min_score:
        return False
    if hit_filter.min_identity is not None and not (
        identity >= hit_filter.min_identity
    ):
        return False
    if hit_filter.min_coverage is not None and not (
        coverage >= hit_filter.min_coverage
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


def read_hit_tsv(
    path: str | Path, engine: str, *, batch_size: int = 1 << 20
) -> Iterator[pa.RecordBatch]:
    """Stream a tool's hit TSV as `similarity` template batches.

    The file is the headerless TSV MMseqs2 or Foldseek writes with
    Pandora's `--format-output` columns (`HIT_COLUMNS[engine]`).
    Columns are renamed to the template (`query` -> `item_id_1`,
    `fident` -> `identity`, `qcov` -> `coverage_1`, `qstart` ->
    `start_1`, ...) and `score` is added: `fident` for MMseqs2,
    `alntmscore` for Foldseek. Other columns keep their names.

    Args:
        path: The hit TSV.
        engine: "MMseqs2" or "Foldseek".
        batch_size: Rows per batch, at most.

    Returns:
        An iterator of record batches.

    Raises:
        ValueError: `engine` is unknown, or a row doesn't fit the
            engine's columns (while iterating).
        FileNotFoundError: `path` does not exist (while iterating).
    """

    if engine not in HIT_COLUMNS:
        raise ValueError(
            f"engine must be one of {sorted(HIT_COLUMNS)}, got {engine!r}"
        )
    return _read_hit_tsv(Path(path), engine, batch_size)


def _read_hit_tsv(
    path: Path, engine: str, batch_size: int
) -> Iterator[pa.RecordBatch]:
    """Generator behind `read_hit_tsv` (engine already checked)."""

    import pyarrow as pa
    from pyarrow import csv

    columns = HIT_COLUMNS[engine]
    names = [_TEMPLATE_NAMES.get(c, c) for c in columns]
    types = {
        name: pa.string()
        if name.startswith("item_id")
        else pa.int64()
        if name in SIMILARITY_INTS or name == "alnlen"
        else pa.float64()
        for name in names
    }
    if not path.exists():
        raise FileNotFoundError(path)
    if path.stat().st_size == 0:
        return
    reader = _open_tsv(
        path,
        csv.ReadOptions(column_names=names, block_size=_block(batch_size)),
        csv.ConvertOptions(column_types=types),
    )
    score = names[columns.index(_SCORE[engine])]
    for batch in _batches(reader, path):
        yield pa.RecordBatch.from_arrays(
            [*batch.columns, batch.column(score)],
            names=[*batch.schema.names, "score"],
        )


def read_similarity_table(
    path: str | Path, *, block_size: int | None = None
) -> Iterator[pa.RecordBatch]:
    """Stream a user's similarity table as record batches.

    Parquet is read batch by batch. A TSV needs a header row; its item
    columns and recognised value columns are read as text, so that
    `match_similarity_batches` decides which values are numbers.

    Args:
        path: A `.parquet` or `.tsv` file.
        block_size: TSV only: bytes read per batch (pyarrow's default
            if None).

    Returns:
        An iterator of record batches.

    Raises:
        ValueError: `path` is neither `.parquet` nor `.tsv`, or a TSV
            column changes type part-way through the file (while
            iterating; the message names the file and the column).
        FileNotFoundError: `path` does not exist.
    """

    path = Path(path)
    if path.suffix not in (".parquet", ".tsv"):
        raise ValueError(
            f"{path}: similarity tables are .parquet or .tsv (with a "
            "header row)"
        )
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq

        return pq.ParquetFile(path).iter_batches()

    import pyarrow as pa
    from pyarrow import csv

    text = ("item_id_1", "item_id_2", *SIMILARITY_FLOATS, *SIMILARITY_INTS)
    reader = _open_tsv(
        path,
        csv.ReadOptions(block_size=block_size),
        csv.ConvertOptions(column_types=dict.fromkeys(text, pa.string())),
    )
    return _batches(reader, path)


def _block(batch_size: int) -> int:
    """TSV bytes per block for about `batch_size` hit rows."""

    return max(1 << 16, min(batch_size * 64, 1 << 30))


def _open_tsv(path: Path, read: Any, convert: Any) -> Any:
    """A pyarrow streaming reader over a tab-separated file."""

    from pyarrow import csv
    from pyarrow.lib import ArrowInvalid

    try:
        return csv.open_csv(
            path,
            read_options=read,
            parse_options=csv.ParseOptions(delimiter="\t"),
            convert_options=convert,
        )
    except ArrowInvalid as exc:
        raise ValueError(f"{path}: {exc}") from exc


def _batches(reader: Any, path: Path) -> Iterator[pa.RecordBatch]:
    """`reader`'s batches; pyarrow parse errors become ValueError naming
    `path`."""

    from pyarrow.lib import ArrowInvalid

    try:
        yield from reader
    except ArrowInvalid as exc:
        raise ValueError(f"{path}: {exc}") from exc
