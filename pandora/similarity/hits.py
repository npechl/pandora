from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator, Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pandora.schemas.similarity import Edge, HitFilter
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

# HitFilter.coverage_of -> the coverage columns it reads.
_COVERAGE_COLUMNS = {
    "both": ("coverage_1", "coverage_2"),
    "side_1": ("coverage_1",),
    "side_2": ("coverage_2",),
    "either": ("coverage_1", "coverage_2"),
}
_RANGES = ("start_1", "end_1", "start_2", "end_2")


def similarity_edges(
    batches: Iterable[pa.RecordBatch],
    hit_filter: HitFilter,
    *,
    interface_residues: Mapping[str, set[int]] | None = None,
) -> Iterator[Edge]:
    """Stream the rows of a similarity layer that pass `hit_filter`.

    Thresholds are applied to whole columns at once, so only passing
    rows reach Python; memory does not grow with the layer. Self-hits
    are skipped, each edge has `source_id < target_id`, and a pair is
    yielded once per passing row (so possibly more than once).

    Args:
        batches: The layer's rows, as record batches of the `similarity`
            template (e.g. from `read_similarity_table()`).
        hit_filter: Which rows become edges.
        interface_residues: Item id -> 1-indexed interface residue
            positions in that item's structure file. When given, each
            row's interface coverage is computed from its alignment
            ranges (`start_1`, `end_1`, `start_2`, `end_2`).

    Returns:
        An iterator of `Edge`s.

    Raises:
        ValueError: `min_interface_coverage` is set without
            `interface_residues` (raised immediately). On the first
            batch, before any edge: a threshold names a column the
            layer doesn't have (`score_column`, `identity`, the
            coverage columns `coverage_of` reads), or
            `interface_residues` is given and the alignment ranges are
            missing.
    """

    if (
        hit_filter.min_interface_coverage is not None
        and interface_residues is None
    ):
        raise ValueError(
            "min_interface_coverage is set but no interface_residues were "
            "given, so every edge would be dropped"
        )
    return _similarity_edges(batches, hit_filter, interface_residues)


def _similarity_edges(
    batches: Iterable[pa.RecordBatch],
    hit_filter: HitFilter,
    interface_residues: Mapping[str, set[int]] | None,
) -> Iterator[Edge]:
    """Generator behind `similarity_edges`."""

    checked = False
    for batch in batches:
        if not checked:
            _check_columns(batch.schema.names, hit_filter, interface_residues)
            checked = True
        yield from _batch_edges(batch, hit_filter, interface_residues)


def _check_columns(
    names: list[str],
    hit_filter: HitFilter,
    interface_residues: Mapping[str, set[int]] | None,
) -> None:
    """Raise ValueError if `hit_filter` needs a column `names` lacks."""

    needed = [hit_filter.score_column]
    if hit_filter.min_identity is not None:
        needed.append("identity")
    if hit_filter.min_coverage is not None:
        needed.extend(_COVERAGE_COLUMNS[hit_filter.coverage_of])
    if interface_residues is not None:
        needed.extend(_RANGES)
    missing = [name for name in needed if name not in names]
    if missing:
        raise ValueError(
            f"the similarity layer has no {', '.join(missing)} column(s), "
            "which the hit filter needs"
        )


def _batch_edges(
    batch: pa.RecordBatch,
    hit_filter: HitFilter,
    interface_residues: Mapping[str, set[int]] | None,
) -> Iterator[Edge]:
    """The edges of one batch."""

    import pyarrow as pa
    import pyarrow.compute as pc

    names = batch.schema.names

    def numbers(name: str) -> pa.Array | None:
        """Column `name` as float64, or None if absent."""

        if name not in names:
            return None
        return pc.cast(batch.column(name), pa.float64())

    score = numbers(hit_filter.score_column)
    identity = numbers("identity")
    coverage = _coverage(numbers, hit_filter.coverage_of)
    keep = pc.not_equal(batch.column("item_id_1"), batch.column("item_id_2"))
    for threshold, values in (
        (hit_filter.min_score, score),
        (hit_filter.min_identity, identity),
        (hit_filter.min_coverage, coverage),
    ):
        if threshold is not None:
            # NaN and null fail: NaN >= x is false, null becomes False.
            passes = pc.fill_null(pc.greater_equal(values, threshold), False)
            keep = pc.and_(keep, passes)
    keep = pc.fill_null(keep, False)

    def kept(values: pa.Array | None) -> list[Any]:
        """The kept rows' values, or Nones for an absent column."""

        if values is None:
            return [None] * n_kept
        return values.filter(keep).to_pylist()

    n_kept = pc.sum(pc.cast(keep, pa.int64())).as_py() or 0
    rows = zip(
        kept(batch.column("item_id_1")),
        kept(batch.column("item_id_2")),
        kept(score),
        kept(identity),
        kept(coverage),
        *(kept(batch.column(n) if n in names else None) for n in _RANGES),
        strict=True,
    )
    for query, target, value, ident, cov, *ranges in rows:
        interface_coverage = None
        if (
            interface_residues is not None
            and query in interface_residues
            and target in interface_residues
        ):
            start_1, end_1, start_2, end_2 = ranges
            interface_coverage = min(
                _interface_coverage(interface_residues[query], start_1, end_1),
                _interface_coverage(interface_residues[target], start_2, end_2),
            )
        if hit_filter.min_interface_coverage is not None and (
            interface_coverage is None
            or not interface_coverage >= hit_filter.min_interface_coverage
        ):
            continue
        source_id, target_id = sorted((query, target))
        yield Edge(source_id, target_id, value, ident, cov, interface_coverage)


def _coverage(
    numbers: Callable[[str], pa.Array | None], coverage_of: str
) -> pa.Array | None:
    """The coverage `coverage_of` selects, or None if a column is absent."""

    import pyarrow.compute as pc

    columns = [numbers(name) for name in _COVERAGE_COLUMNS[coverage_of]]
    if any(column is None for column in columns):
        return None
    # min/max_element_wise skip NaN; as null, a NaN side fails instead.
    columns = [pc.if_else(pc.is_nan(c), None, c) for c in columns]
    if coverage_of == "both":
        return pc.min_element_wise(*columns, skip_nulls=False)
    if coverage_of == "either":
        return pc.max_element_wise(*columns, skip_nulls=False)
    return columns[0]


def _interface_coverage(
    residues: set[int], start: int | None, end: int | None
) -> float:
    """Fraction of residues (1-indexed positions) that fall within the
    inclusive [start, end] alignment range. 0.0 if residues is empty or
    the range is unknown."""

    if not residues or start is None or end is None:
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
