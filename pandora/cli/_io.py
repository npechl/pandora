"""File helpers shared by the CLI subcommands."""

from __future__ import annotations

import itertools
import json
import os
import shutil
import sys
import tempfile
from argparse import ArgumentParser, Namespace
from array import array
from collections import Counter
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any, NamedTuple, get_args

import yaml
from pydantic import BaseModel

from pandora.export import structure_to_mmcif, write_json, write_records
from pandora.ingestion import list_mirror_files, mirror_path
from pandora.parsing import mmcif_to_structure
from pandora.schemas.annotation import (
    ResidueContactRow,
    SimilarityItems,
)
from pandora.schemas.provenance import StepSummary
from pandora.schemas.similarity import HitFilter
from pandora.schemas.structure import Structure
from pandora.similarity.hits import read_similarity_table, similarity_edges

TABLE_FORMATS = ("parquet", "jsonl")

# Errors that fail one entry; the run records them and carries on.
ENTRY_ERRORS = (ValueError, OSError, RuntimeError, KeyError, IndexError)


# What a step writes besides its tables: cleared by --overwrite.
_STEP_FILES = ("failures.jsonl", "mappings.jsonl")
# Step summaries a structures folder can hold, at most one per folder.
STEP_SUMMARIES = ("dedup.json", "curation.json", "canonicalisation.json")
_STEP_TABLES = ("entries", "chain_copies", "chains", "exclusions", "contained")
# Every output folder's summary.
SUMMARY = "summary.json"
# ponytail: legacy summary names, read until Task 9 of the CLI namespace
# plan removes them.
_LEGACY_SUMMARIES = ("layer.json", *STEP_SUMMARIES, "ingestion.json")


def failure_line(path: Path, error_type: str, message: str) -> str:
    """One `failures.jsonl` line for an input that failed."""

    failure = {"path": str(path), "error_type": error_type, "message": message}
    return json.dumps(failure) + "\n"


def duplicate_line(path: Path, entry_id: str) -> str:
    """`failures.jsonl` line for an input whose entry an earlier input
    already had."""

    return failure_line(
        path, "DuplicateEntry", f"another input is also entry {entry_id}"
    )


def prepare_output_dir(
    command: str, out: Path, outputs: Sequence[str], overwrite: bool
) -> None:
    """Create `out`. Refuse if any of `outputs` (names of files or
    folders in `out`) exists unless `overwrite`; then delete only those."""

    existing = [out / name for name in outputs if (out / name).exists()]
    if existing and not overwrite:
        raise SystemExit(
            f"{command}: {out} already holds step outputs "
            f"({existing[0].name}, ...); pass --overwrite to replace them"
        )
    for path in existing:
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
    out.mkdir(parents=True, exist_ok=True)


def prepare_step_dir(
    command: str, out: Path, summary: str, overwrite: bool
) -> None:
    """Create a step's output folder. Refuse one that already holds step
    outputs unless `overwrite`; then delete only those outputs."""

    outputs = ["structures", summary, *_STEP_FILES]
    outputs += [
        f"{table}.{fmt}" for table in _STEP_TABLES for fmt in TABLE_FORMATS
    ]
    prepare_output_dir(command, out, outputs, overwrite)


def input_summary(folder: Path) -> Path | None:
    """The summary file of an input folder, for `previous_summary` (a
    record only, never followed); None if it has none."""

    for name in (SUMMARY, *_LEGACY_SUMMARIES):
        if (folder / name).is_file():
            return folder / name
    return None


def read_summary(folder: Path) -> dict[str, Any] | None:
    """The parsed summary of `folder` (`summary.json`, or a legacy step
    or layer summary name); None if it has none."""

    path = input_summary(folder)
    return None if path is None else json.loads(path.read_text())


def write_summary(out: Path, summary: StepSummary) -> None:
    """Write `out/summary.json`."""

    write_json(summary, out / SUMMARY)


def resolve_file(command: str, path: Path, name: str) -> Path:
    """`path` when it is a file, else `path/name`; exits if neither is."""

    if path.is_file():
        return path
    if (path / name).is_file():
        return path / name
    raise SystemExit(f"{command}: {path} has no {name}")


def resolve_table(command: str, path: Path, name: str) -> Path:
    """`path` when it is a `.parquet`/`.jsonl` file, else the `name`
    table in folder `path`; exits if neither is."""

    if path.is_file() and path.suffix in (".parquet", ".jsonl"):
        return path
    table = find_table(path, name) if path.is_dir() else None
    if table is None:
        raise SystemExit(f"{command}: {path} has no {name} table")
    return table


def require_tables(command: str, folder: Path) -> tuple[Path, Path]:
    """The `entries` and `chains` tables of a structures folder; exits
    if either is missing."""

    tables = [find_table(folder, name) for name in ("entries", "chains")]
    for name, table in zip(("entries", "chains"), tables, strict=True):
        if table is None:
            raise SystemExit(
                f"{command}: {folder} has no {name} table; run "
                "`pandora structures curate` first"
            )
    return tables[0], tables[1]


def add_out(p: ArgumentParser, *, tables: bool = True) -> None:
    """Add the positional `out` folder, `--overwrite` and, when
    `tables`, `--table-format`."""

    p.add_argument("out", help="output folder")
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="replace this command's earlier outputs in OUT",
    )
    if tables:
        p.add_argument(
            "--table-format",
            choices=TABLE_FORMATS,
            default="parquet",
            help="format of the tables written (default parquet)",
        )


def add_similarity_inputs(p: ArgumentParser, *, positional: bool) -> None:
    """Add the similarity inputs (positional, or repeatable
    `--similarity`) with `--similarity-items` and the hit-filter flags."""

    what = "similarity layer folder or .parquet table"
    if positional:
        p.add_argument("similarity", nargs="+", help=f"{what}(s)")
    else:
        p.add_argument(
            "--similarity", action="append", help=f"{what}; repeatable"
        )
    p.add_argument(
        "--similarity-items",
        action="append",
        choices=get_args(SimilarityItems),
        help="item level of the similarity input; needed for a bare table "
        "or a folder without summary.json; once for every input or once "
        "per input, in order (same for the flags below)",
    )
    p.add_argument("--hit-filter", action="append", help="HitFilter YAML")
    p.add_argument(
        "--min-score",
        type=float,
        action="append",
        help="overrides --hit-filter's min_score",
    )
    p.add_argument(
        "--min-coverage",
        type=float,
        action="append",
        help="overrides --hit-filter's min_coverage",
    )


class SimilarityInput(NamedTuple):
    """One similarity input: its path, its table (None: an empty layer),
    item level and hit filter."""

    path: Path
    table: Path | None
    items: SimilarityItems
    hit_filter: HitFilter


def similarity_inputs(command: str, args: Namespace) -> list[SimilarityInput]:
    """One `SimilarityInput` per similarity input of `args`; exits on a
    missing table, a missing or disagreeing item level or a missing
    threshold."""

    paths = [Path(p) for p in args.similarity or []]
    flags = per_layer(
        "--similarity-items", args.similarity_items, len(paths), command
    )
    filters = hit_filters(args, len(paths), command)
    return [
        SimilarityInput(path, table, items, hit_filter)
        for path, flag, hit_filter in zip(paths, flags, filters, strict=True)
        for table, items in [_similarity_source(command, path, flag)]
    ]


def _similarity_source(
    command: str, path: Path, flag: SimilarityItems | None
) -> tuple[Path | None, SimilarityItems]:
    """The table and item level of one similarity input."""

    if path.is_file():
        table, summary = path, None
    else:
        found = path / "layer.parquet"
        table = found if found.is_file() else None
        summary = read_summary(path) if path.is_dir() else None
        if table is None and summary is None:
            raise SystemExit(f"{command}: {path} has no layer.parquet")
    level = None
    if summary is not None:
        policy = summary.get("policy") or {}
        if policy.get("template") != "similarity":
            raise SystemExit(
                f"{command}: {path} is a {policy.get('template')} layer, "
                "not a similarity layer"
            )
        level = policy.get("items")
    if level is not None and flag is not None and flag != level:
        raise SystemExit(
            f"{command}: --similarity-items {flag} disagrees with {path}, "
            f"a {level} layer"
        )
    if level is None and flag is None:
        raise SystemExit(
            f"{command}: {path} does not record its item level; pass "
            "--similarity-items"
        )
    return table, level or flag


def per_layer(
    flag: str, values: list[Any] | None, n: int, command: str = "cluster"
) -> list[Any]:
    """`values` of a repeatable flag, one per --similarity layer: absent,
    given once for all layers, or once per layer in order."""

    if not values:
        return [None] * n
    if len(values) == 1:
        return values * n
    if len(values) != n:
        raise SystemExit(
            f"{command}: give {flag} once (for every layer) or once per "
            f"--similarity layer ({n}), not {len(values)} times"
        )
    return values


def hit_filters(
    args: Namespace, n: int, command: str = "cluster"
) -> list[HitFilter]:
    """One HitFilter per layer, from --hit-filter with --min-score /
    --min-coverage overriding it; exits if a layer has no threshold."""

    filters = []
    for path, min_score, min_coverage in zip(
        per_layer("--hit-filter", args.hit_filter, n, command),
        per_layer("--min-score", args.min_score, n, command),
        per_layer("--min-coverage", args.min_coverage, n, command),
        strict=True,
    ):
        if path is None and min_score is None and min_coverage is None:
            # Without a filter every reported hit would become an edge
            # (the old --threshold defaulted to 0.9), merging clusters
            # silently.
            raise SystemExit(
                f"{command}: set --min-score, --min-coverage or --hit-filter "
                "for every layer; with no threshold every reported hit "
                "becomes an edge"
            )
        hit_filter = (
            load_yaml_model(HitFilter, Path(path)) if path else HitFilter()
        )
        overrides = {
            field: value
            for field, value in (
                ("min_score", min_score),
                ("min_coverage", min_coverage),
            )
            if value is not None
        }
        filters.append(hit_filter.model_copy(update=overrides))
    return filters


def layer_edges(
    layer: Path,
    hit_filter: HitFilter,
    interface_residues: dict[str, set[int]] | None,
    counts: Counter[int],
    index: int,
    command: str = "cluster",
) -> Iterator[Any]:
    """The edges of one layer (a folder's `layer.parquet`, or a bare
    table), counted into `counts[index]`; exits naming the layer on a
    bad table or filter."""

    table = layer / "layer.parquet" if layer.is_dir() else layer
    if not table.exists():
        return
    try:
        for edge in similarity_edges(
            read_similarity_table(table),
            hit_filter,
            interface_residues=interface_residues,
        ):
            counts[index] += 1
            yield edge
    except ValueError as exc:
        raise SystemExit(f"{command}: {layer}: {exc}") from exc


def check_table_format(table_format: str) -> None:
    """Exit early if Parquet tables are asked for without pandas/pyarrow."""

    if table_format != "parquet":
        return
    try:
        import pandas  # noqa: F401
        import pyarrow  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            "Parquet tables need pandas and pyarrow (the `export` extra); "
            "install it or pass --table-format jsonl"
        ) from exc


def structure_files(folder: Path) -> list[Path]:
    """The mmCIF files of a step folder (`folder/structures/` in the
    mirror layout), else of `folder` itself, flat or mirror layout."""

    if (folder / "structures").is_dir():
        return list_mirror_files(folder / "structures")
    return sorted(
        {
            *folder.glob("*.cif"),
            *folder.glob("*.cif.gz"),
            *list_mirror_files(folder),
            *list_mirror_files(folder, ".cif.gz"),
        }
    )


def find_table(folder: Path, name: str) -> Path | None:
    """`folder/<name>.parquet` or `.jsonl`, whichever exists, else None."""

    for suffix in (".parquet", ".jsonl"):
        path = folder / f"{name}{suffix}"
        if path.exists():
            return path
    return None


def check_contact_layer(command: str, folder: Path) -> None:
    """Exit unless `folder` is a `residue_contacts` layer folder."""

    summary = read_summary(folder)
    if summary is None:
        raise SystemExit(f"{command}: {folder} is not a layer folder")
    template = (summary.get("policy") or {}).get("template")
    if template != "residue_contacts":
        raise SystemExit(
            f"{command}: {folder} is a {template} layer, --contacts needs "
            "residue_contacts"
        )


def read_columns(
    path: Path, columns: Sequence[str]
) -> Iterator[dict[str, Any]]:
    """Stream rows of a Parquet or JSON Lines table as dicts holding only
    `columns` (None where the file lacks one)."""

    if path.suffix == ".parquet":
        import pyarrow.parquet as pq

        parquet = pq.ParquetFile(path)
        present = [c for c in columns if c in parquet.schema_arrow.names]
        for batch in parquet.iter_batches(columns=present):
            for row in batch.to_pylist():
                yield {column: row.get(column) for column in columns}
        return
    with path.open() as handle:
        for line in handle:
            if line.strip():
                row = json.loads(line)
                yield {column: row.get(column) for column in columns}


class ParquetSink:
    """Write record batches to one Parquet file, opened on the first
    non-empty batch; no file if nothing is written."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.rows = 0
        self._writer: Any = None

    def write(self, batch: Any) -> None:
        """Append `batch` (a pyarrow RecordBatch) to the file."""

        if batch.num_rows == 0:
            return
        if self._writer is None:
            import pyarrow.parquet as pq

            self._writer = pq.ParquetWriter(self.path, batch.schema)
        self._writer.write_batch(batch)
        self.rows += batch.num_rows

    def close(self) -> None:
        """Finish the file, if one was opened."""

        if self._writer is not None:
            self._writer.close()

    def __enter__(self) -> ParquetSink:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


def copy_source(source: Path, dest: Path) -> None:
    """Copy a tool's original output file or folder into `dest`."""

    if source.is_dir():
        shutil.copytree(source, dest)
    else:
        dest.mkdir(parents=True)
        shutil.copy2(source, dest / source.name)


def write_table(
    rows: list[BaseModel], folder: Path, name: str, table_format: str
) -> None:
    """Write `rows` to `folder/<name>.<table_format>`; nothing if empty."""

    if rows:
        write_records(rows, folder / f"{name}.{table_format}")


def read_structures_dir(input_dir: Path) -> Iterator[Structure]:
    """Parse the structure files in `input_dir` (see `structure_files`)
    one at a time, so only one is held in memory.

    Files that fail to parse are skipped with a warning rather than
    aborting the whole batch.
    """

    for path in structure_files(input_dir):
        structure = read_structure(path)
        if structure is not None:
            yield structure


def read_structure(path: Path) -> Structure | None:
    """Parse one structure file; None, with a warning, if it fails."""

    structure, _, status = mmcif_to_structure(str(path))
    if structure is None:
        print(f"warning: failed to parse {path}: {status}", file=sys.stderr)
    return structure


def structure_index(folder: Path) -> dict[str, Path]:
    """Upper-cased entry ID (the file name up to its first dot) -> the
    structure file in `folder` (see `structure_files`)."""

    return {
        path.name.split(".")[0].upper(): path
        for path in structure_files(folder)
    }


def read_rows(path: Path) -> Iterator[dict[str, Any]]:
    """Stream every row of a `.parquet`, `.jsonl` or `.json` table as a
    dict; raises ValueError for another suffix."""

    suffix = path.suffix.lower()
    if suffix == ".parquet":
        import pyarrow.parquet as pq

        for batch in pq.ParquetFile(path).iter_batches():
            yield from batch.to_pylist()
    elif suffix == ".jsonl":
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    yield json.loads(line)
    elif suffix == ".json":
        yield from json.loads(path.read_text(encoding="utf-8"))
    else:
        raise ValueError(f"unsupported suffix {suffix!r}")


def table_columns(path: Path) -> set[str]:
    """Every column name of a table (see `read_rows`)."""

    if path.suffix.lower() == ".parquet":
        import pyarrow.parquet as pq

        return set(pq.ParquetFile(path).schema_arrow.names)
    return {column for row in read_rows(path) for column in row}


def contact_rows(
    path: Path | None, counts: Counter[str]
) -> Iterator[ResidueContactRow]:
    """Stream a `residue_contacts` table one row at a time (nothing if
    `path` is None), counting rows in `counts["rows_in"]`."""

    if path is None:
        return
    for row in read_rows(path):
        counts["rows_in"] += 1
        yield ResidueContactRow.model_validate(row)


def entry_key(row: dict[str, Any]) -> str:
    """The upper-cased entry ID a table row names."""

    return str(row.get("entry_id", "")).upper()


def rows_by_entry(path: Path) -> Iterator[tuple[str, list[dict[str, Any]]]]:
    """A table's rows, one entry at a time: `(entry_key(row), rows)` in
    sorted entry order, rows in file order within an entry.

    Only one entry's rows are held as Python objects. A Parquet table is
    sorted in columnar form; a JSON Lines file is indexed by byte offset
    and read back an entry at a time.
    """

    suffix = path.suffix.lower()
    if suffix == ".parquet":
        yield from _parquet_by_entry(path)
    elif suffix == ".jsonl":
        yield from _jsonl_by_entry(path)
    else:
        grouped: dict[str, list[dict[str, Any]]] = {}
        for row in read_rows(path):
            grouped.setdefault(entry_key(row), []).append(row)
        yield from sorted(grouped.items())


def _parquet_by_entry(
    path: Path,
) -> Iterator[tuple[str, list[dict[str, Any]]]]:
    """`rows_by_entry` for Parquet: a stable columnar sort by entry."""

    import pyarrow.compute as pc
    import pyarrow.parquet as pq

    table = pq.read_table(path)
    if table.num_rows == 0:
        return
    if "entry_id" not in table.column_names:
        yield "", table.to_pylist()
        return
    ids = pc.cast(table.column("entry_id"), "string")
    keys = pc.utf8_upper(pc.fill_null(ids, "NONE"))  # as entry_key(None)
    order = pc.sort_indices(keys)  # Arrow's sort is stable
    keys, table = keys.take(order), table.take(order)
    entry, rows = None, []
    for start in range(0, table.num_rows, _CHUNK_ROWS):
        chunk = table.slice(start, _CHUNK_ROWS).to_pylist()
        chunk_keys = keys.slice(start, _CHUNK_ROWS).to_pylist()
        for key, row in zip(chunk_keys, chunk, strict=True):
            if key != entry:
                if rows:
                    yield entry, rows
                entry, rows = key, []
            rows.append(row)
    yield entry, rows


def _jsonl_by_entry(
    path: Path,
) -> Iterator[tuple[str, list[dict[str, Any]]]]:
    """`rows_by_entry` for JSON Lines: byte offsets per entry, then each
    entry's lines read back with seek."""

    offsets: dict[str, array[int]] = {}
    with path.open("rb") as handle:
        position = 0
        for line in handle:
            if line.strip():
                key = entry_key(json.loads(line))
                offsets.setdefault(key, array("q")).append(position)
            position += len(line)
        for key in sorted(offsets):
            rows = []
            for offset in offsets.pop(key):
                handle.seek(offset)
                rows.append(json.loads(handle.readline()))
            yield key, rows


_CHUNK_ROWS = 50_000


class TableWriter:
    """Stream dict rows into `folder/<name>.<table_format>`; no file if
    no rows. A column holding mixed kinds (a number in one row, text in
    the next) is written as text, since Parquet needs one type per column.

    Rows are spooled to a JSON Lines file in `folder` first, so memory
    stays at one chunk of rows whatever the table's size.
    """

    def __init__(self, folder: Path, name: str, table_format: str) -> None:
        self.path = folder / f"{name}.{table_format}"
        self.rows = 0
        self._kinds: dict[str, set[str]] = {}
        self._samples: dict[str, Any] = {}
        handle, spool = tempfile.mkstemp(dir=folder, suffix=".spool")
        self._spool = Path(spool)
        self._handle = os.fdopen(handle, "w", encoding="utf-8")

    def write(self, row: dict[str, Any]) -> None:
        """Add one row."""

        for column, value in row.items():
            kinds = self._kinds.setdefault(column, set())
            if value is not None:
                kinds.add(type(value).__name__)
                self._samples.setdefault(column, value)
        self._handle.write(json.dumps(row) + "\n")
        self.rows += 1

    def close(self) -> None:
        """Write the table from the spool and delete the spool."""

        self._handle.close()
        try:
            if self.rows:
                self._finish()
        finally:
            self._spool.unlink()

    def _finish(self) -> None:
        """Copy the spool into the table, mixed columns as text."""

        mixed = {
            column
            for column, kinds in self._kinds.items()
            if len({"number" if k in _NUMBERS else k for k in kinds}) > 1
        }
        with self._spool.open(encoding="utf-8") as handle:
            rows = (
                {
                    column: str(value)
                    if column in mixed and value is not None
                    else value
                    for column, value in json.loads(line).items()
                }
                for line in handle
            )
            if self.path.suffix == ".jsonl":
                with self.path.open("w", encoding="utf-8") as out:
                    for row in rows:
                        out.write(json.dumps(row) + "\n")
                return
            import pyarrow as pa

            schema = pa.schema(
                (column, self._type(column, column in mixed))
                for column in self._kinds
            )
            with ParquetSink(self.path) as sink:
                while chunk := list(itertools.islice(rows, _CHUNK_ROWS)):
                    sink.write(pa.RecordBatch.from_pylist(chunk, schema))

    def _type(self, column: str, mixed: bool) -> Any:
        """The Arrow type of one output column."""

        import pyarrow as pa

        kinds = self._kinds[column]
        if mixed or kinds == {"str"}:
            return pa.string()
        if not kinds:
            return pa.null()
        if kinds == {"bool"}:
            return pa.bool_()
        if kinds == {"int"}:
            return pa.int64()
        if kinds <= {"int", "float"}:
            return pa.float64()
        return pa.infer_type([self._samples[column]])

    def __enter__(self) -> TableWriter:
        return self

    def __exit__(self, exc_type: object, *exc: object) -> None:
        if exc_type is None:
            self.close()
        else:  # failed part-way: leave no table behind
            self._handle.close()
            self._spool.unlink()


_NUMBERS = {"int", "float"}


def write_structures_dir(
    structures: dict[str, Structure], output_dir: Path
) -> None:
    """Write each structure to `<mid>/<entry_id>.cif` (lowercase) in
    output_dir, the PDB mirror layout."""

    for entry_id, structure in structures.items():
        structure_to_mmcif(structure, mirror_path(output_dir, entry_id.lower()))


def write_json_dict(data: dict[str, BaseModel], path: Path) -> None:
    """Write a dict of pydantic models as one pretty-printed JSON object,
    keyed by dict key."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({k: v.model_dump() for k, v in data.items()}, indent=2),
        encoding="utf-8",
    )


def write_json_list(records: Sequence[BaseModel], path: Path) -> None:
    """Write `records` as a JSON array to `path`, `[]` when there are
    none (a small sample can leave a step with nothing to write)."""

    if records:
        write_records(records, path)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("[]", encoding="utf-8")


def load_json_list(model_cls: type[BaseModel], path: Path) -> list:
    """Load and validate a JSON array as a list of model instances."""

    data = json.loads(path.read_text())
    return [model_cls.model_validate(item) for item in data]


def load_yaml_model(model_cls: type[BaseModel], path: Path) -> BaseModel:
    """Load and validate a single model instance from a YAML file."""

    with path.open() as stream:
        data = yaml.safe_load(stream)
    return model_cls.model_validate(data)


def load_pairs(path: Path) -> list[tuple[str, str]]:
    """Load a JSON array of [item_id_1, item_id_2] pairs."""

    data = json.loads(path.read_text())
    return [(item_1, item_2) for item_1, item_2 in data]


def layer_origin(folder: Path) -> Path | None:
    """The structures folder a layer was made from, or None if unknown.

    Follows `previous_summary` back through layer summaries (a filtered
    layer points at its parent layer) to the structures folder's step
    summary. Recorded paths may be relative to wherever the step ran, so
    a path that can't be read from here, or a loop, means unknown.
    """

    path = input_summary(folder)
    seen: set[Path] = set()
    while path is not None and path.resolve() not in seen:
        seen.add(path.resolve())
        summary = _read_json(path)
        if summary is None:
            return None
        previous = summary.get("previous_summary")
        if previous is None:
            if summary.get("step") != "annotation_compute":
                return None
            origin = Path(summary["input"])
            return origin.resolve() if origin.is_dir() else None
        path = Path(previous)
        if not _is_layer_summary(path):
            return path.parent.resolve() if path.exists() else None
    return None


def _read_json(path: Path) -> dict[str, Any] | None:
    """A JSON file's content, or None if it can't be read."""

    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _is_layer_summary(path: Path) -> bool:
    """Whether `path` is a layer's summary (its policy has a template)."""

    if path.name == "layer.json":
        return True
    summary = _read_json(path) if path.name == SUMMARY else None
    return "template" in ((summary or {}).get("policy") or {})
