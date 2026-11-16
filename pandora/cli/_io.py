"""File helpers shared by the CLI subcommands."""

from __future__ import annotations

import itertools
import json
import os
import shutil
import sys
import tempfile
from array import array
from collections import Counter
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel

from pandora.export import structure_to_mmcif, write_records
from pandora.ingestion import list_mirror_files, mirror_path
from pandora.parsing import mmcif_to_structure
from pandora.schemas.annotation import AnnotationLayer, ResidueContactRow
from pandora.schemas.structure import Structure

TABLE_FORMATS = ("parquet", "jsonl")

# Errors that fail one entry; the run records them and carries on.
ENTRY_ERRORS = (ValueError, OSError, RuntimeError, KeyError, IndexError)


# What a step writes besides its tables: cleared by --overwrite.
_STEP_FILES = ("failures.jsonl", "mappings.jsonl")
# Step summaries a structures folder can hold, at most one per folder.
STEP_SUMMARIES = ("dedup.json", "curation.json", "canonicalisation.json")
_STEP_TABLES = ("entries", "chain_copies", "chains", "exclusions", "contained")


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
    out: Path, outputs: Sequence[str], overwrite: bool
) -> None:
    """Create `out`. Refuse if any of `outputs` (names of files or
    folders in `out`) exists unless `overwrite`; then delete only those."""

    existing = [out / name for name in outputs if (out / name).exists()]
    if existing and not overwrite:
        raise SystemExit(
            f"{out} already holds step outputs ({existing[0].name}, ...); "
            "pass --overwrite to replace them"
        )
    for path in existing:
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
    out.mkdir(parents=True, exist_ok=True)


def prepare_step_dir(out: Path, summary: str, overwrite: bool) -> None:
    """Create a step's output folder. Refuse one that already holds step
    outputs unless `overwrite`; then delete only those outputs."""

    outputs = ["structures", summary, *_STEP_FILES]
    outputs += [
        f"{table}.{fmt}" for table in _STEP_TABLES for fmt in TABLE_FORMATS
    ]
    prepare_output_dir(out, outputs, overwrite)


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
    return sorted({*folder.glob("*.cif"), *list_mirror_files(folder)})


def find_table(folder: Path, name: str) -> Path | None:
    """`folder/<name>.parquet` or `.jsonl`, whichever exists, else None."""

    for suffix in (".parquet", ".jsonl"):
        path = folder / f"{name}{suffix}"
        if path.exists():
            return path
    return None


def check_contact_layer(command: str, folder: Path) -> None:
    """Exit unless `folder` is a `residue_contacts` layer folder."""

    summary = folder / "layer.json"
    if not summary.exists():
        raise SystemExit(f"{command}: {folder} is not a layer folder")
    template = (json.loads(summary.read_text()).get("policy") or {}).get(
        "template"
    )
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


def load_json_model(model_cls: type[BaseModel], path: Path) -> BaseModel:
    """Load and validate a single model instance from a JSON file."""

    return model_cls.model_validate(json.loads(path.read_text()))


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


def load_json_dict(model_cls: type[BaseModel], path: Path) -> dict:
    """Load and validate a JSON object as a dict of model instances."""

    data = json.loads(path.read_text())
    return {k: model_cls.model_validate(v) for k, v in data.items()}


def load_annotations(path: Path) -> dict[str, list[AnnotationLayer]]:
    """Load a JSON object of entry_id -> list of AnnotationLayer."""

    data = json.loads(path.read_text())
    return {
        entry_id: [AnnotationLayer.model_validate(item) for item in layers]
        for entry_id, layers in data.items()
    }


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

    path = folder / "layer.json"
    seen: set[Path] = set()
    while path.resolve() not in seen:
        seen.add(path.resolve())
        try:
            summary = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            return None
        previous = summary.get("previous_summary")
        if previous is None:
            if summary.get("step") != "annotation_compute":
                return None
            origin = Path(summary["input"])
            return origin.resolve() if origin.is_dir() else None
        path = Path(previous)
        if path.name != "layer.json":
            return path.parent.resolve() if path.exists() else None
    return None
