"""File helpers shared by the CLI subcommands."""

from __future__ import annotations

import json
import shutil
import sys
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel

from pandora.export import structure_to_mmcif, write_records
from pandora.ingestion import list_mirror_files, mirror_path
from pandora.parsing import mmcif_to_structure
from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.structure import Structure

TABLE_FORMATS = ("parquet", "jsonl")

# Errors that fail one entry; the run records them and carries on.
ENTRY_ERRORS = (ValueError, OSError, RuntimeError, KeyError, IndexError)


# What a step writes besides its tables: cleared by --overwrite.
_STEP_FILES = ("failures.jsonl", "mappings.jsonl")
_STEP_TABLES = ("entries", "chain_copies", "chains", "exclusions")


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


def read_structures_dir(input_dir: Path) -> dict[str, Structure]:
    """Parse every structure file in `input_dir` (see `structure_files`),
    keyed by `Structure.entry_id`.

    Files that fail to parse are skipped with a warning rather than
    aborting the whole batch.
    """

    structures: dict[str, Structure] = {}
    for path in structure_files(input_dir):
        structure, _, status = mmcif_to_structure(str(path))
        if structure is None:
            print(f"warning: failed to parse {path}: {status}", file=sys.stderr)
            continue
        structures[structure.entry_id] = structure
    return structures


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
