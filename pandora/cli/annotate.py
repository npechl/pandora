"""CLI `annotate`: entry-level and pairwise annotation layers."""

from __future__ import annotations

import argparse
import json
import shutil
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

from pandora._util import now_iso
from pandora.annotations import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    annotate_pairwise_sequence_identity,
    annotate_structure_counts,
    check_layer_columns,
    resolve_layer_rows,
)
from pandora.cli._io import (
    TABLE_FORMATS,
    check_table_format,
    prepare_output_dir,
    read_structures_dir,
    write_table,
)
from pandora.export import read_records, write_json
from pandora.ingestion import load_policy
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import (
    AnnotationLayer,
    LayerDescription,
    LayerRow,
)

_ENTRY_ANNOTATORS = {
    "structure_counts": annotate_structure_counts,
    "ligand_contacts": annotate_ligand_contacts,
    "chain_interfaces": annotate_chain_interfaces,
}


def _cmd_compute(args: argparse.Namespace) -> None:
    """Handle `annotate compute`: compute entry-level and/or pairwise
    annotation layers."""

    output_dir = Path(args.output_dir)
    structures = read_structures_dir(Path(args.input_dir))

    layers_by_entry: dict[str, list[AnnotationLayer]] = {
        entry_id: [] for entry_id in structures
    }
    for layer_name in args.layers:
        annotator = _ENTRY_ANNOTATORS[layer_name]
        for entry_id, structure in structures.items():
            layers_by_entry[entry_id].append(annotator(structure))

    if args.pairwise:
        for left_id, right_id in combinations(sorted(structures), 2):
            layer = annotate_pairwise_sequence_identity(
                structures[left_id], structures[right_id]
            )
            layers_by_entry[left_id].append(layer)
            layers_by_entry[right_id].append(layer)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "annotations.json").write_text(
        json.dumps(
            {
                entry_id: [layer.model_dump() for layer in layers]
                for entry_id, layers in layers_by_entry.items()
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"annotated {len(structures)} entries -> {output_dir}")


# What `annotate import` writes: cleared by --overwrite.
_LAYER_OUTPUTS = [
    *(f"{t}.{fmt}" for t in ("layer", "unmatched") for fmt in TABLE_FORMATS),
    "source",
    "layer.json",
]
_PREVIOUS_SUMMARIES = ("curation.json", "canonicalisation.json")


def _cmd_import(args: argparse.Namespace) -> None:
    """Handle `annotate import`: link a tool's table to a step folder's
    structures and write the layer folder."""

    started_at = now_iso()
    check_table_format(args.table_format)
    out = Path(args.output_dir)
    table = Path(args.table)
    try:
        description = load_policy(args.description, model=LayerDescription)
        rows = [row.model_dump() for row in read_records(table, LayerRow)]
        if rows:
            columns = {column for row in rows for column in row}
            check_layer_columns(columns, description.template)
    except ValueError as exc:
        raise SystemExit(f"cannot import {table}: {exc}") from exc

    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    structures_dir = Path(args.structures_dir)
    structures = read_structures_dir(structures_dir)
    matched, unmatched = resolve_layer_rows(
        rows, description.template, structures
    )
    if args.source_files:
        _copy_source(Path(args.source_files), out / "source")
    for name, table_rows in (("layer", matched), ("unmatched", unmatched)):
        write_table(
            [LayerRow(**row) for row in _uniform_columns(table_rows)],
            out,
            name,
            args.table_format,
        )

    counts = Counter(
        rows_in=len(rows),
        rows_matched=len(matched),
        rows_unmatched=len(unmatched),
    )
    counts.update(f"unmatched:{r['reason']}" for r in unmatched)
    previous = next(
        (
            structures_dir / name
            for name in _PREVIOUS_SUMMARIES
            if (structures_dir / name).exists()
        ),
        None,
    )
    write_json(
        build_step_summary(
            "annotation_import",
            description,
            input=table,
            output=out,
            started_at=started_at,
            counts=dict(sorted(counts.items())),
            previous_summary=previous,
        ),
        out / "layer.json",
    )
    print(f"imported: {len(matched)} rows, {len(unmatched)} unmatched -> {out}")


def _uniform_columns(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Rows with every mixed-type column's values turned into strings.

    User tables are often dirty (a number in one row, "n/a" in the
    next); Parquet needs one type per column.
    """

    kinds: dict[str, set[str]] = {}
    for row in rows:
        for column, value in row.items():
            if value is None:
                continue
            number = isinstance(value, int | float) and not isinstance(
                value, bool
            )
            kind = "number" if number else type(value).__name__
            kinds.setdefault(column, set()).add(kind)
    mixed = {column for column, found in kinds.items() if len(found) > 1}
    return [
        {
            column: str(value)
            if column in mixed and value is not None
            else value
            for column, value in row.items()
        }
        for row in rows
    ]


def _copy_source(source: Path, dest: Path) -> None:
    """Copy a tool's original output file or folder into `dest`."""

    if source.is_dir():
        shutil.copytree(source, dest)
    else:
        dest.mkdir(parents=True)
        shutil.copy2(source, dest / source.name)


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    annotate = subparsers.add_parser(
        "annotate", help="Compute or import annotation layers."
    )
    commands = annotate.add_subparsers(dest="annotate_command", required=True)
    p = commands.add_parser(
        "compute", help="Entry-level and/or pairwise annotation layers."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument(
        "--layers",
        nargs="+",
        choices=list(_ENTRY_ANNOTATORS),
        default=list(_ENTRY_ANNOTATORS),
    )
    p.add_argument("--pairwise", action="store_true")
    p.add_argument("--output-dir", required=True)
    p.set_defaults(func=_cmd_compute)

    p = commands.add_parser(
        "import",
        help="Link a table from an external tool to a step folder's "
        "structures, as an annotation layer.",
    )
    p.add_argument("--structures-dir", required=True)
    p.add_argument("--table", required=True)
    p.add_argument("--description", required=True)
    p.add_argument("--source-files")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--table-format", choices=TABLE_FORMATS, default="parquet")
    p.set_defaults(func=_cmd_import)
