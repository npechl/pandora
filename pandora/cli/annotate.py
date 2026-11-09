"""CLI `annotate`: compute Pandora's own annotation layers, or import a
tool's table as one."""

from __future__ import annotations

import argparse
import inspect
import shutil
from collections import Counter
from pathlib import Path
from typing import Any

from pandora import __version__
from pandora._util import now_iso
from pandora.annotations import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    check_layer_columns,
    resolve_layer_rows,
)
from pandora.annotations.entry import CONTACT_ATOM_SETS
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
from pandora.schemas.annotation import LayerDescription, LayerRow

# Layer -> (function, its settings flags: CLI dest -> function argument).
_COMPUTE = {
    "chain_interfaces": (
        annotate_chain_interfaces,
        ("distance_cutoff", "atom_set", "polymer_types"),
    ),
    "ligand_contacts": (
        annotate_ligand_contacts,
        ("distance_cutoff", "include_waters"),
    ),
}
_SETTINGS = {
    "distance_cutoff": "--distance-cutoff",
    "atom_set": "--atom-set",
    "polymer_types": "--polymer-types",
    "include_waters": "--include-waters",
}


def _cmd_compute(args: argparse.Namespace) -> None:
    """Handle `annotate compute`: one built-in contact layer for every
    structure, written as a layer folder."""

    started_at = now_iso()
    check_table_format(args.table_format)
    function, accepted = _COMPUTE[args.layer]
    for dest, flag in _SETTINGS.items():
        if dest not in accepted and getattr(args, dest) is not None:
            raise SystemExit(
                f"annotate compute: {flag} does not apply to --layer "
                f"{args.layer}"
            )
    settings = {
        dest: getattr(args, dest)
        for dest in accepted
        if getattr(args, dest) is not None
    }
    out = Path(args.output_dir)
    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    structures_dir = Path(args.structures_dir)
    structures = read_structures_dir(structures_dir)
    rows = [
        row
        for _, structure in sorted(structures.items())
        for row in function(structure, **settings)
    ]
    write_table(rows, out, "layer", args.table_format)
    defaults = {
        name: param.default
        for name, param in inspect.signature(function).parameters.items()
        if name in accepted
    }
    description = LayerDescription(
        name=args.layer,
        template="residue_contacts",
        tool="pandora",
        tool_version=__version__,
        parameters={**defaults, **settings},
    )
    write_json(
        build_step_summary(
            "annotation_compute",
            description,
            input=structures_dir,
            output=out,
            started_at=started_at,
            counts={"entries": len(structures), "rows": len(rows)},
            previous_summary=_previous_summary(structures_dir),
        ),
        out / "layer.json",
    )
    print(f"computed: {len(rows)} rows for {len(structures)} entries -> {out}")


def _previous_summary(structures_dir: Path) -> Path | None:
    """The step summary in `structures_dir`, if it has one."""

    return next(
        (
            structures_dir / name
            for name in _PREVIOUS_SUMMARIES
            if (structures_dir / name).exists()
        ),
        None,
    )


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
    previous = _previous_summary(structures_dir)
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
        "compute",
        help="Compute one of Pandora's contact layers (residue_contacts) "
        "for a step folder's structures.",
    )
    p.add_argument("--layer", required=True, choices=list(_COMPUTE))
    p.add_argument("--structures-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--distance-cutoff", type=float)
    p.add_argument("--atom-set", choices=sorted(CONTACT_ATOM_SETS))
    p.add_argument(
        "--polymer-types",
        action="append",
        help="chain_interfaces: keep chains of this entity polymer type "
        "(repeatable), e.g. polypeptide(L)",
    )
    p.add_argument(
        "--include-waters", action="store_const", const=True, default=None
    )
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--table-format", choices=TABLE_FORMATS, default="parquet")
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
