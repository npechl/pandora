"""CLI `annotate`: compute Pandora's own annotation layers, or import a
tool's table as one."""

from __future__ import annotations

import argparse
import inspect
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from pandora import __version__
from pandora._util import now_iso
from pandora.annotations import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    check_layer_columns,
    count_chain_contacts,
    count_interface_labels,
    count_residue_contacts,
    drop_labelled_interfaces,
    resolve_layer_rows,
)
from pandora.annotations.entry import CONTACT_ATOM_SETS
from pandora.cli._io import (
    TABLE_FORMATS,
    check_contact_layer,
    check_table_format,
    copy_source,
    find_table,
    prepare_output_dir,
    read_structures_dir,
    write_table,
)
from pandora.export import read_records, write_json
from pandora.ingestion import load_policy
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import (
    LayerDescription,
    LayerRow,
    ResidueContactRow,
)

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
        if description.template == "similarity":
            check_layer_columns((), "similarity")
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
        copy_source(Path(args.source_files), out / "source")
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


# --level -> (count function, template of the layer it writes).
_COUNT = {
    "residue": (count_residue_contacts, "residue_values"),
    "chain": (count_chain_contacts, "chain_values"),
}


def _cmd_count(args: argparse.Namespace) -> None:
    """Handle `annotate count-contacts`: a contact layer's counts per
    residue or chain, written as a value layer folder."""

    started_at = now_iso()
    check_table_format(args.table_format)
    contacts = Path(args.contacts_dir)
    check_contact_layer("annotate count-contacts", contacts)
    source = json.loads((contacts / "layer.json").read_text())["policy"]
    out = Path(args.output_dir)
    if out.resolve() == contacts.resolve():
        raise SystemExit(
            f"annotate count-contacts: {out} is the --contacts folder; "
            "choose another --output-dir"
        )
    table = find_table(contacts, "layer")
    rows = read_records(table, ResidueContactRow) if table else []
    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    function, template = _COUNT[args.level]
    counted = function(rows)
    write_table(counted, out, "layer", args.table_format)
    description = LayerDescription(
        name="contact_count",
        template=template,
        tool="pandora",
        tool_version=__version__,
        parameters={"level": args.level, "contacts_layer": source.get("name")},
    )
    write_json(
        build_step_summary(
            "annotation_count",
            description,
            input=contacts,
            output=out,
            started_at=started_at,
            counts={"rows_in": len(rows), "rows_out": len(counted)},
            previous_summary=contacts / "layer.json",
        ),
        out / "layer.json",
    )
    print(
        f"counted: {len(counted)} {args.level} rows from {len(rows)} "
        f"contacts -> {out}"
    )


def _layer_origin(folder: Path) -> Path | None:
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


def _cmd_filter(args: argparse.Namespace) -> None:
    """Handle `annotate filter-contacts`: a contact layer without the
    interfaces a labels layer marks for dropping, as a new layer."""

    command = "annotate filter-contacts"
    started_at = now_iso()
    check_table_format(args.table_format)
    contacts = Path(args.contacts_dir)
    labels = Path(args.labels_dir)
    out = Path(args.output_dir)
    check_contact_layer(command, contacts)
    if not (labels / "layer.json").exists():
        raise SystemExit(f"{command}: {labels} is not a layer folder")
    labels_policy = (
        json.loads((labels / "layer.json").read_text()).get("policy") or {}
    )
    if labels_policy.get("template") != "interface_labels":
        raise SystemExit(
            f"{command}: {labels} is a {labels_policy.get('template')} "
            "layer, --labels needs interface_labels"
        )
    for flag, folder in (("contacts", contacts), ("labels", labels)):
        if out.resolve() == folder.resolve():
            raise SystemExit(
                f"{command}: {out} is the --{flag} folder; choose another "
                "--output-dir"
            )
    origins = {
        flag: _layer_origin(folder)
        for flag, folder in (("contacts", contacts), ("labels", labels))
    }
    if None in origins.values():
        unknown = [flag for flag, origin in origins.items() if origin is None]
        print(
            f"{command}: warning: cannot tell which structures folder the "
            f"--{' and --'.join(unknown)} layer came from; check that both "
            "layers use the same chain ids",
            file=sys.stderr,
        )
    elif origins["contacts"] != origins["labels"]:
        raise SystemExit(
            f"{command}: --contacts comes from {origins['contacts']}, "
            f"--labels from {origins['labels']}; both layers must come "
            "from the same structures folder"
        )

    table = find_table(contacts, "layer")
    rows = read_records(table, ResidueContactRow) if table else []
    table = find_table(labels, "layer")
    label_rows = read_records(table, LayerRow) if table else []
    drop = set(args.drop)
    try:
        kept = drop_labelled_interfaces(rows, label_rows, drop)
    except ValueError as exc:
        raise SystemExit(f"{command}: {exc}") from exc

    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    write_table(kept, out, "layer", args.table_format)
    source = json.loads((contacts / "layer.json").read_text())["policy"]
    description = LayerDescription(
        name="filtered_contacts",
        template="residue_contacts",
        tool="pandora",
        tool_version=__version__,
        parameters={
            "contacts_layer": source.get("name"),
            "labels_layer": labels_policy.get("name"),
            "labels_dir": str(labels),
            "drop": sorted(drop),
        },
    )
    counts = count_interface_labels(rows, label_rows, drop)
    counts |= {"rows_in": len(rows), "rows_out": len(kept)}
    write_json(
        build_step_summary(
            "annotation_filter",
            description,
            input=contacts,
            output=out,
            started_at=started_at,
            counts=dict(sorted(counts.items())),
            previous_summary=contacts / "layer.json",
        ),
        out / "layer.json",
    )
    kept_interfaces = counts["interfaces_in"] - counts["interfaces_dropped"]
    print(
        f"filtered: {counts['interfaces_dropped']} interfaces dropped, "
        f"{kept_interfaces} kept, {counts['interfaces_unlabelled']} "
        f"unlabelled -> {out}"
    )


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

    p = commands.add_parser(
        "count-contacts",
        help="Count a contact layer's partners per residue, or residues "
        "in contact per chain, as a value layer.",
    )
    p.add_argument("--contacts", dest="contacts_dir", required=True)
    p.add_argument("--level", required=True, choices=list(_COUNT))
    p.add_argument("--output-dir", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--table-format", choices=TABLE_FORMATS, default="parquet")
    p.set_defaults(func=_cmd_count)

    p = commands.add_parser(
        "filter-contacts",
        help="Drop labelled interfaces from a contact layer.",
    )
    p.add_argument("--contacts", dest="contacts_dir", required=True)
    p.add_argument("--labels", dest="labels_dir", required=True)
    p.add_argument(
        "--drop",
        action="append",
        required=True,
        help="an interface label that means drop (repeatable)",
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--table-format", choices=TABLE_FORMATS, default="parquet")
    p.set_defaults(func=_cmd_filter)
