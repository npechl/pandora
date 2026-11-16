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
    label_interface_gaps,
    label_interfaces,
    resolve_layer_rows,
)
from pandora.annotations.entry import CONTACT_ATOM_SETS
from pandora.cli._io import (
    STEP_SUMMARIES,
    TABLE_FORMATS,
    TableWriter,
    check_contact_layer,
    check_table_format,
    contact_rows,
    copy_source,
    entry_key,
    find_table,
    layer_origin,
    prepare_output_dir,
    read_rows,
    read_structure,
    read_structures_dir,
    rows_by_entry,
    structure_index,
    table_columns,
    write_table,
)
from pandora.datasets import missing_positions
from pandora.export import read_records, write_json
from pandora.ingestion import load_policy
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import (
    LayerDescription,
    ResidueContactRow,
)
from pandora.schemas.dataset import DatasetChainRecord

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
    entries = 0
    with TableWriter(out, "layer", args.table_format) as writer:
        for _, path in sorted(structure_index(structures_dir).items()):
            structure = read_structure(path)
            if structure is None:
                continue
            entries += 1
            for row in function(structure, **settings):
                writer.write(row.model_dump())
    rows = writer.rows
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
            counts={"entries": entries, "rows": rows},
            previous_summary=_previous_summary(structures_dir),
        ),
        out / "layer.json",
    )
    print(f"computed: {rows} rows for {entries} entries -> {out}")


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
_PREVIOUS_SUMMARIES = STEP_SUMMARIES


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
        columns = table_columns(table)
        if columns:
            check_layer_columns(columns, description.template)
    except ValueError as exc:
        raise SystemExit(f"cannot import {table}: {exc}") from exc

    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    structures_dir = Path(args.structures_dir)
    paths = structure_index(structures_dir)
    counts: Counter[str] = Counter(rows_in=0)
    with (
        TableWriter(out, "layer", args.table_format) as layer,
        TableWriter(out, "unmatched", args.table_format) as unmatched,
    ):
        for entry, rows in rows_by_entry(table):
            path = paths.get(entry)
            structure = read_structure(path) if path else None
            found, missed = resolve_layer_rows(
                rows,
                description.template,
                {structure.entry_id: structure} if structure else {},
            )
            for row in found:
                layer.write(row)
            for row in missed:
                unmatched.write(row)
                counts[f"unmatched:{row['reason']}"] += 1
            counts["rows_in"] += len(rows)
    counts["rows_matched"] = layer.rows
    counts["rows_unmatched"] = unmatched.rows
    if args.source_files:
        copy_source(Path(args.source_files), out / "source")
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
    print(f"imported: {layer.rows} rows, {unmatched.rows} unmatched -> {out}")


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
    rows_in: Counter[str] = Counter()
    rows = contact_rows(find_table(contacts, "layer"), rows_in)
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
            counts={"rows_in": rows_in["rows_in"], "rows_out": len(counted)},
            previous_summary=contacts / "layer.json",
        ),
        out / "layer.json",
    )
    print(
        f"counted: {len(counted)} {args.level} rows from {rows_in['rows_in']} "
        f"contacts -> {out}"
    )


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
        flag: layer_origin(folder)
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

    label_table = find_table(labels, "layer")
    labels_by_entry: dict[str, list[dict[str, Any]]] = {}
    for row in read_rows(label_table) if label_table else ():
        labels_by_entry.setdefault(entry_key(row), []).append(row)
    drop = set(args.drop)
    # A label the layer can give (`parameters.labels`, recorded by Pandora's
    # labelling commands) is no typo even when no interface has it.
    present = {
        row["label"] for rows in labels_by_entry.values() for row in rows
    }
    known = (labels_policy.get("parameters") or {}).get("labels") or []
    unknown = drop - present - set(known)
    if unknown:
        raise SystemExit(
            f"{command}: no interface has the label(s) {sorted(unknown)}; "
            "check for typos"
        )

    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    counts = _filter_entries(
        find_table(contacts, "layer"),
        labels_by_entry,
        drop,
        TableWriter(out, "layer", args.table_format),
    )
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


def _filter_entries(
    table: Path | None,
    labels_by_entry: dict[str, list[dict[str, Any]]],
    drop: set[str],
    writer: TableWriter,
) -> Counter[str]:
    """Write `table`'s contact rows minus dropped interfaces to `writer`,
    one entry at a time; returns the filter's counts."""

    counts: Counter[str] = Counter(rows_in=0)
    with writer:
        for entry, rows in rows_by_entry(table) if table else ():
            entry_rows = [ResidueContactRow.model_validate(r) for r in rows]
            entry_labels = labels_by_entry.pop(entry, [])
            given = {row["label"] for row in entry_labels}
            for row in drop_labelled_interfaces(
                entry_rows, entry_labels, drop & given
            ):
                writer.write(row.model_dump())
            counts.update(
                count_interface_labels(entry_rows, entry_labels, drop)
            )
            counts["rows_in"] += len(rows)
    # Labels of entries with no contacts only add to labels_unmatched.
    leftover = [row for rows in labels_by_entry.values() for row in rows]
    counts.update(count_interface_labels([], leftover, drop))
    counts["rows_out"] = writer.rows
    return counts


def _cmd_label(args: argparse.Namespace) -> None:
    """Handle `annotate label-interfaces`: label a contact layer's
    interfaces by contacts per side and dimer size, as a layer folder."""

    command = "annotate label-interfaces"
    started_at = now_iso()
    check_table_format(args.table_format)
    contacts = Path(args.contacts_dir)
    check_contact_layer(command, contacts)
    out = Path(args.output_dir)
    if out.resolve() == contacts.resolve():
        raise SystemExit(
            f"{command}: {out} is the --contacts folder; choose another "
            "--output-dir"
        )
    origin = layer_origin(contacts)
    if origin is None:
        raise SystemExit(
            f"{command}: cannot tell which structures folder {contacts} "
            "came from; the paths in its layer.json may be relative to "
            "where `annotate compute` ran, so run from that folder"
        )
    chains_table = find_table(origin, "chains")
    if chains_table is None:
        raise SystemExit(
            f"{command}: {origin} has no chains table; it is written by "
            "`pandora curate`"
        )
    lengths = {
        (chain.entry_id.upper(), chain.chain_id): chain.seqres_length
        for chain in read_records(chains_table, DatasetChainRecord)
        if chain.seqres_length is not None
    }
    rows_in: Counter[str] = Counter()
    rows = contact_rows(find_table(contacts, "layer"), rows_in)
    labelled = label_interfaces(
        rows, lengths, args.min_residues_per_chain, args.max_total_residues
    )

    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    write_table(labelled, out, "layer", args.table_format)
    source = json.loads((contacts / "layer.json").read_text())["policy"]
    description = LayerDescription(
        name="contact_labels",
        template="interface_labels",
        tool="pandora",
        tool_version=__version__,
        parameters={
            "contacts_layer": source.get("name"),
            "structures_dir": str(origin),
            "min_residues_per_chain": args.min_residues_per_chain,
            "max_total_residues": args.max_total_residues,
            "labels": ["ok", "few_contacts", "too_large"],
        },
    )
    labels = Counter(row.label for row in labelled)
    counts = {
        "rows_in": rows_in["rows_in"],
        "interfaces": len(labelled),
        "ok": labels["ok"],
        "few_contacts": labels["few_contacts"],
        "too_large": labels["too_large"],
        "unknown_length": sum(r.total_residues is None for r in labelled),
    }
    write_json(
        build_step_summary(
            "annotation_label",
            description,
            input=contacts,
            output=out,
            started_at=started_at,
            counts=counts,
            previous_summary=contacts / "layer.json",
        ),
        out / "layer.json",
    )
    print(
        f"labelled: {len(labelled)} interfaces ({labels['ok']} ok, "
        f"{labels['few_contacts']} few_contacts, {labels['too_large']} "
        f"too_large) -> {out}"
    )


def _cmd_gaps(args: argparse.Namespace) -> None:
    """Handle `annotate interface-gaps`: label a contact layer's
    interfaces with a missing residue inside them, as a layer folder."""

    command = "annotate interface-gaps"
    started_at = now_iso()
    check_table_format(args.table_format)
    contacts = Path(args.contacts_dir)
    check_contact_layer(command, contacts)
    out = Path(args.output_dir)
    if out.resolve() == contacts.resolve():
        raise SystemExit(
            f"{command}: {out} is the --contacts folder; choose another "
            "--output-dir"
        )
    origin = layer_origin(contacts)
    if origin is None:
        raise SystemExit(
            f"{command}: cannot tell which structures folder {contacts} "
            "came from; the paths in its layer.json may be relative to "
            "where `annotate compute` ran, so run from that folder"
        )
    missing: dict[tuple[str, str], set[int]] = {}
    unmeasured = 0
    for structure in read_structures_dir(origin):
        chains, diagnostics = missing_positions(structure)
        unmeasured += len(diagnostics.warnings)
        for chain_id, positions in chains.items():
            missing[(structure.entry_id.upper(), chain_id)] = positions
    rows_in: Counter[str] = Counter()
    rows = contact_rows(find_table(contacts, "layer"), rows_in)
    labelled = label_interface_gaps(rows, missing)

    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    write_table(labelled, out, "layer", args.table_format)
    source = json.loads((contacts / "layer.json").read_text())["policy"]
    description = LayerDescription(
        name="interface_gaps",
        template="interface_labels",
        tool="pandora",
        tool_version=__version__,
        parameters={
            "contacts_layer": source.get("name"),
            "structures_dir": str(origin),
            "missing_residue_definition": "incomplete_backbone",
            "labels": ["ok", "interface_gap", "unknown"],
        },
    )
    labels = Counter(row.label for row in labelled)
    counts = {
        "rows_in": rows_in["rows_in"],
        "interfaces": len(labelled),
        "ok": labels["ok"],
        "interface_gap": labels["interface_gap"],
        "unknown": labels["unknown"],
        "chains_unmeasured": unmeasured,
    }
    write_json(
        build_step_summary(
            "annotation_label",
            description,
            input=contacts,
            output=out,
            started_at=started_at,
            counts=counts,
            previous_summary=contacts / "layer.json",
        ),
        out / "layer.json",
    )
    if labels["unknown"]:
        print(
            f"{command}: warning: {labels['unknown']} interfaces have a "
            "chain whose missing residues can't be measured (no SEQRES, "
            "or label_seq_id renumbered); they are labelled unknown",
            file=sys.stderr,
        )
    print(
        f"labelled: {len(labelled)} interfaces ({labels['ok']} ok, "
        f"{labels['interface_gap']} interface_gap, {labels['unknown']} "
        f"unknown) -> {out}"
    )


def _register_labels(commands: argparse._SubParsersAction) -> None:
    """Add the subcommands that write interface_labels layers."""

    p = commands.add_parser(
        "label-interfaces",
        help="Label a contact layer's interfaces few_contacts, too_large "
        "or ok, with their counts, as an interface_labels layer.",
    )
    p.add_argument("--contacts", dest="contacts_dir", required=True)
    p.add_argument(
        "--min-residues-per-chain",
        type=int,
        help="label few_contacts when either chain has fewer residues "
        "in contact with the other",
    )
    p.add_argument(
        "--max-total-residues",
        type=int,
        help="label too_large when the two chains' sequence lengths add "
        "up to more",
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--table-format", choices=TABLE_FORMATS, default="parquet")
    p.set_defaults(func=_cmd_label)

    p = commands.add_parser(
        "interface-gaps",
        help="Label a contact layer's interfaces interface_gap when a "
        "residue missing backbone atoms lies inside either side, as an "
        "interface_labels layer.",
    )
    p.add_argument("--contacts", dest="contacts_dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--table-format", choices=TABLE_FORMATS, default="parquet")
    p.set_defaults(func=_cmd_gaps)


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

    _register_labels(commands)

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
