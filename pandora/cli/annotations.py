"""CLI `annotations`: compute Pandora's own annotation layers, label and
filter interfaces, or import a tool's table as a layer."""

from __future__ import annotations

import argparse
import inspect
import sys
from collections import Counter
from collections.abc import Callable
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
    SUMMARY,
    TABLE_FORMATS,
    TableWriter,
    add_out,
    check_out,
    check_table_format,
    contact_rows,
    contact_table,
    copy_source,
    entry_key,
    input_error,
    input_summary,
    labels_table,
    layer_name,
    load_yaml_model,
    prepare_output_dir,
    read_rows,
    read_structure,
    read_summary,
    require_tables,
    rows_by_entry,
    structure_inputs,
    table_columns,
    write_summary,
    write_table,
)
from pandora.datasets import missing_positions
from pandora.export import read_records
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import (
    LayerDescription,
    ResidueContactRow,
)
from pandora.schemas.dataset import DatasetChainRecord

# What a layer folder holds: cleared by --overwrite.
_LAYER_OUTPUTS = [
    *(f"{t}.{fmt}" for t in ("layer", "unmatched") for fmt in TABLE_FORMATS),
    "source",
    SUMMARY,
]
_PANDORA = {"tool": "pandora", "tool_version": __version__}


def _write_layer_summary(
    out: Path,
    step: str,
    description: LayerDescription,
    source: Path,
    started_at: str,
    counts: dict[str, int],
    previous: Path | None,
) -> None:
    """Write the layer folder's `summary.json`."""

    write_summary(
        out,
        build_step_summary(
            step,
            description,
            input=source.resolve(),
            output=out,
            started_at=started_at,
            counts=dict(sorted(counts.items())),
            previous_summary=previous,
        ),
    )


def _compute(
    command: str,
    args: argparse.Namespace,
    function: Callable[..., Any],
    accepted: tuple[str, ...],
    name: str,
) -> None:
    """One built-in contact layer for every structure of `args.structures`,
    written as a layer folder."""

    started_at = now_iso()
    source, out = Path(args.structures), Path(args.out)
    check_out(command, out, source)
    check_table_format(command, args.table_format)
    paths, previous = structure_inputs(command, source)
    prepare_output_dir(command, out, _LAYER_OUTPUTS, args.overwrite)
    settings = {
        dest: getattr(args, dest)
        for dest in accepted
        if getattr(args, dest) is not None
    }
    entries = 0
    with TableWriter(out, "layer", args.table_format) as writer:
        for path in paths:
            structure = read_structure(path)
            if structure is None:
                continue
            entries += 1
            for row in function(structure, **settings):
                writer.write(row.model_dump())
    defaults = {
        key: param.default
        for key, param in inspect.signature(function).parameters.items()
        if key in accepted
    }
    description = LayerDescription(
        name=name,
        template="residue_contacts",
        parameters={**defaults, **settings},
        **_PANDORA,
    )
    _write_layer_summary(
        out,
        "annotation_compute",
        description,
        source,
        started_at,
        {"entries": entries, "rows": writer.rows},
        previous,
    )
    print(f"computed: {writer.rows} rows for {entries} entries -> {out}")


def _cmd_contacts(args: argparse.Namespace) -> None:
    """Handle `annotations contacts`: chain-interface contacts."""

    _compute(
        "annotations contacts",
        args,
        annotate_chain_interfaces,
        ("distance_cutoff", "atom_set", "polymer_types"),
        "chain_interfaces",
    )


def _cmd_ligand_contacts(args: argparse.Namespace) -> None:
    """Handle `annotations ligand-contacts`: ligand contacts."""

    _compute(
        "annotations ligand-contacts",
        args,
        annotate_ligand_contacts,
        ("distance_cutoff", "include_waters"),
        "ligand_contacts",
    )


def _cmd_import(args: argparse.Namespace) -> None:
    """Handle `annotations import`: link a tool's table to structures and
    write the layer folder."""

    command = "annotations import"
    started_at = now_iso()
    out, table = Path(args.out), Path(args.table)
    source = Path(args.structures)
    check_out(command, out, table, source)
    check_table_format(command, args.table_format)
    description = load_yaml_model(
        command, LayerDescription, Path(args.description)
    )
    try:
        if description.template == "similarity":
            check_layer_columns((), "similarity")
        columns = table_columns(table)
        if columns:
            check_layer_columns(columns, description.template)
    except (OSError, ValueError) as exc:
        raise input_error(command, table, exc) from exc

    files, previous = structure_inputs(command, source)
    prepare_output_dir(command, out, _LAYER_OUTPUTS, args.overwrite)
    paths = {path.name.split(".")[0].upper(): path for path in files}
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
    _write_layer_summary(
        out,
        "annotation_import",
        description,
        table,
        started_at,
        dict(counts),
        previous,
    )
    print(f"imported: {layer.rows} rows, {unmatched.rows} unmatched -> {out}")


# --level -> (count function, template of the layer it writes).
_COUNT = {
    "residue": (count_residue_contacts, "residue_values"),
    "chain": (count_chain_contacts, "chain_values"),
}


def _cmd_count(args: argparse.Namespace) -> None:
    """Handle `annotations contact-counts`: a contact layer's counts per
    residue or chain, written as a value layer folder."""

    command = "annotations contact-counts"
    started_at = now_iso()
    contacts, out = Path(args.contacts), Path(args.out)
    check_out(command, out, contacts)
    check_table_format(command, args.table_format)
    table = contact_table(command, contacts)
    source = layer_name(contacts)
    rows_in: Counter[str] = Counter()
    rows = contact_rows(table, rows_in)
    prepare_output_dir(command, out, _LAYER_OUTPUTS, args.overwrite)
    function, template = _COUNT[args.level]
    counted = function(rows)
    write_table(counted, out, "layer", args.table_format)
    description = LayerDescription(
        name="contact_count",
        template=template,
        parameters={"level": args.level, "contacts_layer": source},
        **_PANDORA,
    )
    _write_layer_summary(
        out,
        "annotation_count",
        description,
        contacts,
        started_at,
        {"rows_in": rows_in["rows_in"], "rows_out": len(counted)},
        input_summary(contacts),
    )
    print(
        f"counted: {len(counted)} {args.level} rows from {rows_in['rows_in']} "
        f"contacts -> {out}"
    )


def _write_labels(
    command: str,
    args: argparse.Namespace,
    started_at: str,
    layer: str | None,
    rows_in: Counter[str],
    labelled: list[Any],
    description: tuple[str, dict[str, Any]],
    extra: dict[str, int],
) -> Counter[str]:
    """Write an interface_labels layer folder for the contact layer
    `args.contacts` and print its result line; returns the label counts.
    `description` is the layer's name and parameters (with `labels`)."""

    contacts, out = Path(args.contacts), Path(args.out)
    name, parameters = description
    prepare_output_dir(command, out, _LAYER_OUTPUTS, args.overwrite)
    write_table(labelled, out, "layer", args.table_format)
    layer_description = LayerDescription(
        name=name,
        template="interface_labels",
        parameters={"contacts_layer": layer, **parameters},
        **_PANDORA,
    )
    labels = Counter(row.label for row in labelled)
    counts = {"rows_in": rows_in["rows_in"], "interfaces": len(labelled)}
    counts.update({label: labels[label] for label in parameters["labels"]})
    _write_layer_summary(
        out,
        "annotation_label",
        layer_description,
        contacts,
        started_at,
        counts | extra,
        input_summary(contacts),
    )
    shown = ", ".join(f"{labels[x]} {x}" for x in parameters["labels"])
    print(f"labelled: {len(labelled)} interfaces ({shown}) -> {out}")
    return labels


def _cmd_interface_contacts(args: argparse.Namespace) -> None:
    """Handle `annotations interface-contacts`: label few_contacts."""

    command = "annotations interface-contacts"
    started_at = now_iso()
    contacts = Path(args.contacts)
    check_out(command, Path(args.out), contacts)
    check_table_format(command, args.table_format)
    table = contact_table(command, contacts)
    layer = layer_name(contacts)
    rows_in: Counter[str] = Counter()
    rows = contact_rows(table, rows_in)
    labelled = label_interfaces(rows, {}, args.min_residues_per_chain, None)
    parameters = {
        "min_residues_per_chain": args.min_residues_per_chain,
        "labels": ["ok", "few_contacts"],
    }
    _write_labels(
        command,
        args,
        started_at,
        layer,
        rows_in,
        labelled,
        ("interface_contacts", parameters),
        {},
    )


def _cmd_dimer_length(args: argparse.Namespace) -> None:
    """Handle `annotations dimer-length`: label too_large."""

    command = "annotations dimer-length"
    started_at = now_iso()
    contacts, structures = Path(args.contacts), Path(args.structures)
    check_out(command, Path(args.out), contacts, structures)
    check_table_format(command, args.table_format)
    table = contact_table(command, contacts)
    layer = layer_name(contacts)
    _, chains_table = require_tables(command, structures)
    lengths = {
        (chain.entry_id.upper(), chain.chain_id): chain.seqres_length
        for chain in read_records(chains_table, DatasetChainRecord)
        if chain.seqres_length is not None
    }
    rows_in: Counter[str] = Counter()
    rows = contact_rows(table, rows_in)
    labelled = label_interfaces(rows, lengths, None, args.max_total_residues)
    parameters = {
        "structures_dir": str(structures),
        "max_total_residues": args.max_total_residues,
        "labels": ["ok", "too_large"],
    }
    unknown = sum(row.total_residues is None for row in labelled)
    _write_labels(
        command,
        args,
        started_at,
        layer,
        rows_in,
        labelled,
        ("dimer_length", parameters),
        {"unknown_length": unknown},
    )


def _cmd_gaps(args: argparse.Namespace) -> None:
    """Handle `annotations interface-gaps`: label interface_gap."""

    command = "annotations interface-gaps"
    started_at = now_iso()
    contacts, source = Path(args.contacts), Path(args.structures)
    check_out(command, Path(args.out), contacts, source)
    check_table_format(command, args.table_format)
    table = contact_table(command, contacts)
    layer = layer_name(contacts)
    paths, _ = structure_inputs(command, source)
    missing: dict[tuple[str, str], set[int]] = {}
    unmeasured = 0
    for path in paths:
        structure = read_structure(path)
        if structure is None:
            continue
        chains, diagnostics = missing_positions(structure)
        unmeasured += len(diagnostics.warnings)
        for chain_id, positions in chains.items():
            missing[(structure.entry_id.upper(), chain_id)] = positions
    rows_in: Counter[str] = Counter()
    rows = contact_rows(table, rows_in)
    labelled = label_interface_gaps(rows, missing)
    parameters = {
        "structures_dir": str(source),
        "missing_residue_definition": "incomplete_backbone",
        "labels": ["ok", "interface_gap", "unknown"],
    }
    labels = _write_labels(
        command,
        args,
        started_at,
        layer,
        rows_in,
        labelled,
        ("interface_gaps", parameters),
        {"chains_unmeasured": unmeasured},
    )
    if labels["unknown"]:
        print(
            f"{command}: warning: {labels['unknown']} interfaces have a "
            "chain whose missing residues can't be measured (no SEQRES, "
            "or label_seq_id renumbered); they are labelled unknown",
            file=sys.stderr,
        )


def _label_tables(command: str, paths: list[Path]) -> list[Path | None]:
    """The tables of the `--labels` layers (None: no rows); exits unless
    each is an interface_labels layer."""

    tables = []
    for path in paths:
        table, template = labels_table(command, path)
        if template != "interface_labels":
            raise SystemExit(
                f"{command}: {path} is a {template} layer; --labels needs "
                "interface_labels"
            )
        tables.append(table)
    return tables


def _source_of(
    summary: dict[str, Any], by_folder: dict[Path, dict[str, Any]]
) -> Path | None:
    """The input a layer's summary records; for a layer made from another
    given layer, that layer's input. None if not recorded."""

    seen: set[int] = set()
    while summary.get("input") is not None and id(summary) not in seen:
        seen.add(id(summary))
        recorded = Path(summary["input"])
        if not recorded.is_absolute() and not recorded.exists():
            return None  # relative to where an older step ran
        recorded = recorded.resolve()
        if recorded not in by_folder:
            return recorded
        summary = by_folder[recorded]
    return None


def _warn_if_sources_differ(
    command: str, folders: list[Path], summaries: list[dict[str, Any] | None]
) -> None:
    """Warn when the layers do not record the same input; no check
    unless every layer has a summary with one (a bare table has none)."""

    by_folder = {
        f.resolve(): s
        for f, s in zip(folders, summaries, strict=True)
        if s is not None
    }
    sources = {
        None if s is None else _source_of(s, by_folder) for s in summaries
    }
    if len(sources) > 1 and None not in sources:
        print(
            f"{command}: warning: the layers record different inputs "
            f"({', '.join(sorted(str(s) for s in sources))}); check that "
            "they use the same chain ids",
            file=sys.stderr,
        )


def _cmd_filter(args: argparse.Namespace) -> None:
    """Handle `annotations filter`: a contact layer without the
    interfaces the labels layers mark for dropping, as a new layer."""

    command = "annotations filter"
    started_at = now_iso()
    contacts, out = Path(args.contacts), Path(args.out)
    label_dirs = [Path(p) for p in args.labels]
    check_out(command, out, contacts, *label_dirs)
    check_table_format(command, args.table_format)
    table = contact_table(command, contacts)
    layer = layer_name(contacts)
    label_tables = _label_tables(command, label_dirs)
    summaries = [read_summary(path) for path in label_dirs]
    _warn_if_sources_differ(
        command,
        [contacts, *label_dirs],
        [read_summary(contacts), *summaries],
    )

    labels_by_entry: dict[str, list[dict[str, Any]]] = {}
    for label_table in label_tables:
        for row in read_rows(label_table) if label_table else ():
            labels_by_entry.setdefault(entry_key(row), []).append(row)
    drop = set(args.drop)
    # A label a layer can give (`parameters.labels`, recorded by Pandora's
    # labelling commands) is no typo even when no interface has it.
    present = {
        row["label"] for rows in labels_by_entry.values() for row in rows
    }
    known = {
        label
        for summary in summaries
        for label in (
            ((summary or {}).get("policy") or {}).get("parameters") or {}
        ).get("labels", [])
    }
    unknown = drop - present - known
    if unknown:
        raise SystemExit(
            f"{command}: no interface has the label(s) {sorted(unknown)}; "
            "check for typos"
        )

    prepare_output_dir(command, out, _LAYER_OUTPUTS, args.overwrite)
    counts = _filter_entries(
        table,
        labels_by_entry,
        drop,
        TableWriter(out, "layer", args.table_format),
    )
    description = LayerDescription(
        name="filtered_contacts",
        template="residue_contacts",
        parameters={
            "contacts_layer": layer,
            "labels_layers": [layer_name(path) for path in label_dirs],
            "labels_dirs": [str(p) for p in label_dirs],
            "drop": sorted(drop),
        },
        **_PANDORA,
    )
    _write_layer_summary(
        out,
        "annotation_filter",
        description,
        contacts,
        started_at,
        dict(counts),
        input_summary(contacts),
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


def _add_contacts(p: argparse.ArgumentParser) -> None:
    """Add the positional `contacts` input."""

    p.add_argument(
        "contacts",
        help="residue_contacts layer folder or .parquet/.jsonl table",
    )


def _register_contacts(
    commands: argparse._SubParsersAction, anything: str
) -> None:
    """Add the subcommands that compute contact layers."""

    p = commands.add_parser(
        "contacts",
        help="Compute chain-interface residue contacts (a residue_contacts "
        "layer).",
    )
    p.add_argument("structures", help=anything)
    add_out(p)
    p.add_argument(
        "--distance-cutoff", type=float, help="contact distance in angstrom"
    )
    p.add_argument(
        "--atom-set",
        choices=sorted(CONTACT_ATOM_SETS),
        help="which atoms count as in contact",
    )
    p.add_argument(
        "--polymer-types",
        action="append",
        help="keep chains of this entity polymer type (repeatable), e.g. "
        "polypeptide(L)",
    )
    p.set_defaults(func=_cmd_contacts)

    p = commands.add_parser(
        "ligand-contacts",
        help="Compute residue contacts of ligands (a residue_contacts layer).",
    )
    p.add_argument("structures", help=anything)
    add_out(p)
    p.add_argument(
        "--distance-cutoff", type=float, help="contact distance in angstrom"
    )
    p.add_argument(
        "--include-waters",
        action="store_const",
        const=True,
        default=None,
        help="count water molecules as ligands",
    )
    p.set_defaults(func=_cmd_ligand_contacts)


def _register_labels(commands: argparse._SubParsersAction) -> None:
    """Add the subcommands that write interface_labels layers."""

    p = commands.add_parser(
        "interface-contacts",
        help="Label interfaces few_contacts or ok by the residues in "
        "contact per chain (an interface_labels layer).",
    )
    _add_contacts(p)
    add_out(p)
    p.add_argument(
        "--min-residues-per-chain",
        type=int,
        required=True,
        help="label few_contacts when either chain has fewer residues in "
        "contact with the other",
    )
    p.set_defaults(func=_cmd_interface_contacts)

    p = commands.add_parser(
        "dimer-length",
        help="Label interfaces too_large or ok by the two chains' sequence "
        "lengths (an interface_labels layer).",
    )
    _add_contacts(p)
    add_out(p)
    p.add_argument(
        "--structures",
        required=True,
        help="structures folder with tables (from curate)",
    )
    p.add_argument(
        "--max-total-residues",
        type=int,
        required=True,
        help="label too_large when the two chains' sequence lengths add up "
        "to more",
    )
    p.set_defaults(func=_cmd_dimer_length)

    p = commands.add_parser(
        "interface-gaps",
        help="Label interfaces interface_gap when a residue missing "
        "backbone atoms lies inside either side (an interface_labels "
        "layer).",
    )
    _add_contacts(p)
    add_out(p)
    p.add_argument(
        "--structures",
        required=True,
        help="structures the contacts came from: any structures form",
    )
    p.set_defaults(func=_cmd_gaps)


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add the `annotations` group to the `pandora` parser."""

    group = subparsers.add_parser(
        "annotations",
        help="Compute, label, filter or import annotation layers.",
    )
    commands = group.add_subparsers(dest="annotations_command", required=True)
    _register_contacts(
        commands,
        "structures folder, folder of .cif/.cif.gz files, entries.txt list "
        "or one .cif/.cif.gz file",
    )

    p = commands.add_parser(
        "contact-counts",
        help="Count a contact layer's partners per residue, or residues "
        "in contact per chain (a value layer).",
    )
    _add_contacts(p)
    add_out(p)
    p.add_argument(
        "--level",
        required=True,
        choices=list(_COUNT),
        help="count per residue or per chain",
    )
    p.set_defaults(func=_cmd_count)

    _register_labels(commands)

    p = commands.add_parser(
        "filter",
        help="Drop labelled interfaces from a contact layer.",
    )
    _add_contacts(p)
    add_out(p)
    p.add_argument(
        "--labels",
        action="append",
        required=True,
        help="interface_labels layer folder or .parquet/.jsonl table "
        "(repeatable); an interface is dropped when any layer gives it a "
        "--drop label",
    )
    p.add_argument(
        "--drop",
        action="append",
        required=True,
        help="an interface label that means drop (repeatable)",
    )
    p.set_defaults(func=_cmd_filter)

    p = commands.add_parser(
        "import",
        help="Link a table from an external tool to structures, as an "
        "annotation layer.",
    )
    p.add_argument("table", help="the tool's table (.parquet, .jsonl, .json)")
    add_out(p)
    p.add_argument(
        "--structures",
        required=True,
        help="structures the table's rows refer to: any structures form",
    )
    p.add_argument("--description", required=True, help="LayerDescription YAML")
    p.add_argument(
        "--source-files",
        help="the tool's original output file or folder, copied into OUT",
    )
    p.set_defaults(func=_cmd_import)
