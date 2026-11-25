"""CLI `similarity`: searches and imports."""

from __future__ import annotations

import argparse
import itertools
import tempfile
from collections import Counter
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from pandora import __version__
from pandora._util import now_iso
from pandora.cli._io import (
    SUMMARY,
    ParquetSink,
    add_out,
    check_out,
    contact_rows,
    contact_table,
    copy_source,
    error_message,
    input_summary,
    layer_name,
    load_json_list,
    load_yaml_model,
    prepare_output_dir,
    read_columns,
    read_structures_dir,
    require_tables,
    resolve_file,
    write_summary,
)
from pandora.export import (
    export_chain_mmcif,
    read_records,
    write_records,
)
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import (
    LayerDescription,
    SimilarityItems,
)
from pandora.schemas.dataset import DatasetChainRecord
from pandora.schemas.similarity import (
    RepresentativeCluster,
    SimilarityCluster,
)
from pandora.similarity import (
    compute_sequence_similarity,
    compute_structure_similarity,
    interface_similarity,
)
from pandora.similarity.hits import (
    read_hit_tsv,
    read_similarity_table,
)
from pandora.similarity.layer import (
    ItemSet,
    chain_items,
    check_similarity_columns,
    entry_items,
    interface_items,
    item_id,
    match_similarity_batches,
    parse_item_id,
)

_CLI_ENGINES = {"mmseqs2": "MMseqs2", "foldseek": "Foldseek"}

# What a similarity command writes: cleared by --overwrite.
_LAYER_OUTPUTS = [
    "layer.parquet",
    "unmatched.parquet",
    "source",
    SUMMARY,
    "hits.tsv",
]


def _require_pyarrow(command: str) -> None:
    """Exit unless pyarrow (the `export` extra) is installed."""

    try:
        import pyarrow  # noqa: F401
    except ImportError as exc:
        raise SystemExit(f"{command}: similarity layers need pyarrow (the export extra)") from exc


_CONTACT_COLUMNS = (
    "entry_id",
    "chain_id_1",
    "chain_id_2",
    "seq_id_1",
    "seq_id_2",
)


def _item_set(
    command: str,
    level: SimilarityItems,
    structures: Path | None,
    contacts: Path | None,
) -> tuple[ItemSet, Path | None]:
    """The items a layer at `level` may name, read from a structures
    folder's tables or a contact layer, and the summary of where they
    came from."""

    if level == "interface":
        if contacts is None:
            raise SystemExit(
                f"{command}: interface items come from --contacts "
                "(a residue_contacts layer folder or table)"
            )
        rows_file = contact_table(command, contacts)
        rows = read_columns(rows_file, _CONTACT_COLUMNS) if rows_file is not None else iter(())
        return interface_items(rows), input_summary(contacts)
    if structures is None:
        raise SystemExit(f"{command}: {level} items need --structures")
    entries_table, chains_table = require_tables(command, structures)
    previous = input_summary(structures)
    if level == "entry":
        rows = read_columns(entries_table, ("entry_id", "structure_path"))
        entries = [r["entry_id"] for r in rows if r["structure_path"]]
        return entry_items(entries), previous
    rows = read_columns(chains_table, ("entry_id", "chain_id"))
    return chain_items((r["entry_id"], r["chain_id"]) for r in rows), previous


def _write_layer(batches: Iterable[Any], items: ItemSet, out: Path) -> Counter[str]:
    """Match `batches` against `items` into `out`'s layer and unmatched
    tables; returns the counts."""

    counts: Counter[str] = Counter(rows_in=0, rows_matched=0, rows_unmatched=0)
    with (
        ParquetSink(out / "layer.parquet") as layer,
        ParquetSink(out / "unmatched.parquet") as unmatched,
    ):
        for matched, rest in match_similarity_batches(batches, items):
            layer.write(matched)
            unmatched.write(rest)
            counts["rows_in"] += matched.num_rows + rest.num_rows
            counts.update(f"unmatched:{r}" for r in rest.column("reason").to_pylist())
        counts["rows_matched"] = layer.rows
        counts["rows_unmatched"] = unmatched.rows
    return counts


def _write_summary(
    step: str,
    description: LayerDescription,
    *,
    input: Path,
    out: Path,
    started_at: str,
    counts: Counter[str],
    previous: Path | None,
) -> None:
    """Write the layer folder's `summary.json`."""

    write_summary(
        out,
        build_step_summary(
            step,
            description,
            input=input.resolve(),
            output=out,
            started_at=started_at,
            counts=dict(sorted(counts.items())),
            previous_summary=previous,
        ),
    )


def _search_options(args: argparse.Namespace, engine: str) -> dict[str, Any]:
    """Keyword arguments for the search function from the CLI flags."""

    options: dict[str, Any] = {
        "sensitivity": args.sensitivity,
        "max_seqs": args.max_seqs,
    }
    if engine == "mmseqs2":
        options["mmseqs_bin"] = args.bin
    else:
        options |= {
            "foldseek_bin": args.bin,
            "alignment_type": args.alignment_type,
            "exhaustive_search": args.exhaustive_search,
        }
    return options


def _read_chain_list(path: Path) -> list[str]:
    """Chain item ids from a file, one per line (blank lines and `#`
    comments skipped)."""

    lines = (line.strip() for line in path.read_text().splitlines())
    return [line for line in lines if line and not line.startswith("#")]


def _chains_arg(command: str, path: Path) -> list[str]:
    """Chain ids from `--chains`: a representatives folder or `.json`
    (its representatives), or a text list."""

    if path.is_dir() or path.suffix == ".json":
        file = resolve_file(command, path, "representatives.json")
        clusters = load_json_list(command, RepresentativeCluster, file)
        return [cluster.representative for cluster in clusters]
    if not path.is_file():
        raise SystemExit(f"{command}: {path} does not exist")
    return _read_chain_list(path)


def _search(args: argparse.Namespace, command: str, engine: str) -> None:
    """Run MMseqs2 or Foldseek over chains and write the hits as a
    similarity layer folder."""

    started_at = now_iso()
    _require_pyarrow(command)
    chain_files = getattr(args, "chain_files", None)
    if args.chains is not None and chain_files is not None:
        raise SystemExit(
            f"{command}: --chains needs the chains cut from STRUCTURES; drop --chain-files"
        )
    structures, out = Path(args.structures), Path(args.out)
    check_out(
        command,
        out,
        structures,
        *(Path(p) for p in (args.chains, chain_files) if p is not None),
    )
    require_tables(command, structures)
    items, previous = _item_set(command, "chain", structures, None)
    chains: set[str] | None = None
    if args.chains is not None:
        known = set(items.ids.values())
        listed = _chains_arg(command, Path(args.chains))
        unknown = next((c for c in listed if c not in known), None)
        if unknown is not None:
            raise SystemExit(f"{command}: --chains names {unknown}, not a chain of {structures}")
        chains = set(listed)
    prepare_output_dir(command, out, _LAYER_OUTPUTS, args.overwrite)
    hits = out / "hits.tsv"
    options = _search_options(args, engine)
    counts_extra: dict[str, int] = {}
    if engine == "mmseqs2":
        sequences = _chain_sequences(command, structures)
        if chains is not None:
            sequences = {k: v for k, v in sequences.items() if k in chains}
        counts_extra["chains_searched"] = len(sequences)
        description = compute_sequence_similarity(sequences, hits, **options)
        searched: Path = structures
    elif chain_files is not None:
        searched = Path(chain_files)
        description = compute_structure_similarity(searched, hits, **options)
    else:
        searched = structures
        with tempfile.TemporaryDirectory() as tmp:
            files = _chain_files(command, structures, chains, Path(tmp))
            counts_extra["chains_searched"] = len(files)
            description = compute_structure_similarity(files, hits, **options)
    if args.chains is not None:
        searched = Path(args.chains)
    counts = _write_layer(read_hit_tsv(hits, _CLI_ENGINES[engine]), items, out)
    counts.update(counts_extra)
    hits.unlink()
    _write_summary(
        "similarity_search",
        description,
        input=searched,
        out=out,
        started_at=started_at,
        counts=counts,
        previous=previous,
    )
    print(f"searched: {counts['rows_matched']} hits, {counts['rows_unmatched']} unmatched -> {out}")


def _cmd_mmseqs(args: argparse.Namespace) -> None:
    """Handle `similarity mmseqs`."""

    _search(args, "similarity mmseqs", "mmseqs2")


def _cmd_foldseek(args: argparse.Namespace) -> None:
    """Handle `similarity foldseek`."""

    _search(args, "similarity foldseek", "foldseek")


def _chain_files(
    command: str, structures_dir: Path, chains: set[str] | None, out: Path
) -> dict[str, Path]:
    """Write one mmCIF per non-copy chain of the step folder (or only
    `chains`) into `out`; returns item id -> file."""

    table = require_tables(command, structures_dir)[1]
    wanted = {
        item_id("chain", r["entry_id"], r["chain_id"]): r
        for r in read_columns(table, ("entry_id", "chain_id", "copy_of"))
        if not r["copy_of"]
    }
    if chains is not None:
        wanted = {k: r for k, r in wanted.items() if k in chains}
    by_entry: dict[str, list[tuple[str, str]]] = {}
    for key, row in wanted.items():
        by_entry.setdefault(row["entry_id"], []).append((key, row["chain_id"]))
    # ponytail: parses every entry; limit to the listed entries if slow.
    files: dict[str, Path] = {}
    for structure in read_structures_dir(structures_dir):
        for key, chain_id in by_entry.pop(structure.entry_id, []):
            files[key] = export_chain_mmcif(structure, chain_id, out / f"{key}.cif")
    if by_entry:
        raise SystemExit(f"{command}: {structures_dir} has no structure for {min(by_entry)}")
    return dict(sorted(files.items()))


def _chain_sequences(command: str, structures_dir: Path) -> dict[str, str]:
    """Item ID -> sequence for the step folder's chains, without
    symmetry copies (they would crowd out real homologues)."""

    table = require_tables(command, structures_dir)[1]
    columns = ("entry_id", "chain_id", "sequence", "copy_of")
    return {
        item_id("chain", r["entry_id"], r["chain_id"]): r["sequence"]
        for r in read_columns(table, columns)
        if r["sequence"] and not r["copy_of"]
    }


def _import_error(command: str, table: Path, exc: Exception) -> SystemExit:
    """The import's one-line exit, naming `table` once."""

    message = error_message(exc).removeprefix(f"{table}: ")
    return SystemExit(f"{command}: cannot import {table}: {message}")


def _cmd_import(args: argparse.Namespace) -> None:
    """Handle `similarity import`: link a user's edge list to the items
    of a step folder or contact layer, as a similarity layer folder."""

    command = "similarity import"
    started_at = now_iso()
    _require_pyarrow(command)
    table, out = Path(args.table), Path(args.out)
    structures = Path(args.structures) if args.structures else None
    contacts = Path(args.contacts) if args.contacts else None
    check_out(command, out, table, structures, contacts)
    description = load_yaml_model(command, LayerDescription, Path(args.description))
    try:
        if description.template != "similarity":
            raise ValueError(
                f"template is {description.template!r}; `similarity import` takes similarity layers"
            )
        batches = read_similarity_table(table)
        first = next(batches, None)
        if first is not None:
            check_similarity_columns(first.schema.names)
    except (OSError, ValueError) as exc:
        raise _import_error(command, table, exc) from exc

    if description.items != "interface" and structures is None:
        raise SystemExit(f"{command}: {description.items} items need --structures")
    if structures is not None and description.items != "interface":
        require_tables(command, structures)
    items, previous = _item_set(command, description.items, structures, contacts)
    prepare_output_dir(command, out, _LAYER_OUTPUTS, args.overwrite)
    if args.source_files:
        copy_source(Path(args.source_files), out / "source")
    rest = [] if first is None else itertools.chain([first], batches)
    try:
        counts = _write_layer(rest, items, out)
    except ValueError as exc:
        raise _import_error(command, table, exc) from exc
    _write_summary(
        "similarity_import",
        description,
        input=table,
        out=out,
        started_at=started_at,
        counts=counts,
        previous=previous,
    )
    print(f"imported: {counts['rows_matched']} rows, {counts['rows_unmatched']} unmatched -> {out}")


def _cmd_interface_overlap(args: argparse.Namespace) -> None:
    """Handle `similarity interface-overlap`: score pairs of a contact
    layer's interfaces by per-side Jaccard, as a similarity layer
    folder."""

    command = "similarity interface-overlap"
    started_at = now_iso()
    _require_pyarrow(command)
    contacts, out = Path(args.contacts), Path(args.out)
    structures = Path(args.structures)
    check_out(command, out, contacts, structures, Path(args.chain_clusters))
    table = contact_table(command, contacts)
    _, chains_table = require_tables(command, structures)
    sequences = {
        (chain.entry_id.upper(), chain.chain_id): chain.sequence
        for chain in read_records(chains_table, DatasetChainRecord)
        if chain.sequence
    }
    clusters_file = resolve_file(command, Path(args.chain_clusters), "clusters.json")
    clusters = load_json_list(command, SimilarityCluster, clusters_file)
    try:
        for cluster in clusters:
            for component in cluster.components:
                parse_item_id("chain", component)
    except ValueError as exc:
        raise SystemExit(
            f"{command}: {clusters_file} must cluster chain items (ENTRY_CHAIN): {exc}"
        ) from exc
    rows_in: Counter[str] = Counter()
    scored, counts = interface_similarity(
        contact_rows(table, rows_in), clusters, sequences, args.min_identity
    )

    prepare_output_dir(command, out, _LAYER_OUTPUTS, args.overwrite)
    if scored:
        write_records(scored, out / "layer.parquet")
    description = LayerDescription(
        name="iface_jaccard",
        template="similarity",
        items="interface",
        tool="pandora",
        tool_version=__version__,
        parameters={
            "contacts_layer": layer_name(contacts),
            "structures": str(structures.resolve()),
            "chain_clusters": str(clusters_file.resolve()),
            "min_identity": args.min_identity,
        },
    )
    _write_summary(
        "similarity_interfaces",
        description,
        input=contacts,
        out=out,
        started_at=started_at,
        counts=Counter(rows_in | counts),
        previous=input_summary(contacts),
    )
    print(
        f"scored: {counts['pairs_compared']} interface pairs "
        f"({counts['pairs_zero_score']} zero, "
        f"{counts['pairs_below_identity']} below identity) from "
        f"{counts['interfaces']} interfaces -> {out}"
    )


def _register_search(
    commands: argparse._SubParsersAction,
    name: str,
    engine: str,
    help: str,
    func: Any,
) -> argparse.ArgumentParser:
    """Add `similarity <name>` with the flags both searches share."""

    p = commands.add_parser(name, help=help)
    p.add_argument(
        "structures",
        help="structures folder with entries and chains tables (pandora structures curate)",
    )
    add_out(p, tables=False)
    p.add_argument(
        "--chains",
        help="limit the search to these chains: a representatives folder "
        "or .json, or a text file of chain IDs (one per line, # comments)",
    )
    p.add_argument("--bin", default=name, help=f"the {engine} executable (default {name})")
    p.set_defaults(func=func)
    return p


def _register_interface_overlap(commands: argparse._SubParsersAction) -> None:
    """Add `similarity interface-overlap`."""

    p = commands.add_parser(
        "interface-overlap",
        help="Score pairs of interfaces of the same two proteins by "
        "per-side Jaccard of their interface residues.",
    )
    p.add_argument("contacts", help="contact layer folder or .parquet/.jsonl table")
    add_out(p, tables=False)
    p.add_argument(
        "--chain-clusters",
        required=True,
        help="clusters.json of chains (folder or file), e.g. 90%% sequence "
        "clusters: only interfaces whose chains share clusters are compared",
    )
    p.add_argument(
        "--structures",
        required=True,
        help="structures folder with a chains table: the chain sequences",
    )
    p.add_argument(
        "--min-identity",
        type=float,
        default=0.9,
        help="lowest pairwise identity per chain pair (default 0.9)",
    )
    p.set_defaults(func=_cmd_interface_overlap)


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    similarity = subparsers.add_parser(
        "similarity",
        help="Similarity layers: search with MMseqs2 or Foldseek, score "
        "interface overlap, or import.",
    )
    commands = similarity.add_subparsers(dest="similarity_command", required=True)
    p = _register_search(
        commands,
        "mmseqs",
        "MMseqs2",
        "All-vs-all MMseqs2 sequence search over a structures folder's "
        "chains, written as a similarity layer.",
        _cmd_mmseqs,
    )
    p.add_argument(
        "--sensitivity",
        type=float,
        default=5.7,
        help="MMseqs2 -s (default 5.7)",
    )
    p.add_argument(
        "--max-seqs",
        type=int,
        default=300,
        help="hits kept per query (default 300; raise it for leakage control)",
    )

    p = _register_search(
        commands,
        "foldseek",
        "Foldseek",
        "All-vs-all Foldseek structure search over a structures folder's "
        "chains, written as a similarity layer.",
        _cmd_foldseek,
    )
    p.add_argument(
        "--chain-files",
        help="folder of per-chain structure files to search instead of "
        "cutting them from STRUCTURES (not with --chains)",
    )
    p.add_argument(
        "--sensitivity",
        type=float,
        default=9.5,
        help="Foldseek -s (default 9.5)",
    )
    p.add_argument(
        "--alignment-type",
        type=int,
        default=2,
        help="Foldseek --alignment-type (default 2)",
    )
    p.add_argument(
        "--max-seqs",
        type=int,
        default=1000,
        help="hits kept per query (default 1000; raise it for leakage control)",
    )
    p.add_argument(
        "--exhaustive-search",
        action="store_true",
        help="skip the prefilter",
    )

    _register_interface_overlap(commands)

    p = commands.add_parser(
        "import",
        help="Link an edge list (Parquet, or TSV with a header) to a "
        "structures folder's entries/chains or a contact layer's "
        "interfaces.",
    )
    p.add_argument("table", help="edge list (.parquet, or .tsv with header)")
    add_out(p, tables=False)
    p.add_argument("--description", required=True, help="LayerDescription YAML")
    p.add_argument(
        "--structures",
        help="structures folder with tables, for entry or chain items",
    )
    p.add_argument(
        "--contacts",
        help="residue_contacts layer folder or table, for interface items",
    )
    p.add_argument(
        "--source-files",
        help="the tool's original output file or folder, copied into OUT",
    )
    p.set_defaults(func=_cmd_import)
