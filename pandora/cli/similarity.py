"""CLI `similarity`, `cluster` and `partition`: searches, clusters and leakage-safe splits."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from pandora._util import now_iso
from pandora.cli._io import (
    ParquetSink,
    check_contact_layer,
    copy_source,
    find_table,
    load_json_list,
    load_pairs,
    load_yaml_model,
    prepare_output_dir,
    read_columns,
    read_structures_dir,
)
from pandora.export import read_records, write_json, write_records
from pandora.ingestion import load_policy
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import (
    LayerDescription,
    LayerRow,
    ResidueContactRow,
    SimilarityItems,
)
from pandora.schemas.similarity import (
    HitFilter,
    SimilarityCluster,
)
from pandora.similarity import (
    cluster_edges,
    compute_sequence_similarity,
    compute_structure_similarity,
    entry_edges,
    interface_residues_from_annotation,
    pair_cluster_keys,
    partition_dataset,
)
from pandora.similarity.hits import (
    read_hit_tsv,
    read_similarity_table,
    similarity_edges,
)
from pandora.similarity.layer import (
    ItemSet,
    chain_items,
    check_similarity_columns,
    entry_items,
    interface_items,
    item_id,
    item_labels,
    match_similarity_batches,
)

_CLI_ENGINES = {"mmseqs2": "MMseqs2", "foldseek": "Foldseek"}

# What `similarity search` / `import` write: cleared by --overwrite.
_LAYER_OUTPUTS = [
    "layer.parquet",
    "unmatched.parquet",
    "source",
    "layer.json",
    "hits.tsv",
]
_STEP_SUMMARIES = ("curation.json", "canonicalisation.json")
_CONTACT_COLUMNS = (
    "entry_id",
    "chain_id_1",
    "chain_id_2",
    "seq_id_1",
    "seq_id_2",
)


def _require_pyarrow(command: str) -> None:
    """Exit unless pyarrow (the `export` extra) is installed."""

    try:
        import pyarrow  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            f"{command}: similarity layers need pyarrow (the export extra)"
        ) from exc


def _step_table(command: str, folder: Path, name: str, writer: str) -> Path:
    """`folder`'s `<name>.<fmt>` table, or exit naming the step that
    writes it."""

    path = find_table(folder, name)
    if path is None:
        raise SystemExit(
            f"{command}: {folder} has no {name} table; it is written by "
            f"`{writer}`"
        )
    return path


def _contact_layer(command: str, folder: Path) -> Path:
    """The rows file of a `residue_contacts` layer folder."""

    check_contact_layer(command, folder)
    return find_table(folder, "layer") or folder / "layer.json"


def _item_set(
    command: str,
    level: SimilarityItems,
    structures_dir: Path | None,
    contacts: Path | None,
) -> tuple[ItemSet, Path | None]:
    """The items a layer at `level` may name, read from a step folder's
    tables or a contact layer, and the summary of where they came from."""

    if level == "interface":
        if contacts is None:
            raise SystemExit(
                f"{command}: interface items come from --contacts "
                "(a residue_contacts layer folder)"
            )
        rows_file = _contact_layer(command, contacts)
        rows = (
            read_columns(rows_file, _CONTACT_COLUMNS)
            if rows_file.name != "layer.json"
            else iter(())
        )
        return interface_items(rows), contacts / "layer.json"
    if structures_dir is None:
        raise SystemExit(f"{command}: {level} items need --structures-dir")
    previous = next(
        (
            structures_dir / name
            for name in _STEP_SUMMARIES
            if (structures_dir / name).exists()
        ),
        None,
    )
    if level == "entry":
        table = _step_table(
            command, structures_dir, "entries", "pandora canonicalise"
        )
        rows = read_columns(table, ("entry_id", "structure_path"))
        entries = [r["entry_id"] for r in rows if r["structure_path"]]
        return entry_items(entries), previous
    table = _step_table(command, structures_dir, "chains", "pandora curate")
    rows = read_columns(table, ("entry_id", "chain_id"))
    return chain_items((r["entry_id"], r["chain_id"]) for r in rows), previous


def _write_layer(
    batches: Iterable[Any], items: ItemSet, out: Path
) -> Counter[str]:
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
            counts.update(
                f"unmatched:{r}" for r in rest.column("reason").to_pylist()
            )
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
    """Write the layer folder's `layer.json`."""

    write_json(
        build_step_summary(
            step,
            description,
            input=input,
            output=out,
            started_at=started_at,
            counts=dict(sorted(counts.items())),
            previous_summary=previous,
        ),
        out / "layer.json",
    )


def _search_options(args: argparse.Namespace) -> dict[str, Any]:
    """Keyword arguments for the search function from the CLI flags."""

    if args.engine == "mmseqs2":
        options: dict[str, Any] = {"mmseqs_bin": args.mmseqs_bin}
    else:
        options = {
            "foldseek_bin": args.foldseek_bin,
            "alignment_type": args.alignment_type,
            "exhaustive_search": args.exhaustive_search,
        }
    for name in ("sensitivity", "max_seqs"):
        if getattr(args, name) is not None:
            options[name] = getattr(args, name)
    return options


def _cmd_search(args: argparse.Namespace) -> None:
    """Handle `similarity search`: run MMseqs2 or Foldseek over chains
    and write the hits as a similarity layer folder."""

    started_at = now_iso()
    _require_pyarrow("similarity search")
    if args.engine == "mmseqs2" and args.input_dir is not None:
        raise SystemExit(
            "similarity search: --input-dir is for foldseek; mmseqs2 reads "
            "chain sequences from --structures-dir"
        )
    if args.engine == "foldseek" and args.input_dir is None:
        raise SystemExit(
            "similarity search: foldseek needs --input-dir (the per-chain "
            "structure files)"
        )
    structures_dir = Path(args.structures_dir)
    items, previous = _item_set(
        "similarity search", "chain", structures_dir, None
    )
    out = Path(args.output_dir)
    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    hits = out / "hits.tsv"
    if args.engine == "mmseqs2":
        description = compute_sequence_similarity(
            _chain_sequences(structures_dir), hits, **_search_options(args)
        )
        searched: Path = structures_dir
    else:
        searched = Path(args.input_dir)
        description = compute_structure_similarity(
            searched, hits, **_search_options(args)
        )
    counts = _write_layer(
        read_hit_tsv(hits, _CLI_ENGINES[args.engine]), items, out
    )
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
    print(
        f"searched: {counts['rows_matched']} hits, "
        f"{counts['rows_unmatched']} unmatched -> {out}"
    )


def _chain_sequences(structures_dir: Path) -> dict[str, str]:
    """Item ID -> sequence for the step folder's chains, without
    symmetry copies (they would crowd out real homologues)."""

    table = _step_table(
        "similarity search", structures_dir, "chains", "pandora curate"
    )
    columns = ("entry_id", "chain_id", "sequence", "copy_of")
    return {
        item_id("chain", r["entry_id"], r["chain_id"]): r["sequence"]
        for r in read_columns(table, columns)
        if r["sequence"] and not r["copy_of"]
    }


def _import_error(table: Path, exc: Exception) -> str:
    """The import's error message, naming `table` once."""

    message = str(exc)
    if message.startswith(str(table)):
        return f"cannot import {message}"
    return f"cannot import {table}: {message}"


def _cmd_import(args: argparse.Namespace) -> None:
    """Handle `similarity import`: link a user's edge list to the items
    of a step folder or contact layer, as a similarity layer folder."""

    started_at = now_iso()
    _require_pyarrow("similarity import")
    table = Path(args.table)
    try:
        description = load_policy(args.description, model=LayerDescription)
        if description.template != "similarity":
            raise ValueError(
                f"template is {description.template!r}; `similarity "
                "import` takes similarity layers"
            )
        batches = read_similarity_table(table)
        first = next(batches, None)
        if first is not None:
            check_similarity_columns(first.schema.names)
    except (ValueError, FileNotFoundError) as exc:
        raise SystemExit(_import_error(table, exc)) from exc

    structures_dir = Path(args.structures_dir) if args.structures_dir else None
    contacts = Path(args.contacts) if args.contacts else None
    items, previous = _item_set(
        "similarity import", description.items, structures_dir, contacts
    )
    out = Path(args.output_dir)
    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    if args.source_files:
        copy_source(Path(args.source_files), out / "source")
    rest = [] if first is None else itertools.chain([first], batches)
    try:
        counts = _write_layer(rest, items, out)
    except ValueError as exc:
        raise SystemExit(_import_error(table, exc)) from exc
    _write_summary(
        "similarity_import",
        description,
        input=table,
        out=out,
        started_at=started_at,
        counts=counts,
        previous=previous,
    )
    print(
        f"imported: {counts['rows_matched']} rows, "
        f"{counts['rows_unmatched']} unmatched -> {out}"
    )


def _hit_filter(args: argparse.Namespace) -> HitFilter:
    """The HitFilter from --hit-filter, with --min-score / --min-coverage
    overriding it; exits if no threshold is set."""

    if (
        args.hit_filter is None
        and args.min_score is None
        and (args.min_coverage is None)
    ):
        # Without a filter every reported hit would become an edge (the
        # old --threshold defaulted to 0.9), merging clusters silently.
        raise SystemExit(
            "cluster: set --min-score, --min-coverage or --hit-filter; "
            "with no threshold every reported hit becomes an edge"
        )
    hit_filter = (
        load_yaml_model(HitFilter, Path(args.hit_filter))
        if args.hit_filter
        else HitFilter()
    )
    overrides = {
        field: value
        for field, value in (
            ("min_score", args.min_score),
            ("min_coverage", args.min_coverage),
        )
        if value is not None
    }
    return hit_filter.model_copy(update=overrides)


def _similarity_level(folder: Path) -> SimilarityItems:
    """The item level of a similarity layer folder."""

    summary = folder / "layer.json"
    if not summary.exists():
        raise SystemExit(
            f"cluster: {folder} is not a layer folder (no layer.json)"
        )
    policy = json.loads(summary.read_text()).get("policy") or {}
    if policy.get("template") != "similarity":
        raise SystemExit(
            f"cluster: {folder} is a {policy.get('template')} layer, "
            "--similarity needs a similarity layer"
        )
    return policy["items"]


def _cluster_level(
    args: argparse.Namespace, level: SimilarityItems
) -> SimilarityItems:
    """The level to cluster at: the layer's, or entries from chains."""

    items = args.items or level
    if items != level and not (items == "entry" and level == "chain"):
        raise SystemExit(
            f"cluster: --items {items} does not fit a {level} layer; only "
            "--items entry on a chain layer rolls up"
        )
    if level == "interface" and args.structures_dir is not None:
        raise SystemExit(
            "cluster: interface items come from --contacts; "
            "--structures-dir does not apply"
        )
    if level == "entry" and args.contacts is not None:
        raise SystemExit("cluster: --contacts does not apply to an entry layer")
    return items


def _interface_residues(
    structures_dir: Path, contacts: Path
) -> dict[str, set[int]]:
    """Interface residue positions per chain item, from a contact layer
    and the structures it refers to."""

    rows_file = _contact_layer("cluster", contacts)
    rows = (
        read_records(rows_file, ResidueContactRow)
        if rows_file.name != "layer.json"
        else []
    )
    # ponytail: parses every structure; fine per dataset, not per archive.
    structures = read_structures_dir(structures_dir)
    try:
        return interface_residues_from_annotation(structures, rows)
    except KeyError as exc:
        raise SystemExit(
            f"cluster: {contacts} names an entry missing from "
            f"{structures_dir}: {exc}"
        ) from exc


def _copy_links(structures_dir: Path) -> Iterator[tuple[str, str]]:
    """(copy, source) chain item pairs from the step folder's chains
    table, so symmetry copies join their source chain's cluster."""

    table = _step_table("cluster", structures_dir, "chains", "pandora curate")
    for row in read_columns(table, ("entry_id", "chain_id", "copy_of")):
        if row["copy_of"]:
            yield (
                item_id("chain", row["entry_id"], row["chain_id"]),
                item_id("chain", row["entry_id"], row["copy_of"]),
            )


def _cmd_cluster(args: argparse.Namespace) -> None:
    """Handle the `cluster` subcommand: connected-component clustering
    of a similarity layer."""

    hit_filter = _hit_filter(args)
    _require_pyarrow("cluster")
    layer = Path(args.similarity)
    level = _similarity_level(layer)
    items_level = _cluster_level(args, level)
    structures_dir = Path(args.structures_dir) if args.structures_dir else None
    contacts = Path(args.contacts) if args.contacts else None
    items, _ = _item_set("cluster", items_level, structures_dir, contacts)

    interface_residues = None
    if level == "chain" and contacts is not None:
        interface_residues = _interface_residues(structures_dir, contacts)
    table = layer / "layer.parquet"
    try:
        edges: Iterable[Any] = (
            similarity_edges(
                read_similarity_table(table),
                hit_filter,
                interface_residues=interface_residues,
            )
            if table.exists()
            else iter(())
        )
        if items_level != level:
            edges = entry_edges(edges)
        elif level == "chain":
            edges = itertools.chain(edges, _copy_links(structures_dir))
        clusters, prov = cluster_edges(sorted(set(items.ids.values())), edges)
    except ValueError as exc:
        raise SystemExit(f"cluster: {layer}: {exc}") from exc
    prov = prov.model_copy(
        update={
            "hit_filter": hit_filter,
            "similarity_layer": str(layer / "layer.json"),
        }
    )

    output = Path(args.output)
    write_records(clusters, output)
    write_json(prov, output.parent / "cluster_provenance.json")
    print(f"{len(clusters)} clusters from {prov.n_edges} edges -> {output}")

    if args.pairs:
        pairs = load_pairs(Path(args.pairs))
        keys = pair_cluster_keys(pairs, clusters)
        pairs_output = output.parent / "cluster_pairs.json"
        pairs_output.write_text(
            json.dumps(
                [
                    {
                        "item_id_1": item_1,
                        "item_id_2": item_2,
                        "cluster_id_1": cluster_1,
                        "cluster_id_2": cluster_2,
                    }
                    for (item_1, item_2), (
                        cluster_1,
                        cluster_2,
                    ) in keys.items()
                ],
                indent=2,
            ),
            encoding="utf-8",
        )
        print(f"{len(keys)} paired cluster keys -> {pairs_output}")


def _partition_labels(args: argparse.Namespace) -> dict[str, str] | None:
    """The item labels from `--labels`, or None; exits on a bad layer."""

    if args.labels_dir is None:
        return None
    command = "partition"
    folder = Path(args.labels_dir)
    if args.no_keep_similar_items:
        raise SystemExit(
            f"{command}: --labels needs whole clusters; drop "
            "--no-keep-similar-items"
        )
    if not (folder / "layer.json").exists():
        raise SystemExit(f"{command}: {folder} is not a layer folder")
    policy = json.loads((folder / "layer.json").read_text()).get("policy")
    template = (policy or {}).get("template")
    table = find_table(folder, "layer")
    rows = read_records(table, LayerRow) if table else []
    try:
        return item_labels(rows, template)
    except ValueError as exc:
        raise SystemExit(f"{command}: {folder}: {exc}") from exc


def _cmd_partition(args: argparse.Namespace) -> None:
    """Handle the `partition` subcommand: leakage-safe train/val/test
    split of clusters."""

    labels = _partition_labels(args)
    clusters = load_json_list(SimilarityCluster, Path(args.clusters))
    splits, prov = partition_dataset(
        clusters,
        args.pct_train,
        args.pct_val,
        args.pct_test,
        keep_similar_items=not args.no_keep_similar_items,
        labels=labels,
    )
    if labels is not None:
        prov = prov.model_copy(update={"labels_layer": args.labels_dir})

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(splits, indent=2), encoding="utf-8")
    write_json(prov, output.parent / "partition_provenance.json")
    print(f"split sizes: {prov.split_sizes} -> {output}")
    if labels is not None:
        items = [i for c in clusters for i in c.components]
        unlabelled = sum(i not in labels for i in items)
        print(f"unlabelled items: {unlabelled} of {len(items)}")
        for split, missing in prov.missing_labels.items():
            print(
                f"partition: warning: {split} has no item labelled "
                f"{', '.join(missing)}",
                file=sys.stderr,
            )


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    similarity = subparsers.add_parser(
        "similarity", help="Similarity layers: search or import."
    )
    commands = similarity.add_subparsers(
        dest="similarity_command", required=True
    )
    p = commands.add_parser(
        "search",
        help="All-vs-all MMseqs2 or Foldseek search over a step folder's "
        "chains, written as a similarity layer.",
    )
    p.add_argument("--engine", choices=["mmseqs2", "foldseek"], required=True)
    p.add_argument(
        "--structures-dir",
        required=True,
        help="step folder with a chains table (pandora curate)",
    )
    p.add_argument(
        "--input-dir", help="foldseek only: dir of per-chain structure files"
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--mmseqs-bin", default="mmseqs")
    p.add_argument("--foldseek-bin", default="foldseek")
    p.add_argument("--sensitivity", type=float, default=None)
    p.add_argument("--alignment-type", type=int, default=2)
    p.add_argument(
        "--max-seqs",
        type=int,
        default=None,
        help="hits kept per query (default: the tool's own, MMseqs2 300 "
        "/ Foldseek 1000; raise it for leakage control)",
    )
    p.add_argument(
        "--exhaustive-search",
        action="store_true",
        help="foldseek only: skip the prefilter",
    )
    p.set_defaults(func=_cmd_search)

    p = commands.add_parser(
        "import",
        help="Link an edge list (Parquet, or TSV with a header) to a step "
        "folder's entries/chains or a contact layer's interfaces.",
    )
    p.add_argument("--table", required=True)
    p.add_argument("--description", required=True)
    p.add_argument(
        "--structures-dir", help="step folder, for entry or chain items"
    )
    p.add_argument(
        "--contacts", help="residue_contacts layer folder, for interface items"
    )
    p.add_argument("--source-files")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=_cmd_import)

    p = subparsers.add_parser(
        "cluster",
        help="Connected-component clustering of a similarity layer.",
    )
    p.add_argument(
        "--similarity",
        required=True,
        help="a similarity layer folder (similarity search/import)",
    )
    p.add_argument(
        "--items",
        choices=["entry", "chain", "interface"],
        help="cluster at this level (default: the layer's); entry on a "
        "chain layer rolls chains up to their entries",
    )
    p.add_argument("--hit-filter", help="HitFilter YAML (thresholds)")
    p.add_argument(
        "--min-score",
        type=float,
        default=None,
        help="overrides --hit-filter's min_score",
    )
    p.add_argument(
        "--min-coverage",
        type=float,
        default=None,
        help="overrides --hit-filter's min_coverage",
    )
    p.add_argument(
        "--contacts",
        help=(
            "a residue_contacts layer folder: the items of an interface "
            "layer, or interface residues for min_interface_coverage on a "
            "chain layer"
        ),
    )
    p.add_argument(
        "--structures-dir",
        help="the step folder whose entries/chains are the items",
    )
    p.add_argument(
        "--pairs",
        help=(
            "JSON file of [item_id_1, item_id_2] PPI pairs to also derive "
            "per-side paired cluster keys for -> cluster_pairs.json"
        ),
    )
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_cluster)

    p = subparsers.add_parser(
        "partition", help="Leakage-safe train/val/test split of clusters."
    )
    p.add_argument("--clusters", required=True)
    p.add_argument("--pct-train", type=float, default=0.6)
    p.add_argument("--pct-val", type=float, default=0.2)
    p.add_argument("--pct-test", type=float, default=0.2)
    p.add_argument("--no-keep-similar-items", action="store_true")
    p.add_argument(
        "--labels",
        dest="labels_dir",
        help="chain_labels or interface_labels layer folder: place clusters "
        "so every label reaches every split",
    )
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_partition)
