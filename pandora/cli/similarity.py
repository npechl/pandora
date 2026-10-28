"""CLI `similarity`, `cluster` and `partition`: searches, clusters and leakage-safe splits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pandora.cli._io import (
    load_interface_residues,
    load_json_list,
    load_json_model,
    load_pairs,
    load_yaml_model,
    read_structures_dir,
)
from pandora.datasets import (
    entry_sequences,
)
from pandora.export import write_json, write_records
from pandora.schemas.similarity import (
    HitFilter,
    SimilarityCluster,
    SimilaritySearch,
)
from pandora.similarity import (
    cluster_similar_items,
    compute_sequence_similarity,
    compute_structure_similarity,
    load_similarity_search,
    pair_cluster_keys,
    partition_dataset,
)

_CLI_ENGINES = {"mmseqs2": "MMseqs2", "foldseek": "Foldseek"}


def _cmd_similarity(args: argparse.Namespace) -> None:
    """Handle the `similarity` subcommand: run (or wrap a precomputed)
    all-vs-all sequence or structure search, keeping its hit file."""

    output = Path(args.output)

    if args.precomputed_hits:
        ignored = [
            flag
            for flag, value in (
                ("--max-seqs", args.max_seqs),
                ("--sensitivity", args.sensitivity),
            )
            if value is not None
        ]
        if ignored:
            raise SystemExit(
                f"similarity: {', '.join(ignored)} has no effect with "
                "--precomputed-hits (no search is run)"
            )
        search = load_similarity_search(
            args.precomputed_hits, _CLI_ENGINES[args.engine]
        )
    elif args.input_dir is None:
        raise SystemExit(
            "similarity: --input-dir is required unless --precomputed-hits "
            "is given"
        )
    elif args.engine == "mmseqs2":
        structures = read_structures_dir(Path(args.input_dir))
        sequences = entry_sequences(structures)
        kwargs = {"mmseqs_bin": args.mmseqs_bin}
        if args.sensitivity is not None:
            kwargs["sensitivity"] = args.sensitivity
        if args.max_seqs is not None:
            kwargs["max_seqs"] = args.max_seqs
        search = compute_sequence_similarity(sequences, output, **kwargs)
    else:
        kwargs = {
            "foldseek_bin": args.foldseek_bin,
            "alignment_type": args.alignment_type,
            "exhaustive_search": args.exhaustive_search,
        }
        if args.sensitivity is not None:
            kwargs["sensitivity"] = args.sensitivity
        if args.max_seqs is not None:
            kwargs["max_seqs"] = args.max_seqs
        search = compute_structure_similarity(
            Path(args.input_dir), output, **kwargs
        )

    record = output.with_name(output.name + ".search.json")
    write_json(search, record)
    print(f"hits -> {search.hits_path}; search record -> {record}")


def _cmd_cluster(args: argparse.Namespace) -> None:
    """Handle the `cluster` subcommand: connected-component clustering
    of a relationship network."""

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
    input_dir = Path(args.input_dir)
    search = load_json_model(SimilaritySearch, Path(args.search))
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
    if overrides:
        hit_filter = hit_filter.model_copy(update=overrides)
    interface_residues = (
        load_interface_residues(Path(args.interface_residues))
        if args.interface_residues
        else None
    )
    item_ids = sorted(path.stem.upper() for path in input_dir.glob("*.cif"))

    clusters, prov = cluster_similar_items(
        item_ids, search, hit_filter, interface_residues=interface_residues
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


def _cmd_partition(args: argparse.Namespace) -> None:
    """Handle the `partition` subcommand: leakage-safe train/val/test
    split of clusters."""

    clusters = load_json_list(SimilarityCluster, Path(args.clusters))
    splits, prov = partition_dataset(
        clusters,
        args.pct_train,
        args.pct_val,
        args.pct_test,
        keep_similar_items=not args.no_keep_similar_items,
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(splits, indent=2), encoding="utf-8")
    write_json(prov, output.parent / "partition_provenance.json")
    print(f"split sizes: {prov.split_sizes} -> {output}")


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    p = subparsers.add_parser(
        "similarity", help="All-vs-all sequence or structure similarity."
    )
    p.add_argument("--engine", choices=["mmseqs2", "foldseek"], required=True)
    p.add_argument("--mmseqs-bin", default="mmseqs")
    p.add_argument("--foldseek-bin", default="foldseek")
    p.add_argument("--sensitivity", type=float, default=None)
    p.add_argument("--alignment-type", type=int, default=2)
    p.add_argument(
        "--input-dir", help="required unless --precomputed-hits is given"
    )
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
    p.add_argument(
        "--precomputed-hits",
        help="wrap an existing hit TSV (Pandora's --format-output "
        "columns) instead of running a search",
    )
    p.add_argument(
        "--output",
        required=True,
        help="hit TSV path; the search record goes to <output>.search.json",
    )
    p.set_defaults(func=_cmd_similarity)

    p = subparsers.add_parser(
        "cluster",
        help="Connected-component clustering of a relationship network.",
    )
    p.add_argument(
        "--input-dir", required=True, help="dir of the clustered *.cif"
    )
    p.add_argument(
        "--search",
        required=True,
        help="<hits>.search.json written by `pandora similarity`",
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
        "--interface-residues",
        help=(
            "JSON file of {item_id: [residue position, ...]} for "
            "min_interface_coverage (Foldseek searches only)"
        ),
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
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_partition)
