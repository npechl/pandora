"""CLI `rebalance` and `negatives`: pick the positive dimers of each
dimer group, and sample negative chain pairs for them (PPI v2.1 Steps 7
and 8)."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import Any

from pandora._util import now_iso, rank_sort_key
from pandora.cli._io import (
    TABLE_FORMATS,
    check_table_format,
    load_json_list,
    load_yaml_model,
    prepare_output_dir,
    write_table,
)
from pandora.cli.similarity import (
    _chain_map,
    _hit_filters,
    _interface_chains,
    _interface_values,
    _item_set,
    _layer_edges,
    _require_pyarrow,
    _similarity_level,
)
from pandora.export import write_json
from pandora.provenance import build_step_summary
from pandora.schemas.similarity import (
    ClusteredLayer,
    HitFilter,
    RebalanceSettings,
    RepresentativePolicy,
    SimilarityCluster,
)
from pandora.similarity import (
    dimer_adjacency,
    item_id,
    parse_item_id,
    rebalance,
)

_REBALANCE_OUTPUTS = [
    f"{name}.{fmt}"
    for name in ("positives", "dropped")
    for fmt in TABLE_FORMATS
] + ["rebalance.json"]
_REBALANCE_RANK_KEYS = {"interface_residues", "resolution", "missing_residues"}


def _layer_flags(p: argparse.ArgumentParser, help_layer: str) -> None:
    """Add --similarity and its per-layer hit filter flags."""

    p.add_argument(
        "--similarity", required=True, action="append", help=help_layer
    )
    p.add_argument(
        "--hit-filter",
        action="append",
        help="HitFilter YAML; once for every layer or once per "
        "--similarity, in order (same for the two flags below)",
    )
    p.add_argument(
        "--min-score",
        type=float,
        action="append",
        help="overrides --hit-filter's min_score",
    )
    p.add_argument(
        "--min-coverage",
        type=float,
        action="append",
        help="overrides --hit-filter's min_coverage",
    )


def _layers(
    args: argparse.Namespace, command: str
) -> tuple[list[Path], list[HitFilter], list[str]]:
    """The --similarity layers, their hit filters and item levels."""

    _require_pyarrow(command)
    layers = [Path(folder) for folder in args.similarity]
    filters = _hit_filters(args, len(layers), command)
    levels = [_similarity_level(layer) for layer in layers]
    return layers, filters, levels


def _homomers(
    path: Path, dimers: list[str], command: str = "rebalance"
) -> set[str]:
    """The dimers whose two chains share a cluster of `path`; exits
    naming a chain in no cluster."""

    cluster_of = {
        member: index
        for index, cluster in enumerate(load_json_list(SimilarityCluster, path))
        for member in cluster.components
    }
    missing = next(
        (c for c in sorted(_interface_chains(dimers)) if c not in cluster_of),
        None,
    )
    if missing is not None:
        raise SystemExit(
            f"{command}: chain {missing} is in no cluster of {path}"
        )
    homomers = set()
    for dimer in dimers:
        entry, chain_1, chain_2 = parse_item_id("interface", dimer)
        if (
            cluster_of[item_id("chain", entry, chain_1)]
            == cluster_of[item_id("chain", entry, chain_2)]
        ):
            homomers.add(dimer)
    return homomers


def _rank_policy(args: argparse.Namespace) -> RepresentativePolicy:
    """The tie-break ranking: --policy, or the interface default."""

    loaded = (
        load_yaml_model(RepresentativePolicy, Path(args.policy))
        if args.policy
        else RepresentativePolicy()
    )
    policy = RepresentativePolicy.model_validate(
        loaded.model_dump(exclude_unset=True) | {"items": "interface"}
    )
    unknown = sorted({k.key for k in policy.rank} - _REBALANCE_RANK_KEYS)
    if unknown:
        raise SystemExit(
            f"rebalance: unknown rank key {unknown[0]!r}; use "
            + ", ".join(sorted(_REBALANCE_RANK_KEYS))
        )
    return policy


def _cmd_rebalance(args: argparse.Namespace) -> None:
    """Handle `rebalance`: keep the positive dimers of each dimer group
    and write them as a step folder."""

    command = "rebalance"
    started_at = now_iso()
    check_table_format(args.table_format)
    layers, filters, levels = _layers(args, command)
    policy = _rank_policy(args)
    contacts = Path(args.contacts)
    items, previous = _item_set(command, "interface", None, contacts)
    known = set(items.ids.values())
    groups = load_json_list(SimilarityCluster, Path(args.clusters))
    dimers = sorted(d for g in groups for d in g.components)
    missing = next((d for d in dimers if d not in known), None)
    if missing is not None:
        raise SystemExit(f"{command}: dimer {missing} is not in {contacts}")
    values, _, _ = _interface_values(contacts, command)
    homomers = _homomers(Path(args.chain_clusters), dimers)
    chain_map = (
        _chain_map(
            command,
            Path(args.chain_representatives),
            _interface_chains(dimers),
        )
        if args.chain_representatives
        else None
    )

    counts: Counter[int] = Counter()
    chain_edges: list[Any] = []
    dimer_edges: list[Any] = []
    for index, (layer, level, hit_filter) in enumerate(
        zip(layers, levels, filters, strict=True)
    ):
        if level not in ("chain", "interface"):
            raise SystemExit(f"{command}: {layer} is a {level} layer")
        target = chain_edges if level == "chain" else dimer_edges
        target.extend(_layer_edges(layer, hit_filter, None, counts, index))
    adjacency = dimer_adjacency(dimers, chain_edges, dimer_edges, chain_map)
    try:
        positives, dropped = rebalance(
            groups,
            adjacency,
            homomers,
            lambda d: rank_sort_key(d, values[d], policy.rank),
        )
    except ValueError as exc:
        raise SystemExit(f"{command}: {exc}") from exc

    out = Path(args.output_dir)
    prepare_output_dir(out, _REBALANCE_OUTPUTS, args.overwrite)
    write_table(positives, out, "positives", args.table_format)
    write_table(dropped, out, "dropped", args.table_format)
    settings = RebalanceSettings(
        rank=policy,
        layers=[
            ClusteredLayer(
                similarity_layer=str(layer / "layer.json"),
                items=level,
                hit_filter=hit_filter,
                n_edges=counts[index],
            )
            for index, (layer, level, hit_filter) in enumerate(
                zip(layers, levels, filters, strict=True)
            )
        ],
        clusters=str(args.clusters),
        chain_clusters=str(args.chain_clusters),
        chain_representatives=args.chain_representatives,
    )
    tallies = Counter(f"positives:{p.kind}" for p in positives)
    tallies.update(f"dropped:{d.reason}" for d in dropped)
    write_json(
        build_step_summary(
            command,
            settings,
            input=args.clusters,
            output=out,
            started_at=started_at,
            counts={
                "groups": len(groups),
                "dimers_in": len(dimers),
                **dict(sorted(tallies.items())),
            },
            previous_summary=previous,
        ),
        out / "rebalance.json",
    )
    print(f"{len(positives)} positives from {len(dimers)} dimers -> {out}")


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    p = subparsers.add_parser(
        "rebalance",
        help="Pick the positive dimers of each dimer group: one homomer "
        "per homomer-only group, else mutually dissimilar heteromers.",
    )
    p.add_argument(
        "--clusters", required=True, help="dimer groups (pandora cluster)"
    )
    p.add_argument(
        "--contacts",
        required=True,
        help="the residue_contacts layer folder the groups were built from",
    )
    p.add_argument(
        "--chain-clusters",
        required=True,
        help="90%% chain clusters (pandora cluster): a dimer is a homomer "
        "when both chains share a cluster",
    )
    p.add_argument(
        "--chain-representatives",
        help="representatives.json: map chains to their representative "
        "(as for the Step 5 cluster run)",
    )
    _layer_flags(p, "similarity layer folder of the dimer graph; repeatable")
    p.add_argument(
        "--policy",
        help="RepresentativePolicy YAML (items: interface): tie-break rank",
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument("--table-format", choices=TABLE_FORMATS, default="parquet")
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=_cmd_rebalance)
