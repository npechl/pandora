"""CLI `rebalance` and `negatives`: pick the positive dimers of each
dimer group, and sample negative chain pairs for them (PPI v2.1 Steps 7
and 8)."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from typing import Any

from pandora import __version__
from pandora._util import now_iso, rank_sort_key
from pandora.cli._io import (
    TABLE_FORMATS,
    check_table_format,
    find_table,
    load_json_list,
    load_yaml_model,
    prepare_output_dir,
    read_columns,
    write_table,
)
from pandora.cli.similarity import (
    _CONTACT_COLUMNS,
    _chain_map,
    _contact_layer,
    _hit_filters,
    _interface_chains,
    _interface_residues,
    _interface_values,
    _item_set,
    _layer_edges,
    _require_pyarrow,
    _similarity_level,
)
from pandora.export import read_records, write_json, write_records
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import LayerDescription, LayerRow
from pandora.schemas.similarity import (
    ClusteredLayer,
    HitFilter,
    NegativeSamplingSettings,
    PositiveDimer,
    RebalanceSettings,
    RepresentativeCluster,
    RepresentativePolicy,
    SimilarityCluster,
)
from pandora.similarity import (
    dimer_adjacency,
    item_id,
    member_representatives,
    parse_item_id,
    rebalance,
    restrict_clusters,
    sample_negatives,
)

_REBALANCE_OUTPUTS = [
    f"{name}.{fmt}"
    for name in ("positives", "dropped")
    for fmt in TABLE_FORMATS
] + ["rebalance.json", "groups.json", "kinds"]
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
    values, _, origin = _interface_values(contacts, command)
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

    residues = None
    if any(
        level == "chain" and f.min_interface_coverage is not None
        for level, f in zip(levels, filters, strict=True)
    ):
        residues = _interface_residues(origin, contacts)
    counts: Counter[int] = Counter()
    chain_edges: list[Any] = []
    dimer_edges: list[Any] = []
    for index, (layer, level, hit_filter) in enumerate(
        zip(layers, levels, filters, strict=True)
    ):
        if level not in ("chain", "interface"):
            raise SystemExit(f"{command}: {layer} is a {level} layer")
        target = chain_edges if level == "chain" else dimer_edges
        scoped = residues if level == "chain" else None
        target.extend(
            _layer_edges(layer, hit_filter, scoped, counts, index, command)
        )
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
    _write_groups_and_kinds(
        groups, positives, out, args.table_format, started_at
    )
    print(f"{len(positives)} positives from {len(dimers)} dimers -> {out}")


def _write_groups_and_kinds(
    groups: list[SimilarityCluster],
    positives: list[PositiveDimer],
    out: Path,
    table_format: str,
    started_at: str,
) -> None:
    """Write `groups.json` (the groups limited to the positives) and
    `kinds/`, an interface_labels layer of homomer/heteromer."""

    kept = restrict_clusters(groups, {p.dimer for p in positives})
    write_records(kept, out / "groups.json")
    kinds = out / "kinds"
    kinds.mkdir()
    rows = [
        LayerRow(
            entry_id=p.entry_id,
            chain_id_1=p.chain_id_1,
            chain_id_2=p.chain_id_2,
            label=p.kind,
        )
        for p in positives
    ]
    write_table(rows, kinds, "layer", table_format)
    description = LayerDescription(
        name="positive_kinds",
        template="interface_labels",
        tool="pandora",
        tool_version=__version__,
        parameters={},
    )
    labels = Counter(p.kind for p in positives)
    write_json(
        build_step_summary(
            "rebalance_kinds",
            description,
            input=out / "positives",
            output=kinds,
            started_at=started_at,
            counts={"interfaces": len(rows), **dict(sorted(labels.items()))},
            previous_summary=out / "rebalance.json",
        ),
        kinds / "layer.json",
    )


_NEGATIVES_OUTPUTS = [f"negatives.{fmt}" for fmt in TABLE_FORMATS] + [
    "negatives.json"
]


def _positives(command: str, folder: Path) -> tuple[list[PositiveDimer], Path]:
    """The positives of a `rebalance` step folder and its summary."""

    summary = folder / "rebalance.json"
    if not summary.exists():
        raise SystemExit(
            f"{command}: {folder} has no rebalance.json; it is written by "
            "`pandora rebalance`"
        )
    table = find_table(folder, "positives")
    return (read_records(table, PositiveDimer) if table else []), summary


def _binding_pairs(command: str, contacts: Path) -> list[tuple[str, str]]:
    """The chain item pairs in contact in a residue_contacts layer."""

    rows_file = _contact_layer(command, contacts)
    if rows_file.name == "layer.json":
        return []
    return sorted(
        {
            (
                item_id("chain", r["entry_id"], r["chain_id_1"]),
                item_id("chain", r["entry_id"], r["chain_id_2"]),
            )
            for r in read_columns(rows_file, _CONTACT_COLUMNS)
        }
    )


def _cmd_negatives(args: argparse.Namespace) -> None:
    """Handle `negatives`: draw one negative chain pair per positive and
    write them as a step folder."""

    command = "negatives"
    started_at = now_iso()
    check_table_format(args.table_format)
    layers, filters, levels = _layers(args, command)
    if len(layers) != 1 or levels[0] != "chain":
        raise SystemExit(
            f"{command}: --similarity takes one chain layer (structure "
            "similarity, e.g. Foldseek)"
        )
    positives, previous = _positives(command, Path(args.positives))
    representative_of = member_representatives(
        load_json_list(RepresentativeCluster, Path(args.representatives))
    )
    counts: Counter[int] = Counter()
    tm_edges = list(
        _layer_edges(layers[0], filters[0], None, counts, 0, command)
    )
    binding = _binding_pairs(command, Path(args.contacts))
    try:
        negatives, unpaired = sample_negatives(
            positives, representative_of, tm_edges, binding, args.seed
        )
    except ValueError as exc:
        raise SystemExit(f"{command}: {exc} in {args.representatives}") from exc

    out = Path(args.output_dir)
    prepare_output_dir(out, _NEGATIVES_OUTPUTS, args.overwrite)
    write_table(negatives, out, "negatives", args.table_format)
    settings = NegativeSamplingSettings(
        seed=args.seed,
        hit_filter=filters[0],
        similarity_layer=str(layers[0] / "layer.json"),
        contacts=str(args.contacts),
        representatives=str(args.representatives),
        unpaired=unpaired,
    )
    pool = set(representative_of.values()) - {
        representative_of[item_id("chain", p.entry_id, chain)]
        for p in positives
        for chain in (p.chain_id_1, p.chain_id_2)
    }
    write_json(
        build_step_summary(
            command,
            settings,
            input=args.positives,
            output=out,
            started_at=started_at,
            counts={
                "pool": len(pool),
                "positives": len(positives),
                "negatives": len(negatives),
                "unpaired": len(unpaired),
                "similarity_edges": counts[0],
            },
            previous_summary=previous,
        ),
        out / "negatives.json",
    )
    print(
        f"{len(negatives)} negatives for {len(positives)} positives "
        f"({len(unpaired)} unpaired) -> {out}"
    )


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

    p = subparsers.add_parser(
        "negatives",
        help="Draw one negative chain pair per positive dimer, avoiding "
        "partners similar to anything the chain might bind.",
    )
    p.add_argument(
        "--positives", required=True, help="pandora rebalance step folder"
    )
    p.add_argument(
        "--representatives",
        required=True,
        help="representatives.json (90%% chains): the pool, and each "
        "chain's representative",
    )
    _layer_flags(p, "the chain structure-similarity layer (one)")
    p.add_argument(
        "--contacts",
        required=True,
        help="the 'any contact' residue_contacts layer folder",
    )
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--table-format", choices=TABLE_FORMATS, default="parquet")
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=_cmd_negatives)
