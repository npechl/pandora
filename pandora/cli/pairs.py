"""CLI `pairs positives`, `pairs negatives` and `pairs assemble`: pick the
positive dimers of each dimer group, sample negative chain pairs for
them, and write the final pairs table (PPI v2.1 Steps 7 and 8)."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from pandora import __version__
from pandora._util import now_iso, rank_sort_key
from pandora.cli._io import (
    SUMMARY,
    TABLE_FORMATS,
    SimilarityInput,
    add_out,
    add_similarity_inputs,
    check_out,
    check_table_format,
    contact_table,
    input_summary,
    layer_edges,
    load_json,
    load_json_list,
    prepare_output_dir,
    read_columns,
    read_summary,
    require_tables,
    resolve_file,
    resolve_table,
    similarity_inputs,
    write_json_list,
    write_summary,
    write_table,
)
from pandora.cli.clusters import (
    _chain_map,
    _interface_chains,
    _interface_residues,
    _interface_values,
    _rank_policy,
)
from pandora.cli.similarity import (
    _CONTACT_COLUMNS,
    _item_set,
    _require_pyarrow,
)
from pandora.export import read_records
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import LayerDescription, LayerRow
from pandora.schemas.similarity import (
    AssemblySettings,
    ClusteredLayer,
    NegativePair,
    NegativeSamplingSettings,
    PositiveDimer,
    RebalanceSettings,
    RepresentativeCluster,
    RepresentativePolicy,
    SimilarityCluster,
)
from pandora.similarity import (
    assemble_pairs,
    cluster_edges,
    dimer_adjacency,
    item_id,
    member_representatives,
    parse_item_id,
    rebalance,
    restrict_clusters,
    sample_negatives,
)

_POSITIVES_OUTPUTS = [
    f"{name}.{fmt}" for name in ("positives", "dropped") for fmt in TABLE_FORMATS
] + ["clusters.json", "kinds", SUMMARY]


def _positives_policy(command: str, args: argparse.Namespace) -> RepresentativePolicy:
    """The tie-break ranking: --policy, or the interface default; n_hits
    does not apply."""

    policy = _rank_policy(command, args, "interface")
    if any(key.key == "n_hits" for key in policy.rank):
        raise SystemExit(
            f"{command}: rank key 'n_hits' does not apply; use "
            "interface_residues, missing_residues, resolution"
        )
    return policy


def _homomers(command: str, path: Path, dimers: list[str]) -> set[str]:
    """The dimers whose two chains share a cluster of `path`; exits
    naming a chain in no cluster."""

    cluster_of = {
        member: index
        for index, cluster in enumerate(load_json_list(command, SimilarityCluster, path))
        for member in cluster.components
    }
    missing = next(
        (c for c in sorted(_interface_chains(dimers)) if c not in cluster_of),
        None,
    )
    if missing is not None:
        raise SystemExit(f"{command}: chain {missing} is in no cluster of {path}")
    homomers = set()
    for dimer in dimers:
        entry, chain_1, chain_2 = parse_item_id("interface", dimer)
        if (
            cluster_of[item_id("chain", entry, chain_1)]
            == cluster_of[item_id("chain", entry, chain_2)]
        ):
            homomers.add(dimer)
    return homomers


def _dimer_edges(
    command: str,
    inputs: list[SimilarityInput],
    structures: Path,
    contacts: Path,
    counts: Counter[int],
) -> tuple[list[Any], list[Any]]:
    """The chain and interface edges of the similarity inputs; exits on
    an input at another level."""

    residues = None
    if any(i.items == "chain" and i.hit_filter.min_interface_coverage is not None for i in inputs):
        residues = _interface_residues(command, structures, contacts)
    chain_edges: list[Any] = []
    dimer_edges: list[Any] = []
    for index, source in enumerate(inputs):
        if source.items not in ("chain", "interface"):
            raise SystemExit(
                f"{command}: {source.path} is a {source.items} layer; the "
                "dimer graph takes chain and interface layers"
            )
        chain = source.items == "chain"
        (chain_edges if chain else dimer_edges).extend(
            layer_edges(
                source.table or source.path,
                source.hit_filter,
                residues if chain else None,
                counts,
                index,
                command,
            )
        )
    return chain_edges, dimer_edges


def _cmd_positives(args: argparse.Namespace) -> None:
    """Handle `pairs positives`: group the dimers by the similarity
    inputs, keep the positives of each group, and write them."""

    command = "pairs positives"
    started_at = now_iso()
    out, structures = Path(args.out), Path(args.structures)
    contacts = Path(args.contacts)
    reps = args.chain_representatives
    check_out(
        command,
        out,
        *(Path(p) for p in args.similarity),
        contacts,
        structures,
        Path(args.chain_clusters),
        Path(reps) if reps else None,
    )
    check_table_format(command, args.table_format)
    _require_pyarrow(command)
    inputs = similarity_inputs(command, args)
    policy = _positives_policy(command, args)
    _, chains_table = require_tables(command, structures)
    items, _ = _item_set(command, "interface", None, contacts)
    dimers = sorted(set(items.ids.values()))
    values, _ = _interface_values(command, contacts, chains_table)
    chain_clusters = resolve_file(command, Path(args.chain_clusters), "clusters.json")
    homomers = _homomers(command, chain_clusters, dimers)
    reps_file = (
        resolve_file(command, Path(args.chain_representatives), "representatives.json")
        if args.chain_representatives
        else None
    )
    chain_map = _chain_map(command, reps_file, _interface_chains(dimers)) if reps_file else None
    counts: Counter[int] = Counter()
    chain_edges, dimer_edges = _dimer_edges(command, inputs, structures, contacts, counts)
    try:
        adjacency = dimer_adjacency(dimers, chain_edges, dimer_edges, chain_map)
        groups, _ = cluster_edges(dimers, ((a, b) for a, ns in adjacency.items() for b in ns))
        positives, dropped = rebalance(
            groups,
            adjacency,
            homomers,
            lambda d: rank_sort_key(d, values[d], policy.rank),
        )
    except ValueError as exc:
        raise SystemExit(f"{command}: {exc}") from exc

    prepare_output_dir(command, out, _POSITIVES_OUTPUTS, args.overwrite)
    write_table(positives, out, "positives", args.table_format)
    write_table(dropped, out, "dropped", args.table_format)
    settings = RebalanceSettings(
        rank=policy,
        layers=[
            ClusteredLayer(
                similarity_layer=str(source.path.resolve()),
                items=source.items,
                hit_filter=source.hit_filter,
                n_edges=counts[index],
            )
            for index, source in enumerate(inputs)
        ],
        contacts=str(contacts.resolve()),
        structures_dir=str(structures.resolve()),
        chain_clusters=str(chain_clusters.resolve()),
        chain_representatives=None if reps_file is None else str(reps_file.resolve()),
    )
    tallies = Counter(f"positives:{p.kind}" for p in positives)
    tallies.update(f"dropped:{d.reason}" for d in dropped)
    first = inputs[0].path
    write_summary(
        out,
        build_step_summary(
            "rebalance",
            settings,
            input=first.resolve(),
            output=out,
            started_at=started_at,
            counts={
                "groups": len(groups),
                "dimers_in": len(dimers),
                **dict(sorted(tallies.items())),
            },
            previous_summary=input_summary(first) if first.is_dir() else None,
        ),
    )
    _write_clusters_and_kinds(groups, positives, out, args.table_format, started_at)
    print(f"{len(positives)} positives from {len(dimers)} dimers -> {out}")


def _write_clusters_and_kinds(
    groups: list[SimilarityCluster],
    positives: list[PositiveDimer],
    out: Path,
    table_format: str,
    started_at: str,
) -> None:
    """Write `clusters.json` (the groups limited to the positives) and
    `kinds/`, an interface_labels layer of homomer/heteromer."""

    kept = restrict_clusters(groups, {p.dimer for p in positives})
    write_json_list(kept, out / "clusters.json")
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
    write_summary(
        kinds,
        build_step_summary(
            "rebalance_kinds",
            description,
            input=(out / f"positives.{table_format}").resolve(),
            output=kinds,
            started_at=started_at,
            counts={"interfaces": len(rows), **dict(sorted(labels.items()))},
            previous_summary=out / SUMMARY,
        ),
    )


def _step_rows(
    command: str, path: Path, name: str, step: str, model: type
) -> tuple[list, Path, Path]:
    """The rows of a `name` table given as a folder or the table itself,
    that table (or the folder, when `step` wrote no rows) and its
    folder."""

    if path.is_dir():
        summary = read_summary(path)
        empty = (summary or {}).get("step") == step
        if empty and not any((path / f"{name}.{fmt}").exists() for fmt in TABLE_FORMATS):
            return [], path, path
    table = resolve_table(command, path, name)
    return read_records(table, model), table, table.parent


def _one_chain_input(command: str, args: argparse.Namespace, what: str) -> SimilarityInput:
    """The one chain similarity input of `args`; exits otherwise."""

    _require_pyarrow(command)
    inputs = similarity_inputs(command, args)
    if len(inputs) != 1 or inputs[0].items != "chain":
        raise SystemExit(f"{command}: --similarity takes one chain layer ({what})")
    return inputs[0]


def _representatives(command: str, path: str) -> tuple[dict[str, str], Path]:
    """Chain -> representative from a representatives folder or file, and
    the file."""

    file = resolve_file(command, Path(path), "representatives.json")
    clusters = load_json_list(command, RepresentativeCluster, file)
    return member_representatives(clusters), file


_NEGATIVES_OUTPUTS = [f"negatives.{fmt}" for fmt in TABLE_FORMATS] + [SUMMARY]


def _binding_pairs(command: str, contacts: Path) -> list[tuple[str, str]]:
    """The chain item pairs in contact in a residue_contacts layer."""

    rows_file = contact_table(command, contacts)
    if rows_file is None:
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
    """Handle `pairs negatives`: draw one negative chain pair per
    positive and write them."""

    command = "pairs negatives"
    started_at = now_iso()
    out, contacts = Path(args.out), Path(args.contacts)
    check_out(
        command,
        out,
        Path(args.positives),
        Path(args.representatives),
        contacts,
        *(Path(p) for p in args.similarity or []),
    )
    check_table_format(command, args.table_format)
    source = _one_chain_input(command, args, "structure similarity, e.g. Foldseek")
    positives, pos_path, pos_folder = _step_rows(
        command, Path(args.positives), "positives", "rebalance", PositiveDimer
    )
    representative_of, reps_file = _representatives(command, args.representatives)
    counts: Counter[int] = Counter()
    tm_edges = list(
        layer_edges(
            source.table or source.path,
            source.hit_filter,
            None,
            counts,
            0,
            command,
        )
    )
    binding = _binding_pairs(command, contacts)
    try:
        negatives, unpaired = sample_negatives(
            positives, representative_of, tm_edges, binding, args.seed
        )
    except ValueError as exc:
        raise SystemExit(f"{command}: {exc} in {reps_file}") from exc

    prepare_output_dir(command, out, _NEGATIVES_OUTPUTS, args.overwrite)
    write_table(negatives, out, "negatives", args.table_format)
    settings = NegativeSamplingSettings(
        seed=args.seed,
        hit_filter=source.hit_filter,
        similarity_layer=str(source.path.resolve()),
        contacts=str(contacts.resolve()),
        representatives=str(reps_file.resolve()),
        unpaired=unpaired,
    )
    pool = set(representative_of.values()) - {
        representative_of[item_id("chain", p.entry_id, chain)]
        for p in positives
        for chain in (p.chain_id_1, p.chain_id_2)
    }
    write_summary(
        out,
        build_step_summary(
            "negatives",
            settings,
            input=pos_path.resolve(),
            output=out,
            started_at=started_at,
            counts={
                "pool": len(pool),
                "positives": len(positives),
                "negatives": len(negatives),
                "unpaired": len(unpaired),
                "similarity_edges": counts[0],
            },
            previous_summary=input_summary(pos_folder),
        ),
    )
    print(
        f"{len(negatives)} negatives for {len(positives)} positives "
        f"({len(unpaired)} unpaired) -> {out}"
    )


_ASSEMBLE_OUTPUTS = [f"pairs.{fmt}" for fmt in TABLE_FORMATS] + [SUMMARY]


def _same_path(recorded: str | None, given: Path) -> bool:
    """Whether a path recorded in a summary names the same file as
    `given`."""

    return recorded is not None and Path(recorded).resolve() == given.resolve()


def _layer_file(path: Path) -> Path:
    """A similarity input as its table: a folder's `layer.parquet`, or
    the bare table itself."""

    return path / "layer.parquet" if path.is_dir() else path


def _warn_other_inputs(command: str, folder: Path, reps_file: Path, similarity: Path) -> None:
    """Warn when the representatives or similarity input differ from
    what the `pairs negatives` run in `folder` recorded; skip without
    its summary."""

    summary = read_summary(folder)
    if summary is None or summary.get("step") != "negatives":
        return
    used = summary.get("policy") or {}
    recorded = used.get("similarity_layer")
    if not _same_path(used.get("representatives"), reps_file) or not _same_path(
        recorded and str(_layer_file(Path(recorded))), _layer_file(similarity)
    ):
        print(
            f"{command}: warning: --representatives or --similarity differs "
            f"from what {input_summary(folder)} used; leak counts are not "
            "comparable",
            file=sys.stderr,
        )


def _cmd_assemble(args: argparse.Namespace) -> None:
    """Handle `pairs assemble`: write the final pairs table with splits
    and C-chain leak flags."""

    command = "pairs assemble"
    started_at = now_iso()
    out = Path(args.out)
    check_out(
        command,
        out,
        *(
            Path(p)
            for p in (
                args.positives,
                args.negatives,
                args.splits,
                args.representatives,
                *(args.similarity or []),
            )
        ),
    )
    check_table_format(command, args.table_format)
    source = _one_chain_input(command, args, "the one `pairs negatives` used")
    positives, pos_path, pos_folder = _step_rows(
        command, Path(args.positives), "positives", "rebalance", PositiveDimer
    )
    negatives, neg_path, neg_folder = _step_rows(
        command, Path(args.negatives), "negatives", "negatives", NegativePair
    )
    representative_of, reps_file = _representatives(command, args.representatives)
    _warn_other_inputs(command, neg_folder, reps_file, source.path)
    splits_file = resolve_file(command, Path(args.splits), "splits.json")
    splits = load_json(command, splits_file)
    counts: Counter[int] = Counter()
    tm_edges = list(
        layer_edges(
            source.table or source.path,
            source.hit_filter,
            None,
            counts,
            0,
            command,
        )
    )
    try:
        rows = assemble_pairs(positives, negatives, splits, representative_of, tm_edges)
    except ValueError as exc:
        raise SystemExit(f"{command}: {exc}") from exc

    prepare_output_dir(command, out, _ASSEMBLE_OUTPUTS, args.overwrite)
    write_table(rows, out, "pairs", args.table_format)
    tallies = Counter({f"pairs:{n}:{lab}": 0 for n in splits for lab in (0, 1)})
    tallies.update(f"pairs:{r.split}:{r.label}" for r in rows)
    paired = {r.positive for r in rows if r.label == 0}
    settings = AssemblySettings(
        positives=str(pos_path.resolve()),
        negatives=str(neg_path.resolve()),
        splits=str(splits_file.resolve()),
        representatives=str(reps_file.resolve()),
        similarity_layer=str(source.path.resolve()),
        hit_filter=source.hit_filter,
    )
    write_summary(
        out,
        build_step_summary(
            "assemble",
            settings,
            input=pos_path.resolve(),
            output=out,
            started_at=started_at,
            counts={
                **dict(sorted(tallies.items())),
                "negatives_c_leak_positive": sum(bool(r.c_leak_positive) for r in rows),
                "negatives_c_leak_negative": sum(bool(r.c_leak_negative) for r in rows),
                "positives_unpaired": len(positives) - len(paired),
            },
            previous_summary=input_summary(pos_folder),
        ),
    )
    print(f"{len(rows)} pairs -> {out}")


_POSITIVES_HELP = "positives folder or table (pandora pairs positives)"


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add the `pairs` group to the `pandora` parser."""

    group = subparsers.add_parser(
        "pairs",
        help="Pick positive dimers, sample negative chain pairs, and "
        "assemble the final pairs table.",
    )
    commands = group.add_subparsers(dest="pairs_command", required=True)

    p = commands.add_parser(
        "positives",
        help="Group the dimers of --contacts by the similarity inputs, then "
        "keep one homomer per homomer-only group, else mutually dissimilar "
        "heteromers.",
        description="The dimer groups are the connected components of the "
        "dimer graph of the similarity inputs. They match `pandora clusters "
        "build --items interface` run on the same similarity inputs without "
        "--structures: symmetry copies are not linked to their source chain.",
    )
    add_similarity_inputs(p, positional=True)
    add_out(p)
    p.add_argument(
        "--contacts",
        required=True,
        help="residue_contacts layer folder or table: its interfaces are the dimers",
    )
    p.add_argument(
        "--structures",
        required=True,
        help="structures folder with entries and chains tables "
        "(pandora structures curate): rank values, and the structures for "
        "min_interface_coverage",
    )
    p.add_argument(
        "--chain-clusters",
        required=True,
        help="90%% chain clusters folder or clusters.json: a dimer is a "
        "homomer when both chains share a cluster",
    )
    p.add_argument(
        "--chain-representatives",
        help="representatives folder or representatives.json: map each "
        "chain to its representative, so dimers of the same two proteins "
        "are similar",
    )
    p.add_argument(
        "--policy",
        help="RepresentativePolicy YAML: tie-break rank (interface_residues, "
        "missing_residues, resolution)",
    )
    p.set_defaults(func=_cmd_positives)

    p = commands.add_parser(
        "negatives",
        help="Draw one negative chain pair per positive dimer, avoiding "
        "partners similar to anything the chain might bind.",
    )
    p.add_argument("positives", help=_POSITIVES_HELP)
    add_out(p)
    p.add_argument(
        "--representatives",
        required=True,
        help="representatives folder or representatives.json (90%% "
        "chains): the pool, and each chain's representative",
    )
    add_similarity_inputs(p, positional=False)
    p.add_argument(
        "--contacts",
        required=True,
        help="the 'any contact' residue_contacts layer folder or table",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=0,
        help="seed of the random draws (default 0)",
    )
    p.set_defaults(func=_cmd_negatives)

    p = commands.add_parser(
        "assemble",
        help="Write the final pairs table: positives and negatives with "
        "their split, kind, group and C-chain leak flags.",
    )
    p.add_argument("positives", help=_POSITIVES_HELP)
    add_out(p)
    p.add_argument(
        "--negatives",
        required=True,
        help="negatives folder or table (pandora pairs negatives)",
    )
    p.add_argument(
        "--splits",
        required=True,
        help="splits folder or splits.json (pandora splits assign)",
    )
    p.add_argument(
        "--representatives",
        required=True,
        help="representatives folder or representatives.json (the one `pairs negatives` used)",
    )
    add_similarity_inputs(p, positional=False)
    p.set_defaults(func=_cmd_assemble)
