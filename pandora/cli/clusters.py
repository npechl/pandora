"""CLI `clusters build` and `clusters representatives`: connected-component
clusters of similarity layers, and one representative item per cluster."""

from __future__ import annotations

import argparse
import itertools
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from pandora._util import now_iso
from pandora.annotations import label_interfaces
from pandora.cli._io import (
    SUMMARY,
    add_out,
    add_similarity_inputs,
    check_out,
    contact_rows,
    contact_table,
    input_summary,
    layer_edges,
    load_json_list,
    load_yaml_model,
    prepare_output_dir,
    read_columns,
    read_structure,
    require_tables,
    resolve_file,
    rows_by_entry,
    similarity_inputs,
    structure_index,
    write_json_list,
    write_summary,
)
from pandora.cli.similarity import _item_set, _require_pyarrow
from pandora.datasets import interface_residue_pairs
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import ResidueContactRow, SimilarityItems
from pandora.schemas.similarity import (
    ClusteredLayer,
    RepresentativeCluster,
    RepresentativePolicy,
    RepresentativesProvenance,
    SimilarityCluster,
)
from pandora.similarity import (
    cluster_edges,
    cluster_hit_counts,
    entry_edges,
    interface_edges,
    interface_residues_from_annotation,
    member_representatives,
    pick_representatives,
    restrict_clusters,
)
from pandora.similarity.layer import (
    interface_items,
    item_id,
    parse_item_id,
)

_CHAIN_RANK_COLUMNS = (
    "residue_count",
    "atom_count",
    "seqres_length",
    "missing_n_term",
    "missing_c_term",
    "missing_middle",
    "resolution",
)
_MISSING_COLUMNS = ("missing_n_term", "missing_c_term", "missing_middle")
_INTERFACE_RANK_KEYS = (
    "interface_residues",
    "resolution",
    "missing_residues",
    "n_hits",
)


def _cluster_level(args: argparse.Namespace, levels: list[SimilarityItems]) -> SimilarityItems:
    """The level to cluster at: --items, or the layers' shared level.
    Chain layers roll up to entries or interfaces."""

    command = "clusters build"
    items = args.items or levels[0]
    for level in levels:
        if level != items and not (level == "chain" and items != "chain"):
            raise SystemExit(
                f"{command}: --items {items} does not fit a {level} layer; "
                "only chain layers roll up (to entry or interface)"
            )
    if items == "interface" and "chain" not in levels and args.structures is not None:
        raise SystemExit(
            f"{command}: interface items come from --contacts; --structures does not apply"
        )
    if "entry" in levels and args.contacts is not None:
        raise SystemExit(f"{command}: --contacts does not apply to an entry layer")
    return items


def _interface_residues(
    command: str, structures: Path | None, contacts: Path | None
) -> dict[str, set[int]] | None:
    """Interface residue positions per chain item, from a contact layer
    and the structures it refers to; None without --contacts."""

    if contacts is None:
        return None
    if structures is None:
        raise SystemExit(
            f"{command}: min_interface_coverage needs --structures "
            "(the structures the --contacts layer was computed from)"
        )
    rows_file = contact_table(command, contacts)
    groups = rows_by_entry(rows_file) if rows_file is not None else ()
    paths = structure_index(structures)
    result: dict[str, set[int]] = {}
    for entry, rows in groups:
        entry_rows = [ResidueContactRow.model_validate(r) for r in rows]
        if not interface_residue_pairs(entry_rows):
            continue  # no inter-chain polymer contact: nothing to place
        path = paths.get(entry)
        structure = read_structure(path) if path else None
        if structure is None:
            raise SystemExit(
                f"{command}: {contacts} names an entry missing from {structures}: {entry!r}"
            )
        by_id = {row.entry_id: structure for row in entry_rows}
        result.update(interface_residues_from_annotation(by_id, entry_rows))
    return result


def _copy_links(chains_table: Path) -> Iterator[tuple[str, str]]:
    """(copy, source) chain item pairs from a chains table, so symmetry
    copies join their source chain's cluster."""

    for row in read_columns(chains_table, ("entry_id", "chain_id", "copy_of")):
        if row["copy_of"]:
            yield (
                item_id("chain", row["entry_id"], row["chain_id"]),
                item_id("chain", row["entry_id"], row["copy_of"]),
            )


def _chain_map(command: str, path: Path, chains: Iterable[str]) -> dict[str, str]:
    """Chain item id -> 90% representative from a representatives JSON;
    exits naming the first of `chains` in no members list."""

    mapping = member_representatives(load_json_list(command, RepresentativeCluster, path))
    missing = next((c for c in sorted(chains) if c not in mapping), None)
    if missing is not None:
        raise SystemExit(f"{command}: chain {missing} is in no members list of {path}")
    return mapping


def _interface_chains(interface_ids: Iterable[str]) -> set[str]:
    """The chain item ids of `interface_ids`."""

    chains = set()
    for interface in interface_ids:
        entry, *sides = parse_item_id("interface", interface)
        chains.update(item_id("chain", entry, side) for side in sides)
    return chains


def _build_chain_map(
    args: argparse.Namespace, items_level: str, item_ids: list[str]
) -> tuple[dict[str, str] | None, Path | None]:
    """`build`'s --chain-representatives file and its chain map, or
    (None, None) without it."""

    if args.chain_representatives is None:
        return None, None
    command = "clusters build"
    if items_level != "interface":
        raise SystemExit(f"{command}: --chain-representatives needs --items interface")
    file = resolve_file(command, Path(args.chain_representatives), "representatives.json")
    return _chain_map(command, file, _interface_chains(item_ids)), file


def _cmd_build(args: argparse.Namespace) -> None:
    """Handle `clusters build`: connected-component clustering of the
    union of one or more similarity inputs."""

    command = "clusters build"
    started_at = now_iso()
    out = Path(args.out)
    structures = Path(args.structures) if args.structures else None
    contacts = Path(args.contacts) if args.contacts else None
    reps = args.chain_representatives
    check_out(
        command,
        out,
        *(Path(p) for p in args.similarity),
        structures,
        contacts,
        Path(reps) if reps else None,
    )
    _require_pyarrow(command)
    inputs = similarity_inputs(command, args)
    levels = [i.items for i in inputs]
    items_level = _cluster_level(args, levels)
    if structures is not None:
        require_tables(command, structures)
    if contacts is not None:
        contact_table(command, contacts)
    items, _ = _item_set(command, items_level, structures, contacts)
    item_ids = sorted(set(items.ids.values()))
    chain_map, map_file = _build_chain_map(args, items_level, item_ids)

    residues = None
    if any(i.items == "chain" and i.hit_filter.min_interface_coverage is not None for i in inputs):
        residues = _interface_residues(command, structures, contacts)
    counts: Counter[int] = Counter()
    same_level: list[Iterator[Any]] = []
    chain_streams: list[Iterator[Any]] = []
    for index, source in enumerate(inputs):
        chain = source.items == "chain"
        (chain_streams if chain else same_level).append(
            layer_edges(
                source.table or source.path,
                source.hit_filter,
                residues if chain else None,
                counts,
                index,
                command,
            )
        )
    chain_edges: Iterable[Any] = itertools.chain(*chain_streams)
    if chain_streams and structures is not None:
        chain_edges = itertools.chain(
            chain_edges, _copy_links(require_tables(command, structures)[1])
        )
    if chain_streams and items_level == "entry":
        chain_edges = entry_edges(chain_edges)
    elif (chain_streams or chain_map) and items_level == "interface":
        # With a chain map, dimers sharing a representative link even
        # without a chain layer.
        chain_edges = interface_edges(chain_edges, item_ids, chain_map)
    edges = itertools.chain(*same_level, chain_edges)
    try:
        clusters, prov = cluster_edges(item_ids, edges)
    except ValueError as exc:
        raise SystemExit(f"{command}: {exc}") from exc
    prov = prov.model_copy(
        update={
            "layers": [
                ClusteredLayer(
                    similarity_layer=str(source.path.resolve()),
                    items=source.items,
                    hit_filter=source.hit_filter,
                    n_edges=counts[index],
                )
                for index, source in enumerate(inputs)
            ],
            "chain_representatives": (None if map_file is None else str(map_file.resolve())),
        }
    )

    prepare_output_dir(command, out, ["clusters.json", SUMMARY], args.overwrite)
    write_json_list(clusters, out / "clusters.json")
    first = inputs[0].path
    write_summary(
        out,
        build_step_summary(
            "cluster",
            prov,
            input=first.resolve(),
            output=out,
            started_at=started_at,
            counts={"clusters": len(clusters), "edges": prov.n_edges},
            previous_summary=input_summary(first) if first.is_dir() else None,
        ),
    )
    print(f"{len(clusters)} clusters from {prov.n_edges} edges -> {out}")


def _chain_values(table: Path, keys: set[str]) -> tuple[dict[str, dict[str, Any]], set[str]]:
    """Rank values per chain item from a chains table (missing_residues
    summed; n_hits left for later), and the non-copy chain items."""

    columns = ("entry_id", "chain_id", "copy_of", *_CHAIN_RANK_COLUMNS)
    values: dict[str, dict[str, Any]] = {}
    eligible: set[str] = set()
    for row in read_columns(table, columns):
        key = item_id("chain", row["entry_id"], row["chain_id"])
        parts = [row[c] for c in _MISSING_COLUMNS]
        missing = None if None in parts else sum(parts)
        values[key] = {k: row[k] for k in keys & set(_CHAIN_RANK_COLUMNS)}
        values[key]["missing_residues"] = missing
        if not row["copy_of"]:
            eligible.add(key)
    return values, eligible


def _interface_values(
    command: str, contacts: Path, chains_table: Path
) -> tuple[dict[str, dict[str, Any]], set[str]]:
    """Rank values per interface item of a contact layer (n_hits left for
    later; chain values from `chains_table`) and the interfaces with no
    symmetry-copy chain."""

    rows_file = contact_table(command, contacts)
    chains, non_copies = _chain_values(chains_table, {"resolution"})
    rows = (row.model_dump() for row in contact_rows(rows_file, Counter()))
    interfaces = set(interface_items(rows).ids.values())
    values: dict[str, dict[str, Any]] = {}
    eligible: set[str] = set()
    for row in label_interfaces(contact_rows(rows_file, Counter()), {}):
        key = item_id("interface", row.entry_id, row.chain_id_1, row.chain_id_2)
        if key not in interfaces:
            continue
        sides = [
            item_id("chain", row.entry_id, chain) for chain in (row.chain_id_1, row.chain_id_2)
        ]
        found = [chains.get(side, {}) for side in sides]
        missing = [f.get("missing_residues") for f in found]
        resolution = [f.get("resolution") for f in found]
        values[key] = {
            "interface_residues": row.residues_1 + row.residues_2,
            "resolution": next((r for r in resolution if r is not None), None),
            "missing_residues": None if None in missing else sum(missing),
        }
        if all(side in non_copies for side in sides):
            eligible.add(key)
    return values, eligible


def _rank_policy(command: str, args: argparse.Namespace, level: str) -> RepresentativePolicy:
    """The ranking policy for `level`: --policy (or the default), with
    its rank keys checked."""

    loaded = (
        load_yaml_model(command, RepresentativePolicy, Path(args.policy))
        if args.policy
        else RepresentativePolicy()
    )
    # Revalidate so a policy without `rank` gets the default for `items`.
    policy = RepresentativePolicy.model_validate(
        loaded.model_dump(exclude_unset=True) | {"items": level}
    )
    allowed = (
        set(_INTERFACE_RANK_KEYS)
        if level == "interface"
        else {"missing_residues", "n_hits", *_CHAIN_RANK_COLUMNS}
    )
    unknown = sorted({k.key for k in policy.rank} - allowed)
    if unknown:
        raise SystemExit(
            f"{command}: unknown rank key {unknown[0]!r}; use " + ", ".join(sorted(allowed))
        )
    return policy


def _cmd_representatives(args: argparse.Namespace) -> None:
    """Handle `clusters representatives`: one representative item per
    cluster."""

    command = "clusters representatives"
    started_at = now_iso()
    out = Path(args.out)
    structures = Path(args.structures)
    contacts = Path(args.contacts) if args.contacts else None
    check_out(
        command,
        out,
        Path(args.clusters),
        structures,
        contacts,
        *(Path(p) for p in args.similarity or []),
    )
    level = "interface" if args.contacts else "chain"
    policy = _rank_policy(command, args, level)
    keys = {rank_key.key for rank_key in policy.rank}
    inputs = similarity_inputs(command, args)
    if "n_hits" in keys and not inputs:
        raise SystemExit(f"{command}: rank key n_hits needs --similarity (and a hit filter)")
    clusters_file = resolve_file(command, Path(args.clusters), "clusters.json")
    _, chains_table = require_tables(command, structures)
    if contacts is None:
        values, eligible = _chain_values(chains_table, keys)
    else:
        values, eligible = _interface_values(command, contacts, chains_table)
    clusters_in = load_json_list(command, SimilarityCluster, clusters_file)
    clusters = restrict_clusters(clusters_in, values.keys())

    counts: Counter[int] = Counter()
    if "n_hits" in keys:
        _require_pyarrow(command)
        for source in inputs:
            if source.items != level:
                raise SystemExit(
                    f"{command}: {source.path} is a {source.items} layer, "
                    f"not {level} (--contacts picks interface items)"
                )
        edges = itertools.chain(
            *(
                layer_edges(
                    source.table or source.path,
                    source.hit_filter,
                    None,
                    counts,
                    index,
                    command,
                )
                for index, source in enumerate(inputs)
            )
        )
        for item, n in cluster_hit_counts(clusters, edges).items():
            values[item]["n_hits"] = n
    reps = pick_representatives(clusters, values, policy, eligible)

    n_in = sum(c.n_components for c in clusters_in)
    n_members = sum(c.n_components for c in clusters)
    prov = RepresentativesProvenance(
        selected_at=now_iso(),
        items=policy.items,
        policy=policy,
        clusters=str(clusters_file.resolve()),
        structures_dir=str(structures.resolve()),
        contacts=None if contacts is None else str(contacts.resolve()),
        layers=[
            ClusteredLayer(
                similarity_layer=str(source.path.resolve()),
                items=source.items,
                hit_filter=source.hit_filter,
                n_edges=counts[index],
            )
            for index, source in enumerate(inputs)
        ],
        n_clusters=len(reps),
        n_members=n_members,
        n_members_dropped=n_in - n_members,
        n_clusters_dropped=len(clusters_in) - len(clusters),
    )
    prepare_output_dir(command, out, ["representatives.json", SUMMARY], args.overwrite)
    write_json_list(reps, out / "representatives.json")
    write_summary(
        out,
        build_step_summary(
            "representatives",
            prov,
            input=clusters_file.resolve(),
            output=out,
            started_at=started_at,
            counts={"clusters": len(reps), "members": n_members},
            previous_summary=input_summary(clusters_file.parent),
        ),
    )
    print(f"{len(reps)} representatives from {n_members} {level}s -> {out}")


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add the `clusters` group to the `pandora` parser."""

    group = subparsers.add_parser(
        "clusters",
        help="Cluster similarity layers and pick one representative per cluster.",
    )
    commands = group.add_subparsers(dest="clusters_command", required=True)

    p = commands.add_parser(
        "build",
        help="Connected-component clustering of one or more similarity "
        "inputs (an edge from any input links).",
    )
    add_similarity_inputs(p, positional=True)
    add_out(p, tables=False)
    p.add_argument(
        "--items",
        choices=["entry", "chain", "interface"],
        help="cluster at this level (default: the inputs'); chain inputs "
        "roll up to entries, or to the interfaces of --contacts",
    )
    p.add_argument(
        "--structures",
        help="structures folder with entries and chains tables "
        "(pandora structures curate): the entry or chain items, and the "
        "structures for min_interface_coverage",
    )
    p.add_argument(
        "--contacts",
        help="residue_contacts layer folder or table: the interface items, "
        "or the interface residues for min_interface_coverage on a chain "
        "input",
    )
    p.add_argument(
        "--chain-representatives",
        help="representatives folder or representatives.json: map each "
        "chain to its representative before rolling chain inputs up to "
        "interfaces",
    )
    p.set_defaults(func=_cmd_build)

    p = commands.add_parser(
        "representatives",
        help="Pick one representative item per cluster, by a ranking policy.",
    )
    p.add_argument("clusters", help="clusters folder or clusters.json")
    add_out(p, tables=False)
    p.add_argument(
        "--structures",
        required=True,
        help="structures folder with entries and chains tables "
        "(pandora structures curate): the chains the clusters are limited "
        "to, their rank values and symmetry copies",
    )
    p.add_argument(
        "--contacts",
        help="residue_contacts layer folder or table: pick interface items "
        "(the clusters are limited to its interfaces); without it, chain "
        "items",
    )
    add_similarity_inputs(p, positional=False)
    p.add_argument(
        "--policy",
        help="RepresentativePolicy YAML (rank keys); n_hits needs "
        "--similarity at the clusters' level",
    )
    p.set_defaults(func=_cmd_representatives)
