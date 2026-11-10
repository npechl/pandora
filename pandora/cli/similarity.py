"""CLI `similarity`, `cluster` and `partition`: searches, clusters and leakage-safe splits."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import tempfile
from collections import Counter
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from pandora import __version__
from pandora._util import now_iso
from pandora.annotations import label_interfaces
from pandora.cli._io import (
    STEP_SUMMARIES,
    ParquetSink,
    check_contact_layer,
    copy_source,
    find_table,
    layer_origin,
    load_json_list,
    load_pairs,
    load_yaml_model,
    prepare_output_dir,
    read_columns,
    read_structures_dir,
)
from pandora.export import (
    export_chain_mmcif,
    read_records,
    write_json,
    write_records,
)
from pandora.ingestion import load_policy
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import (
    LayerDescription,
    LayerRow,
    ResidueContactRow,
    SimilarityItems,
)
from pandora.schemas.dataset import DatasetChainRecord
from pandora.schemas.similarity import (
    ClusteredLayer,
    HitFilter,
    RepresentativeCluster,
    RepresentativePolicy,
    RepresentativesProvenance,
    SimilarityCluster,
)
from pandora.similarity import (
    cluster_edges,
    cluster_hit_counts,
    compute_sequence_similarity,
    compute_structure_similarity,
    entry_edges,
    interface_edges,
    interface_residues_from_annotation,
    interface_similarity,
    member_representatives,
    pair_cluster_keys,
    partition_dataset,
    pick_representatives,
    restrict_clusters,
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
    parse_item_id,
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
_STEP_SUMMARIES = STEP_SUMMARIES
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


def _read_chain_list(path: Path) -> list[str]:
    """Chain item ids from a file, one per line (blank lines and `#`
    comments skipped)."""

    lines = (line.strip() for line in path.read_text().splitlines())
    return [line for line in lines if line and not line.startswith("#")]


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
    if args.chains is not None and args.input_dir is not None:
        raise SystemExit(
            "similarity search: --chains needs the chains cut from "
            "--structures-dir; drop --input-dir"
        )
    structures_dir = Path(args.structures_dir)
    items, previous = _item_set(
        "similarity search", "chain", structures_dir, None
    )
    chains: set[str] | None = None
    if args.chains is not None:
        known = set(items.ids.values())
        listed = _read_chain_list(Path(args.chains))
        unknown = next((c for c in listed if c not in known), None)
        if unknown is not None:
            raise SystemExit(
                f"similarity search: --chains names {unknown}, not a chain "
                f"of {structures_dir}"
            )
        chains = set(listed)
    out = Path(args.output_dir)
    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    hits = out / "hits.tsv"
    counts_extra: dict[str, int] = {}
    if args.engine == "mmseqs2":
        sequences = _chain_sequences(structures_dir)
        if chains is not None:
            sequences = {k: v for k, v in sequences.items() if k in chains}
        counts_extra["chains_searched"] = len(sequences)
        description = compute_sequence_similarity(
            sequences, hits, **_search_options(args)
        )
        searched: Path = structures_dir
    elif args.input_dir is not None:
        searched = Path(args.input_dir)
        description = compute_structure_similarity(
            searched, hits, **_search_options(args)
        )
    else:
        searched = structures_dir
        with tempfile.TemporaryDirectory() as tmp:
            files = _chain_files(structures_dir, chains, Path(tmp))
            counts_extra["chains_searched"] = len(files)
            description = compute_structure_similarity(
                files, hits, **_search_options(args)
            )
    if args.chains is not None:
        searched = Path(args.chains)
    counts = _write_layer(
        read_hit_tsv(hits, _CLI_ENGINES[args.engine]), items, out
    )
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
    print(
        f"searched: {counts['rows_matched']} hits, "
        f"{counts['rows_unmatched']} unmatched -> {out}"
    )


def _chain_files(
    structures_dir: Path, chains: set[str] | None, out: Path
) -> dict[str, Path]:
    """Write one mmCIF per non-copy chain of the step folder (or only
    `chains`) into `out`; returns item id -> file."""

    table = _step_table(
        "similarity search", structures_dir, "chains", "pandora curate"
    )
    wanted = {
        item_id("chain", r["entry_id"], r["chain_id"]): r
        for r in read_columns(table, ("entry_id", "chain_id", "copy_of"))
        if not r["copy_of"]
    }
    if chains is not None:
        wanted = {k: r for k, r in wanted.items() if k in chains}
    # ponytail: parses every entry; limit to the listed entries if slow.
    structures = read_structures_dir(structures_dir)
    files: dict[str, Path] = {}
    for key, row in sorted(wanted.items()):
        structure = structures.get(row["entry_id"])
        if structure is None:
            raise SystemExit(
                f"similarity search: {structures_dir} has no structure for "
                f"{row['entry_id']}"
            )
        files[key] = export_chain_mmcif(
            structure, row["chain_id"], out / f"{key}.cif"
        )
    return files


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


def _per_layer(
    flag: str, values: list[Any] | None, n: int, command: str = "cluster"
) -> list[Any]:
    """`values` of a repeatable flag, one per --similarity layer: absent,
    given once for all layers, or once per layer in order."""

    if not values:
        return [None] * n
    if len(values) == 1:
        return values * n
    if len(values) != n:
        raise SystemExit(
            f"{command}: give {flag} once (for every layer) or once per "
            f"--similarity layer ({n}), not {len(values)} times"
        )
    return values


def _hit_filters(
    args: argparse.Namespace, n: int, command: str = "cluster"
) -> list[HitFilter]:
    """One HitFilter per layer, from --hit-filter with --min-score /
    --min-coverage overriding it; exits if a layer has no threshold."""

    filters = []
    for path, min_score, min_coverage in zip(
        _per_layer("--hit-filter", args.hit_filter, n, command),
        _per_layer("--min-score", args.min_score, n, command),
        _per_layer("--min-coverage", args.min_coverage, n, command),
        strict=True,
    ):
        if path is None and min_score is None and min_coverage is None:
            # Without a filter every reported hit would become an edge
            # (the old --threshold defaulted to 0.9), merging clusters
            # silently.
            raise SystemExit(
                f"{command}: set --min-score, --min-coverage or --hit-filter "
                "for every layer; with no threshold every reported hit "
                "becomes an edge"
            )
        hit_filter = (
            load_yaml_model(HitFilter, Path(path)) if path else HitFilter()
        )
        overrides = {
            field: value
            for field, value in (
                ("min_score", min_score),
                ("min_coverage", min_coverage),
            )
            if value is not None
        }
        filters.append(hit_filter.model_copy(update=overrides))
    return filters


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
    args: argparse.Namespace, levels: list[SimilarityItems]
) -> SimilarityItems:
    """The level to cluster at: --items, or the layers' shared level.
    Chain layers roll up to entries or interfaces."""

    items = args.items or levels[0]
    for level in levels:
        if level != items and not (level == "chain" and items != "chain"):
            raise SystemExit(
                f"cluster: --items {items} does not fit a {level} layer; "
                "only chain layers roll up (to entry or interface)"
            )
    if (
        items == "interface"
        and "chain" not in levels
        and args.structures_dir is not None
    ):
        raise SystemExit(
            "cluster: interface items come from --contacts; "
            "--structures-dir does not apply"
        )
    if "entry" in levels and args.contacts is not None:
        raise SystemExit("cluster: --contacts does not apply to an entry layer")
    return items


def _interface_residues(
    structures_dir: Path | None, contacts: Path | None
) -> dict[str, set[int]] | None:
    """Interface residue positions per chain item, from a contact layer
    and the structures it refers to; None without --contacts."""

    if contacts is None:
        return None
    if structures_dir is None:
        raise SystemExit(
            "cluster: min_interface_coverage needs --structures-dir "
            "(the structures the --contacts layer was computed from)"
        )
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


def _layer_edges(
    layer: Path,
    hit_filter: HitFilter,
    interface_residues: dict[str, set[int]] | None,
    counts: Counter[int],
    index: int,
) -> Iterator[Any]:
    """The edges of one layer, counted into `counts[index]`; exits naming
    the layer on a bad table or filter."""

    table = layer / "layer.parquet"
    if not table.exists():
        return
    try:
        for edge in similarity_edges(
            read_similarity_table(table),
            hit_filter,
            interface_residues=interface_residues,
        ):
            counts[index] += 1
            yield edge
    except ValueError as exc:
        raise SystemExit(f"cluster: {layer}: {exc}") from exc


def _chain_map(
    command: str, path: Path, chains: Iterable[str]
) -> dict[str, str]:
    """Chain item id -> 90% representative from a representatives JSON;
    exits naming the first of `chains` in no members list."""

    mapping = member_representatives(
        load_json_list(RepresentativeCluster, path)
    )
    missing = next((c for c in sorted(chains) if c not in mapping), None)
    if missing is not None:
        raise SystemExit(
            f"{command}: chain {missing} is in no members list of {path}"
        )
    return mapping


def _interface_chains(interface_ids: Iterable[str]) -> set[str]:
    """The chain item ids of `interface_ids`."""

    chains = set()
    for interface in interface_ids:
        entry, *sides = parse_item_id("interface", interface)
        chains.update(item_id("chain", entry, side) for side in sides)
    return chains


def _cluster_chain_map(
    args: argparse.Namespace, items_level: str, item_ids: list[str]
) -> dict[str, str] | None:
    """`cluster`'s --chain-representatives map, or None without it."""

    if args.chain_representatives is None:
        return None
    if items_level != "interface":
        raise SystemExit(
            "cluster: --chain-representatives needs --items interface"
        )
    return _chain_map(
        "cluster", Path(args.chain_representatives), _interface_chains(item_ids)
    )


def _cmd_cluster(args: argparse.Namespace) -> None:
    """Handle the `cluster` subcommand: connected-component clustering
    of the union of one or more similarity layers."""

    layers = [Path(folder) for folder in args.similarity]
    hit_filters = _hit_filters(args, len(layers))
    _require_pyarrow("cluster")
    levels = [_similarity_level(layer) for layer in layers]
    items_level = _cluster_level(args, levels)
    structures_dir = Path(args.structures_dir) if args.structures_dir else None
    contacts = Path(args.contacts) if args.contacts else None
    if contacts is not None:
        check_contact_layer("cluster", contacts)
    items, _ = _item_set("cluster", items_level, structures_dir, contacts)
    item_ids = sorted(set(items.ids.values()))
    chain_map = _cluster_chain_map(args, items_level, item_ids)

    residues = None
    if any(
        level == "chain" and f.min_interface_coverage is not None
        for level, f in zip(levels, hit_filters, strict=True)
    ):
        residues = _interface_residues(structures_dir, contacts)
    counts: Counter[int] = Counter()
    same_level: list[Iterator[Any]] = []
    chain_streams: list[Iterator[Any]] = []
    for index, (layer, level, hit_filter) in enumerate(
        zip(layers, levels, hit_filters, strict=True)
    ):
        if level == "chain":
            chain_streams.append(
                _layer_edges(layer, hit_filter, residues, counts, index)
            )
        else:
            same_level.append(
                _layer_edges(layer, hit_filter, None, counts, index)
            )
    chain_edges: Iterable[Any] = itertools.chain(*chain_streams)
    if chain_streams and structures_dir is not None:
        chain_edges = itertools.chain(chain_edges, _copy_links(structures_dir))
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
        raise SystemExit(f"cluster: {exc}") from exc
    records = [
        ClusteredLayer(
            similarity_layer=str(layer / "layer.json"),
            items=level,
            hit_filter=hit_filter,
            n_edges=counts[index],
        )
        for index, (layer, level, hit_filter) in enumerate(
            zip(layers, levels, hit_filters, strict=True)
        )
    ]
    prov = prov.model_copy(
        update={
            "layers": records,
            "chain_representatives": args.chain_representatives,
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


def _chain_values(
    table: Path, keys: set[str]
) -> tuple[dict[str, dict[str, Any]], set[str]]:
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


_INTERFACE_RANK_KEYS = (
    "interface_residues",
    "resolution",
    "missing_residues",
    "n_hits",
)


def _interface_values(
    contacts: Path, command: str = "representatives"
) -> tuple[dict[str, dict[str, Any]], set[str], Path]:
    """Rank values per interface item of a contact layer (n_hits left for
    later), the interfaces with no symmetry-copy chain, and the
    structures folder the layer traces to."""

    rows_file = _contact_layer(command, contacts)
    origin = layer_origin(contacts)
    if origin is None:
        raise SystemExit(
            f"{command}: cannot tell which structures folder {contacts} "
            "came from; the paths in its layer.json may be relative to "
            "where `annotate compute` ran, so run from that folder"
        )
    table = _step_table(command, origin, "chains", "pandora curate")
    chains, non_copies = _chain_values(table, {"resolution"})
    rows = (
        read_records(rows_file, ResidueContactRow)
        if rows_file.name != "layer.json"
        else []
    )
    interfaces = set(interface_items(r.model_dump() for r in rows).ids.values())
    values: dict[str, dict[str, Any]] = {}
    eligible: set[str] = set()
    for row in label_interfaces(rows, {}):
        key = item_id("interface", row.entry_id, row.chain_id_1, row.chain_id_2)
        if key not in interfaces:
            continue
        sides = [
            item_id("chain", row.entry_id, chain)
            for chain in (row.chain_id_1, row.chain_id_2)
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
    return values, eligible, origin


def _cmd_representatives(args: argparse.Namespace) -> None:
    """Handle the `representatives` subcommand: one representative item
    per cluster."""

    command = "representatives"
    items = args.items
    if items == "chain" and (args.structures_dir is None or args.contacts):
        raise SystemExit(
            f"{command}: chain items need --structures-dir, not --contacts"
        )
    if items == "interface" and (args.contacts is None or args.structures_dir):
        raise SystemExit(
            f"{command}: interface items need --contacts (a "
            "residue_contacts layer folder), not --structures-dir"
        )
    loaded = (
        load_yaml_model(RepresentativePolicy, Path(args.policy))
        if args.policy
        else RepresentativePolicy()
    )
    # Revalidate so a policy without `rank` gets the default for `items`.
    policy = RepresentativePolicy.model_validate(
        loaded.model_dump(exclude_unset=True) | {"items": items}
    )
    keys = {rank_key.key for rank_key in policy.rank}
    allowed = (
        set(_INTERFACE_RANK_KEYS)
        if items == "interface"
        else {"missing_residues", "n_hits", *_CHAIN_RANK_COLUMNS}
    )
    unknown = sorted(keys - allowed)
    if unknown:
        raise SystemExit(
            f"representatives: unknown rank key {unknown[0]!r}; use "
            + ", ".join(sorted(allowed))
        )
    layers = [Path(folder) for folder in args.similarity or []]
    if "n_hits" in keys and not layers:
        raise SystemExit(
            "representatives: rank key n_hits needs --similarity (and a "
            "hit filter)"
        )
    hit_filters = _hit_filters(args, len(layers), command) if layers else []
    if items == "interface":
        values, eligible, structures_dir = _interface_values(
            Path(args.contacts)
        )
    else:
        structures_dir = Path(args.structures_dir)
        table = _step_table(command, structures_dir, "chains", "pandora curate")
        values, eligible = _chain_values(table, keys)
    clusters_in = load_json_list(SimilarityCluster, Path(args.clusters))
    clusters = restrict_clusters(clusters_in, values.keys())

    counts: Counter[int] = Counter()
    if "n_hits" in keys:
        _require_pyarrow("representatives")
        for layer in layers:
            if _similarity_level(layer) != items:
                raise SystemExit(
                    f"representatives: {layer} is not a {items} layer"
                )
        edges = itertools.chain(
            *(
                _layer_edges(layer, hit_filter, None, counts, index)
                for index, (layer, hit_filter) in enumerate(
                    zip(layers, hit_filters, strict=True)
                )
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
        clusters=str(args.clusters),
        structures_dir=str(structures_dir),
        contacts=args.contacts,
        layers=[
            ClusteredLayer(
                similarity_layer=str(layer / "layer.json"),
                items=items,
                hit_filter=hit_filter,
                n_edges=counts[index],
            )
            for index, (layer, hit_filter) in enumerate(
                zip(layers, hit_filters, strict=True)
            )
        ],
        n_clusters=len(reps),
        n_members=n_members,
        n_members_dropped=n_in - n_members,
        n_clusters_dropped=len(clusters_in) - len(clusters),
    )
    output = Path(args.output)
    write_records(reps, output)
    write_json(prov, output.parent / "representatives_provenance.json")
    print(f"{len(reps)} representatives from {n_members} {items}s -> {output}")


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
        prov = prov.model_copy(
            update={"labels_layer": str(Path(args.labels_dir) / "layer.json")}
        )

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


def _cmd_interfaces(args: argparse.Namespace) -> None:
    """Handle `similarity interfaces`: score pairs of a contact layer's
    interfaces by per-side Jaccard, as a similarity layer folder."""

    command = "similarity interfaces"
    started_at = now_iso()
    _require_pyarrow(command)
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
    sequences = {
        (chain.entry_id.upper(), chain.chain_id): chain.sequence
        for chain in read_records(chains_table, DatasetChainRecord)
        if chain.sequence
    }
    clusters_file = Path(args.chain_clusters)
    clusters = load_json_list(SimilarityCluster, clusters_file)
    try:
        for cluster in clusters:
            for component in cluster.components:
                parse_item_id("chain", component)
    except ValueError as exc:
        raise SystemExit(
            f"{command}: {clusters_file} must cluster chain items "
            f"(ENTRY_CHAIN): {exc}"
        ) from exc
    table = find_table(contacts, "layer")
    rows = read_records(table, ResidueContactRow) if table else []
    scored, counts = interface_similarity(
        rows, clusters, sequences, args.min_identity
    )

    prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)
    if scored:
        write_records(scored, out / "layer.parquet")
    source = json.loads((contacts / "layer.json").read_text())["policy"]
    description = LayerDescription(
        name="iface_jaccard",
        template="similarity",
        items="interface",
        tool="pandora",
        tool_version=__version__,
        parameters={
            "contacts_layer": source.get("name"),
            "structures_dir": str(origin),
            "chain_clusters": str(clusters_file),
            "min_identity": args.min_identity,
        },
    )
    _write_summary(
        "similarity_interfaces",
        description,
        input=contacts,
        out=out,
        started_at=started_at,
        counts=Counter({"rows_in": len(rows)} | counts),
        previous=contacts / "layer.json",
    )
    print(
        f"scored: {counts['pairs_compared']} interface pairs "
        f"({counts['pairs_zero_score']} zero, "
        f"{counts['pairs_below_identity']} below identity) from "
        f"{counts['interfaces']} interfaces -> {out}"
    )


def _register_interfaces(commands: argparse._SubParsersAction) -> None:
    """Add `similarity interfaces` to the `similarity` subcommands."""

    p = commands.add_parser(
        "interfaces",
        help="Score pairs of interfaces of the same two proteins by "
        "per-side Jaccard of their interface residues.",
    )
    p.add_argument(
        "--contacts",
        dest="contacts_dir",
        required=True,
        help="residue_contacts layer folder: the interfaces",
    )
    p.add_argument(
        "--chain-clusters",
        required=True,
        help="clusters.json of chains (pandora cluster), e.g. 90%% sequence "
        "clusters: only interfaces whose chains share clusters are compared",
    )
    p.add_argument(
        "--min-identity",
        type=float,
        default=0.9,
        help="lowest pairwise identity per chain pair (default 0.9)",
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument("--overwrite", action="store_true")
    p.set_defaults(func=_cmd_interfaces)


def _register_representatives(subparsers: argparse._SubParsersAction) -> None:
    """Add the `representatives` subcommand."""

    p = subparsers.add_parser(
        "representatives",
        help="Pick one representative item per cluster, by a ranking policy.",
    )
    p.add_argument(
        "--clusters", required=True, help="clusters.json (pandora cluster)"
    )
    p.add_argument(
        "--structures-dir",
        help="chain items: the step folder whose chains the clusters are "
        "limited to",
    )
    p.add_argument(
        "--contacts",
        help="interface items: the residue_contacts layer folder whose "
        "interfaces the clusters are limited to",
    )
    p.add_argument(
        "--items",
        choices=["chain", "interface"],
        default="chain",
        help="the clusters' item level (default chain)",
    )
    p.add_argument(
        "--similarity",
        action="append",
        help="similarity layer folder (at --items level) for the n_hits "
        "rank key; repeatable",
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
    p.add_argument("--policy", help="RepresentativePolicy YAML (rank keys)")
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_representatives)


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
        "--input-dir",
        help="foldseek only: dir of per-chain structure files (default: "
        "Pandora writes them from --structures-dir)",
    )
    p.add_argument(
        "--chains",
        help="file of chain item ids (one per line, # comments) to limit "
        "the search to",
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

    _register_interfaces(commands)

    p = subparsers.add_parser(
        "cluster",
        help="Connected-component clustering of one or more similarity "
        "layers (an edge from any layer links).",
    )
    p.add_argument(
        "--similarity",
        required=True,
        action="append",
        help="a similarity layer folder (similarity search/import); "
        "repeat to join several layers",
    )
    p.add_argument(
        "--items",
        choices=["entry", "chain", "interface"],
        help="cluster at this level (default: the layers'); chain layers "
        "roll up to entries, or to the interfaces of --contacts",
    )
    p.add_argument(
        "--hit-filter",
        action="append",
        help="HitFilter YAML (thresholds); once for every layer or once "
        "per --similarity, in order (same for the two flags below)",
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
    p.add_argument(
        "--contacts",
        help=(
            "a residue_contacts layer folder: the interface items, or "
            "interface residues for min_interface_coverage on a chain "
            "layer"
        ),
    )
    p.add_argument(
        "--structures-dir",
        help="the step folder whose entries/chains are the items",
    )
    p.add_argument(
        "--chain-representatives",
        help="representatives.json (pandora representatives): map each "
        "chain to its representative before rolling chain layers up to "
        "interfaces",
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

    _register_representatives(subparsers)

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
