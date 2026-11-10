"""CLI `curate` and `dedup`: filter and deduplicate a directory of structures."""

from __future__ import annotations

import argparse
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

from pandora._util import now_iso
from pandora.cli._io import (
    ENTRY_ERRORS,
    TABLE_FORMATS,
    check_table_format,
    duplicate_line,
    failure_line,
    find_table,
    load_json_list,
    load_yaml_model,
    prepare_step_dir,
    read_structures_dir,
    structure_files,
    write_structures_dir,
    write_table,
)
from pandora.datasets import (
    curate_structure,
    dataset_chain_records,
    deduplicate_structures,
    find_contained,
)
from pandora.export import (
    read_records,
    structure_to_mmcif,
    write_json,
    write_records,
)
from pandora.ingestion import mirror_path
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.provenance import build_step_summary, curation_entry_record
from pandora.schemas.canonicalisation import ChainCopyRecord
from pandora.schemas.dataset import (
    DatasetChainRecord,
    DatasetCurationPolicy,
    DeduplicationRules,
    ExclusionRecord,
)
from pandora.schemas.provenance import CurationEntryRecord
from pandora.schemas.similarity import SimilarityCluster
from pandora.schemas.structure import Structure
from pandora.similarity import item_id


def _read_copies(folder: Path) -> dict[str, dict[str, str]]:
    """Entry id -> {copy chain id: source chain id}, from a canonicalise
    step folder's `chain_copies` table (empty if there is none)."""

    copies: dict[str, dict[str, str]] = defaultdict(dict)
    table = find_table(folder, "chain_copies")
    for row in read_records(table, ChainCopyRecord) if table else []:
        copies[row.entry_id][row.chain_id] = row.source_chain_id
    return copies


def _curate_one(
    structure: Structure,
    policy: DatasetCurationPolicy,
    copies: dict[str, dict[str, str]],
    out: Path,
) -> tuple[
    list[ExclusionRecord], list[DatasetChainRecord], CurationEntryRecord | None
]:
    """Curate one structure; write it under `out` if kept. Returns its
    exclusions, chain rows and entry row (None if not kept)."""

    metadata = collect_metadata(structure)
    curated, dropped, prov = curate_structure(structure, metadata, policy)
    if curated is None:
        return dropped, [], None
    written = mirror_path(Path("structures"), curated.entry_id.lower())
    structure_to_mmcif(curated, out / written)
    rows = dataset_chain_records(
        curated,
        metadata,
        policy.quality_rules.missing_residue_definition,
        copy_of=copies.get(curated.entry_id),
        structure_path=written,
    )
    return (
        dropped,
        rows,
        curation_entry_record(curated, metadata, prov, written),
    )


def _cmd_curate(args: argparse.Namespace) -> None:
    """Handle the `curate` subcommand: curate each structure and write
    the step folder."""

    started_at = now_iso()
    in_dir = Path(args.input_dir)
    out = Path(args.output_dir)
    policy = load_yaml_model(DatasetCurationPolicy, Path(args.policy))
    copies = _read_copies(in_dir)

    check_table_format(args.table_format)
    prepare_step_dir(out, "curation.json", args.overwrite)
    counts: Counter[str] = Counter(entries_failed=0)
    entries: list[CurationEntryRecord] = []
    chains: list[DatasetChainRecord] = []
    exclusions: list[ExclusionRecord] = []
    seen: set[str] = set()
    with (out / "failures.jsonl").open("w") as failures:
        for path in structure_files(in_dir):
            counts["entries_in"] += 1
            try:
                structure, _, _ = mmcif_to_structure(str(path))
                if structure is None:
                    raise ValueError(f"could not parse {path}")
                key = structure.entry_id.lower()
                if key in seen:
                    counts["entries_failed"] += 1
                    failures.write(duplicate_line(path, structure.entry_id))
                    continue
                seen.add(key)
                dropped, rows, entry = _curate_one(
                    structure, policy, copies, out
                )
            except ENTRY_ERRORS as exc:
                counts["entries_failed"] += 1
                failures.write(failure_line(path, type(exc).__name__, str(exc)))
                continue
            exclusions += dropped
            chains += rows
            if entry is not None:
                entries.append(entry)

    for e in exclusions:
        level = "entry" if e.chain_id is None else "chain"
        counts[f"excluded_{level}:{e.reason_code}"] += 1
    counts["entries_kept"] = len(entries)
    counts["chains_kept"] = len(chains)
    write_table(entries, out, "entries", args.table_format)
    write_table(chains, out, "chains", args.table_format)
    write_table(exclusions, out, "exclusions", args.table_format)
    previous = in_dir / "canonicalisation.json"
    write_json(
        build_step_summary(
            "curation",
            policy,
            input=in_dir,
            output=out,
            started_at=started_at,
            counts=dict(sorted(counts.items())),
            previous_summary=previous if previous.exists() else None,
        ),
        out / "curation.json",
    )
    failed = counts["entries_failed"]
    excluded = counts["entries_in"] - len(entries) - failed
    print(
        f"curated: {len(entries)} retained, {excluded} excluded -> {out}"
        + (f", {failed} failed" if failed else "")
    )


_CONTAINED_RANK_KEYS = ("resolution", "chain_count")


def _read_step_tables(
    folder: Path,
) -> tuple[list[CurationEntryRecord], list[DatasetChainRecord]]:
    """The entries and chains tables of a curate step folder."""

    entries, chains = (
        find_table(folder, "entries"),
        find_table(folder, "chains"),
    )
    if entries is None or chains is None:
        raise SystemExit(
            "dedup: --contained needs a curate step folder (entries and "
            "chains tables)"
        )
    return (
        read_records(entries, CurationEntryRecord),
        read_records(chains, DatasetChainRecord),
    )


def _cluster_counts(
    chains: list[DatasetChainRecord], clusters: list[SimilarityCluster]
) -> tuple[dict[str, Counter[str]], int]:
    """Entry id -> chain count per cluster, and how many chains were in
    no cluster (each counts as its own cluster)."""

    cluster_of = {
        member: str(index)
        for index, cluster in enumerate(clusters)
        for member in cluster.components
    }
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    unclustered = 0
    for chain in chains:
        key = item_id("chain", chain.entry_id, chain.chain_id)
        if key not in cluster_of:
            unclustered += 1
        counts[chain.entry_id][cluster_of.get(key, key)] += 1
    return counts, unclustered


def _cmd_dedup_contained(args: argparse.Namespace) -> None:
    """Handle `dedup --contained`: drop entries contained in another and
    write the step folder."""

    if not args.clusters:
        raise SystemExit("dedup: --contained needs --clusters")
    started_at = now_iso()
    in_dir, out = Path(args.input_dir), Path(args.output_dir)
    rules = (
        load_yaml_model(DeduplicationRules, Path(args.policy))
        if args.policy
        else DeduplicationRules()
    )
    for rank_key in rules.contained_rank:
        if rank_key.key not in _CONTAINED_RANK_KEYS:
            raise SystemExit(
                f"dedup: unknown contained_rank key {rank_key.key!r}; use "
                + ", ".join(_CONTAINED_RANK_KEYS)
            )
    check_table_format(args.table_format)
    entries, chains = _read_step_tables(in_dir)
    clusters = load_json_list(SimilarityCluster, Path(args.clusters))
    prepare_step_dir(out, "dedup.json", args.overwrite)

    counts, unclustered = _cluster_counts(chains, clusters)
    for entry in entries:
        counts.setdefault(entry.entry_id, Counter())
    values = {
        e.entry_id: {k: getattr(e, k) for k in _CONTAINED_RANK_KEYS}
        for e in entries
    }
    contained = find_contained(counts, values, rules.contained_rank)
    removed = set() if args.disable else {r.entry_id for r in contained}
    kept = [e for e in entries if e.entry_id not in removed]
    kept_chains = [c for c in chains if c.entry_id not in removed]
    for entry in kept:
        if entry.structure_path:
            dest = out / entry.structure_path
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(in_dir / entry.structure_path, dest)

    write_table(kept, out, "entries", args.table_format)
    write_table(kept_chains, out, "chains", args.table_format)
    write_table(contained, out, "contained", args.table_format)
    (out / "failures.jsonl").write_text("")
    previous = in_dir / "curation.json"
    write_json(
        build_step_summary(
            "dedup",
            rules,
            input=in_dir,
            output=out,
            started_at=started_at,
            counts={
                "chains_kept": len(kept_chains),
                "chains_unclustered": unclustered,
                "entries_contained": len(contained),
                "entries_in": len(entries),
                "entries_kept": len(kept),
            },
            previous_summary=previous if previous.exists() else None,
        ),
        out / "dedup.json",
    )
    if unclustered:
        print(
            f"dedup: warning: {unclustered} chains are in no cluster of "
            f"{args.clusters} (each counted as its own cluster); were the "
            "clusters built from this folder?",
            file=sys.stderr,
        )
    print(f"dedup: {len(kept)} retained, {len(contained)} contained -> {out}")


def _cmd_dedup(args: argparse.Namespace) -> None:
    """Handle the `dedup` subcommand: remove structures sharing the same
    entry_id, or with --contained, entries contained in another."""

    if args.contained:
        _cmd_dedup_contained(args)
        return
    for flag, value in (
        ("--clusters", args.clusters),
        ("--policy", args.policy),
        ("--overwrite", args.overwrite),
    ):
        if value:
            raise SystemExit(f"dedup: {flag} only applies with --contained")
    output_dir = Path(args.output_dir)
    structures = read_structures_dir(Path(args.input_dir))
    rules = DeduplicationRules(enabled=not args.disable)

    retained, removed, prov = deduplicate_structures(
        list(structures.values()), rules
    )

    write_structures_dir({s.entry_id: s for s in retained}, output_dir)
    write_json(prov, output_dir / "dedup_provenance.json")
    if removed:
        write_records(removed, output_dir / "dedup_exclusions.json")
    print(
        f"dedup: {len(retained)} retained, {len(removed)} removed "
        f"-> {output_dir}"
    )


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    p = subparsers.add_parser(
        "curate",
        help="Apply quality/organism/content rules; write a step folder.",
    )
    p.add_argument(
        "--input-dir",
        required=True,
        help="a canonicalise step folder, or a folder of mmCIF files",
    )
    p.add_argument("--policy", required=True, help="DatasetCurationPolicy YAML")
    p.add_argument("--output-dir", required=True)
    p.add_argument(
        "--table-format",
        choices=TABLE_FORMATS,
        default="parquet",
        help="format of the entries/chains/exclusions tables",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="replace step outputs already in --output-dir",
    )
    p.set_defaults(func=_cmd_curate)

    p = subparsers.add_parser(
        "dedup", help="Remove structures sharing the same entry_id."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument(
        "--disable", action="store_true", help="record but keep dupes"
    )
    p.add_argument(
        "--contained",
        action="store_true",
        help="remove entries contained in another (by chain counts per "
        "cluster of --clusters); writes a step folder",
    )
    p.add_argument(
        "--clusters", help="clusters.json from `pandora cluster` (chains)"
    )
    p.add_argument("--policy", help="DeduplicationRules YAML (contained_rank)")
    p.add_argument(
        "--table-format",
        choices=TABLE_FORMATS,
        default="parquet",
        help="format of the step folder tables (--contained)",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="replace step outputs already in --output-dir (--contained)",
    )
    p.set_defaults(func=_cmd_dedup)
