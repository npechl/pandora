"""CLI `curate` and `dedup`: filter and deduplicate a directory of structures."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

from pandora._util import now_iso
from pandora.cli._io import (
    ENTRY_ERRORS,
    TABLE_FORMATS,
    check_table_format,
    failure_line,
    find_table,
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


def _read_copies(folder: Path) -> dict[str, dict[str, str]]:
    """Entry id -> {copy chain id: source chain id}, from a canonicalise
    step folder's `chain_copies` table (empty if there is none)."""

    copies: dict[str, dict[str, str]] = defaultdict(dict)
    table = find_table(folder, "chain_copies")
    for row in read_records(table, ChainCopyRecord) if table else []:
        copies[row.entry_id][row.chain_id] = row.source_chain_id
    return copies


def _curate_one(
    path: Path,
    policy: DatasetCurationPolicy,
    copies: dict[str, dict[str, str]],
    out: Path,
) -> tuple[
    list[ExclusionRecord], list[DatasetChainRecord], CurationEntryRecord | None
]:
    """Curate one structure file; write it under `out` if kept. Returns
    its exclusions, chain rows and entry row (None if not kept)."""

    structure, _, _ = mmcif_to_structure(str(path))
    if structure is None:
        raise ValueError(f"could not parse {path}")
    metadata = collect_metadata(structure)
    curated, dropped, prov = curate_structure(structure, metadata, policy)
    if curated is None:
        return dropped, [], None
    written = mirror_path(Path("structures"), path.stem)
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
            if path.stem in seen:
                counts["entries_failed"] += 1
                failures.write(
                    failure_line(
                        path,
                        "DuplicateStem",
                        f"another input is also named {path.stem}",
                    )
                )
                continue
            seen.add(path.stem)
            try:
                dropped, rows, entry = _curate_one(path, policy, copies, out)
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


def _cmd_dedup(args: argparse.Namespace) -> None:
    """Handle the `dedup` subcommand: remove structures sharing the same
    entry_id."""

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
    p.set_defaults(func=_cmd_dedup)
