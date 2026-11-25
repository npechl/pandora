"""CLI `structures canonicalise`, `curate`, `drop-contained` and `convert`:
transform structures into step folders, or convert one file."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO

from pandora._util import now_iso
from pandora.canonicalisation import (
    canonicalise_structure,
    chain_copy_records,
    estimate_atom_count,
)
from pandora.cli._io import (
    ENTRY_ERRORS,
    add_out,
    check_out,
    check_table_format,
    duplicate_line,
    failure_line,
    find_table,
    input_summary,
    load_json_list,
    load_yaml_model,
    prepare_step_dir,
    require_tables,
    resolve_file,
    structure_inputs,
    write_summary,
    write_table,
)
from pandora.datasets import (
    curate_structure,
    dataset_chain_records,
    find_contained,
)
from pandora.export import read_records, structure_to_mmcif, write_json
from pandora.ingestion import mirror_path
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.provenance import (
    build_step_summary,
    canonicalisation_entry_record,
    curation_entry_record,
)
from pandora.schemas.canonicalisation import (
    ChainCopyRecord,
    canonicalisationPolicy,
)
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


def _skip_too_large(
    path: Path,
    policy: canonicalisationPolicy,
    counts: Counter[str],
    failures: TextIO,
) -> bool:
    """Record a `TooLarge` failure and return True when the entry's
    estimated atom count is above `policy.max_atoms`."""

    if policy.max_atoms is None:
        return False
    n_atoms = estimate_atom_count(path, policy)
    if n_atoms <= policy.max_atoms:
        return False
    counts["entries_failed"] += 1
    counts["entries_too_large"] += 1
    failures.write(
        failure_line(
            path,
            "TooLarge",
            f"estimated {n_atoms} atoms > max_atoms {policy.max_atoms}",
        )
    )
    return True


def _cmd_canonicalise(args: argparse.Namespace) -> None:
    """Handle `structures canonicalise`: canonicalise each input file and
    write the step folder."""

    command = "structures canonicalise"
    started_at = now_iso()
    source, out = Path(args.structures), Path(args.out)
    check_out(command, out, source)
    policy = load_yaml_model(command, canonicalisationPolicy, Path(args.policy))
    paths, previous = structure_inputs(command, source)

    check_table_format(command, args.table_format)
    prepare_step_dir(command, out, args.overwrite)
    counts: Counter[str] = Counter(
        entries_failed=0,
        **({} if policy.max_atoms is None else {"entries_too_large": 0}),
    )
    entries, copies = [], []
    seen: set[str] = set()
    with (
        (out / "failures.jsonl").open("w") as failures,
        (out / "mappings.jsonl").open("w") as mappings_file,
    ):
        for path in paths:
            counts["entries_in"] += 1
            try:
                if _skip_too_large(path, policy, counts, failures):
                    continue
                parsed, _, _ = mmcif_to_structure(str(path))
                if parsed is None:
                    raise ValueError(f"could not parse {path}")
                stem = parsed.entry_id.lower()
                if stem in seen:
                    counts["entries_failed"] += 1
                    failures.write(duplicate_line(path, parsed.entry_id))
                    continue
                seen.add(stem)
                canon, mappings, prov = canonicalise_structure(parsed, policy)
                written = mirror_path(Path("structures"), stem)
                structure_to_mmcif(canon, out / written)
            except ENTRY_ERRORS as exc:
                counts["entries_failed"] += 1
                failures.write(failure_line(path, type(exc).__name__, str(exc)))
                continue
            entries.append(canonicalisation_entry_record(canon, prov, written))
            copies += chain_copy_records(canon, mappings)
            line = {"entry_id": canon.entry_id}
            line.update(mappings.model_dump(mode="json"))
            mappings_file.write(json.dumps(line) + "\n")

    counts["entries_kept"] = len(entries)
    write_table(entries, out, "entries", args.table_format)
    write_table(copies, out, "chain_copies", args.table_format)
    write_summary(
        out,
        build_step_summary(
            "canonicalisation",
            policy,
            input=source.resolve(),
            output=out,
            started_at=started_at,
            counts=dict(sorted(counts.items())),
            previous_summary=previous,
        ),
    )
    failed = counts["entries_failed"]
    print(
        f"canonicalised {len(entries)} structures -> {out}"
        + (f", {failed} failed" if failed else "")
    )


def _read_copies(folder: Path) -> dict[str, dict[str, str]]:
    """Entry id -> {copy chain id: source chain id}, from a structures
    folder's `chain_copies` table (empty if it has none)."""

    copies: dict[str, dict[str, str]] = defaultdict(dict)
    table = find_table(folder, "chain_copies") if folder.is_dir() else None
    for row in read_records(table, ChainCopyRecord) if table else []:
        copies[row.entry_id][row.chain_id] = row.source_chain_id
    return copies


def _curate_one(
    structure: Structure,
    policy: DatasetCurationPolicy,
    copies: dict[str, dict[str, str]],
    out: Path,
) -> tuple[list[ExclusionRecord], list[DatasetChainRecord], CurationEntryRecord | None]:
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
    """Handle `structures curate`: curate each structure and write the
    step folder."""

    command = "structures curate"
    started_at = now_iso()
    source, out = Path(args.structures), Path(args.out)
    check_out(command, out, source)
    policy = load_yaml_model(command, DatasetCurationPolicy, Path(args.policy))
    paths, previous = structure_inputs(command, source)
    copies = _read_copies(source)

    check_table_format(command, args.table_format)
    prepare_step_dir(command, out, args.overwrite)
    counts: Counter[str] = Counter(entries_failed=0)
    entries: list[CurationEntryRecord] = []
    chains: list[DatasetChainRecord] = []
    exclusions: list[ExclusionRecord] = []
    seen: set[str] = set()
    with (out / "failures.jsonl").open("w") as failures:
        for path in paths:
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
                dropped, rows, entry = _curate_one(structure, policy, copies, out)
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
    write_summary(
        out,
        build_step_summary(
            "curation",
            policy,
            input=source.resolve(),
            output=out,
            started_at=started_at,
            counts=dict(sorted(counts.items())),
            previous_summary=previous,
        ),
    )
    failed = counts["entries_failed"]
    excluded = counts["entries_in"] - len(entries) - failed
    print(
        f"curated: {len(entries)} retained, {excluded} excluded -> {out}"
        + (f", {failed} failed" if failed else "")
    )


_CONTAINED_RANK_KEYS = ("resolution", "chain_count")


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


def _cmd_drop_contained(args: argparse.Namespace) -> None:
    """Handle `structures drop-contained`: drop entries contained in
    another and write the step folder."""

    command = "structures drop-contained"
    started_at = now_iso()
    source, out = Path(args.structures), Path(args.out)
    check_out(command, out, source, Path(args.clusters))
    rules = (
        load_yaml_model(command, DeduplicationRules, Path(args.policy))
        if args.policy
        else DeduplicationRules()
    )
    for rank_key in rules.contained_rank:
        if rank_key.key not in _CONTAINED_RANK_KEYS:
            raise SystemExit(
                f"{command}: unknown contained_rank key {rank_key.key!r}; "
                "use " + ", ".join(_CONTAINED_RANK_KEYS)
            )
    check_table_format(command, args.table_format)
    entries_table, chains_table = require_tables(command, source)
    entries = read_records(entries_table, CurationEntryRecord)
    chains = read_records(chains_table, DatasetChainRecord)
    clusters_path = resolve_file(command, Path(args.clusters), "clusters.json")
    clusters = load_json_list(command, SimilarityCluster, clusters_path)

    counts, unclustered = _cluster_counts(chains, clusters)
    for entry in entries:
        counts.setdefault(entry.entry_id, Counter())
    values = {e.entry_id: {k: getattr(e, k) for k in _CONTAINED_RANK_KEYS} for e in entries}
    contained = find_contained(counts, values, rules.contained_rank)
    removed = set() if args.disable else {r.entry_id for r in contained}
    kept = [e for e in entries if e.entry_id not in removed]
    kept_chains = [c for c in chains if c.entry_id not in removed]
    files = [e.structure_path for e in kept if e.structure_path]
    missing = next((f for f in files if not (source / f).is_file()), None)
    if missing is not None:
        raise SystemExit(f"{command}: {source} has no structure file {missing}")
    prepare_step_dir(command, out, args.overwrite)
    for file in files:
        dest = out / file
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / file, dest)

    write_table(kept, out, "entries", args.table_format)
    write_table(kept_chains, out, "chains", args.table_format)
    write_table(contained, out, "contained", args.table_format)
    (out / "failures.jsonl").write_text("")
    write_summary(
        out,
        build_step_summary(
            "dedup",
            rules,
            input=source.resolve(),
            output=out,
            started_at=started_at,
            counts={
                "chains_kept": len(kept_chains),
                "chains_unclustered": unclustered,
                "entries_contained": len(contained),
                "entries_in": len(entries),
                "entries_kept": len(kept),
            },
            previous_summary=input_summary(source),
        ),
    )
    if unclustered:
        print(
            f"{command}: warning: {unclustered} chains are in no cluster of "
            f"{clusters_path} (each counted as its own cluster); were the "
            "clusters built from this folder?",
            file=sys.stderr,
        )
    print(f"{len(kept)} retained, {len(contained)} contained -> {out}")


def _cmd_convert(args: argparse.Namespace) -> None:
    """Handle `structures convert`: write one mmCIF file as mmCIF or
    JSON, by the output's suffix."""

    structure, _, status = mmcif_to_structure(args.input)
    if structure is None:
        raise SystemExit(f"structures convert: failed to parse {args.input}: {status}")
    output = Path(args.output)
    if output.suffix == ".json":
        write_json(structure, output)
    else:
        structure_to_mmcif(structure, output)
    print(f"converted -> {output}")


def _add_structures(p: argparse.ArgumentParser, what: str) -> None:
    """Add the positional `structures` input."""

    p.add_argument("structures", help=what)


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add the `structures` group to the `pandora` parser."""

    group = subparsers.add_parser(
        "structures",
        help="Transform structures: canonicalise, curate, drop contained entries, convert a file.",
    )
    commands = group.add_subparsers(dest="structures_command", required=True)
    anything = (
        "structures folder, folder of .cif/.cif.gz files, entries.txt list or one .cif/.cif.gz file"
    )

    p = commands.add_parser(
        "canonicalise",
        help="Normalise chain IDs, residues, altlocs and more with a "
        "policy; write structures and tables.",
    )
    _add_structures(p, anything)
    add_out(p)
    p.add_argument("--policy", required=True, help="canonicalisationPolicy YAML")
    p.set_defaults(func=_cmd_canonicalise)

    p = commands.add_parser(
        "curate",
        help="Apply quality, organism and content rules; write entry and chain tables.",
    )
    _add_structures(p, anything)
    add_out(p)
    p.add_argument("--policy", required=True, help="DatasetCurationPolicy YAML")
    p.set_defaults(func=_cmd_curate)

    p = commands.add_parser(
        "drop-contained",
        help="Drop entries whose chains are all contained in another entry's.",
    )
    _add_structures(p, "structures folder with tables (from curate)")
    add_out(p)
    p.add_argument(
        "--clusters",
        required=True,
        help="chain clusters: a clusters.json or the folder holding it",
    )
    p.add_argument("--policy", help="DeduplicationRules YAML (contained_rank)")
    p.add_argument(
        "--disable",
        action="store_true",
        help="record contained entries but keep them",
    )
    p.set_defaults(func=_cmd_drop_contained)

    p = commands.add_parser(
        "convert",
        help="Convert one mmCIF file to mmCIF or JSON, by the output's suffix.",
    )
    p.add_argument("input", help="mmCIF file to read")
    p.add_argument("output", help="file to write (.json for JSON)")
    p.set_defaults(func=_cmd_convert)
