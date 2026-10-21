"""CLI `curate` and `dedup`: filter and deduplicate a directory of structures."""

from __future__ import annotations

import argparse
from pathlib import Path

from pandora.cli._io import (
    load_yaml_model,
    read_structures_dir,
    write_json_dict,
    write_structures_dir,
)
from pandora.datasets import (
    curate_structure,
    deduplicate_structures,
)
from pandora.export import write_json, write_records
from pandora.metadata import collect_metadata
from pandora.schemas.dataset import (
    DatasetCurationPolicy,
    DeduplicationRules,
    ExclusionRecord,
)
from pandora.schemas.structure import Structure


def _cmd_curate(args: argparse.Namespace) -> None:
    """Handle the `curate` subcommand: apply quality/organism/content
    rules to a directory of structures."""

    output_dir = Path(args.output_dir)
    policy = load_yaml_model(DatasetCurationPolicy, Path(args.policy))
    structures = read_structures_dir(Path(args.input_dir))

    retained: dict[str, Structure] = {}
    provenance = {}
    exclusions: list[ExclusionRecord] = []
    for entry_id, structure in structures.items():
        metadata = collect_metadata(structure)
        curated, records, prov = curate_structure(structure, metadata, policy)
        provenance[entry_id] = prov
        exclusions.extend(records)
        if curated is not None:
            retained[entry_id] = curated

    write_structures_dir(retained, output_dir)
    write_json_dict(provenance, output_dir / "curation_provenance.json")
    if exclusions:
        write_records(exclusions, output_dir / "curation_exclusions.json")
    print(
        f"curated: {len(retained)} retained, "
        f"{len(structures) - len(retained)} excluded -> {output_dir}"
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
        "curate", help="Apply quality/organism/content rules to a directory."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument("--policy", required=True, help="DatasetCurationPolicy YAML")
    p.add_argument("--output-dir", required=True)
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
