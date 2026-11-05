"""CLI `canonicalise`: parse and canonicalise mmCIF files into a step
folder."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from pandora._util import now_iso
from pandora.canonicalisation import canonicalise_structure, chain_copy_records
from pandora.cli._io import (
    ENTRY_ERRORS,
    TABLE_FORMATS,
    check_table_format,
    duplicate_line,
    failure_line,
    prepare_step_dir,
    structure_files,
    write_table,
)
from pandora.export import structure_to_mmcif, write_json
from pandora.ingestion import load_policy, mirror_path, read_entry_paths
from pandora.parsing import mmcif_to_structure
from pandora.provenance import build_step_summary, canonicalisation_entry_record


def _cmd_canonicalise(args: argparse.Namespace) -> None:
    """Handle the `canonicalise` subcommand: canonicalise each input file
    and write the step folder."""

    started_at = now_iso()
    out = Path(args.output_dir)
    policy = load_policy(args.policy)
    previous = None
    if args.entries:
        source = Path(args.entries)
        paths = read_entry_paths(source)
        if (source.parent / "ingestion.json").exists():
            previous = source.parent / "ingestion.json"
    else:
        source = Path(args.input_dir)
        paths = structure_files(source)

    check_table_format(args.table_format)
    prepare_step_dir(out, "canonicalisation.json", args.overwrite)
    counts: Counter[str] = Counter(entries_failed=0)
    entries, copies = [], []
    seen: set[str] = set()
    with (
        (out / "failures.jsonl").open("w") as failures,
        (out / "mappings.jsonl").open("w") as mappings_file,
    ):
        for path in paths:
            counts["entries_in"] += 1
            try:
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
    write_json(
        build_step_summary(
            "canonicalisation",
            policy,
            input=source,
            output=out,
            started_at=started_at,
            counts=dict(sorted(counts.items())),
            previous_summary=previous,
        ),
        out / "canonicalisation.json",
    )
    failed = counts["entries_failed"]
    print(
        f"canonicalised {len(entries)} structures -> {out}"
        + (f", {failed} failed" if failed else "")
    )


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    p = subparsers.add_parser(
        "canonicalise",
        help="Canonicalise mmCIF files into a step folder.",
    )
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--input-dir", help="folder of mmCIF files (flat or mirror layout)"
    )
    source.add_argument(
        "--entries", help="text file listing mmCIF paths, one per line"
    )
    p.add_argument(
        "--policy", required=True, help="canonicalisationPolicy YAML"
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument(
        "--table-format",
        choices=TABLE_FORMATS,
        default="parquet",
        help="format of the entries/chain_copies tables",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="replace step outputs already in --output-dir",
    )
    p.set_defaults(func=_cmd_canonicalise)
