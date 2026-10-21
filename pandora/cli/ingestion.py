"""CLI `fetch` and `ingest`: get mmCIF files and record where they came from."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pandora.cli._io import (
    write_json_dict,
)
from pandora.schemas.ingestion import FetchOptions, IngestionProvenance


def _cmd_fetch(args: argparse.Namespace) -> None:
    """Handle the `fetch` subcommand: download raw mmCIF files and write
    ingestion provenance."""

    from pandora.ingestion.mmcif import fetch_mmcif

    output_dir = Path(args.output_dir)
    fetch_options = FetchOptions(use_cache=not args.no_cache)

    provenance: dict[str, IngestionProvenance] = {}
    for entry_id in args.entry_ids:
        try:
            provenance[entry_id.upper()] = fetch_mmcif(
                entry_id, args.provider, None, output_dir, fetch_options
            )
        except (RuntimeError, ValueError) as exc:
            if not args.allow_partial:
                raise
            print(
                f"warning: failed to fetch {entry_id}: {exc}", file=sys.stderr
            )

    write_json_dict(provenance, output_dir / "ingestion_provenance.json")
    print(
        f"fetched {len(provenance)}/{len(args.entry_ids)} entries -> {output_dir}"
    )


def _cmd_ingest(args: argparse.Namespace) -> None:
    """Handle the `ingest` subcommand: record ingestion provenance for a
    directory of already-downloaded mmCIF files (e.g. a local bulk
    mirror/snapshot), without fetching or copying anything."""

    from pandora.ingestion.mmcif import ingest_local_mmcif

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    provenance: dict[str, IngestionProvenance] = {}
    paths = sorted(input_dir.glob("*.cif"))
    for path in paths:
        try:
            provenance[path.stem.upper()] = ingest_local_mmcif(
                path, args.source_uri
            )
        except ValueError as exc:
            if not args.allow_partial:
                raise
            print(f"warning: failed to ingest {path}: {exc}", file=sys.stderr)

    write_json_dict(provenance, output_dir / "ingestion_provenance.json")
    print(f"ingested {len(provenance)}/{len(paths)} entries -> {output_dir}")


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    p = subparsers.add_parser("fetch", help="Fetch raw mmCIF files.")
    p.add_argument("entry_ids", nargs="+")
    p.add_argument("--provider", choices=["pdbe", "pdb"], default="pdbe")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--no-cache", action="store_true")
    p.add_argument("--allow-partial", action="store_true")
    p.set_defaults(func=_cmd_fetch)

    p = subparsers.add_parser(
        "ingest",
        help="Record ingestion provenance for an already-downloaded "
        "directory of mmCIF files (e.g. a local bulk mirror/snapshot).",
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument(
        "--source-uri",
        help="Label for where these files came from (e.g. "
        "'pdb_snapshot_2024-01-29'), applied to every entry. Defaults to "
        "each file's own path.",
    )
    p.add_argument("--allow-partial", action="store_true")
    p.set_defaults(func=_cmd_ingest)
