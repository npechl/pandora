"""CLI `fetch` and `ingest`: get mmCIF files and record where they came from."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pandora._util import now_iso
from pandora.cli._io import write_json_dict
from pandora.export import write_json
from pandora.ingestion import template_entry_paths, write_entry_paths
from pandora.provenance import build_step_summary
from pandora.schemas.ingestion import (
    FetchOptions,
    IngestionDiscovery,
    IngestionProvenance,
)


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
    root = output_dir.resolve()
    paths = []
    for entry_id in provenance:
        gz = root / f"{entry_id.lower()}.cif.gz"
        paths.append(gz if gz.exists() else root / f"{entry_id.lower()}.cif")
    write_entry_paths(sorted(paths), output_dir / "entries.txt")
    print(
        f"fetched {len(provenance)}/{len(args.entry_ids)} entries -> {output_dir}"
    )


def _cmd_ingest(args: argparse.Namespace) -> None:
    """Handle the `ingest` subcommand: list an archive's mmCIF files into
    `entries.txt` and write the `ingestion.json` step summary."""

    started_at = now_iso()
    if args.ids and not args.path_template:
        raise SystemExit("--ids needs --path-template")
    root = Path(args.input_dir).resolve()
    if not root.is_dir():
        raise SystemExit(f"{root} is not a folder")
    out = Path(args.output_dir)
    written = [out / "entries.txt", out / "ingestion.json"]
    existing = [p for p in written if p.exists()]
    if existing and not args.overwrite:
        raise SystemExit(
            f"{out} already holds {existing[0].name}; "
            "pass --overwrite to replace it"
        )

    if args.glob:
        paths = sorted(root.glob(args.glob))
    else:
        ids = None
        if args.ids:
            lines = Path(args.ids).read_text().splitlines()
            ids = [
                line.strip()
                for line in lines
                if line.strip() and not line.strip().startswith("#")
            ]
        try:
            paths = template_entry_paths(root, args.path_template, ids)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc

    write_entry_paths(paths, out / "entries.txt")
    discovery = IngestionDiscovery(
        input_dir=str(root),
        glob=args.glob,
        path_template=args.path_template,
        ids_file=args.ids,
        source_label=args.source_label,
        release_date=args.release_date,
    )
    write_json(
        build_step_summary(
            "ingestion",
            discovery,
            input=root,
            output=out,
            started_at=started_at,
            counts={"entries_found": len(paths)},
        ),
        out / "ingestion.json",
    )
    print(f"ingested {len(paths)} entries -> {out}")


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
        help="List an archive's mmCIF files into entries.txt, without "
        "opening them.",
    )
    p.add_argument("--input-dir", required=True, help="the archive's root")
    find = p.add_mutually_exclusive_group(required=True)
    find.add_argument(
        "--glob", help="glob pattern relative to --input-dir, e.g. '*/*.cif'"
    )
    find.add_argument(
        "--path-template",
        help="path relative to --input-dir with {id} and {mid}, e.g. "
        "'{mid}/{id}/clean_mmcif/{id}_updated.cif.gz'; lists only the "
        "folders above {id}",
    )
    p.add_argument(
        "--ids",
        help="file of entry IDs, one per line, to fill --path-template "
        "with instead of listing",
    )
    p.add_argument("--source-label", help="free-text label for the source")
    p.add_argument("--release-date", help="the source's release date")
    p.add_argument("--output-dir", required=True)
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="replace entries.txt/ingestion.json in --output-dir",
    )
    p.set_defaults(func=_cmd_ingest)
