"""CLI `entries fetch` and `entries scan`: get mmCIF files and record
where they came from."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pandora._util import now_iso
from pandora.cli._io import (
    SUMMARY,
    add_out,
    prepare_output_dir,
    write_summary,
)
from pandora.ingestion import (
    read_entry_paths,
    template_entry_paths,
    write_entry_paths,
)
from pandora.provenance import build_step_summary
from pandora.schemas.ingestion import (
    FetchOptions,
    IngestionDiscovery,
    IngestionProvenance,
)

_PROVENANCE = "ingestion_provenance.json"


def _read_ids(command: str, path: Path) -> list[str]:
    """Entry IDs of `path`, one per line; blanks and `#` lines skipped."""

    if not path.is_file():
        raise SystemExit(f"{command}: {path} is not a file")
    lines = (line.strip() for line in path.read_text().splitlines())
    return [line for line in lines if line and not line.startswith("#")]


def _entry_key(path: Path) -> str:
    """The upper-case entry ID of an mmCIF file path."""

    return path.name.split(".")[0].upper()


def _previous(out: Path) -> tuple[dict[str, Path], dict[str, object]]:
    """The entry paths and provenance a fetch folder already holds."""

    listed = out / "entries.txt"
    paths = (
        {_entry_key(p): p for p in read_entry_paths(listed)}
        if listed.is_file()
        else {}
    )
    recorded = out / _PROVENANCE
    try:
        provenance = (
            json.loads(recorded.read_text()) if recorded.is_file() else {}
        )
    except ValueError as exc:
        raise SystemExit(f"entries fetch: {recorded}: {exc}") from exc
    return paths, provenance


def _cmd_fetch(args: argparse.Namespace) -> None:
    """Handle `entries fetch`: download mmCIF files into OUT, adding to
    what it already holds."""

    from pandora.ingestion.mmcif import fetch_mmcif

    started_at = now_iso()
    ids = list(args.ids)
    if args.ids_file:
        ids += _read_ids("entries fetch", Path(args.ids_file))
    ids = list(dict.fromkeys(i.upper() for i in ids))
    if not ids:
        raise SystemExit("entries fetch: give entry IDs or --ids FILE")

    out = Path(args.out)
    options = FetchOptions(
        use_cache=not args.no_cache, allow_partial=args.allow_partial
    )
    fetched: dict[str, IngestionProvenance] = {}
    for entry_id in ids:
        try:
            fetched[entry_id] = fetch_mmcif(
                entry_id, args.provider, None, out, options
            )
        except (RuntimeError, ValueError) as exc:
            if not args.allow_partial:
                raise SystemExit(f"entries fetch: {entry_id}: {exc}") from exc
            print(
                f"warning: failed to fetch {entry_id}: {exc}", file=sys.stderr
            )

    paths, provenance = _previous(out)
    root = out.resolve()
    for entry_id in fetched:
        gz = root / f"{entry_id.lower()}.cif.gz"
        paths[entry_id] = (
            gz if gz.exists() else root / f"{entry_id.lower()}.cif"
        )
    provenance.update({k: v.model_dump() for k, v in fetched.items()})

    out.mkdir(parents=True, exist_ok=True)
    (out / _PROVENANCE).write_text(json.dumps(provenance, indent=2))
    write_entry_paths(sorted(paths.values()), out / "entries.txt")
    write_summary(
        out,
        build_step_summary(
            "fetch",
            options,
            input=f"provider:{args.provider}",
            output=out,
            started_at=started_at,
            counts={
                "fetched": len(fetched),
                "failed": len(ids) - len(fetched),
                "entries_total": len(paths),
            },
        ),
    )
    print(f"fetched {len(fetched)}/{len(ids)} entries -> {out}")


def _scan_paths(args: argparse.Namespace, root: Path) -> list[Path]:
    """The mmCIF files `entries scan` lists, from `--glob` or
    `--path-template`, restricted to `--ids` when given."""

    ids = _read_ids("entries scan", Path(args.ids)) if args.ids else None
    if args.glob:
        paths = sorted(root.glob(args.glob))
        if ids is not None:
            wanted = {i.upper() for i in ids}
            paths = [p for p in paths if _entry_key(p) in wanted]
        return paths
    try:
        return template_entry_paths(root, args.path_template, ids)
    except ValueError as exc:
        raise SystemExit(f"entries scan: {exc}") from exc


def _cmd_scan(args: argparse.Namespace) -> None:
    """Handle `entries scan`: list an archive's mmCIF files into
    `entries.txt` and write `summary.json`."""

    started_at = now_iso()
    root = Path(args.root).resolve()
    if not root.is_dir():
        raise SystemExit(f"entries scan: {root} is not a folder")
    out = Path(args.out)
    paths = _scan_paths(args, root)
    prepare_output_dir(
        "entries scan", out, ["entries.txt", SUMMARY], args.overwrite
    )

    write_entry_paths(paths, out / "entries.txt")
    discovery = IngestionDiscovery(
        input_dir=str(root),
        glob=args.glob,
        path_template=args.path_template,
        ids_file=args.ids,
        source_label=args.source_label,
        release_date=args.release_date,
    )
    write_summary(
        out,
        build_step_summary(
            "ingestion",
            discovery,
            input=root,
            output=out,
            started_at=started_at,
            counts={"entries_found": len(paths)},
        ),
    )
    print(f"scanned {len(paths)} entries -> {out}")


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add the `entries` group to the `pandora` parser."""

    group = subparsers.add_parser(
        "entries", help="Get mmCIF files: download them or list an archive."
    )
    commands = group.add_subparsers(dest="entries_command", required=True)

    p = commands.add_parser(
        "fetch",
        help="Download mmCIF files into OUT, adding to what it holds.",
    )
    p.add_argument("ids", nargs="*", help="entry IDs to fetch")
    p.add_argument("out", help="output folder; also the download cache")
    p.add_argument(
        "--ids",
        dest="ids_file",
        help="file of entry IDs, one per line, # comments",
    )
    p.add_argument(
        "--provider",
        choices=["pdbe", "pdb"],
        default="pdbe",
        help="where to download from (default pdbe)",
    )
    p.add_argument(
        "--no-cache",
        action="store_true",
        help="download again even if the file is already in OUT",
    )
    p.add_argument(
        "--allow-partial",
        action="store_true",
        help="warn about entries that fail and carry on",
    )
    p.set_defaults(func=_cmd_fetch)

    p = commands.add_parser(
        "scan",
        help="List an archive's mmCIF files into entries.txt, without "
        "opening them.",
    )
    p.add_argument("root", help="the archive's root folder")
    add_out(p, tables=False)
    find = p.add_mutually_exclusive_group(required=True)
    find.add_argument(
        "--glob", help="glob pattern relative to ROOT, e.g. '*/*.cif'"
    )
    find.add_argument(
        "--path-template",
        help="path relative to ROOT with {id} and {mid}, e.g. "
        "'{mid}/{id}/clean_mmcif/{id}_updated.cif.gz'; lists only the "
        "folders above {id}",
    )
    p.add_argument(
        "--ids",
        help="file of entry IDs, one per line: fills --path-template "
        "instead of listing, or keeps only these entries of --glob",
    )
    p.add_argument("--source-label", help="free-text label for the source")
    p.add_argument("--release-date", help="the source's release date")
    p.set_defaults(func=_cmd_scan)
