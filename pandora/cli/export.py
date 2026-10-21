"""CLI `export`: convert one mmCIF file to mmCIF or JSON."""

from __future__ import annotations

import argparse
from pathlib import Path

from pandora.export import structure_to_mmcif, write_json
from pandora.parsing import mmcif_to_structure


def _cmd_export(args: argparse.Namespace) -> None:
    """Handle the `export` subcommand: convert one mmCIF file to mmCIF
    or JSON."""

    structure, _, status = mmcif_to_structure(args.input)
    if structure is None:
        raise SystemExit(f"failed to parse {args.input}: {status}")

    output = Path(args.output)
    if output.suffix == ".json":
        write_json(structure, output)
    else:
        structure_to_mmcif(structure, output)
    print(f"exported -> {output}")


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    p = subparsers.add_parser(
        "export", help="Convert one mmCIF file to mmCIF or JSON."
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_export)
