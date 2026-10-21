"""The `pandora` command: assembles the subcommands from one module per
library area (see `pandora/cli/README.md` for the I/O convention)."""

from __future__ import annotations

import argparse

from pandora.cli import (
    annotate,
    canonicalisation,
    datasets,
    export,
    ingestion,
    provenance,
    similarity,
)

# Order of the subcommands in `pandora --help`.
_MODULES = (
    ingestion,
    canonicalisation,
    datasets,
    similarity,
    annotate,
    provenance,
    export,
)


def _build_parser() -> argparse.ArgumentParser:
    """Build the `pandora` CLI's argparse parser with all subcommands
    registered."""

    parser = argparse.ArgumentParser(
        prog="pandora",
        description="Pandora command-line tools for building ML-ready structure datasets.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    for module in _MODULES:
        module.register(subparsers)
    return parser


def main(argv: list[str] | None = None) -> None:
    """Parse CLI arguments and dispatch to the selected subcommand's
    handler."""

    parser = _build_parser()
    args = parser.parse_args(argv)
    args.func(args)


def app() -> None:
    """Pandora CLI entry point (the `pandora` console script).

    Parses arguments from `sys.argv` and dispatches to the selected
    subcommand's handler (`fetch`, `canonicalise`, `curate`, `dedup`,
    `similarity`, `cluster`, `partition`, `annotate`, `manifest`,
    `reproduce`, or `export`); see `pandora/cli/README.md` for the
    I/O convention between subcommands. Takes no arguments and returns
    nothing — each subcommand handler writes its output to disk and
    prints a summary.
    """

    main()
