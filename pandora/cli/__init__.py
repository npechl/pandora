"""Pandora CLI — entry point."""
from __future__ import annotations

import typer

from pandora.cli.c01 import ingest

app = typer.Typer(
    name="pandora",
    help="A toolkit for building AI-ready structural biology datasets.",
    no_args_is_help=True,
)


@app.callback()
def _root() -> None:  # noqa: D401 — prevents Typer from collapsing a single command into the root
    """A toolkit for building AI-ready structural biology datasets."""


app.command("ingest")(ingest)
