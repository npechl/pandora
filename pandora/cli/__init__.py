"""Pandora CLI — entry point."""
from __future__ import annotations

import typer

from pandora.cli import c01

app = typer.Typer(
    name="pandora",
    help="Pandora — AI-ready structural biology dataset pipeline.",
    no_args_is_help=True,
)

app.add_typer(c01.app, name="c01")
