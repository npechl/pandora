"""CLI subcommands for C01 — mmCIF Ingestion."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.table import Table

from pandora.c01_ingestion import ingest_list_mmcif, ingest_mmcif
from pandora.cli._writers import OutputFormat, write_ingestion_result
from pandora.schemas.c01_ingestion import (
    FetchOptions,
    MmCIFBatchInput,
    MmCIFIngestionInput,
)

app = typer.Typer(help="C01 — mmCIF ingestion from PDBe, RCSB, or local files.")
console = Console()


def _provider_callback(value: str) -> str:
    allowed = {"pdbe", "pdb", "local", "custom"}
    if value not in allowed:
        raise typer.BadParameter(f"must be one of {sorted(allowed)}")
    return value


def _format_callback(value: str) -> str:
    allowed = {"parquet", "json", "jsonl"}
    if value not in allowed:
        raise typer.BadParameter(f"must be one of {sorted(allowed)}")
    return value


@app.command("ingest")
def ingest(
    entry_ids: Annotated[list[str], typer.Argument(help="One or more PDB entry IDs (e.g. 1cbs 4hhb).")],
    provider: Annotated[str, typer.Option("--provider", "-p",
        help="Data provider: pdbe | pdb | local | custom.",
        callback=_provider_callback)] = "pdbe",
    source_uri: Annotated[Optional[str], typer.Option("--source-uri",
        help="Explicit URL or file path (required for --provider local/custom).")] = None,
    output: Annotated[Path, typer.Option("--output", "-o",
        help="Root output directory. Component results go into <output>/c01/.")] = Path("pandora_output"),
    fmt: Annotated[str, typer.Option("--format", "-f",
        help="Tabular output format: parquet (default) | json | jsonl.",
        callback=_format_callback)] = "parquet",
    no_cache: Annotated[bool, typer.Option("--no-cache",
        help="Bypass the on-disk mmCIF cache.")] = False,
) -> None:
    """Fetch, parse, and validate one or more mmCIF entries and write results to disk."""
    fetch_opts = FetchOptions(use_cache=not no_cache)

    inputs = [
        MmCIFIngestionInput(
            entry_id=eid,
            provider=provider,
            source_uri=source_uri,
            fetch_options=fetch_opts,
        )
        for eid in entry_ids
    ]

    batch = MmCIFBatchInput(entries=inputs)

    with console.status(f"Ingesting {len(entry_ids)} entr{'y' if len(entry_ids) == 1 else 'ies'}…"):
        batch_result = ingest_list_mmcif(batch)

    # ── per-entry results table ───────────────────────────────────────────────
    table = Table(title="C01 Ingestion Results", show_lines=False)
    table.add_column("Entry",    style="bold cyan",  no_wrap=True)
    table.add_column("Status",   no_wrap=True)
    table.add_column("Chains",   justify="right")
    table.add_column("Residues", justify="right")
    table.add_column("Atoms",    justify="right")
    table.add_column("Cached",   justify="center")
    table.add_column("Output",   style="dim")

    for result in batch_result.results:
        status_style = {
            "success": "green",
            "warning": "yellow",
            "failed":  "red",
        }.get(result.status, "white")

        if result.parsed_structure:
            ps = result.parsed_structure
            entry_dir = write_ingestion_result(result, output, fmt)  # type: ignore[arg-type]
            table.add_row(
                result.entry_id,
                f"[{status_style}]{result.status}[/{status_style}]",
                str(len(ps.chains)),
                str(len(ps.residues)),
                str(len(ps.atoms)),
                "✓" if result.provenance.from_cache else "–",
                str(entry_dir),
            )
        else:
            err = result.diagnostics.errors[0].message if result.diagnostics.errors else "–"
            table.add_row(
                result.entry_id,
                f"[{status_style}]{result.status}[/{status_style}]",
                "–", "–", "–", "–",
                err[:60],
            )

    console.print(table)

    # ── batch summary JSON ────────────────────────────────────────────────────
    summary_dir = output / "c01"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_path = summary_dir / "summary.json"
    summary_path.write_text(
        json.dumps({
            "total":   batch_result.summary.total,
            "success": batch_result.summary.success,
            "warning": batch_result.summary.warning,
            "failed":  batch_result.summary.failed,
            "format":  fmt,
            "entries": [
                {
                    "entry_id": r.entry_id,
                    "status":   r.status,
                    "from_cache": r.provenance.from_cache,
                    "atoms":    len(r.parsed_structure.atoms)    if r.parsed_structure else None,
                    "residues": len(r.parsed_structure.residues) if r.parsed_structure else None,
                    "chains":   len(r.parsed_structure.chains)   if r.parsed_structure else None,
                }
                for r in batch_result.results
            ],
        }, indent=2),
        encoding="utf-8",
    )

    s = batch_result.summary
    console.print(
        f"\nSummary: [green]{s.success} success[/green]  "
        f"[yellow]{s.warning} warning[/yellow]  "
        f"[red]{s.failed} failed[/red]  "
        f"— written to [bold]{output}/c01/[/bold]"
    )
