"""Output-writing utilities shared across CLI subcommands."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

from pandora.schemas.c01_ingestion import MmCIFIngestionResult

OutputFormat = Literal["parquet", "json", "jsonl"]


def _write_json(path: Path, data: object) -> None:
    path.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")


def write_ingestion_result(
    result: MmCIFIngestionResult,
    base_dir: Path,
    fmt: OutputFormat,
) -> Path:
    """Write one MmCIFIngestionResult under base_dir/c01/<entry_id>/."""
    entry_dir = base_dir / "c01" / result.entry_id.lower()
    entry_dir.mkdir(parents=True, exist_ok=True)

    # ── metadata (always JSON) ────────────────────────────────────────────────
    _write_json(entry_dir / "meta.json", {
        "entry_id": result.entry_id,
        "status": result.status,
        "provenance": result.provenance.model_dump(),
        "diagnostics": result.diagnostics.model_dump(),
    })

    if result.parsed_structure is None:
        return entry_dir

    ps = result.parsed_structure

    # Entities and assemblies have nested structures — always JSON
    _write_json(entry_dir / "entities.json",   [e.model_dump() for e in ps.entities])
    _write_json(entry_dir / "assemblies.json", [a.model_dump() for a in ps.assemblies])

    # Flatten residues and chains (strip nested lists before tabular write)
    atoms_rows    = [a.model_dump() for a in ps.atoms]
    residues_rows = [r.model_dump(exclude={"atoms"}) for r in ps.residues]
    chains_rows   = [c.model_dump(exclude={"residues"}) for c in ps.chains]

    if fmt == "parquet":
        import pandas as pd
        pd.DataFrame(atoms_rows).to_parquet(entry_dir / "atoms.parquet",    index=False)
        pd.DataFrame(residues_rows).to_parquet(entry_dir / "residues.parquet", index=False)
        pd.DataFrame(chains_rows).to_parquet(entry_dir / "chains.parquet",  index=False)

    elif fmt == "json":
        _write_json(entry_dir / "atoms.json",    atoms_rows)
        _write_json(entry_dir / "residues.json", residues_rows)
        _write_json(entry_dir / "chains.json",   chains_rows)

    elif fmt == "jsonl":
        for fname, rows in [
            ("atoms.jsonl",    atoms_rows),
            ("residues.jsonl", residues_rows),
            ("chains.jsonl",   chains_rows),
        ]:
            with open(entry_dir / fname, "w", encoding="utf-8") as fh:
                for row in rows:
                    fh.write(json.dumps(row, default=str) + "\n")

    return entry_dir
