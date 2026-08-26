from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from pydantic import BaseModel


def write_records(records: Sequence[BaseModel], path: str | Path) -> Path:
    """Write a list of record models (Chain/Residue/InterfaceRecord, ...).

    Dispatches on the file suffix:
        .parquet — columnar, via pandas/pyarrow (the `export` extra)
        .json    — a JSON array of the records
        .jsonl   — one JSON object per line

    Args:
        records: The record models to write; must be non-empty.
        path: Destination file path; its suffix (`.json`/`.jsonl`/
            `.parquet`) selects the output format. Parent directories
            are created if missing.

    Returns:
        The resolved `Path` the records were written to.

    Raises:
        ValueError: `records` is empty, or the suffix is unsupported.
        RuntimeError: `.parquet` requested without pandas/pyarrow installed.
    """

    if not records:
        raise ValueError("no records to write")

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = out_path.suffix.lower()

    if suffix == ".jsonl":
        lines = (record.model_dump_json() for record in records)
        out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return out_path

    if suffix == ".json":
        body = ",".join(record.model_dump_json() for record in records)
        out_path.write_text(f"[{body}]", encoding="utf-8")
        return out_path

    if suffix == ".parquet":
        try:
            import pandas as pd
        except ImportError as exc:
            raise RuntimeError(
                "writing .parquet requires pandas/pyarrow "
                "(install the `export` extra)"
            ) from exc
        pd.DataFrame([record.model_dump() for record in records]).to_parquet(
            out_path
        )
        return out_path

    raise ValueError(f"unsupported record file suffix: {suffix!r}")
