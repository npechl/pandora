from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any, TypeVar

from pydantic import BaseModel

M = TypeVar("M", bound=BaseModel)


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


def read_records(path: str | Path, model: type[M]) -> list[M]:
    """Read records written by `write_records` back as models.

    Dispatches on the file suffix, as `write_records` does. Parquet null
    values come back as None (not NaN) and list columns as lists.

    Args:
        path: The `.jsonl`, `.json` or `.parquet` file to read.
        model: The record model to validate each row against.

    Returns:
        One `model` per row, in file order.

    Raises:
        ValueError: The suffix is unsupported.
        RuntimeError: `.parquet` requested without pandas/pyarrow installed.
    """

    in_path = Path(path)
    suffix = in_path.suffix.lower()

    if suffix == ".jsonl":
        lines = in_path.read_text(encoding="utf-8").splitlines()
        return [model.model_validate_json(line) for line in lines if line]

    if suffix == ".json":
        rows = json.loads(in_path.read_text(encoding="utf-8"))
        return [model.model_validate(row) for row in rows]

    if suffix == ".parquet":
        try:
            import pandas as pd
        except ImportError as exc:
            raise RuntimeError(
                "reading .parquet requires pandas/pyarrow "
                "(install the `export` extra)"
            ) from exc
        table = pd.read_parquet(in_path)
        return [
            model.model_validate(
                {k: _plain(v, pd.isna) for k, v in row.items()}
            )
            for row in table.to_dict("records")
        ]

    raise ValueError(f"unsupported suffix {suffix!r}")


def _plain(value: Any, isna: Any) -> Any:
    """A Parquet cell as plain Python: numpy values to Python, NaN to None."""

    if hasattr(value, "tolist"):  # numpy scalars and arrays (list columns)
        value = value.tolist()
    if isinstance(value, list):
        return value
    return None if isna(value) else value
