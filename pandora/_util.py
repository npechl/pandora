from __future__ import annotations

import importlib
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from types import ModuleType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pandora.schemas.common import RankKey


def now_iso() -> str:
    """Current UTC time as an ISO 8601 string."""

    return datetime.now(timezone.utc).isoformat()


def require(module: str, extra: str) -> ModuleType:
    """Import an optional dependency, or raise a RuntimeError naming the
    pyproject extra that installs it."""

    try:
        return importlib.import_module(module)
    except ImportError as exc:
        raise RuntimeError(
            f"{module} is required here (install the `{extra}` extra)"
        ) from exc


def rank_sort_key(
    item_id: str,
    values: Mapping[str, float | int | None],
    rank: Sequence[RankKey],
) -> tuple:
    """Sort key ranking an item by `rank` (smaller is better): nulls
    last, then the lowest `item_id`. KeyError if a key is missing."""

    parts: list[tuple[bool, float]] = []
    for rank_key in rank:
        value = values[rank_key.key]
        if value is None:
            parts.append((True, 0))
        else:
            parts.append((False, value if rank_key.order == "min" else -value))
    return (*parts, item_id)
