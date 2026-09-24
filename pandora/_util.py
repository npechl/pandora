from __future__ import annotations

import importlib
from datetime import datetime, timezone
from types import ModuleType


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
