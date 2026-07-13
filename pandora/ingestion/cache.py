from __future__ import annotations

import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

from pandora.schemas.ingestion import StaleBehavior

_CACHE_EXTENSIONS = (".cif", ".cif.gz")


class CacheStaleError(RuntimeError):
    """Raised when a cached file is stale and stale_behavior='fail'."""


def cache_path_candidates(entry_id: str, output_dir: Path) -> list[Path]:
    """Possible on-disk cache paths for entry_id (compressed and
    uncompressed), in lookup order."""
    stem = entry_id.lower()
    return [output_dir / f"{stem}{ext}" for ext in _CACHE_EXTENSIONS]


def find_cached(entry_id: str, output_dir: Path) -> Path | None:
    """First existing cached file for entry_id, if any."""
    for path in cache_path_candidates(entry_id, output_dir):
        if path.exists():
            return path
    return None


def is_stale(path: Path, max_age_seconds: int | None) -> bool:
    """True if path's mtime is older than max_age_seconds. A None
    max_age_seconds means cached files never go stale."""
    if max_age_seconds is None:
        return False
    age = time.time() - path.stat().st_mtime
    return age > max_age_seconds


def mtime_iso(path: Path) -> str:
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()


def resolve_cache_hit(
    entry_id: str,
    output_dir: Path,
    max_age_seconds: int | None,
    stale_behavior: StaleBehavior,
) -> Path | None:
    """Look up a usable cached file for entry_id, honoring max_age_seconds
    and stale_behavior.

    Returns the cache path to reuse, or None if a fresh fetch is required
    (no cached file, or stale with stale_behavior='warn'). Raises
    CacheStaleError if stale_behavior='fail' and the cached file is stale.
    """
    cached = find_cached(entry_id, output_dir)
    if cached is None:
        return None

    if not is_stale(cached, max_age_seconds):
        return cached

    if stale_behavior == "use_stale":
        return cached
    if stale_behavior == "warn":
        warnings.warn(
            f"cached file {cached} is older than max_age_seconds="
            f"{max_age_seconds}; refetching",
            stacklevel=2,
        )
        return None
    if stale_behavior == "fail":
        raise CacheStaleError(
            f"cached file {cached} is stale (max_age_seconds="
            f"{max_age_seconds}) and stale_behavior='fail'"
        )
    raise ValueError(f"unknown stale_behavior={stale_behavior!r}")
