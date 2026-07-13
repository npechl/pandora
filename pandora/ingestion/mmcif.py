from pathlib import Path
from datetime import datetime, timezone

import gzip
import httpx

from pandora.ingestion.cache import mtime_iso, resolve_cache_hit
from pandora.schemas.ingestion import (
    FetchOptions,
    IngestionProvenance,
)

_PROVIDER_URLS: dict[str, str] = {
    "pdbe": "https://www.ebi.ac.uk/pdbe/entry-files/download/{id}_updated.cif",
    "pdb": "https://files.rcsb.org/download/{id}.cif",
}


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def fetch_mmcif(
    entry_id: str,
    provider: str,
    source_uri: str | None,
    output_dir: Path,
    fetch_options: FetchOptions | None = None,
) -> IngestionProvenance:
    """Fetch a raw mmCIF file from a provider URL, write it to
    output_dir, return provenance."""

    fetch_options = fetch_options or FetchOptions()

    if source_uri:
        url = source_uri
    elif provider in _PROVIDER_URLS:
        fmt_id = entry_id.lower() if provider == "pdbe" else entry_id.upper()
        url = _PROVIDER_URLS[provider].format(id=fmt_id)
    else:
        raise ValueError(
            f"provider={provider!r} requires an explicit source_uri "
            f"or one of 'pdbe' or 'pdb'"
        )

    if fetch_options.use_cache:
        cached = resolve_cache_hit(
            entry_id=entry_id,
            output_dir=output_dir,
            max_age_seconds=fetch_options.max_age_seconds,
            stale_behavior=fetch_options.stale_behavior,
        )
        if cached is not None:
            return IngestionProvenance(
                provider=provider,
                source_uri=url,
                retrieved_at=mtime_iso(cached),
                from_cache=True,
            )

    try:
        resp = httpx.get(url, follow_redirects=True, timeout=60.0)
        resp.raise_for_status()
    except httpx.HTTPStatusError as exc:
        raise RuntimeError(
            f"HTTP {exc.response.status_code} fetching {url}"
        ) from exc
    except httpx.RequestError as exc:
        raise RuntimeError(f"Network error fetching {url}: {exc}") from exc

    raw_bytes = resp.content
    is_gzipped = url.endswith(".gz") or raw_bytes[:2] == b"\x1f\x8b"

    if is_gzipped and fetch_options.decompress:
        try:
            raw_bytes = gzip.decompress(raw_bytes)
        except gzip.BadGzipFile as exc:
            raise RuntimeError(
                f"Failed to decompress response from {url}"
            ) from exc
        is_gzipped = False

    output_dir.mkdir(parents=True, exist_ok=True)

    if is_gzipped:
        out_path = output_dir / f"{entry_id.lower()}.cif.gz"
        out_path.write_bytes(raw_bytes)
    else:
        try:
            content = raw_bytes.decode("utf-8")
        except UnicodeDecodeError as exc:
            decompressed = "decompressed " if fetch_options.decompress else ""
            raise RuntimeError(
                f"Failed to decode {decompressed}response from {url} as UTF-8"
            ) from exc
        out_path = output_dir / f"{entry_id.lower()}.cif"
        out_path.write_text(content, encoding="utf-8")

    return IngestionProvenance(
        provider=provider,
        source_uri=url,
        retrieved_at=_now_iso(),
        from_cache=False,
    )


def fetch_list_mmcif(
    entry_ids: list[str],
    provider: str,
    source_uri: str | None,
    output_dir: Path,
    fetch_options: FetchOptions | None = None,
) -> list[IngestionProvenance]:
    """Fetch a list of raw mmCIF files from a provider URL, write them to
    output_dir, return provenance for each.

    If fetch_options.allow_partial is True, entries that fail to fetch are
    skipped instead of aborting the whole batch.
    """

    fetch_options = fetch_options or FetchOptions()

    provenance_list = []
    for entry_id in entry_ids:
        try:
            provenance = fetch_mmcif(
                entry_id=entry_id,
                provider=provider,
                source_uri=source_uri,
                output_dir=output_dir,
                fetch_options=fetch_options,
            )
        except (RuntimeError, ValueError):
            if fetch_options.allow_partial:
                continue
            raise
        provenance_list.append(provenance)

    return provenance_list
