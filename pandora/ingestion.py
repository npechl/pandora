from pathlib import Path
from datetime import datetime, timezone

import gzip
import httpx

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


def fetch_pdb():
    """Fetch a raw mmCIF file from a provider URL or local path, with optional disk cache."""
    raise NotImplementedError(
        "fetch_pdb() is not supported; use fetch_mmcif() instead."
    )


def fetch_mmcif(
    entry_id: str,
    provider: str,
    source_uri: str | None,
    output_dir: Path,
    fetch_options: FetchOptions = FetchOptions(),
) -> IngestionProvenance:
    """Fetch a raw mmCIF file from a provider URL, write it to output_dir, return provenance."""

    if source_uri:
        url = source_uri
    elif provider in _PROVIDER_URLS:
        fmt_id = entry_id.lower() if provider == "pdbe" else entry_id.upper()
        url = _PROVIDER_URLS[provider].format(id=fmt_id)
    else:
        raise ValueError(
            f"provider={provider!r} requires an explicit source_uri or one of 'pdbe' or 'pdb'"
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
            raise RuntimeError(
                f"Failed to decode {'decompressed ' if fetch_options.decompress else ''}response from {url} as UTF-8"
            ) from exc
        out_path = output_dir / f"{entry_id.lower()}.cif"
        out_path.write_text(content, encoding="utf-8")

    return IngestionProvenance(
        provider=provider,
        source_uri=url,
        retrieved_at=_now_iso(),
        from_cache=False,
    )
