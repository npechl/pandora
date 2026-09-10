from pathlib import Path

import gzip

import gemmi
import httpx

from pandora._util import now_iso
from pandora.ingestion.cache import mtime_iso, resolve_cache_hit
from pandora.schemas.ingestion import (
    FetchOptions,
    IngestionProvenance,
)

_PROVIDER_URLS: dict[str, str] = {
    "pdbe": "https://www.ebi.ac.uk/pdbe/entry-files/download/{id}_updated.cif",
    "pdb": "https://files.rcsb.org/download/{id}.cif",
}

_PDBE_SUMMARY_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{id}"

_NULL_CIF = frozenset({".", "?"})


def _revision_date_from_cif_file(path: Path) -> str | None:
    """Latest `_pdbx_audit_revision_history.revision_date` in an mmCIF
    file (RCSB files carry this category; PDBe's don't), or None if the
    category is absent or the file can't be parsed. Reads gzip-compressed
    files transparently.
    """

    try:
        block = gemmi.cif.read(str(path)).sole_block()
        dates = list(
            block.find_loop("_pdbx_audit_revision_history.revision_date")
        )
    except (RuntimeError, ValueError):
        return None
    return dates[-1] if dates and dates[-1] not in _NULL_CIF else None


def _pdbe_revision_date(entry_id: str) -> str | None:
    """Revision date for entry_id via PDBe's summary API (its mmCIF
    download has no revision-history category of its own), or None on
    any failure — this enrichment must never fail the fetch itself.
    """

    try:
        resp = httpx.get(
            _PDBE_SUMMARY_URL.format(id=entry_id.lower()), timeout=30.0
        )
        resp.raise_for_status()
        entry = resp.json()[entry_id.lower()][0]
        raw_date = entry["revision_date"]
        return f"{raw_date[0:4]}-{raw_date[4:6]}-{raw_date[6:8]}"
    except (httpx.HTTPError, KeyError, IndexError, ValueError, TypeError):
        return None


def fetch_mmcif(
    entry_id: str,
    provider: str,
    source_uri: str | None,
    output_dir: Path,
    fetch_options: FetchOptions | None = None,
) -> IngestionProvenance:
    """Fetch a raw mmCIF file from a provider URL, write it to
    output_dir, return provenance.

    Resolves the download URL from `source_uri` or from `provider`
    ("pdbe"/"pdb"), reuses a cached file when one is fresh (per
    `fetch_options`), otherwise downloads it, optionally decompresses
    it, and writes it to `output_dir`.

    Args:
        entry_id: The PDB/PDBe entry identifier (e.g. "1abc").
        provider: One of "pdbe" or "pdb", used to build the download
            URL when `source_uri` is not given.
        source_uri: An explicit URL to fetch from, overriding
            `provider`'s default URL template.
        output_dir: Directory to write the downloaded (or cached) file
            into; created if missing.
        fetch_options: Caching/decompression options. Defaults
            to `FetchOptions()` when not given.

    Returns:
        `IngestionProvenance` describing the provider, source URL,
        retrieval timestamp, whether the result came from cache, and
        the entry's own revision date where the provider exposes one.

    Raises:
        ValueError: If `provider` is not "pdbe"/"pdb" and no
            `source_uri` was given.
        RuntimeError: If the HTTP request fails, or the response
            cannot be decompressed/decoded.
    """

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
                revision_date=(
                    _revision_date_from_cif_file(cached)
                    if provider == "pdb"
                    else None
                ),
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

    if provider == "pdb":
        revision_date = _revision_date_from_cif_file(out_path)
    elif provider == "pdbe":
        revision_date = _pdbe_revision_date(entry_id)
    else:
        revision_date = None

    return IngestionProvenance(
        provider=provider,
        source_uri=url,
        retrieved_at=now_iso(),
        from_cache=False,
        revision_date=revision_date,
    )


def fetch_list_mmcif(
    entry_ids: list[str],
    provider: str,
    output_dir: Path,
    source_uri: str | None = None,
    fetch_options: FetchOptions | None = None,
) -> list[IngestionProvenance]:
    """Fetch a list of raw mmCIF files from a provider URL, write them to
    output_dir, return provenance for each.

    Calls `fetch_mmcif()` once per entry ID. If
    `fetch_options.allow_partial` is True, entries that fail to fetch
    (`RuntimeError` or `ValueError`) are skipped instead of aborting
    the whole batch.

    Args:
        entry_ids: The PDB/PDBe entry identifiers to fetch.
        provider: One of "pdbe" or "pdb", used to build the download
            URL when `source_uri` is not given.
        source_uri: A URL/path template overriding `provider`'s default
            URL template, applied to every entry. Must contain an
            `{id}` placeholder (same convention as `_PROVIDER_URLS`),
            formatted per entry_id — e.g.
            `"https://example.org/mirror/{id}.cif"`.
        output_dir: Directory to write the downloaded (or cached)
            files into; created if missing.
        fetch_options: Caching/decompression options, including
            `allow_partial`. Defaults to `FetchOptions()` when not
            given.

    Returns:
        A list of `IngestionProvenance`, one per successfully fetched
        entry, in the same order as `entry_ids` (with failed entries
        omitted when `allow_partial` is True).

    Raises:
        ValueError: If `provider` is invalid and no `source_uri` was
            given, and `allow_partial` is False.
        RuntimeError: If a fetch fails and `allow_partial` is False.
    """

    fetch_options = fetch_options or FetchOptions()

    provenance_list = []
    for entry_id in entry_ids:
        try:
            provenance = fetch_mmcif(
                entry_id=entry_id,
                provider=provider,
                source_uri=(
                    source_uri.format(id=entry_id) if source_uri else None
                ),
                output_dir=output_dir,
                fetch_options=fetch_options,
            )
        except (RuntimeError, ValueError):
            if fetch_options.allow_partial:
                continue
            raise
        provenance_list.append(provenance)

    return provenance_list


def ingest_local_mmcif(
    path: Path,
    source_uri: str | None = None,
) -> IngestionProvenance:
    """Describe an already-downloaded mmCIF file as provenance, without
    fetching or copying anything.

    For a file that arrived via some external bulk source — a local
    rsync mirror of PDB, an archived dated snapshot, a lab file server —
    rather than `fetch_mmcif()`. Pandora has no bulk-download step of its
    own; this just lets such a file's origin be recorded the same way a
    live fetch's is, so it can still be attached to a `ProvenanceBundle`.

    Args:
        path: Path to the existing mmCIF file on disk. Not modified or
            copied.
        source_uri: A label for where this file came from (e.g.
            "pdb_snapshot_2024-01-29"), stored verbatim. Defaults to
            `str(path)` when not given.

    Returns:
        `IngestionProvenance` with `provider="local"`, `from_cache=False`
        (this isn't Pandora's own fetch cache), the file's own mtime as
        `retrieved_at`, and `revision_date` read from the file's
        `_pdbx_audit_revision_history` category when present.

    Raises:
        ValueError: If `path` doesn't exist.
    """

    if not path.exists():
        raise ValueError(f"path={path} does not exist")

    return IngestionProvenance(
        provider="local",
        source_uri=source_uri or str(path),
        retrieved_at=mtime_iso(path),
        from_cache=False,
        revision_date=_revision_date_from_cif_file(path),
    )
