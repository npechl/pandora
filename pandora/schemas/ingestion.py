from typing import Literal

from pydantic import BaseModel

StaleBehavior = Literal["use_stale", "warn", "fail"]


class FetchOptions(BaseModel):
    """Caching and decompression options for fetching one mmCIF file.

    Attributes:
        allow_partial: Whether a batch fetch continues past individual
            failures instead of raising.
        use_cache: Whether an existing cached file is reused instead
            of refetching.
        decompress: Whether a gzip-compressed response is decompressed
            before writing.
        max_age_seconds: How old a cached file may be before it's
            considered stale; `None` means never stale.
        stale_behavior: What to do with a stale cached file: use it
            anyway, warn and refetch, or raise.
    """

    allow_partial: bool = False
    use_cache: bool = True
    decompress: bool = True
    max_age_seconds: int | None = None
    stale_behavior: StaleBehavior = "use_stale"


class IngestionProvenance(BaseModel):
    """Record of where one structure's raw mmCIF file came from and
    when it was retrieved.

    Attributes:
        provider: Which provider the file was fetched from ("pdbe" or
            "pdb").
        source_uri: The exact URL the file was fetched from.
        retrieved_at: When the file was retrieved, as an ISO 8601
            timestamp.
        from_cache: Whether the file came from the local cache rather
            than a fresh download.
        revision_date: The entry's own most recent revision date (ISO
            "YYYY-MM-DD"), when the provider exposes one. Lets a caller
            detect drift on a later re-fetch of the same entry_id
            without Pandora computing a content checksum. `None` when
            the provider doesn't expose revision info or it couldn't
            be determined.
    """

    provider: str
    source_uri: str | None = None
    retrieved_at: str | None = None
    from_cache: bool = False
    revision_date: str | None = None
