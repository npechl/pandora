from typing import Literal

from pydantic import BaseModel

StaleBehavior = Literal["use_stale", "warn", "fail"]


class FetchOptions(BaseModel):
    allow_partial: bool = False
    use_cache: bool = True
    decompress: bool = True
    max_age_seconds: int | None = None
    stale_behavior: StaleBehavior = "use_stale"


class IngestionProvenance(BaseModel):
    provider: str
    source_uri: str | None = None
    retrieved_at: str | None = None
    from_cache: bool = False
