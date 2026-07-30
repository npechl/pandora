from __future__ import annotations

import hashlib
import json

from pydantic import BaseModel


def compute_checksum(model: BaseModel) -> str:
    """SHA-256 hex digest of a pydantic model's canonical JSON form.

    Serialises with sorted keys and no whitespace so the digest is
    stable regardless of field-construction order.

    Args:
        model: Any pydantic model (e.g. a `Structure`).

    Returns:
        A lowercase 64-character hex digest.
    """

    payload = json.dumps(
        model.model_dump(mode="json"), sort_keys=True, separators=(",", ":")
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
