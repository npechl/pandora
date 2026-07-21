from __future__ import annotations

import yaml

from pandora.schemas.canonicalisation import canonicalisationPolicy


def load_policy(path: str) -> canonicalisationPolicy:
    with open(path) as stream:
        data = yaml.safe_load(stream)

    return canonicalisationPolicy.model_validate(data)
