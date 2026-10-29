from __future__ import annotations

from pathlib import Path
from typing import TypeVar

import yaml
from pydantic import BaseModel

from pandora.schemas.canonicalisation import canonicalisationPolicy

PolicyT = TypeVar("PolicyT", bound=BaseModel)


def load_policy(
    path: str | Path,
    model: type[PolicyT] = canonicalisationPolicy,
) -> PolicyT:
    """Load a policy from a YAML file.

    Args:
        path: Path to the policy YAML file.
        model: The policy model to validate against, e.g.
            `DatasetCurationPolicy`. Defaults to `canonicalisationPolicy`.

    Returns:
        The parsed, validated policy.
    """

    with open(path) as stream:
        data = yaml.safe_load(stream)

    return model.model_validate(data)
