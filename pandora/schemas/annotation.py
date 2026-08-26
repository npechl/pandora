from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field

AnnotationScope = Literal[
    "entry",
    "entity",
    "chain",
    "residue",
    "ligand",
    "interface",
    "pair",
    "cluster",
    "dataset",
]


class AnnotationLayer(BaseModel):
    """One derived annotation result (counts, contacts, identity, etc.)
    for an entry, pair, or other scope.

    Attributes:
        layer_name: A human-readable name for the layer.
        layer_type: A machine-readable type identifying the layer
            (e.g. "structure_counts").
        scope: What the layer describes: an entry, a pair, an
            interface, etc.
        method: The identifier of the method/algorithm used to
            compute the layer.
        target_ids: The entry id(s) the layer applies to.
        parameters: The parameters the layer was computed with.
        data: The layer's actual computed output.
        provenance: Additional provenance details about how the layer
            was computed.
    """

    layer_name: str
    layer_type: str
    scope: AnnotationScope
    method: str
    target_ids: list[str] = Field(default_factory=list)
    parameters: dict[str, Any] = Field(default_factory=dict)
    data: dict[str, Any] = Field(default_factory=dict)
    provenance: dict[str, Any] = Field(default_factory=dict)
