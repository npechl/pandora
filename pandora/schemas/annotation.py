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
    layer_name: str
    layer_type: str
    scope: AnnotationScope
    method: str
    target_ids: list[str] = Field(default_factory=list)
    data: dict[str, Any] = Field(default_factory=dict)
    provenance: dict[str, Any] = Field(default_factory=dict)
