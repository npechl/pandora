from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel


def write_json(model: BaseModel, path: str | Path) -> Path:
    """Write any Pandora pydantic model to disk as pretty-printed JSON.

    Covers every stage that returns a single model — canonicalisation
    provenance/mappings, MetadataRecord, AnnotationLayer,
    ProvenanceBundle, etc.
    """

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(model.model_dump_json(indent=2), encoding="utf-8")
    return out_path
