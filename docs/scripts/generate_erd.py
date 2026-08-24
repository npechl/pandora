"""Regenerate the entity-relationship diagrams embedded in
docs/architecture.md, straight from the current pandora.schemas models.

Run after changing any pydantic model in pandora/schemas/:

    uv run --extra docs python docs/scripts/generate_erd.py

Requires Graphviz (the `dot` binary) on PATH, e.g. `brew install graphviz`.
"""

from __future__ import annotations

from pathlib import Path

import erdantic as erd

from pandora.schemas.canonicalisation import canonicalisationPolicy
from pandora.schemas.dataset import (
    ChainRecord,
    DatasetCurationPolicy,
    InterfaceRecord,
    ResidueRecord,
)
from pandora.schemas.provenance import DatasetManifest
from pandora.schemas.similarity import SimilarityCluster, SimilarityRelationship
from pandora.schemas.structure import Structure

OUTPUT_DIR = Path(__file__).parent.parent / "assets" / "diagrams"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    erd.create(Structure, ChainRecord, ResidueRecord, InterfaceRecord).draw(
        OUTPUT_DIR / "structure.svg"
    )
    erd.create(canonicalisationPolicy).draw(
        OUTPUT_DIR / "canonicalisation-policy.svg"
    )
    erd.create(DatasetCurationPolicy).draw(
        OUTPUT_DIR / "curation-policy.svg"
    )
    erd.create(
        DatasetManifest,
        SimilarityRelationship,
        SimilarityCluster,
        terminal_models=[canonicalisationPolicy, DatasetCurationPolicy],
    ).draw(OUTPUT_DIR / "provenance-manifest.svg")

    print(f"wrote 4 diagrams to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
