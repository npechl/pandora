from __future__ import annotations

from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.structure import Structure

MISSING_VALUES = frozenset({"", ".", "?"})


def annotate_pairwise_sequence_identity(
    left: Structure,
    right: Structure,
) -> AnnotationLayer:
    """Compute a simple ungapped sequence identity between two entries.

    Compares every polymer entity in `left` against every polymer
    entity in `right` with a position-wise (ungapped, no alignment)
    identity, then reports the best-scoring pair.

    Args:
        left: The first structure to compare.
        right: The second structure to compare.

    Returns:
        An `AnnotationLayer` of type "pairwise_sequence_identity" whose
        `data` holds the best identity score, the best-matching entity
        pair, and every entity-pair comparison.
    """

    comparisons = []
    for left_entity_id, left_sequence in _entity_sequences(left):
        for right_entity_id, right_sequence in _entity_sequences(right):
            comparisons.append(
                {
                    "left_entity_id": left_entity_id,
                    "right_entity_id": right_entity_id,
                    **_ungapped_identity(left_sequence, right_sequence),
                }
            )

    best = max(
        comparisons,
        key=lambda item: item["identity"] or 0.0,
        default=None,
    )

    return AnnotationLayer(
        layer_name="Pairwise sequence identity",
        layer_type="pairwise_sequence_identity",
        scope="pair",
        method="pandora.basic.ungapped_entity_identity.v1",
        target_ids=[left.entry_id, right.entry_id],
        data={
            "best_identity": best["identity"] if best else None,
            "best_match": best,
            "comparisons": comparisons,
        },
        provenance={"inputs": ["EntityRecord.poly"]},
    )


def _entity_sequences(structure: Structure) -> list[tuple[str, str]]:
    sequences: list[tuple[str, str]] = []
    for entity in structure.entities:
        if entity.poly is None:
            continue
        sequence = _clean_sequence(
            entity.poly.pdbx_seq_one_letter_code_can
            or entity.poly.pdbx_seq_one_letter_code
        )
        if sequence:
            sequences.append((entity.id, sequence.upper()))
    return sequences


def _clean_sequence(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    if stripped in MISSING_VALUES:
        return None
    return "".join(stripped.split())


def _ungapped_identity(left: str, right: str) -> dict[str, float | int]:
    if not left or not right:
        return {
            "identity": 0.0,
            "aligned_length": 0,
            "left_length": len(left),
            "right_length": len(right),
            "coverage": 0.0,
        }

    aligned_length = min(len(left), len(right))
    denominator = max(len(left), len(right))
    matches = sum(
        1
        for left_char, right_char in zip(left, right)
        if left_char == right_char
    )
    return {
        "identity": matches / denominator,
        "aligned_length": aligned_length,
        "left_length": len(left),
        "right_length": len(right),
        "coverage": aligned_length / denominator,
    }
