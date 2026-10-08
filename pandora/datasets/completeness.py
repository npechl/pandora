from __future__ import annotations

from collections import defaultdict

from pandora.canonicalisation.missing_data import (
    _PROTEIN_BACKBONE,
    _backbone_by_entity,
)
from pandora.canonicalisation.modified_residues import STANDARD_RESIDUES
from pandora.schemas.common import Diagnostic, DiagnosticBundle
from pandora.schemas.dataset import ChainCompleteness, MissingResidueDefinition
from pandora.schemas.structure import Structure


def _seqres_by_entity(structure: Structure) -> dict[str, dict[int, set[str]]]:
    """mon_ids per SEQRES position, per entity, from `_entity_poly_seq`."""

    seqres: dict[str, dict[int, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for row in structure.raw.get("_entity_poly_seq", []):
        num, entity_id = row.get("num"), row.get("entity_id")
        if entity_id and num and num.isdigit():
            seqres[entity_id][int(num)].add(row.get("mon_id") or "")
    return seqres


def _matches_seqres(
    residues: dict[int, tuple[str, set[str]]],
    positions: dict[int, set[str]],
) -> bool:
    """Whether every residue sits on a SEQRES position with its comp_id."""

    for seq_id, (comp_id, _) in residues.items():
        mon_ids = positions.get(seq_id)
        if mon_ids is None:
            return False
        if comp_id not in mon_ids and mon_ids <= STANDARD_RESIDUES:
            return False
    return True


def chain_completeness(
    structure: Structure,
    definition: MissingResidueDefinition = "incomplete_backbone",
) -> tuple[list[ChainCompleteness], DiagnosticBundle]:
    """Count the missing residues of each polymer chain against SEQRES.

    SEQRES is `_entity_poly_seq`; a residue lines up with it by
    `label_seq_id`. A chain that has no SEQRES, or whose residues don't
    line up (for example after renumbering), gets no record and a
    `NO_SEQRES` or `SEQRES_MISMATCH` warning instead.

    Args:
        structure: The structure to measure. Its `label_seq_id`s must
            still index SEQRES.
        definition: `"unobserved"` counts a residue as missing when it
            has no atoms; `"incomplete_backbone"` also when any backbone
            atom (N, CA, C, O for protein) is missing.

    Returns:
        `(records, diagnostics)`: one `ChainCompleteness` per measured
        polymer chain, and a warning per chain that could not be
        measured.
    """

    diagnostics = DiagnosticBundle()
    seqres = _seqres_by_entity(structure)
    backbone_by_entity = _backbone_by_entity(structure.entities)
    polymer_entities = {e.id for e in structure.entities if e.type == "polymer"}

    # chain -> seq_id -> (comp_id, atom names)
    residues: dict[str, dict[int, tuple[str, set[str]]]] = defaultdict(dict)
    for a in structure.atoms:
        if a.label_seq_id is None:
            continue
        chain = residues[a.label_asym_id]
        if a.label_seq_id not in chain:
            chain[a.label_seq_id] = (a.label_comp_id, set())
        chain[a.label_seq_id][1].add(a.label_atom_id)

    records: list[ChainCompleteness] = []
    for asym in structure.asym_units:
        if asym.entity_id not in polymer_entities:
            continue
        positions = seqres.get(asym.entity_id, {})
        chain_residues = residues.get(asym.id, {})
        code = None
        if not positions:
            code = "NO_SEQRES"
        elif not _matches_seqres(chain_residues, positions):
            code = "SEQRES_MISMATCH"
        if code is not None:
            diagnostics.warnings.append(
                Diagnostic(
                    code=code,
                    severity="warning",
                    message=f"Chain {asym.id}: completeness not measured "
                    f"({code})",
                    entry_id=structure.entry_id,
                    context={"chain": asym.id},
                )
            )
            continue

        length = len(positions)
        backbone = backbone_by_entity.get(asym.entity_id, _PROTEIN_BACKBONE)
        present = sorted(
            seq_id
            for seq_id, (_, names) in chain_residues.items()
            if definition == "unobserved" or backbone <= names
        )
        if present:
            first, last = present[0], present[-1]
            n_term, c_term = first - 1, length - last
            middle = (last - first + 1) - len(present)
        else:
            n_term, c_term, middle = length, 0, 0
        records.append(
            ChainCompleteness(
                entry_id=structure.entry_id,
                chain_id=asym.id,
                seqres_length=length,
                missing_n_term=n_term,
                missing_c_term=c_term,
                missing_middle=middle,
            )
        )
    return records, diagnostics
