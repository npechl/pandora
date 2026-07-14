from __future__ import annotations

from datetime import datetime, timezone

from pandora.schemas.structure import (
    AsymRecord,
    AssemblyRecord,
    AtomSiteRecord,
    EntityRecord,
    Structure,
)
from pandora.schemas.canonicalisation import (
    canonicalisationPolicy,
    CanonicalMappings,
    canonicalisationProvenance,
)

from pandora.schemas.common import DiagnosticBundle

from pandora.canonicalisation.altlocs import _resolve_altlocs
from pandora.canonicalisation.assemblies import _normalize_assemblies
from pandora.canonicalisation.chain_ids import (
    _normalize_chain_ids,
    _apply_chain_map,
)
from pandora.canonicalisation.entities import _normalize_entities
from pandora.canonicalisation.residues import _normalize_residue_numbering
from pandora.canonicalisation.missing_data import (
    _handle_missing_atoms,
    _handle_missing_residues,
    _handle_incomplete_chains,
)
from pandora.canonicalisation.ligands import _filter_ligands
from pandora.canonicalisation.validation import _validate

def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonicalise_structure(
    structure: Structure,
    policy: canonicalisationPolicy,
) -> tuple[Structure, CanonicalMappings, canonicalisationProvenance]:
    """Convert a parsed structure into a canonical structure."""
    ir = policy.identifier_rules
    mdr = policy.missing_data_rules
    ar = policy.altloc_rules
    asmr = policy.assembly_rules
    er = policy.entity_rules
    lr = policy.ligand_rules
    vr = policy.validation_rules
    pr = policy.provenance_rules

    record = pr.record_original_mappings
    transforms: list[str] = []
    diagnostics = DiagnosticBundle()

    atoms: list[AtomSiteRecord] = list(structure.atoms)
    asym_units: list[AsymRecord] = list(structure.asym_units)
    entities: list[EntityRecord] = list(structure.entities)
    assemblies: list[AssemblyRecord] = list(structure.assemblies)

    # normalize_chain_ids --------------------------------
    chain_map, chain_id_mapping = _normalize_chain_ids(
        asym_units, ir.chain_id.strategy, record
    )
    atoms, asym_units, structure = _apply_chain_map(
        atoms, asym_units, structure, chain_map
    )
    assemblies = list(structure.assemblies)
    if ir.chain_id.strategy != "preserve":
        transforms.append(f"chain_id:{ir.chain_id.strategy}")

    # normalize_residue_numbering ---------------------------------
    atoms, residue_number_mapping = _normalize_residue_numbering(
        atoms,
        chain_map,
        ir.residue_numbering.strategy,
        ir.residue_numbering.preserve_insertion_codes,
        record,
    )
    if ir.residue_numbering.strategy != "preserve":
        transforms.append(f"residue_numbering:{ir.residue_numbering.strategy}")

    # normalize_assemblies --------------------------------
    assemblies, assembly_mapping = _normalize_assemblies(
        assemblies, asmr, ir.assembly_id.strategy, record
    )
    if (
        ir.assembly_id.strategy != "preserve"
        or asmr.strategy != "preserve_as_reported"
    ):
        transforms.append(
            f"assembly:{asmr.strategy}/id:{ir.assembly_id.strategy}"
        )

    # handle_missing_atoms ---------------------------------------
    atoms = _handle_missing_atoms(
        atoms, mdr.missing_atoms, diagnostics, structure.entry_id
    )
    if mdr.missing_atoms.strategy not in ("preserve",):
        transforms.append(f"missing_atoms:{mdr.missing_atoms.strategy}")

    # handle_missing_residues ---------------------------------------
    atoms = _handle_missing_residues(
        atoms, mdr.missing_residues, diagnostics, structure.entry_id
    )
    if mdr.missing_residues.strategy not in ("preserve",):
        transforms.append(f"missing_residues:{mdr.missing_residues.strategy}")

    # handle_incomplete_chains -------------------------------------
    atoms, asym_units = _handle_incomplete_chains(
        atoms, asym_units, mdr.incomplete_chains
    )
    if mdr.incomplete_chains.strategy != "preserve":
        transforms.append(f"incomplete_chains:{mdr.incomplete_chains.strategy}")

    # esolve_altlocs ------------------------------
    atoms, altloc_selection_mapping = _resolve_altlocs(atoms, ar)
    if ar.strategy != "preserve":
        transforms.append(f"altloc:{ar.strategy}")

    # normalize_entities ---------------------------------------
    entities, asym_units, atoms, entity_mapping = _normalize_entities(
        entities, asym_units, atoms, er, record
    )
    if er.strategy != "preserve":
        transforms.append(f"entity:{er.strategy}")

    # filter_ligands ---------------------------------------
    atoms, asym_units = _filter_ligands(
        atoms, asym_units, entities, lr, diagnostics, structure.entry_id
    )
    if lr.strategy != "preserve":
        transforms.append(f"ligands:{lr.strategy}")

    # validate_canonical_structure ---------------------------------------
    _validate(atoms, asym_units, vr, diagnostics, structure.entry_id)

    canonical = structure.model_copy(
        update={
            "atoms": atoms,
            "asym_units": asym_units,
            "entities": entities,
            "assemblies": assemblies,
        }
    )

    mappings = CanonicalMappings(
        chain_id_mapping=chain_id_mapping,
        residue_number_mapping=residue_number_mapping,
        assembly_mapping=assembly_mapping,
        entity_mapping=entity_mapping,
        altloc_selection_mapping=altloc_selection_mapping,
    )

    provenance = canonicalisationProvenance(
        canonicalised_at=_now_iso(),
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
        transforms=transforms,
        report={
            "warnings": len(diagnostics.warnings),
            "errors": len(diagnostics.errors),
        }
        if pr.emit_canonicalisation_report
        else {},
    )

    return canonical, mappings, provenance
