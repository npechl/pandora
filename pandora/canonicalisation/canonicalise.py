from __future__ import annotations

from pandora._util import now_iso
from pandora.schemas.structure import (
    AsymRecord,
    AssemblyRecord,
    AtomSiteRecord,
    ConnRecord,
    EntityRecord,
    SSRecord,
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
from pandora.canonicalisation.entities import (
    _normalize_entities,
    _remap_raw_entity_ids,
)
from pandora.canonicalisation.residues import _normalize_residue_numbering
from pandora.canonicalisation.missing_data import (
    _handle_missing_atoms,
    _handle_missing_residues,
    _handle_incomplete_chains,
)
from pandora.canonicalisation.ligands import filter_ligands
from pandora.canonicalisation.modified_residues import (
    _map_connections,
    _map_entity_poly_seq,
    _map_modified_residues,
)
from pandora.canonicalisation.raw_categories import _drop_stale_raw_categories
from pandora.canonicalisation.validation import _validate


def _transforms(policy: canonicalisationPolicy) -> list[str]:
    """Labels of the steps whose strategy changes the structure, in the
    order canonicalise_structure runs them."""

    ir = policy.identifier_rules
    mdr = policy.missing_data_rules
    asm = policy.assembly_rules.strategy
    asm_id = ir.assembly_id.strategy
    assembly = (
        "preserve"
        if (asm, asm_id) == ("preserve_as_reported", "preserve")
        else f"{asm}/id:{asm_id}"
    )
    steps = [
        ("modified_residues", policy.modified_residue_rules.strategy),
        ("chain_id", ir.chain_id.strategy),
        ("assembly", assembly),
        ("missing_atoms", mdr.missing_atoms.strategy),
        ("missing_residues", mdr.missing_residues.strategy),
        ("incomplete_chains", mdr.incomplete_chains.strategy),
        ("altloc", policy.altloc_rules.strategy),
        ("residue_numbering", ir.residue_numbering.strategy),
        ("entity", policy.entity_rules.strategy),
        ("ligands", policy.ligand_rules.strategy),
    ]
    return [
        f"{step}:{strategy}"
        for step, strategy in steps
        if strategy != "preserve"
    ]


def canonicalise_structure(
    structure: Structure,
    policy: canonicalisationPolicy,
) -> tuple[Structure, CanonicalMappings, canonicalisationProvenance]:
    """Convert a parsed structure into a canonical structure.

    Applies modified-residue mapping, chain ID normalization, assembly
    normalization, missing-data
    handling, altloc resolution, residue renumbering, entity
    normalization, and ligand filtering, in that order, according to the
    given policy. Each
    step records its mapping (for reversibility) and, if the step's
    strategy deviates from "preserve", appends a transform label to the
    returned provenance.

    Args:
        structure: The parsed input structure to canonicalise.
        policy: The canonicalisation policy whose rules govern each
            normalization step.

    Returns:
        A tuple of:
            - The canonicalised `Structure`.
            - `CanonicalMappings` recording original-to-canonical
              identifier mappings for chains, residues, assemblies,
              entities, altloc selections, and modified residues.
            - `canonicalisationProvenance` describing which policy was
              applied, which transforms ran, and diagnostic counts.
    """
    ir = policy.identifier_rules
    mdr = policy.missing_data_rules
    ar = policy.altloc_rules
    asmr = policy.assembly_rules
    er = policy.entity_rules
    lr = policy.ligand_rules
    mrr = policy.modified_residue_rules
    vr = policy.validation_rules
    pr = policy.provenance_rules

    record = pr.record_original_mappings
    diagnostics = DiagnosticBundle()

    original_atoms = structure.atoms
    original_raw = structure.raw
    atoms: list[AtomSiteRecord] = list(structure.atoms)
    asym_units: list[AsymRecord] = list(structure.asym_units)
    entities: list[EntityRecord] = list(structure.entities)
    assemblies: list[AssemblyRecord] = list(structure.assemblies)
    connections: list[ConnRecord] = list(structure.connections)
    secondary_structure: SSRecord = structure.secondary_structure

    # map_modified_residues --------------------------------
    # First, so chain ids and seq ids still match
    # _pdbx_struct_mod_residue.
    atoms, modified_residue_mapping = _map_modified_residues(
        atoms,
        structure.raw.get("_pdbx_struct_mod_residue", []),
        mrr,
        diagnostics,
        structure.entry_id,
    )
    # Bonds and SEQRES name the same residues, so they follow. The bonds
    # are read back from `structure` after the chain map below.
    structure = structure.model_copy(
        update={
            "connections": _map_connections(
                structure.connections, modified_residue_mapping
            )
        }
    )
    original_raw = _map_entity_poly_seq(
        original_raw, modified_residue_mapping, asym_units
    )

    # normalize_chain_ids --------------------------------
    chain_map, chain_id_mapping = _normalize_chain_ids(
        asym_units, ir.chain_id.strategy, record
    )
    atoms, asym_units, structure = _apply_chain_map(
        atoms, asym_units, structure, chain_map
    )
    assemblies = list(structure.assemblies)
    connections = list(structure.connections)
    secondary_structure = structure.secondary_structure

    # normalize_assemblies --------------------------------
    (
        assemblies,
        atoms,
        asym_units,
        connections,
        secondary_structure,
        assembly_mapping,
    ) = _normalize_assemblies(
        assemblies,
        atoms,
        asym_units,
        connections,
        secondary_structure,
        asmr,
        ir.assembly_id.strategy,
        record and asmr.record_original_assembly_mapping,
        diagnostics,
        structure.entry_id,
    )

    # handle_missing_atoms ---------------------------------------
    atoms = _handle_missing_atoms(
        atoms, mdr.missing_atoms, diagnostics, structure.entry_id, entities
    )

    # handle_missing_residues ---------------------------------------
    atoms = _handle_missing_residues(
        atoms, mdr.missing_residues, diagnostics, structure.entry_id
    )

    # handle_incomplete_chains -------------------------------------
    atoms, asym_units = _handle_incomplete_chains(
        atoms, asym_units, mdr.incomplete_chains
    )

    # resolve_altlocs ------------------------------
    atoms, altloc_selection_mapping = _resolve_altlocs(atoms, ar)

    # normalize_residue_numbering ---------------------------------
    # After the gap checks (renumbering makes every chain contiguous)
    # and altloc resolution (so compositional disorder gets one number).
    atoms, residue_number_mapping = _normalize_residue_numbering(
        atoms,
        chain_map,
        ir.residue_numbering.strategy,
        ir.residue_numbering.preserve_insertion_codes,
        record,
    )

    # normalize_entities ---------------------------------------
    entities, asym_units, atoms, entity_mapping, entity_id_map = (
        _normalize_entities(entities, asym_units, atoms, er, record)
    )

    # filter_ligands ---------------------------------------
    atoms, asym_units = filter_ligands(
        atoms, asym_units, entities, lr, diagnostics, structure.entry_id
    )

    # validate_canonical_structure ---------------------------------------
    status = _validate(atoms, asym_units, vr, diagnostics, structure.entry_id)
    if status == "failed":
        raise ValueError(
            f"Validation of canonicalised structure failed for entry "
            f"{structure.entry_id!r}"
        )

    # Verbatim categories that index chains or residues would contradict
    # the new atoms, so they go (and are listed in the provenance).
    raw, dropped_raw_categories = _drop_stale_raw_categories(
        original_raw, original_atoms, atoms
    )
    # Verbatim `entity_id` columns follow the canonical entity ids.
    raw = _remap_raw_entity_ids(raw, entity_id_map)

    canonical = structure.model_copy(
        update={
            "raw": raw,
            "atoms": atoms,
            "asym_units": asym_units,
            "entities": entities,
            "assemblies": assemblies,
            "connections": connections,
            "secondary_structure": secondary_structure,
        }
    )

    mappings = CanonicalMappings(
        chain_id_mapping=chain_id_mapping,
        residue_number_mapping=residue_number_mapping,
        assembly_mapping=assembly_mapping,
        entity_mapping=entity_mapping,
        altloc_selection_mapping=altloc_selection_mapping,
        modified_residue_mapping=modified_residue_mapping,
    )

    provenance = canonicalisationProvenance(
        canonicalised_at=now_iso(),
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
        transforms=_transforms(policy),
        dropped_raw_categories=dropped_raw_categories,
        report=(
            {
                "warnings": len(diagnostics.warnings),
                "errors": len(diagnostics.errors),
            }
            if pr.emit_canonicalisation_report
            else {}
        ),
    )

    return canonical, mappings, provenance
