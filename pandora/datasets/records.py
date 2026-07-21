from __future__ import annotations

from collections import defaultdict

from pandora.annotations.entry import (
    annotate_chain_interfaces,
    atoms_by_asym_id,
    polymer_asym_ids,
)
from pandora.schemas.dataset import ChainRecord, InterfaceRecord, ResidueRecord
from pandora.schemas.structure import AtomSiteRecord, Structure

DEFAULT_INTERFACE_CUTOFF = 4.0


def extract_chain_records(structure: Structure) -> list[ChainRecord]:
    """Build one ChainRecord per polymer chain (asym unit) in `structure`."""

    entities_by_id = {entity.id: entity for entity in structure.entities}
    atoms_by_chain = atoms_by_asym_id(structure)

    records = []
    for asym in structure.asym_units:
        entity = entities_by_id.get(asym.entity_id)
        if entity is None or entity.type != "polymer":
            continue

        chain_atoms = atoms_by_chain.get(asym.id, [])
        if not chain_atoms:
            continue

        residue_count = len(
            {(atom.label_seq_id, atom.auth_seq_id) for atom in chain_atoms}
        )
        sequence = None
        if entity.poly is not None:
            sequence = (
                entity.poly.pdbx_seq_one_letter_code_can
                or entity.poly.pdbx_seq_one_letter_code
            )

        records.append(
            ChainRecord(
                entry_id=structure.entry_id,
                chain_id=asym.id,
                auth_chain_id=asym.auth_id,
                entity_id=asym.entity_id,
                sequence=sequence,
                residue_count=residue_count,
                atom_count=len(chain_atoms),
            )
        )
    return records


def extract_residue_records(structure: Structure) -> list[ResidueRecord]:
    """Build one ResidueRecord per polymer residue in `structure`.

    Each record carries its full atom list (coordinates, B-factor, etc.)
    rather than a bare count, so residue-level ML use cases (contact
    maps, solvent exposure) don't need to go back to the structure.
    """

    polymer_ids = polymer_asym_ids(structure)

    residues: dict[
        tuple[str, int | None, str | None, str], list[AtomSiteRecord]
    ] = defaultdict(list)
    for atom in structure.atoms:
        if atom.label_asym_id not in polymer_ids:
            continue
        key = (
            atom.label_asym_id,
            atom.label_seq_id,
            atom.auth_seq_id,
            atom.label_comp_id,
        )
        residues[key].append(atom)

    return [
        ResidueRecord(
            entry_id=structure.entry_id,
            chain_id=chain_id,
            seq_id=seq_id,
            auth_seq_id=auth_seq_id,
            comp_id=comp_id,
            atoms=atoms,
        )
        for (chain_id, seq_id, auth_seq_id, comp_id), atoms in sorted(
            residues.items(), key=lambda item: (item[0][0], item[0][1] or 0)
        )
    ]


def extract_interface_records(
    structure: Structure,
    distance_cutoff: float = DEFAULT_INTERFACE_CUTOFF,
) -> list[InterfaceRecord]:
    """Build one InterfaceRecord per polymer chain pair in contact.

    Contact detection itself lives in
    `annotations.entry.annotate_chain_interfaces` — this function only
    reshapes that annotation layer's output into records.
    """

    layer = annotate_chain_interfaces(structure, distance_cutoff)
    return [
        InterfaceRecord(
            entry_id=structure.entry_id,
            chain_id_1=interface["chain_id_1"],
            chain_id_2=interface["chain_id_2"],
            distance_cutoff=distance_cutoff,
            interface_residues_chain_1=interface["interface_residues_chain_1"],
            interface_residues_chain_2=interface["interface_residues_chain_2"],
            contact_count=interface["contact_count"],
        )
        for interface in layer.data["interfaces"]
    ]
