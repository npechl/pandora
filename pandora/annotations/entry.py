from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any

from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.structure import AtomSiteRecord, Structure


WATER_COMP_IDS = frozenset({"HOH", "WAT", "DOD"})


def annotate_structure_counts(structure: Structure) -> AnnotationLayer:
    """Compute simple per-entry counts from the canonical structure.

    Args:
        structure: The canonical structure to summarize.

    Returns:
        An `AnnotationLayer` of type "structure_counts" whose `data`
        holds atom/residue/asym/entity/assembly/connection counts,
        secondary-structure record counts, entity-type and atom-group
        breakdowns, and the altloc atom count.
    """

    residue_keys = {
        (
            atom.label_asym_id,
            atom.label_seq_id,
            atom.auth_seq_id,
            atom.label_comp_id,
        )
        for atom in structure.atoms
    }
    entity_type_counts = Counter(entity.type for entity in structure.entities)
    atom_group_counts = Counter(atom.group_PDB for atom in structure.atoms)
    altloc_atom_count = sum(1 for atom in structure.atoms if atom.label_alt_id)

    return AnnotationLayer(
        layer_name="Structure counts",
        layer_type="structure_counts",
        scope="entry",
        method="pandora.basic.structure_counts.v1",
        target_ids=[structure.entry_id],
        data={
            "atom_count": len(structure.atoms),
            "residue_count": len(residue_keys),
            "asym_unit_count": len(structure.asym_units),
            "entity_count": len(structure.entities),
            "assembly_count": len(structure.assemblies),
            "connection_count": len(structure.connections),
            "secondary_structure": {
                "conf_records": len(structure.secondary_structure.conf_records),
                "sheet_strands": len(
                    structure.secondary_structure.sheet_strands
                ),
            },
            "entity_type_counts": dict(entity_type_counts),
            "atom_group_counts": dict(atom_group_counts),
            "altloc_atom_count": altloc_atom_count,
        },
        provenance={"inputs": ["Structure.atoms", "Structure.entities"]},
    )


def annotate_ligand_contacts(
    structure: Structure,
    distance_cutoff: float = 4.0,
    include_waters: bool = False,
) -> AnnotationLayer:
    """Compute polymer residues near non-polymer atoms within a cutoff.

    Groups HETATM records into ligands (by asym/auth_seq/comp id) and,
    for each ligand, finds polymer residues with at least one atom
    within `distance_cutoff` angstroms of a ligand atom.

    Args:
        structure: The structure to scan for ligand contacts.
        distance_cutoff: Contact distance in angstroms.
        include_waters: If True, treat water residues as ligands too.

    Returns:
        An `AnnotationLayer` of type "ligand_contacts" whose `data`
        holds the cutoff/flags used and, per ligand, the contacting
        polymer residues with their nearest distance.
    """

    polymer_atoms = [
        atom for atom in structure.atoms if atom.group_PDB == "ATOM"
    ]
    ligand_atoms = [
        atom
        for atom in structure.atoms
        if _is_ligand_atom(atom, include_waters)
    ]

    ligand_groups: dict[tuple[str, str | None, str], list[AtomSiteRecord]]
    ligand_groups = defaultdict(list)
    for atom in ligand_atoms:
        ligand_groups[
            (atom.label_asym_id, atom.auth_seq_id, atom.label_comp_id)
        ].append(atom)

    cutoff_sq = distance_cutoff * distance_cutoff
    contacts = []
    for ligand_key, atoms in ligand_groups.items():
        residue_contacts = _residues_within_cutoff(
            atoms, polymer_atoms, cutoff_sq
        )
        contacts.append(
            {
                "ligand_asym_id": ligand_key[0],
                "ligand_auth_seq_id": ligand_key[1],
                "ligand_comp_id": ligand_key[2],
                "contact_count": len(residue_contacts),
                "contacts": residue_contacts,
            }
        )

    return AnnotationLayer(
        layer_name="Ligand contacts",
        layer_type="ligand_contacts",
        scope="entry",
        method="pandora.basic.distance_cutoff_contacts.v1",
        target_ids=[structure.entry_id],
        data={
            "distance_cutoff": distance_cutoff,
            "include_waters": include_waters,
            "ligands": contacts,
        },
        provenance={"inputs": ["Structure.atoms"]},
    )


def _is_ligand_atom(atom: AtomSiteRecord, include_waters: bool) -> bool:
    if atom.group_PDB != "HETATM":
        return False
    if include_waters:
        return True
    return atom.label_comp_id.upper() not in WATER_COMP_IDS


def _residues_within_cutoff(
    ligand_atoms: list[AtomSiteRecord],
    polymer_atoms: list[AtomSiteRecord],
    cutoff_sq: float,
) -> list[dict[str, Any]]:
    contacts: dict[tuple[str, int | None, str], float] = {}
    for ligand_atom in ligand_atoms:
        for polymer_atom in polymer_atoms:
            dist_sq = _squared_distance(ligand_atom, polymer_atom)
            if dist_sq > cutoff_sq:
                continue

            key = (
                polymer_atom.label_asym_id,
                polymer_atom.label_seq_id,
                polymer_atom.label_comp_id,
            )
            current = contacts.get(key)
            if current is None or dist_sq < current:
                contacts[key] = dist_sq

    return [
        {
            "label_asym_id": asym_id,
            "label_seq_id": seq_id,
            "label_comp_id": comp_id,
            "distance": round(dist_sq**0.5, 3),
        }
        for (asym_id, seq_id, comp_id), dist_sq in sorted(contacts.items())
    ]


def _squared_distance(left: AtomSiteRecord, right: AtomSiteRecord) -> float:
    return (
        (left.Cartn_x - right.Cartn_x) ** 2
        + (left.Cartn_y - right.Cartn_y) ** 2
        + (left.Cartn_z - right.Cartn_z) ** 2
    )
