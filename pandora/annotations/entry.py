from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations
from typing import Any, get_args

from pandora.canonicalisation.missing_data import (
    _PROTEIN_BACKBONE,
    _backbone_by_entity,
)
from pandora.schemas.annotation import AnnotationLayer, ContactAtomSet
from pandora.schemas.structure import AtomSiteRecord, Structure

WATER_COMP_IDS = frozenset({"HOH", "WAT", "DOD"})
HYDROGEN_ELEMENTS = frozenset({"H", "D"})
CONTACT_ATOM_SETS = frozenset(get_args(ContactAtomSet))


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
    polymer_grid = _spatial_grid(polymer_atoms, distance_cutoff)
    contacts = []
    for ligand_key, atoms in ligand_groups.items():
        residue_contacts = _residues_within_cutoff(
            atoms, polymer_grid, distance_cutoff, cutoff_sq
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
        parameters={
            "distance_cutoff": distance_cutoff,
            "include_waters": include_waters,
        },
        data={
            "distance_cutoff": distance_cutoff,
            "include_waters": include_waters,
            "ligands": contacts,
        },
        provenance={"inputs": ["Structure.atoms"]},
    )


def polymer_asym_ids(structure: Structure) -> set[str]:
    """label_asym_ids of asym units whose entity is a polymer."""

    entities_by_id = {entity.id: entity for entity in structure.entities}
    return {
        asym.id
        for asym in structure.asym_units
        if (entity := entities_by_id.get(asym.entity_id)) is not None
        and entity.type == "polymer"
    }


def atoms_by_asym_id(structure: Structure) -> dict[str, list[AtomSiteRecord]]:
    """Structure.atoms grouped by label_asym_id."""

    grouped: dict[str, list[AtomSiteRecord]] = defaultdict(list)
    for atom in structure.atoms:
        grouped[atom.label_asym_id].append(atom)
    return grouped


def annotate_chain_interfaces(
    structure: Structure,
    distance_cutoff: float = 4.0,
    atom_set: ContactAtomSet = "heavy",
    polymer_types: list[str] | None = None,
) -> AnnotationLayer:
    """Compute polymer chain pairs with residues in contact within a cutoff.

    For each pair of polymer chains, finds residue pairs (one per chain)
    with at least one pair of selected atoms within `distance_cutoff`
    angstroms. Chain pairs whose bounding boxes are further apart than
    the cutoff are skipped without comparing atoms.

    Args:
        structure: The structure to scan for chain-chain contacts.
        distance_cutoff: Contact distance in angstroms.
        atom_set: Which atoms count: "all", "heavy" (every atom except
            H/D), "backbone" (protein N/CA/C/O, nucleic-acid
            O5'/C5'/C4'/C3'/O3'), or "ca".
        polymer_types: Keep only chains whose entity polymer type
            (e.g. "polypeptide(L)") is listed. None keeps every
            polymer chain; an empty list keeps none.

    Returns:
        An `AnnotationLayer` of type "chain_interfaces" whose `data`
        holds the settings used and, per chain pair in contact, the
        contacting residue ids on each side and the residue pairs.

    Raises:
        ValueError: `atom_set` is not one of the supported values.
        TypeError: `polymer_types` is a single string, not a list.
    """

    if isinstance(polymer_types, str):
        # A bare string would match by substring here and be recorded
        # character by character, so the recipe would rebuild differently.
        raise TypeError(
            "polymer_types must be a list of polymer types, e.g. "
            f"[{polymer_types!r}], not a single string"
        )
    if atom_set not in CONTACT_ATOM_SETS:
        raise ValueError(
            f"atom_set must be one of {sorted(CONTACT_ATOM_SETS)}, "
            f"got {atom_set!r}"
        )

    selected = _select_contact_atoms(structure, atom_set, polymer_types)
    boxes = {
        chain_id: _bounding_box(atoms) for chain_id, atoms in selected.items()
    }
    grids = {
        chain_id: _spatial_grid(atoms, distance_cutoff)
        for chain_id, atoms in selected.items()
    }
    cutoff_sq = distance_cutoff * distance_cutoff

    interfaces = []
    for chain_a, chain_b in combinations(sorted(selected), 2):
        if _boxes_apart(boxes[chain_a], boxes[chain_b], distance_cutoff):
            continue
        pairs = _chain_pair_contacts(
            selected[chain_a],
            boxes[chain_b],
            grids[chain_b],
            distance_cutoff,
            cutoff_sq,
        )
        if not pairs:
            continue
        residues_a = sorted({residue_a for residue_a, _ in pairs})
        residues_b = sorted({residue_b for _, residue_b in pairs})
        interfaces.append(
            {
                "chain_id_1": chain_a,
                "chain_id_2": chain_b,
                "interface_residues_chain_1": residues_a,
                "interface_residues_chain_2": residues_b,
                "residue_pairs": [list(pair) for pair in sorted(pairs)],
                "contact_count": len(residues_a) + len(residues_b),
            }
        )

    return AnnotationLayer(
        layer_name="Chain-chain interfaces",
        layer_type="chain_interfaces",
        scope="interface",
        method="pandora.basic.distance_cutoff_contacts.v2",
        target_ids=[structure.entry_id],
        parameters=_contact_settings(distance_cutoff, atom_set, polymer_types),
        data={
            **_contact_settings(distance_cutoff, atom_set, polymer_types),
            "interfaces": interfaces,
        },
        provenance={"inputs": ["Structure.atoms", "Structure.entities"]},
    )


def _contact_settings(
    distance_cutoff: float,
    atom_set: ContactAtomSet,
    polymer_types: list[str] | None,
) -> dict[str, Any]:
    """Settings dict with its own copy of polymer_types (never aliased)."""

    return {
        "distance_cutoff": distance_cutoff,
        "atom_set": atom_set,
        "polymer_types": (
            list(polymer_types) if polymer_types is not None else None
        ),
    }


def _select_contact_atoms(
    structure: Structure,
    atom_set: ContactAtomSet,
    polymer_types: list[str] | None,
) -> dict[str, list[AtomSiteRecord]]:
    """Selected atoms per kept polymer chain; chains with none are dropped."""

    entities = {entity.id: entity for entity in structure.entities}
    chain_entity: dict[str, str] = {}
    for asym in structure.asym_units:
        entity = entities.get(asym.entity_id)
        if entity is None or entity.type != "polymer":
            continue
        if polymer_types is not None and (
            entity.poly is None or entity.poly.type not in polymer_types
        ):
            continue
        chain_entity[asym.id] = entity.id

    nucleic_backbone = _backbone_by_entity(structure.entities)
    selected: dict[str, list[AtomSiteRecord]] = defaultdict(list)
    for atom in structure.atoms:
        entity_id = chain_entity.get(atom.label_asym_id)
        if entity_id is None:
            continue
        backbone = nucleic_backbone.get(entity_id, _PROTEIN_BACKBONE)
        if _in_atom_set(atom, atom_set, backbone):
            selected[atom.label_asym_id].append(atom)
    return dict(selected)


def _in_atom_set(
    atom: AtomSiteRecord, atom_set: ContactAtomSet, backbone: frozenset[str]
) -> bool:
    """Whether atom belongs to atom_set (backbone: its chain's names)."""

    if atom_set == "all":
        return True
    if atom_set == "heavy":
        return atom.type_symbol.upper() not in HYDROGEN_ELEMENTS
    if atom_set == "backbone":
        return atom.label_atom_id in backbone
    return atom.label_atom_id == "CA"


def _chain_pair_contacts(
    atoms_a: list[AtomSiteRecord],
    box_b: tuple[float, float, float, float, float, float],
    grid_b: dict[tuple[int, int, int], list[AtomSiteRecord]],
    cutoff: float,
    cutoff_sq: float,
) -> set[tuple[str, str]]:
    """Residue-id pairs with an atom of atoms_a within cutoff of grid_b."""

    min_x, min_y, min_z, max_x, max_y, max_z = box_b
    cell_size = _cell_size(cutoff)
    pairs: set[tuple[str, str]] = set()
    for atom_a in atoms_a:
        # Only atoms inside chain b's box grown by the cutoff can touch it.
        if not (
            min_x - cutoff <= atom_a.Cartn_x <= max_x + cutoff
            and min_y - cutoff <= atom_a.Cartn_y <= max_y + cutoff
            and min_z - cutoff <= atom_a.Cartn_z <= max_z + cutoff
        ):
            continue
        for atom_b in _neighbor_atoms(grid_b, _cell_key(atom_a, cell_size)):
            if _squared_distance(atom_a, atom_b) <= cutoff_sq:
                pairs.add((_residue_id(atom_a), _residue_id(atom_b)))
    return pairs


def _bounding_box(
    atoms: list[AtomSiteRecord],
) -> tuple[float, float, float, float, float, float]:
    """Axis-aligned box (min x, min y, min z, max x, max y, max z)."""

    xs = [atom.Cartn_x for atom in atoms]
    ys = [atom.Cartn_y for atom in atoms]
    zs = [atom.Cartn_z for atom in atoms]
    return min(xs), min(ys), min(zs), max(xs), max(ys), max(zs)


def _boxes_apart(
    box_a: tuple[float, float, float, float, float, float],
    box_b: tuple[float, float, float, float, float, float],
    cutoff: float,
) -> bool:
    """Whether the two boxes are more than cutoff apart on any axis."""

    return any(
        box_a[axis] - box_b[axis + 3] > cutoff
        or box_b[axis] - box_a[axis + 3] > cutoff
        for axis in range(3)
    )


def _residue_id(atom: AtomSiteRecord) -> str:
    """The `label_asym_id:label_seq_id` id of atom's residue."""

    return f"{atom.label_asym_id}:{atom.label_seq_id}"


def _is_ligand_atom(atom: AtomSiteRecord, include_waters: bool) -> bool:
    """Whether atom is a HETATM ligand atom (optionally including waters)."""

    if atom.group_PDB != "HETATM":
        return False
    if include_waters:
        return True
    return atom.label_comp_id.upper() not in WATER_COMP_IDS


def _residues_within_cutoff(
    ligand_atoms: list[AtomSiteRecord],
    polymer_grid: dict[tuple[int, int, int], list[AtomSiteRecord]],
    cutoff: float,
    cutoff_sq: float,
) -> list[dict[str, Any]]:
    """Polymer residues in polymer_grid within cutoff of any ligand_atoms
    atom, with nearest distance."""

    contacts: dict[tuple[str, int | None, str], float] = {}
    cell_size = _cell_size(cutoff)
    for ligand_atom in ligand_atoms:
        key = _cell_key(ligand_atom, cell_size)
        for polymer_atom in _neighbor_atoms(polymer_grid, key):
            dist_sq = _squared_distance(ligand_atom, polymer_atom)
            if dist_sq > cutoff_sq:
                continue

            rk = (
                polymer_atom.label_asym_id,
                polymer_atom.label_seq_id,
                polymer_atom.label_comp_id,
            )
            current = contacts.get(rk)
            if current is None or dist_sq < current:
                contacts[rk] = dist_sq

    return [
        {
            "label_asym_id": asym_id,
            "label_seq_id": seq_id,
            "label_comp_id": comp_id,
            "distance": round(dist_sq**0.5, 3),
        }
        for (asym_id, seq_id, comp_id), dist_sq in sorted(contacts.items())
    ]


# Cell-list spatial index: bucket atoms into cutoff-sized grid cells so a
# distance query only has to check the 27 cells around a point instead of
# every atom. Cell size == cutoff guarantees any atom within cutoff of a
# point falls in that point's cell or an immediate neighbor (standard
# Verlet cell-list argument), so this is exact, not approximate.


def _cell_size(cutoff: float) -> float:
    """Grid cell edge length for a distance cutoff (guards against a
    zero/negative cutoff)."""

    return cutoff if cutoff > 0 else 1e-6


def _cell_key(atom: AtomSiteRecord, cell_size: float) -> tuple[int, int, int]:
    """Grid cell coordinates containing atom, for a given cell_size."""

    return (
        int(atom.Cartn_x // cell_size),
        int(atom.Cartn_y // cell_size),
        int(atom.Cartn_z // cell_size),
    )


def _spatial_grid(
    atoms: list[AtomSiteRecord], cutoff: float
) -> dict[tuple[int, int, int], list[AtomSiteRecord]]:
    """Bucket atoms into cutoff-sized grid cells for neighbor queries."""

    cell_size = _cell_size(cutoff)
    grid: dict[tuple[int, int, int], list[AtomSiteRecord]] = defaultdict(list)
    for atom in atoms:
        grid[_cell_key(atom, cell_size)].append(atom)
    return grid


def _neighbor_atoms(
    grid: dict[tuple[int, int, int], list[AtomSiteRecord]],
    key: tuple[int, int, int],
):
    """Atoms in the 27 grid cells around key (the cell plus its
    immediate neighbors)."""

    x, y, z = key
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                yield from grid.get((x + dx, y + dy, z + dz), ())


def _squared_distance(left: AtomSiteRecord, right: AtomSiteRecord) -> float:
    """Squared Euclidean distance between two atoms' Cartesian coordinates."""

    return (
        (left.Cartn_x - right.Cartn_x) ** 2
        + (left.Cartn_y - right.Cartn_y) ** 2
        + (left.Cartn_z - right.Cartn_z) ** 2
    )
