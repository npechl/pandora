from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from typing import Any, get_args

from pandora.canonicalisation.missing_data import (
    _PROTEIN_BACKBONE,
    _backbone_by_entity,
)
from pandora.schemas.annotation import (
    ContactAtomSet,
    ResidueContactRow,
)
from pandora.schemas.structure import AtomSiteRecord, Structure

WATER_COMP_IDS = frozenset({"HOH", "WAT", "DOD"})
HYDROGEN_ELEMENTS = frozenset({"H", "D"})
CONTACT_ATOM_SETS = frozenset(get_args(ContactAtomSet))


def annotate_ligand_contacts(
    structure: Structure,
    distance_cutoff: float = 4.0,
    include_waters: bool = False,
) -> list[ResidueContactRow]:
    """Compute polymer residues near ligands within a cutoff.

    A ligand is a residue of a non-polymer chain (waters only with
    `include_waters`). A modified residue inside a polymer chain is
    not a ligand, even when written as HETATM: it counts on the
    polymer side.

    Args:
        structure: The structure to scan for ligand contacts.
        distance_cutoff: Contact distance in angstroms.
        include_waters: If True, treat water residues as ligands too.

    Returns:
        A `residue_contacts` layer: one `ResidueContactRow` per
        (polymer residue, ligand residue) in contact, side 1 the
        polymer residue, side 2 the ligand (`seq_id_2` None, named by
        its author ids), `distance` the nearest atom distance (rounded
        to 3 decimals). Sorted by ligand, then polymer residue.
    """

    polymers = polymer_asym_ids(structure)
    polymer_atoms = [
        atom for atom in structure.atoms if atom.label_asym_id in polymers
    ]
    ligand_atoms = [
        atom
        for atom in structure.atoms
        if atom.label_asym_id not in polymers
        and (include_waters or atom.label_comp_id.upper() not in WATER_COMP_IDS)
    ]

    cutoff_sq = distance_cutoff * distance_cutoff
    cell_size = _cell_size(distance_cutoff)
    polymer_grid = _spatial_grid(polymer_atoms, distance_cutoff)
    nearest: dict[tuple[tuple, tuple], float] = {}
    atoms_of: dict[tuple, AtomSiteRecord] = {}
    for ligand_atom in ligand_atoms:
        ligand_key = _residue_key(ligand_atom)
        atoms_of.setdefault(ligand_key, ligand_atom)
        cell = _cell_key(ligand_atom, cell_size)
        for polymer_atom in _neighbor_atoms(polymer_grid, cell):
            dist_sq = _squared_distance(ligand_atom, polymer_atom)
            if dist_sq > cutoff_sq:
                continue
            polymer_key = _residue_key(polymer_atom)
            atoms_of.setdefault(polymer_key, polymer_atom)
            key = (ligand_key, polymer_key)
            if dist_sq < nearest.get(key, cutoff_sq + 1.0):
                nearest[key] = dist_sq

    return [
        _contact_row(
            structure.entry_id,
            atoms_of[polymer_key],
            atoms_of[ligand_key],
            dist_sq,
        )
        for (ligand_key, polymer_key), dist_sq in sorted(
            nearest.items(), key=lambda item: _sort_key(item[0])
        )
    ]


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
) -> list[ResidueContactRow]:
    """Compute residue pairs in contact across polymer chains.

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
        A `residue_contacts` layer: one `ResidueContactRow` per
        contacting residue pair, side 1 on the chain that sorts first,
        both namings filled in, `distance` the closest selected-atom
        distance of the pair (rounded to 3 decimals). Sorted by chain
        pair, then residue pair.

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

    residues = {
        _residue_key(atom): atom
        for atoms in selected.values()
        for atom in reversed(atoms)
    }
    rows = []
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
        for (key_a, key_b), dist_sq in sorted(
            pairs.items(), key=lambda item: _sort_key(item[0])
        ):
            rows.append(
                _contact_row(
                    structure.entry_id,
                    residues[key_a],
                    residues[key_b],
                    dist_sq,
                )
            )
    return rows


def _residue_key(
    atom: AtomSiteRecord,
) -> tuple[str, int | None, str | None, str]:
    """An atom's residue: (label_asym_id, label_seq_id, auth_seq_id,
    insertion code or "")."""

    return (
        atom.label_asym_id,
        atom.label_seq_id,
        atom.auth_seq_id,
        _ins_code(atom),
    )


def _ins_code(atom: AtomSiteRecord) -> str:
    """An atom's insertion code, "" for none ("?" and "." included)."""

    ins = (atom.pdbx_PDB_ins_code or "").strip()
    return "" if ins in {"?", "."} else ins


def _sort_key(pair: tuple[tuple, tuple]) -> tuple:
    """Sort residue pairs by chain, then residue number (None first)."""

    return tuple(
        (chain, seq is not None, seq or 0, auth or "", ins)
        for chain, seq, auth, ins in pair
    )


def _side(atom: AtomSiteRecord, side: str) -> dict[str, Any]:
    """Both namings of atom's residue, as `ResidueContactRow` side fields."""

    ins = _ins_code(atom)
    return {
        f"chain_id{side}": atom.label_asym_id,
        f"seq_id{side}": atom.label_seq_id,
        f"auth_chain_id{side}": atom.auth_asym_id,
        f"auth_seq_id{side}": atom.auth_seq_id,
        f"ins_code{side}": ins or None,
        f"comp_id{side}": atom.label_comp_id,
    }


def _contact_row(
    entry_id: str,
    atom_1: AtomSiteRecord,
    atom_2: AtomSiteRecord,
    dist_sq: float,
) -> ResidueContactRow:
    """A contact row between the residues of atom_1 and atom_2."""

    return ResidueContactRow(
        entry_id=entry_id,
        **_side(atom_1, "_1"),
        **_side(atom_2, "_2"),
        distance=round(dist_sq**0.5, 3),
    )


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
) -> dict[tuple[tuple, tuple], float]:
    """Residue pairs with an atom of atoms_a within cutoff of grid_b,
    with their closest squared distance."""

    min_x, min_y, min_z, max_x, max_y, max_z = box_b
    cell_size = _cell_size(cutoff)
    pairs: dict[tuple[tuple, tuple], float] = {}
    for atom_a in atoms_a:
        # Only atoms inside chain b's box grown by the cutoff can touch it.
        if not (
            min_x - cutoff <= atom_a.Cartn_x <= max_x + cutoff
            and min_y - cutoff <= atom_a.Cartn_y <= max_y + cutoff
            and min_z - cutoff <= atom_a.Cartn_z <= max_z + cutoff
        ):
            continue
        for atom_b in _neighbor_atoms(grid_b, _cell_key(atom_a, cell_size)):
            dist_sq = _squared_distance(atom_a, atom_b)
            if dist_sq > cutoff_sq:
                continue
            key = (_residue_key(atom_a), _residue_key(atom_b))
            if dist_sq < pairs.get(key, cutoff_sq + 1.0):
                pairs[key] = dist_sq
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
