from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path

from pandora.annotations.entry import atoms_by_asym_id, polymer_asym_ids
from pandora.datasets.completeness import chain_completeness
from pandora.schemas.annotation import ContactAtomSet, ResidueContactRow
from pandora.schemas.dataset import (
    ChainRecord,
    DatasetChainRecord,
    InterfaceRecord,
    MissingResidueDefinition,
    ResidueRecord,
)
from pandora.schemas.metadata import MetadataRecord
from pandora.schemas.structure import AtomSiteRecord, Structure


def extract_chain_records(structure: Structure) -> list[ChainRecord]:
    """Build one ChainRecord per polymer chain (asym unit) in `structure`.

    Args:
        structure: The canonical structure to extract chain records
            from. Non-polymer asym units, and polymer asym units with
            no atoms, are skipped.

    Returns:
        A `ChainRecord` per polymer chain, with its sequence (from the
        entity's polymer record, if any), residue count, and atom count.
    """

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

    Args:
        structure: The canonical structure to extract residue records
            from. Only atoms belonging to polymer asym units are
            considered.

    Returns:
        A `ResidueRecord` per distinct polymer residue, sorted by chain
        then sequence id, each carrying every atom belonging to it.
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
    rows: Iterable[ResidueContactRow],
    distance_cutoff: float | None = None,
    atom_set: ContactAtomSet | None = None,
) -> list[InterfaceRecord]:
    """Build one InterfaceRecord per polymer chain pair in contact.

    Works on any `residue_contacts` layer: Pandora's own
    (`annotate_chain_interfaces`) or one brought in with
    `annotate import`. Only rows between polymer residues (both
    `seq_id`s present) on different chains count; ligand and
    same-chain rows are skipped. Pairs are unordered: (A, B) and
    (B, A) rows make one interface, and a residue pair on several rows
    (e.g. several contact types) counts once.

    Args:
        rows: `residue_contacts` rows, for any number of entries.
        distance_cutoff: The contact cutoff the rows were computed
            with, recorded on every record; None if unknown.
        atom_set: The atoms that counted towards contacts, recorded on
            every record; None if unknown.

    Returns:
        An `InterfaceRecord` per entry and polymer chain pair with at
        least one contact, side 1 on the chain that sorts first,
        ordered by entry and chain pair.
    """

    pairs = interface_residue_pairs(rows)
    records = []
    for (entry_id, chain_1, chain_2), residue_pairs in sorted(pairs.items()):
        ordered = sorted(residue_pairs)
        side_1 = sorted({a for a, _ in ordered})
        side_2 = sorted({b for _, b in ordered})
        records.append(
            InterfaceRecord(
                entry_id=entry_id,
                chain_id_1=chain_1,
                chain_id_2=chain_2,
                distance_cutoff=distance_cutoff,
                atom_set=atom_set,
                interface_residues_chain_1=side_1,
                interface_residues_chain_2=side_2,
                contact_count=len(side_1) + len(side_2),
                residue_pairs=ordered,
            )
        )
    return records


def interface_residue_pairs(
    rows: Iterable[ResidueContactRow],
) -> dict[tuple[str, str, str], set[tuple[str, str]]]:
    """Residue pairs across polymer chains, per entry and chain pair.

    Keeps rows with both `seq_id`s present on different chains,
    oriented so side 1 is the chain that sorts first.

    Args:
        rows: `residue_contacts` rows.

    Returns:
        `{(entry_id, chain_1, chain_2): {("A:5", "B:7"), ...}}`, residue
        ids as `label_asym_id:label_seq_id`.
    """

    pairs: dict[tuple[str, str, str], set[tuple[str, str]]] = defaultdict(set)
    for row in rows:
        if row.seq_id_1 is None or row.seq_id_2 is None:
            continue
        if row.chain_id_1 == row.chain_id_2:
            continue
        side_1 = (row.chain_id_1, f"{row.chain_id_1}:{row.seq_id_1}")
        side_2 = (row.chain_id_2, f"{row.chain_id_2}:{row.seq_id_2}")
        (chain_1, res_1), (chain_2, res_2) = sorted((side_1, side_2))
        pairs[(row.entry_id, chain_1, chain_2)].add((res_1, res_2))
    return dict(pairs)


def dataset_chain_records(
    structure: Structure,
    metadata: MetadataRecord | None,
    definition: MissingResidueDefinition = "incomplete_backbone",
    copy_of: Mapping[str, str] | None = None,
    structure_path: str | Path | None = None,
) -> list[DatasetChainRecord]:
    """Build one dataset table row per polymer chain of `structure`.

    Joins `extract_chain_records` with `chain_completeness` and the
    entry's quality metadata. A chain whose completeness can't be
    measured (no SEQRES, or residues that don't line up with it) gets
    None in the completeness fields.

    Args:
        structure: The (curated) structure.
        metadata: The entry's metadata, for method and resolution.
        definition: What counts as a missing residue; pass the curation
            policy's `quality_rules.missing_residue_definition`.
        copy_of: Symmetry-copy chain id to its source chain id, e.g.
            from `chain_copy_records`. Chains not in it get None.
        structure_path: The file the structure was written to, stored
            on every row.

    Returns:
        One `DatasetChainRecord` per row of
        `extract_chain_records(structure)`, in the same order.
    """

    completeness = {
        c.chain_id: c for c in chain_completeness(structure, definition)[0]
    }
    quality = metadata.quality if metadata else None
    copies = copy_of or {}
    rows = []
    for record in extract_chain_records(structure):
        c = completeness.get(record.chain_id)
        rows.append(
            DatasetChainRecord(
                **record.model_dump(),
                seqres_length=c.seqres_length if c else None,
                missing_n_term=c.missing_n_term if c else None,
                missing_c_term=c.missing_c_term if c else None,
                missing_middle=c.missing_middle if c else None,
                experimental_method=(
                    quality.experimental_method if quality else None
                ),
                resolution=quality.resolution if quality else None,
                structure_path=(
                    None if structure_path is None else str(structure_path)
                ),
                copy_of=copies.get(record.chain_id),
            )
        )
    return rows
