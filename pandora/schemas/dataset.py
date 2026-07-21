from __future__ import annotations

from pydantic import BaseModel

from pandora.schemas.structure import AtomSiteRecord


class ChainRecord(BaseModel):
    entry_id: str
    chain_id: str  # label_asym_id
    auth_chain_id: str | None = None
    entity_id: str
    sequence: str | None = None
    residue_count: int
    atom_count: int


class ResidueRecord(BaseModel):
    entry_id: str
    chain_id: str  # label_asym_id
    seq_id: int | None = None
    auth_seq_id: str | None = None
    comp_id: str
    atoms: list[AtomSiteRecord]
    # Coordinates/B-factor live on each AtomSiteRecord — this record is a
    # residue-scoped slice of Structure.atoms, not a re-derived summary.


class InterfaceRecord(BaseModel):
    entry_id: str
    chain_id_1: str
    chain_id_2: str
    distance_cutoff: float
    interface_residues_chain_1: list[str]
    interface_residues_chain_2: list[str]
    contact_count: int
    # Populated from annotations.entry.annotate_chain_interfaces() — this
    # record never computes contacts itself, only reshapes that layer's
    # output per chain pair.


# TODO: Add Dataset, ChainDataset, InterfaceDataset, ResidueDataset, and
# DatasetStoreRef (materialized-mode collection types, C04/C05).
