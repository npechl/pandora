from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

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


# Curation ---------------------------------


ExclusionReason = Literal[
    "RESOLUTION_THRESHOLD",
    "NULL_RESOLUTION",
    "CHAIN_TOO_SHORT",
    "METHOD_EXCLUDED",
    "ORGANISM_EXCLUDED",
    "MISSING_TAXONOMY",
    "DUPLICATE",
]


class ExclusionRecord(BaseModel):
    entry_id: str
    reason_code: ExclusionReason
    message: str


class QualityRules(BaseModel):
    max_resolution: float | None = None
    null_resolution_behavior: Literal["exclude", "include"] = "exclude"
    min_chain_length: int | None = None
    include_experimental_methods: list[str] = Field(default_factory=list)
    exclude_experimental_methods: list[str] = Field(default_factory=list)


class OrganismRules(BaseModel):
    include_taxa: list[str] = Field(default_factory=list)
    exclude_taxa: list[str] = Field(default_factory=list)


class ContentRules(BaseModel):
    keep_ligands: bool = True
    keep_waters: bool = True
    keep_ions: bool = True


class DeduplicationRules(BaseModel):
    enabled: bool = False
    strategy: Literal["entry_id", "exact_hash"] = "entry_id"


class DatasetCurationPolicy(BaseModel):
    policy_id: str
    policy_name: str
    policy_version: str
    description: str = ""
    quality_rules: QualityRules = Field(default_factory=QualityRules)
    organism_rules: OrganismRules = Field(default_factory=OrganismRules)
    content_rules: ContentRules = Field(default_factory=ContentRules)
    deduplication_rules: DeduplicationRules = Field(
        default_factory=DeduplicationRules
    )


class CurationProvenance(BaseModel):
    curated_at: str
    policy_id: str
    policy_name: str
    policy_version: str


class DeduplicationProvenance(BaseModel):
    deduplicated_at: str
    enabled: bool
    strategy: str | None = None
    duplicates_found: int = 0
