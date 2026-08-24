from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from pandora.schemas.structure import AtomSiteRecord


class ChainRecord(BaseModel):
    """One polymer chain, reshaped from a canonical `Structure` for
    dataset export.

    Attributes:
        entry_id: The structure's entry id.
        chain_id: The chain's `label_asym_id`.
        auth_chain_id: The chain's `auth_asym_id`, if different.
        entity_id: The id of the entity this chain belongs to.
        sequence: The chain's polymer sequence, if any.
        residue_count: The number of distinct residues in the chain.
        atom_count: The number of atoms in the chain.
    """

    entry_id: str
    chain_id: str  # label_asym_id
    auth_chain_id: str | None = None
    entity_id: str
    sequence: str | None = None
    residue_count: int
    atom_count: int


class ResidueRecord(BaseModel):
    """One polymer residue and its atoms, reshaped from a canonical
    `Structure` for dataset export.

    Attributes:
        entry_id: The structure's entry id.
        chain_id: The residue's `label_asym_id`.
        seq_id: The residue's `label_seq_id`, if polymeric.
        auth_seq_id: The residue's `auth_seq_id`.
        comp_id: The residue's component id (e.g. "ALA").
        atoms: Every atom belonging to this residue.
    """

    entry_id: str
    chain_id: str  # label_asym_id
    seq_id: int | None = None
    auth_seq_id: str | None = None
    comp_id: str
    atoms: list[AtomSiteRecord]
    # Coordinates/B-factor live on each AtomSiteRecord — this record is a
    # residue-scoped slice of Structure.atoms, not a re-derived summary.


class InterfaceRecord(BaseModel):
    """One polymer chain pair in contact, reshaped from a
    chain-interfaces annotation layer.

    Attributes:
        entry_id: The structure's entry id.
        chain_id_1: The first chain in the interface.
        chain_id_2: The second chain in the interface.
        distance_cutoff: The contact distance cutoff used, in
            angstroms.
        interface_residues_chain_1: Contacting residue ids on
            `chain_id_1`'s side.
        interface_residues_chain_2: Contacting residue ids on
            `chain_id_2`'s side.
        contact_count: The total number of contacting residues on
            both sides.
    """

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
    """Record of why one structure was excluded during curation or
    deduplication.

    Attributes:
        entry_id: The excluded structure's entry id.
        reason_code: A machine-readable reason the structure was
            excluded.
        message: A human-readable description of the exclusion.
    """

    entry_id: str
    reason_code: ExclusionReason
    message: str


class QualityRules(BaseModel):
    """Curation policy for resolution, experimental method, and minimum
    chain length.

    Attributes:
        max_resolution: The maximum resolution (in angstroms) allowed;
            structures coarser than this are excluded.
        null_resolution_behavior: Whether a structure with no
            resolution is excluded or included when `max_resolution`
            is set.
        min_chain_length: The minimum residue count a structure's
            longest chain must reach.
        include_experimental_methods: If non-empty, only these
            experimental methods are allowed.
        exclude_experimental_methods: Experimental methods that are
            always excluded.
    """

    max_resolution: float | None = None
    null_resolution_behavior: Literal["exclude", "include"] = "exclude"
    min_chain_length: int | None = None
    include_experimental_methods: list[str] = Field(default_factory=list)
    exclude_experimental_methods: list[str] = Field(default_factory=list)


class OrganismRules(BaseModel):
    """Curation policy for including/excluding structures by source
    organism taxon.

    Attributes:
        include_taxa: If non-empty, only these NCBI taxon ids are
            allowed.
        exclude_taxa: NCBI taxon ids that are always excluded.
    """

    include_taxa: list[str] = Field(default_factory=list)
    exclude_taxa: list[str] = Field(default_factory=list)


class ContentRules(BaseModel):
    """Curation policy for which ligand/water/ion content is kept.

    Attributes:
        keep_ligands: Whether non-polymer ligands are kept in curated
            structures.
        keep_waters: Whether water molecules are kept in curated
            structures.
        keep_ions: Whether ions are kept in curated structures.
    """

    keep_ligands: bool = True
    keep_waters: bool = True
    keep_ions: bool = True


class DeduplicationRules(BaseModel):
    """Policy for whether duplicate entry_ids are removed.

    Attributes:
        enabled: Whether structures sharing the same `entry_id` are
            removed.
    """

    enabled: bool = False


class DatasetCurationPolicy(BaseModel):
    """The full, named policy governing dataset curation: quality,
    organism, and content rules.

    Attributes:
        policy_id: A unique identifier for this policy.
        policy_name: A human-readable name for this policy.
        policy_version: The policy's version string.
        description: A free-text description of the policy.
        quality_rules: The resolution/method/chain-length rules.
        organism_rules: The organism taxon inclusion/exclusion rules.
        content_rules: The ligand/water/ion content rules.
    """

    policy_id: str
    policy_name: str
    policy_version: str
    description: str = ""
    quality_rules: QualityRules = Field(default_factory=QualityRules)
    organism_rules: OrganismRules = Field(default_factory=OrganismRules)
    content_rules: ContentRules = Field(default_factory=ContentRules)


class CurationProvenance(BaseModel):
    """Record of which curation policy ran on one structure, and when.

    Attributes:
        curated_at: When this curation run completed, as an ISO 8601
            timestamp.
        policy_id: The id of the policy that ran.
        policy_name: The name of the policy that ran.
        policy_version: The version of the policy that ran.
    """

    curated_at: str
    policy_id: str
    policy_name: str
    policy_version: str


class DeduplicationProvenance(BaseModel):
    """Record of one deduplication run: whether it was enabled and how
    many duplicates it found.

    Attributes:
        deduplicated_at: When this deduplication run completed, as an
            ISO 8601 timestamp.
        enabled: Whether deduplication was actually applied.
        duplicates_found: How many structures were removed as
            duplicates.
    """

    deduplicated_at: str
    enabled: bool
    duplicates_found: int = 0
