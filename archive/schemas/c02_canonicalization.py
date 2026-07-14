from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from .common import AppliedPolicyRef, DiagnosticBundle, ResultStatus
from .c01_ingestion import Atom, Assembly, Chain, Entity, Ligand, Residue


# ── Canonical structure ────────────────────────────────────────────────────────


class CanonicalStructure(BaseModel):
    """Reuses C01 structural types; all identifier fields reflect canonical values."""

    atoms: list[Atom] = Field(default_factory=list)
    residues: list[Residue] = Field(default_factory=list)
    chains: list[Chain] = Field(default_factory=list)
    entities: list[Entity] = Field(default_factory=list)
    assemblies: list[Assembly] = Field(default_factory=list)
    ligands: list[Ligand] = Field(default_factory=list)


# ── Mapping schemas ────────────────────────────────────────────────────────────


class ChainIdMappingItem(BaseModel):
    canonical_chain_id: str
    original_chain_id: str
    original_auth_chain_id: str


class ChainIdMapping(BaseModel):
    items: list[ChainIdMappingItem] = Field(default_factory=list)


class ResidueNumberMappingItem(BaseModel):
    canonical_chain_id: str
    canonical_seq_id: int
    original_chain_id: str
    original_seq_id: int | None = None
    original_auth_seq_id: str
    original_insertion_code: str | None = None


class ResidueNumberMapping(BaseModel):
    items: list[ResidueNumberMappingItem] = Field(default_factory=list)


class AssemblyMappingItem(BaseModel):
    canonical_assembly_id: str
    original_assembly_id: str


class AssemblyMapping(BaseModel):
    items: list[AssemblyMappingItem] = Field(default_factory=list)


class EntityMappingItem(BaseModel):
    canonical_entity_id: str
    original_entity_ids: list[str]


class EntityMapping(BaseModel):
    items: list[EntityMappingItem] = Field(default_factory=list)


class AltlocSelectionMappingItem(BaseModel):
    canonical_chain_id: str
    residue_id: str
    selected_altloc: str
    available_altlocs: list[str]
    selection_reason: Literal[
        "best_occupancy",
        "first_alphabetical",
        "last_alphabetical",
        "user_defined",
    ]


class AltlocSelectionMapping(BaseModel):
    items: list[AltlocSelectionMappingItem] = Field(default_factory=list)


class CanonicalMappings(BaseModel):
    chain_id_mapping: ChainIdMapping = Field(default_factory=ChainIdMapping)
    residue_number_mapping: ResidueNumberMapping = Field(
        default_factory=ResidueNumberMapping
    )
    assembly_mapping: AssemblyMapping = Field(default_factory=AssemblyMapping)
    entity_mapping: EntityMapping = Field(default_factory=EntityMapping)
    altloc_selection_mapping: AltlocSelectionMapping = Field(
        default_factory=AltlocSelectionMapping
    )


# ── Output schemas ─────────────────────────────────────────────────────────────


class CanonicalStructureProvenance(BaseModel):
    provider: str
    source_uri: str | None = None
    retrieved_at: str | None = None
    canonicalised_at: str | None = None


class CanonicalStructureResult(BaseModel):
    entry_id: str
    status: ResultStatus
    canonical_structure: CanonicalStructure | None = None
    canonical_mappings: CanonicalMappings | None = None
    applied_policy: AppliedPolicyRef
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: CanonicalStructureProvenance


class canonicalisationBatchSummary(BaseModel):
    total: int
    success: int
    warning: int
    failed: int


class canonicalisationBatchResult(BaseModel):
    mode: Literal["sequential", "parallel"]
    summary: canonicalisationBatchSummary
    results: list[CanonicalStructureResult]


# ── Policy schema ──────────────────────────────────────────────────────────────

ChainIdStrategy = Literal["preserve", "remap", "use_auth_chain_id"]
ResidueNumberingStrategy = Literal["preserve", "use_auth_seq", "renumber"]
AssemblyIdStrategy = Literal["preserve", "remap", "standardize"]
MissingAtomsStrategy = Literal[
    "preserve", "annotate", "drop_partial_residue", "impute"
]
MissingResiduesStrategy = Literal[
    "preserve", "annotate", "drop_chain_segment", "impute"
]
IncompleteChainStrategy = Literal[
    "preserve", "exclude", "truncate_to_complete_regions"
]
AltlocStrategy = Literal[
    "preserve", "select_best_occupancy", "select_first", "select_user_defined"
]
AltlocTieBreaker = Literal[
    "alphabetical_first",
    "alphabetical_last",
    "lowest_b_factor",
    "highest_b_factor",
]
AssemblyStrategy = Literal[
    "preserve_as_reported",
    "standardize_biological_assembly",
    "select_first_assembly",
]
PreferredAssemblySource = Literal["author", "pdbe", "pdbx", "first"]
EntityStrategy = Literal["preserve", "standardize", "merge_equivalent_entities"]
LigandStrategy = Literal["preserve", "filter", "annotate_only"]
ValidationStrictness = Literal["strict", "moderate", "permissive"]


class ChainIdRules(BaseModel):
    strategy: ChainIdStrategy = "preserve"


class ResidueNumberingRules(BaseModel):
    strategy: ResidueNumberingStrategy = "preserve"
    preserve_insertion_codes: bool = True


class AssemblyIdRules(BaseModel):
    strategy: AssemblyIdStrategy = "preserve"


class IdentifierRules(BaseModel):
    chain_id: ChainIdRules = Field(default_factory=ChainIdRules)
    residue_numbering: ResidueNumberingRules = Field(
        default_factory=ResidueNumberingRules
    )
    assembly_id: AssemblyIdRules = Field(default_factory=AssemblyIdRules)


class MissingAtomsRules(BaseModel):
    strategy: MissingAtomsStrategy = "annotate"
    allow_imputation: bool = False
    record_missingness: bool = True


class MissingResiduesRules(BaseModel):
    strategy: MissingResiduesStrategy = "annotate"
    record_gaps: bool = True


class IncompleteChainRules(BaseModel):
    strategy: IncompleteChainStrategy = "preserve"


class MissingDataRules(BaseModel):
    missing_atoms: MissingAtomsRules = Field(default_factory=MissingAtomsRules)
    missing_residues: MissingResiduesRules = Field(
        default_factory=MissingResiduesRules
    )
    incomplete_chains: IncompleteChainRules = Field(
        default_factory=IncompleteChainRules
    )


class AltlocRules(BaseModel):
    strategy: AltlocStrategy = "select_best_occupancy"
    tie_breaker: AltlocTieBreaker = "alphabetical_first"
    user_defined_altloc: str | None = None
    record_selection: bool = True


class AssemblyRules(BaseModel):
    strategy: AssemblyStrategy = "preserve_as_reported"
    preferred_assembly_source: PreferredAssemblySource = "author"
    record_original_assembly_mapping: bool = True


class EntityRules(BaseModel):
    strategy: EntityStrategy = "preserve"
    preserve_original_entity_ids: bool = True


class LigandRules(BaseModel):
    strategy: LigandStrategy = "preserve"
    keep_waters: bool = True
    keep_ions: bool = True
    keep_nonpolymer_ligands: bool = True


class ValidationRules(BaseModel):
    strictness: ValidationStrictness = "moderate"
    fail_on_unresolved_issues: bool = True
    warnings_as_errors: bool = False


class canonicalisationProvenanceRules(BaseModel):
    record_original_mappings: bool = True
    record_transforms: bool = True
    record_policy_application: bool = True
    emit_canonicalisation_report: bool = False


class canonicalisationPolicy(BaseModel):
    policy_id: str
    policy_name: str
    policy_version: str
    description: str = ""
    identifier_rules: IdentifierRules = Field(default_factory=IdentifierRules)
    missing_data_rules: MissingDataRules = Field(
        default_factory=MissingDataRules
    )
    altloc_rules: AltlocRules = Field(default_factory=AltlocRules)
    assembly_rules: AssemblyRules = Field(default_factory=AssemblyRules)
    entity_rules: EntityRules = Field(default_factory=EntityRules)
    ligand_rules: LigandRules = Field(default_factory=LigandRules)
    validation_rules: ValidationRules = Field(default_factory=ValidationRules)
    provenance_rules: canonicalisationProvenanceRules = Field(
        default_factory=canonicalisationProvenanceRules
    )
