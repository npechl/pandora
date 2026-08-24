from __future__ import annotations
from typing import Any, Literal
from pydantic import BaseModel, Field

# Mapping classes ----------------------


class AltlocSelectionMappingItem(BaseModel):
    """Record of which altloc was selected for one residue, and why.

    Attributes:
        canonical_chain_id: The canonical chain id the residue belongs
            to.
        residue_id: The residue's identifier (`seq_id:comp_id`).
        selected_altloc: The altloc code that was kept.
        available_altlocs: Every altloc code that was present before
            selection.
        selection_reason: Why this altloc was selected.
    """

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
    """Every altloc selection made while canonicalising a structure.

    Attributes:
        items: One `AltlocSelectionMappingItem` per residue with
            resolved altlocs.
    """

    items: list[AltlocSelectionMappingItem] = Field(default_factory=list)


class AssemblyMappingItem(BaseModel):
    """Original-to-canonical id mapping for one assembly.

    Attributes:
        canonical_assembly_id: The canonicalised assembly id.
        original_assembly_id: The assembly id as originally deposited.
    """

    canonical_assembly_id: str
    original_assembly_id: str


class AssemblyMapping(BaseModel):
    """Original-to-canonical id mapping for every assembly.

    Attributes:
        items: One `AssemblyMappingItem` per assembly.
    """

    items: list[AssemblyMappingItem] = Field(default_factory=list)


class EntityMappingItem(BaseModel):
    """Original-to-canonical id mapping for one entity (possibly merging
    several original ids into one).

    Attributes:
        canonical_entity_id: The canonicalised entity id.
        original_entity_ids: The original entity id(s) merged into
            this canonical id.
    """

    canonical_entity_id: str
    original_entity_ids: list[str]


class EntityMapping(BaseModel):
    """Original-to-canonical id mapping for every entity.

    Attributes:
        items: One `EntityMappingItem` per canonical entity.
    """

    items: list[EntityMappingItem] = Field(default_factory=list)


class ChainIdMappingItem(BaseModel):
    """Original-to-canonical id mapping for one chain.

    Attributes:
        canonical_chain_id: The canonicalised chain id.
        original_chain_id: The `label_asym_id` as originally deposited.
        original_auth_chain_id: The `auth_asym_id` as originally
            deposited.
    """

    canonical_chain_id: str
    original_chain_id: str
    original_auth_chain_id: str


class ChainIdMapping(BaseModel):
    """Original-to-canonical id mapping for every chain.

    Attributes:
        items: One `ChainIdMappingItem` per chain.
    """

    items: list[ChainIdMappingItem] = Field(default_factory=list)


class ResidueNumberMappingItem(BaseModel):
    """Original-to-canonical residue numbering for one residue.

    Attributes:
        canonical_chain_id: The canonical chain id the residue belongs
            to.
        canonical_seq_id: The residue's canonicalised sequence number.
        original_chain_id: The chain id as originally deposited.
        original_seq_id: The `label_seq_id` as originally deposited,
            if any.
        original_auth_seq_id: The `auth_seq_id` as originally
            deposited.
        original_insertion_code: The insertion code as originally
            deposited, if any.
    """

    canonical_chain_id: str
    canonical_seq_id: int
    original_chain_id: str
    original_seq_id: int | None = None
    original_auth_seq_id: str
    original_insertion_code: str | None = None


class ResidueNumberMapping(BaseModel):
    """Original-to-canonical residue numbering for every residue.

    Attributes:
        items: One `ResidueNumberMappingItem` per residue.
    """

    items: list[ResidueNumberMappingItem] = Field(default_factory=list)


# canonicalisation rules ----------------------

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
    """Policy for how chain ids are canonicalised.

    Attributes:
        strategy: Kept as-is, remapped sequentially, or replaced by
            the author chain id.
    """

    strategy: ChainIdStrategy = "preserve"


class ResidueNumberingRules(BaseModel):
    """Policy for how residue numbering is canonicalised.

    Attributes:
        strategy: Kept as-is, replaced by the author numbering, or
            renumbered sequentially.
        preserve_insertion_codes: Whether insertion codes are kept
            when `strategy="preserve"`.
    """

    strategy: ResidueNumberingStrategy = "preserve"
    preserve_insertion_codes: bool = True


class AssemblyIdRules(BaseModel):
    """Policy for how assembly ids are canonicalised.

    Attributes:
        strategy: Kept as-is, remapped sequentially, or standardized.
    """

    strategy: AssemblyIdStrategy = "preserve"


class IdentifierRules(BaseModel):
    """Bundles the chain id, residue numbering, and assembly id policies.

    Attributes:
        chain_id: The chain id canonicalisation policy.
        residue_numbering: The residue numbering canonicalisation
            policy.
        assembly_id: The assembly id canonicalisation policy.
    """

    chain_id: ChainIdRules = Field(default_factory=ChainIdRules)
    residue_numbering: ResidueNumberingRules = Field(
        default_factory=ResidueNumberingRules
    )
    assembly_id: AssemblyIdRules = Field(default_factory=AssemblyIdRules)


class MissingAtomsRules(BaseModel):
    """Policy for handling polymer residues missing backbone atoms.

    Attributes:
        strategy: Kept as-is, flagged, dropped, or imputed.
        allow_imputation: Whether missing atoms may be imputed (not
            yet implemented).
        record_missingness: Whether missing-atom diagnostics are
            recorded.
    """

    strategy: MissingAtomsStrategy = "annotate"
    allow_imputation: bool = False
    record_missingness: bool = True


class MissingResiduesRules(BaseModel):
    """Policy for handling sequence gaps within a chain.

    Attributes:
        strategy: Kept as-is, flagged, dropped, or imputed.
        record_gaps: Whether sequence-gap diagnostics are recorded.
    """

    strategy: MissingResiduesStrategy = "annotate"
    record_gaps: bool = True


class IncompleteChainRules(BaseModel):
    """Policy for handling chains with non-contiguous residue numbering.

    Attributes:
        strategy: Kept as-is, excluded, or truncated to the longest
            contiguous region.
    """

    strategy: IncompleteChainStrategy = "preserve"


class MissingDataRules(BaseModel):
    """Bundles the missing-atoms, missing-residues, and incomplete-chain
    policies.

    Attributes:
        missing_atoms: The missing-atoms policy.
        missing_residues: The missing-residues (sequence gap) policy.
        incomplete_chains: The incomplete-chain policy.
    """

    missing_atoms: MissingAtomsRules = Field(default_factory=MissingAtomsRules)
    missing_residues: MissingResiduesRules = Field(
        default_factory=MissingResiduesRules
    )
    incomplete_chains: IncompleteChainRules = Field(
        default_factory=IncompleteChainRules
    )


class AltlocRules(BaseModel):
    """Policy for resolving altloc-disordered atoms to a single conformer.

    Attributes:
        strategy: Kept as-is, or resolved to the best-occupancy/first/
            user-defined conformer.
        tie_breaker: How ties in `select_best_occupancy` selection are
            broken.
        user_defined_altloc: The altloc code to select when
            `strategy="select_user_defined"`.
        record_selection: Whether altloc selection diagnostics are
            recorded.
    """

    strategy: AltlocStrategy = "select_best_occupancy"
    tie_breaker: AltlocTieBreaker = "alphabetical_first"
    user_defined_altloc: str | None = None
    record_selection: bool = True


class AssemblyRules(BaseModel):
    """Policy for selecting/standardizing which assemblies are kept.

    Attributes:
        strategy: Kept as reported, standardized to the biological
            assembly, or only the first assembly.
        preferred_assembly_source: Which source's "preferred assembly"
            flag to trust when standardizing.
        record_original_assembly_mapping: Whether assembly id mappings
            are recorded.
    """

    strategy: AssemblyStrategy = "preserve_as_reported"
    preferred_assembly_source: PreferredAssemblySource = "author"
    record_original_assembly_mapping: bool = True


class EntityRules(BaseModel):
    """Policy for standardizing or merging entity ids.

    Attributes:
        strategy: Kept as-is, standardized to sequential ids, or
            merged when sequence-identical.
        preserve_original_entity_ids: Whether original entity ids are
            kept in the mapping when merging.
    """

    strategy: EntityStrategy = "preserve"
    preserve_original_entity_ids: bool = True


class LigandRules(BaseModel):
    """Policy for filtering non-polymer ligands, waters, and ions.

    Attributes:
        strategy: Kept as-is, filtered by `keep_waters`/`keep_ions`/
            `keep_nonpolymer_ligands`, or excluded from the canonical
            structure but annotated.
        keep_waters: Whether water molecules are kept.
        keep_ions: Whether ions are kept.
        keep_nonpolymer_ligands: Whether other non-polymer ligands are
            kept.
    """

    strategy: LigandStrategy = "preserve"
    keep_waters: bool = True
    keep_ions: bool = True
    keep_nonpolymer_ligands: bool = True


class ValidationRules(BaseModel):
    """Policy for how canonicalisation validation issues are treated.

    Attributes:
        strictness: How strictly validation issues are treated
            ("permissive" is not yet fully implemented).
        fail_on_unresolved_issues: Whether an unresolved validation
            error raises instead of just being recorded.
        warnings_as_errors: Whether warnings are promoted to errors.
    """

    strictness: ValidationStrictness = "moderate"
    fail_on_unresolved_issues: bool = True
    warnings_as_errors: bool = False


# canonicalisation classes --------------------


class canonicalisationProvenanceRules(BaseModel):
    """Policy for which canonicalisation provenance is recorded.

    Attributes:
        record_original_mappings: Whether id/numbering mappings are
            recorded.
        record_transforms: Whether the list of applied transforms is
            recorded.
        record_policy_application: Whether the policy id/name/version
            is recorded.
        emit_canonicalisation_report: Whether a warning/error count
            summary is recorded in provenance.
    """

    record_original_mappings: bool = True
    record_transforms: bool = True
    record_policy_application: bool = True
    emit_canonicalisation_report: bool = False


class CanonicalMappings(BaseModel):
    """Every id/numbering mapping produced by one canonicalisation run.

    Attributes:
        chain_id_mapping: The chain id mapping produced by this run.
        residue_number_mapping: The residue number mapping produced by
            this run.
        assembly_mapping: The assembly id mapping produced by this
            run.
        entity_mapping: The entity id mapping produced by this run.
        altloc_selection_mapping: The altloc selections made during
            this run.
    """

    chain_id_mapping: ChainIdMapping = Field(default_factory=ChainIdMapping)
    residue_number_mapping: ResidueNumberMapping = Field(
        default_factory=ResidueNumberMapping
    )
    assembly_mapping: AssemblyMapping = Field(default_factory=AssemblyMapping)
    entity_mapping: EntityMapping = Field(default_factory=EntityMapping)
    altloc_selection_mapping: AltlocSelectionMapping = Field(
        default_factory=AltlocSelectionMapping
    )


class canonicalisationPolicy(BaseModel):
    """The full, named policy governing one canonicalisation run: every
    rule group plus provenance settings.

    Attributes:
        policy_id: A unique identifier for this policy.
        policy_name: A human-readable name for this policy.
        policy_version: The policy's version string.
        description: A free-text description of the policy.
        identifier_rules: The chain id / residue numbering / assembly
            id rules.
        missing_data_rules: The missing-atoms / missing-residues /
            incomplete-chain rules.
        altloc_rules: The altloc resolution rules.
        assembly_rules: The assembly selection rules.
        entity_rules: The entity id rules.
        ligand_rules: The ligand/water/ion filtering rules.
        validation_rules: The validation strictness rules.
        provenance_rules: Which provenance is recorded for this run.
    """

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


class canonicalisationProvenance(BaseModel):
    """Record of one canonicalisation run: which policy ran, what it
    changed, and (optionally) a diagnostics summary.

    Attributes:
        canonicalised_at: When this canonicalisation run completed, as
            an ISO 8601 timestamp.
        policy_id: The id of the policy that ran.
        policy_name: The name of the policy that ran.
        policy_version: The version of the policy that ran.
        transforms: The transform labels applied (e.g.
            "chain_id:remap"), one per rule group that deviated from
            "preserve".
        report: A summary of warning/error counts, if
            `provenance_rules.emit_canonicalisation_report` was True.
    """

    canonicalised_at: str
    policy_id: str
    policy_name: str
    policy_version: str
    transforms: list[str] = Field(default_factory=list)
    report: dict[str, Any] = Field(default_factory=dict)
