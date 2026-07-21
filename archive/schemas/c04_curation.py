from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field

from .common import AppliedPolicyRef, DiagnosticBundle, ResultStatus
from .c01_ingestion import Atom, Residue
from .c03_metadata import (
    AnnotatedStructureWithPlugins,
    AnnotationLayer,
    ArchiveMetadata,
    SIFTSMapping,
    TaxonomyRecord,
    UniProtMapping,
)

# ── Exclusion & deduplication records ─────────────────────────────────────────


class ExclusionRecord(BaseModel):
    unit_id: str
    granularity: Literal["structure", "chain", "interface", "residue"]
    reason_code: str
    reason_message: str
    applied_rule: str | None = None


class DeduplicationReport(BaseModel):
    enabled: bool
    strategy: Literal["entry_id", "exact_hash"]
    duplicates_found: int = 0
    removed_items: list[ExclusionRecord] = Field(default_factory=list)


class AppliedFilterRecord(BaseModel):
    filter_name: str
    filter_category: Literal[
        "quality",
        "selection",
        "content",
        "organism",
        "deduplication",
        "extraction",
    ]
    filter_value: dict[str, Any] = Field(default_factory=dict)
    structures_excluded: int


class SelectionSummary(BaseModel):
    applied_filters: list[AppliedFilterRecord] = Field(default_factory=list)


# ── Chain, interface, residue records ─────────────────────────────────────────


class ChainRecordMetadata(BaseModel):
    archive_metadata: ArchiveMetadata | None = None
    uniprot_mappings: list[UniProtMapping] = Field(default_factory=list)
    sifts_mappings: list[SIFTSMapping] = Field(default_factory=list)
    taxonomy: TaxonomyRecord | None = None


class ChainRecord(BaseModel):
    chain_id: str
    entry_id: str
    chain_type: Literal["protein", "rna", "dna"]
    entity_id: str
    sequence: str
    chain_length: int
    residues: list[Residue] = Field(default_factory=list)
    metadata: ChainRecordMetadata = Field(default_factory=ChainRecordMetadata)
    derived_annotations: list[AnnotationLayer] = Field(default_factory=list)
    parent_entry_id: str
    applied_policy: AppliedPolicyRef


class InterfaceRecordMetadata(BaseModel):
    archive_metadata: ArchiveMetadata | None = None


class InterfaceRecord(BaseModel):
    interface_id: str
    entry_id: str
    chain_id_1: str
    chain_id_2: str
    interface_type: Literal[
        "protein_protein", "protein_rna", "protein_dna", "protein_ligand"
    ]
    chain_record_1: ChainRecord
    chain_record_2: ChainRecord
    interface_residues_chain_1: list[str] = Field(default_factory=list)
    interface_residues_chain_2: list[str] = Field(default_factory=list)
    contact_count: int
    interface_area: float | None = None
    source_annotation_layer: AnnotationLayer | None = None
    metadata: InterfaceRecordMetadata = Field(
        default_factory=InterfaceRecordMetadata
    )
    parent_entry_id: str
    applied_policy: AppliedPolicyRef


class NeighborReference(BaseModel):
    residue_id: str
    chain_id: str
    distance: float


class ResidueRecordMetadata(BaseModel):
    archive_metadata: ArchiveMetadata | None = None


class ResidueRecord(BaseModel):
    residue_id: str
    entry_id: str
    chain_id: str
    seq_id: int | None = None
    comp_id: str
    residue_type: Literal["amino_acid", "nucleotide", "non_standard"]
    atoms: list[Atom] = Field(default_factory=list)
    neighboring_residues: list[NeighborReference] | None = None
    metadata: ResidueRecordMetadata = Field(
        default_factory=ResidueRecordMetadata
    )
    derived_annotations: list[AnnotationLayer] = Field(default_factory=list)
    parent_entry_id: str
    parent_chain_id: str
    applied_policy: AppliedPolicyRef


# ── Dataset output schemas ─────────────────────────────────────────────────────


class DatasetCounts(BaseModel):
    total_input: int
    total_selected: int
    total_excluded: int
    total_duplicates_removed: int


class DatasetProvenance(BaseModel):
    created_at: str | None = None
    source_count: int
    input_sources: list[str] = Field(default_factory=list)


class Dataset(BaseModel):
    dataset_id: str
    dataset_name: str
    dataset_version: str
    granularity: Literal["structure"] = "structure"
    structures: list[AnnotatedStructureWithPlugins] = Field(
        default_factory=list
    )
    counts: DatasetCounts
    selection_summary: SelectionSummary = Field(
        default_factory=SelectionSummary
    )
    excluded_items: list[ExclusionRecord] = Field(default_factory=list)
    deduplication_report: DeduplicationReport
    applied_policy: AppliedPolicyRef
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: DatasetProvenance


class ChainDatasetCounts(BaseModel):
    total_structures_input: int
    total_chains_extracted: int
    total_chains_excluded: int


class ChainDatasetProvenance(BaseModel):
    created_at: str
    source_dataset_id: str


class ChainDataset(BaseModel):
    dataset_id: str
    dataset_name: str
    dataset_version: str
    granularity: Literal["chain"] = "chain"
    chains: list[ChainRecord] = Field(default_factory=list)
    source_dataset_id: str
    counts: ChainDatasetCounts
    excluded_items: list[ExclusionRecord] = Field(default_factory=list)
    applied_policy: AppliedPolicyRef
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: ChainDatasetProvenance


class InterfaceDatasetCounts(BaseModel):
    total_structures_input: int
    total_interfaces_extracted: int
    total_interfaces_excluded: int


class InterfaceDatasetProvenance(BaseModel):
    created_at: str
    source_dataset_id: str


class InterfaceDataset(BaseModel):
    dataset_id: str
    dataset_name: str
    dataset_version: str
    granularity: Literal["interface"] = "interface"
    interfaces: list[InterfaceRecord] = Field(default_factory=list)
    source_dataset_id: str
    counts: InterfaceDatasetCounts
    excluded_items: list[ExclusionRecord] = Field(default_factory=list)
    applied_policy: AppliedPolicyRef
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: InterfaceDatasetProvenance


class ResidueDatasetCounts(BaseModel):
    total_source_units_input: int
    total_residues_extracted: int
    total_residues_excluded: int


class ResidueDatasetProvenance(BaseModel):
    created_at: str
    source_dataset_id: str


class ResidueDataset(BaseModel):
    dataset_id: str
    dataset_name: str
    dataset_version: str
    granularity: Literal["residue"] = "residue"
    residues: list[ResidueRecord] = Field(default_factory=list)
    source_dataset_id: str
    source_granularity: Literal["structure", "chain"]
    counts: ResidueDatasetCounts
    excluded_items: list[ExclusionRecord] = Field(default_factory=list)
    applied_policy: AppliedPolicyRef
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: ResidueDatasetProvenance


class DatasetConstructionBatchSummary(BaseModel):
    total: int
    success: int
    warning: int
    failed: int


class DatasetConstructionBatchResultItem(BaseModel):
    dataset_id: str
    status: ResultStatus
    dataset: Dataset | None = None
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)


class DatasetConstructionBatchResult(BaseModel):
    mode: Literal["sequential", "parallel"]
    summary: DatasetConstructionBatchSummary
    results: list[DatasetConstructionBatchResultItem]


# ── Policy schema ──────────────────────────────────────────────────────────────


class IncludeBiomolecules(BaseModel):
    proteins: bool = True
    rna: bool = False
    dna: bool = False
    complexes: bool = True


class SelectionRules(BaseModel):
    include_sources: list[str] = Field(default_factory=list)
    exclude_sources: list[str] = Field(default_factory=list)
    include_biomolecules: IncludeBiomolecules = Field(
        default_factory=IncludeBiomolecules
    )
    include_experimental_methods: list[str] = Field(default_factory=list)
    exclude_experimental_methods: list[str] = Field(default_factory=list)


class QualityRules(BaseModel):
    max_resolution: float | None = None
    null_resolution_behavior: Literal["exclude", "include"] = "exclude"
    min_chain_length: int | None = None
    allow_incomplete_chains: bool = True
    allow_missing_residues: bool = True
    allow_missing_atoms: bool = True


class ContentRules(BaseModel):
    keep_ligands: bool = True
    keep_waters: bool = False
    keep_ions: bool = True
    keep_nonpolymer_entities: bool = True


class OrganismRules(BaseModel):
    include_taxa: list[str] = Field(default_factory=list)
    exclude_taxa: list[str] = Field(default_factory=list)


class DeduplicationRules(BaseModel):
    enabled: bool = True
    strategy: Literal["entry_id", "exact_hash"] = "entry_id"


class IncludeChainTypes(BaseModel):
    proteins: bool = True
    rna: bool = False
    dna: bool = False


class ChainExtractionRules(BaseModel):
    include_chain_types: IncludeChainTypes = Field(
        default_factory=IncludeChainTypes
    )
    min_chain_length: int | None = None


class InterfaceTypes(BaseModel):
    protein_protein: bool = True
    protein_rna: bool = False
    protein_dna: bool = False
    protein_ligand: bool = False


class InterfaceExtractionRules(BaseModel):
    interface_types: InterfaceTypes = Field(default_factory=InterfaceTypes)
    min_contact_residues: int | None = None
    min_interface_area: float | None = None


class IncludeResidueTypes(BaseModel):
    standard_amino_acids: bool = True
    non_standard: bool = False
    nucleotides: bool = False


class ResidueExtractionRules(BaseModel):
    source_granularity: Literal["structure", "chain"] = "chain"
    include_residue_types: IncludeResidueTypes = Field(
        default_factory=IncludeResidueTypes
    )
    require_full_backbone: bool = True
    context_radius: float | None = None


class ExtractionRules(BaseModel):
    granularity: Literal["structure", "chain", "interface", "residue"] = (
        "structure"
    )
    chain_extraction_rules: ChainExtractionRules = Field(
        default_factory=ChainExtractionRules
    )
    interface_extraction_rules: InterfaceExtractionRules = Field(
        default_factory=InterfaceExtractionRules
    )
    residue_extraction_rules: ResidueExtractionRules = Field(
        default_factory=ResidueExtractionRules
    )


class CurationProvenanceRules(BaseModel):
    record_filters: bool = True
    record_exclusions: bool = True
    record_dataset_version: bool = True
    record_input_sources: bool = True


class DatasetCurationPolicy(BaseModel):
    policy_id: str
    policy_name: str
    policy_version: str
    description: str = ""
    selection_rules: SelectionRules = Field(default_factory=SelectionRules)
    quality_rules: QualityRules = Field(default_factory=QualityRules)
    content_rules: ContentRules = Field(default_factory=ContentRules)
    organism_rules: OrganismRules = Field(default_factory=OrganismRules)
    deduplication_rules: DeduplicationRules = Field(
        default_factory=DeduplicationRules
    )
    extraction_rules: ExtractionRules = Field(default_factory=ExtractionRules)
    provenance_rules: CurationProvenanceRules = Field(
        default_factory=CurationProvenanceRules
    )
