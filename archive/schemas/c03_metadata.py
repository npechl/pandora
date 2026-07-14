from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field

from .common import AppliedPolicyRef, DiagnosticBundle, ResultStatus
from .c02_canonicalisation import CanonicalStructureResult


# ── Metadata category schemas ──────────────────────────────────────────────────


class UnitCell(BaseModel):
    a: float | None = None
    b: float | None = None
    c: float | None = None
    alpha: float | None = None
    beta: float | None = None
    gamma: float | None = None


class ArchiveMetadata(BaseModel):
    entry_id: str
    title: str | None = None
    experimental_method: str | None = None
    resolution: float | None = None
    r_factor: float | None = None
    r_free: float | None = None
    deposition_date: str | None = None
    release_date: str | None = None
    revision_date: str | None = None
    keywords: list[str] = Field(default_factory=list)
    space_group: str | None = None
    unit_cell: UnitCell = Field(default_factory=UnitCell)


class UniProtMapping(BaseModel):
    canonical_chain_id: str
    uniprot_accession: str
    uniprot_id: str | None = None
    segment_start: int
    segment_end: int
    identity: float | None = None


class SIFTSResidueLevelMapping(BaseModel):
    canonical_chain_id: str
    canonical_seq_id: int
    uniprot_position: int | None = None
    uniprot_accession: str | None = None
    pfam_id: str | None = None


class SIFTSMapping(BaseModel):
    canonical_chain_id: str
    uniprot_accession: str | None = None
    pfam_id: str | None = None
    cath_id: str | None = None
    scop_id: str | None = None
    interpro_id: str | None = None
    residue_level_mappings: list[SIFTSResidueLevelMapping] | None = None


class TaxonomyRecord(BaseModel):
    ncbi_taxon_id: int | None = None
    organism_scientific: str | None = None
    organism_common: str | None = None
    lineage: list[str] = Field(default_factory=list)


class BiologicalMappings(BaseModel):
    uniprot_mappings: list[UniProtMapping] = Field(default_factory=list)
    sifts_mappings: list[SIFTSMapping] = Field(default_factory=list)
    taxonomy: TaxonomyRecord | None = None


class EntityAnnotation(BaseModel):
    canonical_entity_id: str
    entity_type: Literal["polymer", "non-polymer", "water", "branched"]
    description: str | None = None


class ChainAnnotation(BaseModel):
    canonical_chain_id: str
    chain_length: int | None = None
    is_polymer: bool
    uniprot_accession: str | None = None


class StructuralAnnotations(BaseModel):
    preferred_assembly_id: str | None = None
    entity_annotations: list[EntityAnnotation] = Field(default_factory=list)
    chain_annotations: list[ChainAnnotation] = Field(default_factory=list)


RetrievalStatus = Literal["success", "partial", "failed", "disabled"]


class MetadataSourceRecord(BaseModel):
    source_name: Literal["pdbe", "sifts", "uniprot", "taxonomy"]
    retrieval_status: RetrievalStatus
    source_version: str | None = None
    retrieved_at: str | None = None
    base_url: str | None = None


class ProvenanceMetadata(BaseModel):
    sources: list[MetadataSourceRecord] = Field(default_factory=list)


class MetadataRetrievalStatus(BaseModel):
    archive_metadata: RetrievalStatus | None = None
    biological_mappings: RetrievalStatus | None = None
    structural_annotations: RetrievalStatus | None = None


class MetadataRecord(BaseModel):
    entry_id: str
    archive_metadata: ArchiveMetadata | None = None
    biological_mappings: BiologicalMappings | None = None
    structural_annotations: StructuralAnnotations | None = None
    provenance_metadata: ProvenanceMetadata = Field(
        default_factory=ProvenanceMetadata
    )
    retrieval_status: MetadataRetrievalStatus = Field(
        default_factory=MetadataRetrievalStatus
    )
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)


# ── Annotation layer & plugin ──────────────────────────────────────────────────


class AnnotationLayer(BaseModel):
    layer_name: str
    layer_type: str
    schema_version: str | None = None
    data: dict[str, Any] = Field(default_factory=dict)


class Plugin(BaseModel):
    plugin_id: str
    plugin_name: str
    plugin_version: str
    plugin_type: str
    input_type: Literal["canonical_structure", "metadata_annotated_structure"]
    output_layer_type: str
    description: str | None = None
    config: dict[str, Any] | None = None


class AnnotationPluginResult(BaseModel):
    plugin_id: str
    plugin_name: str
    plugin_version: str
    status: ResultStatus
    annotation_layer: AnnotationLayer | None = None
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: dict[str, str | None] = Field(default_factory=dict)


# ── Output schemas ─────────────────────────────────────────────────────────────


class MetadataAnnotations(BaseModel):
    archive_metadata: ArchiveMetadata | None = None
    biological_mappings: BiologicalMappings | None = None
    structural_annotations: StructuralAnnotations | None = None
    provenance_metadata: ProvenanceMetadata = Field(
        default_factory=ProvenanceMetadata
    )


class MetadataAnnotatedStructureProvenance(BaseModel):
    metadata_sources: list[str] = Field(default_factory=list)
    retrieved_at: str | None = None


class MetadataAnnotatedStructure(BaseModel):
    canonical_structure_result: CanonicalStructureResult
    metadata_annotations: MetadataAnnotations = Field(
        default_factory=MetadataAnnotations
    )
    applied_metadata_policy: AppliedPolicyRef
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: MetadataAnnotatedStructureProvenance = Field(
        default_factory=MetadataAnnotatedStructureProvenance
    )


class AnnotatedStructureProvenance(BaseModel):
    metadata_sources: list[str] = Field(default_factory=list)
    retrieved_at: str | None = None
    annotation_history: list[str] = Field(default_factory=list)


class AnnotatedStructureWithPlugins(BaseModel):
    canonical_structure_result: CanonicalStructureResult
    metadata_annotations: MetadataAnnotations = Field(
        default_factory=MetadataAnnotations
    )
    derived_annotations: list[AnnotationLayer] = Field(default_factory=list)
    applied_metadata_policy: AppliedPolicyRef
    applied_plugins: list[Plugin] = Field(default_factory=list)
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: AnnotatedStructureProvenance = Field(
        default_factory=AnnotatedStructureProvenance
    )


class MetadataAndAnnotationBatchSummary(BaseModel):
    total: int
    success: int
    warning: int
    failed: int


class MetadataAndAnnotationBatchResultItem(BaseModel):
    entry_id: str
    status: ResultStatus
    annotated_structure: AnnotatedStructureWithPlugins | None = None
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)


class MetadataAndAnnotationBatchResult(BaseModel):
    mode: Literal["sequential", "parallel"]
    summary: MetadataAndAnnotationBatchSummary
    results: list[MetadataAndAnnotationBatchResultItem]


# ── Policy schemas ─────────────────────────────────────────────────────────────


class MetadataIncludeSources(BaseModel):
    pdbe: bool = True
    sifts: bool = True
    uniprot: bool = True
    taxonomy: bool = True


class MetadataIncludeCategories(BaseModel):
    archive_metadata: bool = True
    biological_mappings: bool = True
    structural_annotations: bool = True
    provenance_metadata: bool = True


class MetadataRetrievalRules(BaseModel):
    eager: bool = True
    use_cache: bool = True
    fail_on_missing: bool = False
    allow_partial: bool = True


class MetadataProvenanceRules(BaseModel):
    record_source_versions: bool = True
    record_retrieval_time: bool = True
    record_mapping_history: bool = True


class MetadataIntegrationPolicy(BaseModel):
    policy_id: str
    policy_name: str
    policy_version: str
    description: str = ""
    include_sources: MetadataIncludeSources = Field(
        default_factory=MetadataIncludeSources
    )
    include_categories: MetadataIncludeCategories = Field(
        default_factory=MetadataIncludeCategories
    )
    retrieval_rules: MetadataRetrievalRules = Field(
        default_factory=MetadataRetrievalRules
    )
    provenance_rules: MetadataProvenanceRules = Field(
        default_factory=MetadataProvenanceRules
    )


class PluginExecutionParallelOptions(BaseModel):
    max_workers: int | None = None
    fail_fast: bool = False


class PluginExecutionRules(BaseModel):
    mode: Literal["sequential", "parallel"] = "sequential"
    parallel_options: PluginExecutionParallelOptions = Field(
        default_factory=PluginExecutionParallelOptions
    )
    fail_on_plugin_error: bool = False
    allow_partial_annotation: bool = True


class PluginProvenanceRules(BaseModel):
    record_plugin_versions: bool = True
    record_plugin_configuration: bool = True
    record_execution_history: bool = True


class AnnotationPluginPolicy(BaseModel):
    policy_id: str
    policy_name: str
    policy_version: str
    description: str = ""
    execution_rules: PluginExecutionRules = Field(
        default_factory=PluginExecutionRules
    )
    provenance_rules: PluginProvenanceRules = Field(
        default_factory=PluginProvenanceRules
    )
