from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field

from .common import AppliedPolicyRef, Diagnostic, DiagnosticBundle, ResultStatus
from .c04_curation import AppliedFilterRecord, ExclusionRecord
from .c05_splitting import LeakageSafeDataset


# ── Per-stage provenance sub-schemas ──────────────────────────────────────────


class TimestampRange(BaseModel):
    earliest: str | None = None
    latest: str | None = None


class PipelineIngestionProvenance(BaseModel):
    providers: list[str] = Field(default_factory=list)
    source_uris: list[str | None] = Field(default_factory=list)
    retrieved_at_range: TimestampRange = Field(default_factory=TimestampRange)
    from_cache: bool | None = None


class PipelineCanonicalizationProvenance(BaseModel):
    policy_id: str
    policy_name: str
    policy_version: str
    canonicalized_at_range: TimestampRange = Field(
        default_factory=TimestampRange
    )


class PipelineMetadataProvenance(BaseModel):
    metadata_policy_id: str
    metadata_policy_name: str
    metadata_policy_version: str
    sources_used: list[str] = Field(default_factory=list)
    retrieved_at_range: TimestampRange = Field(default_factory=TimestampRange)


class PluginProvenanceRecord(BaseModel):
    plugin_id: str
    plugin_name: str
    plugin_version: str
    output_layer_type: str
    config: dict[str, Any] | None = None


class PipelineAnnotationProvenance(BaseModel):
    annotation_plugin_policy_id: str | None = None
    annotation_plugin_policy_version: str | None = None
    plugin_records: list[PluginProvenanceRecord] = Field(default_factory=list)


class SplitRatios(BaseModel):
    train: float
    validation: float
    test: float


class PipelineCurationProvenance(BaseModel):
    curation_policy_id: str
    curation_policy_name: str
    curation_policy_version: str
    granularity: Literal["structure", "chain", "interface", "residue"]
    total_input: int
    total_selected: int
    total_excluded: int
    total_duplicates_removed: int
    filters_applied: list[AppliedFilterRecord] = Field(default_factory=list)
    exclusions: list[ExclusionRecord] = Field(default_factory=list)
    deduplication_strategy: str | None = None


class PipelineSplittingProvenance(BaseModel):
    leakage_policy_id: str
    leakage_policy_name: str
    leakage_policy_version: str
    similarity_engines: list[str] = Field(default_factory=list)
    similarity_engine_versions: list[str | None] = Field(default_factory=list)
    sequence_similarity_threshold: float | None = None
    structure_similarity_threshold: float | None = None
    clustering_strategy: str | None = None
    split_ratios_requested: SplitRatios
    split_ratios_achieved: SplitRatios
    leakage_detected: bool = False
    split_at: str | None = None


class PipelineProvenance(BaseModel):
    """
    Per-stage provenance. For chain/interface/residue granularities,
    ingestion, canonicalization, metadata, and annotation fields are null
    because upstream AnnotatedStructureWithPlugins objects are not embedded
    (UPSTREAM_PROVENANCE_NOT_EMBEDDED).
    """

    ingestion: PipelineIngestionProvenance | None = None
    canonicalization: PipelineCanonicalizationProvenance | None = None
    metadata_integration: PipelineMetadataProvenance | None = None
    annotation_plugins: PipelineAnnotationProvenance | None = None
    dataset_curation: PipelineCurationProvenance | None = None
    leakage_splitting: PipelineSplittingProvenance | None = None


# ── Source releases & policy versions ─────────────────────────────────────────


class ExternalSourceRelease(BaseModel):
    source_name: str
    release_version: str | None = None
    release_date: str | None = None


class SourceReleaseProvenanceRecord(BaseModel):
    pdbe_release: str | None = None
    sifts_release: str | None = None
    uniprot_release: str | None = None
    other_sources: list[ExternalSourceRelease] = Field(default_factory=list)


class PolicyProvenanceRecord(BaseModel):
    canonicalization_policy_id: str | None = None
    canonicalization_policy_version: str | None = None
    metadata_policy_id: str | None = None
    metadata_policy_version: str | None = None
    annotation_plugin_policy_id: str | None = None
    annotation_plugin_policy_version: str | None = None
    curation_policy_id: str | None = None
    curation_policy_version: str | None = None
    leakage_policy_id: str | None = None
    leakage_policy_version: str | None = None
    provenance_policy_id: str | None = None
    provenance_policy_version: str | None = None


# ── Provenance bundle ──────────────────────────────────────────────────────────


class ProvenanceBundle(BaseModel):
    pipeline_provenance: PipelineProvenance = Field(
        default_factory=PipelineProvenance
    )
    source_release_provenance: SourceReleaseProvenanceRecord = Field(
        default_factory=SourceReleaseProvenanceRecord
    )
    policy_provenance: PolicyProvenanceRecord = Field(
        default_factory=PolicyProvenanceRecord
    )


# ── Manifest ───────────────────────────────────────────────────────────────────


class DatasetSummary(BaseModel):
    dataset_id: str
    dataset_version: str
    granularity: Literal["structure", "chain", "interface", "residue"]
    total_items: int
    train_count: int
    validation_count: int
    test_count: int
    train_fraction_achieved: float
    validation_fraction_achieved: float
    test_fraction_achieved: float


class ManifestChecksums(BaseModel):
    artifact_checksum: str | None = None
    manifest_checksum: str | None = None
    split_checksum: str | None = None
    checksum_algorithm: str = "SHA-256"


class PandoraManifest(BaseModel):
    manifest_id: str
    manifest_format: str
    pandora_version: str
    generated_at: str
    artifact_id: str
    dataset_summary: DatasetSummary
    source_releases: SourceReleaseProvenanceRecord = Field(
        default_factory=SourceReleaseProvenanceRecord
    )
    policies: PolicyProvenanceRecord = Field(
        default_factory=PolicyProvenanceRecord
    )
    pipeline_steps: list[str] = Field(default_factory=list)
    checksums: ManifestChecksums = Field(default_factory=ManifestChecksums)


# ── Reproducibility report ─────────────────────────────────────────────────────


class ReproducibilityReportSummary(BaseModel):
    pipeline_steps: int
    source_count: int
    policy_count: int
    plugin_count: int


class ReproducibilityReport(BaseModel):
    report_id: str
    artifact_id: str
    summary: ReproducibilityReportSummary
    lineage: list[str] = Field(default_factory=list)
    reproducibility_risks: list[Diagnostic] = Field(default_factory=list)


# ── PandoraArtifact ────────────────────────────────────────────────────────────


class ArtifactChecksums(BaseModel):
    artifact_checksum: str | None = None
    manifest_checksum: str | None = None
    split_checksum: str | None = None
    checksum_algorithm: str = "SHA-256"


class ArtifactProvenance(BaseModel):
    generated_at: str
    pandora_version: str


class PandoraArtifact(BaseModel):
    artifact_id: str
    artifact_name: str | None = None
    leakage_safe_dataset: LeakageSafeDataset
    provenance_bundle: ProvenanceBundle
    manifest: PandoraManifest
    checksums: ArtifactChecksums = Field(default_factory=ArtifactChecksums)
    reproducibility_report: ReproducibilityReport | None = None
    applied_policy: AppliedPolicyRef
    provenance: ArtifactProvenance


class ReproducibilityBatchSummary(BaseModel):
    total: int
    success: int
    warning: int
    failed: int


class ReproducibilityBatchResultItem(BaseModel):
    artifact_id: str
    status: ResultStatus
    artifact: PandoraArtifact | None = None
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)


class ReproducibilityBatchResult(BaseModel):
    mode: Literal["sequential", "parallel"]
    summary: ReproducibilityBatchSummary
    results: list[ReproducibilityBatchResultItem]


# ── Policy schemas ─────────────────────────────────────────────────────────────


class ProvenancePolicy(BaseModel):
    policy_id: str
    policy_name: str
    policy_version: str
    description: str = ""
    record_software_versions: bool = True
    record_policy_versions: bool = True
    record_source_releases: bool = True
    record_annotation_plugin_versions: bool = True
    record_curation_history: bool = True
    record_split_history: bool = True
    record_checksums: bool = True
    record_execution_timestamps: bool = True


class ExportPolicy(BaseModel):
    emit_manifest_yaml: bool = True
    emit_manifest_json: bool = False
    emit_provenance_report: bool = True
    emit_checksum_bundle: bool = True
    emit_lineage_graph: bool = True
