# Component 06 — Provenance & Reproducibility Layer

## Purpose

The Provenance & Reproducibility Layer aggregates the provenance generated
by all previous Pandora components and turns a leakage-safe dataset into a
fully traceable, reproducible Pandora artifact.

This component is responsible for:

* aggregating provenance across the pipeline,
* generating reproducibility manifests,
* recording policy and software versions,
* recording source release versions,
* computing integrity checksums,
* tracking lineage from raw archive data to final split artifacts,
* and exporting reproducibility bundles.

This component does **not** perform ingestion, canonicalization, metadata
retrieval, annotation, curation, or splitting. It only records, assembles,
and exports the provenance of those steps.

---

# 1. Architectural Role

```text
MmCIFIngestionResult
  → CanonicalStructureResult
  → MetadataAnnotatedStructure
  → AnnotatedStructureWithPlugins
  → Dataset
  → LeakageSafeDataset
  → PandoraArtifact
```

---

# 2. Core Design Principles

## Provenance is cross-cutting

Every preceding component contributes provenance data.

Component 6 centralises that information into one reproducible artifact by
traversing the `LeakageSafeDataset` object hierarchy.

---

## Policies and source releases are first-class

Pandora records:

* software version,
* policy IDs and versions (all six components),
* metadata source versions,
* archive release versions,
* annotation plugin versions,
* and split policy versions.

---

## Checksums are integrity guards

Checksums allow users to verify that a manifest or exported artifact has
not been modified after generation.

---

## Dataset versioning is optional

Pandora prefers lineage and policy provenance over manual dataset version
numbering. The identity of a generated artifact is primarily defined by
source releases, policies, software version, and transformations applied.

---

## Artifact embedding contract

In V1, `PandoraArtifact.leakage_safe_dataset` is embedded **by value** —
the full nested object is stored, including `source_dataset.structures`.
This makes the artifact self-contained and shareable without external
dependencies, but may produce large files for datasets with many structures.

Users who need lightweight artifacts may serialise structure data separately
and retain only the partition entry_id lists. This is a user-side concern;
Pandora does not implement a by-reference artifact mode in V1.

---

# 3. Input Schemas

## 3.1 Reproducibility input

```yaml
ReproducibilityInput:
  leakage_safe_dataset: LeakageSafeDataset

  provenance_policy: ProvenancePolicy
  # See Section 10.1 for the full policy schema.

  export_policy: ExportPolicy
  # See Section 10.2 for the full policy schema.
```

---

## 3.2 Batch reproducibility input

```yaml
ReproducibilityBatchInput:
  leakage_safe_datasets:
    - LeakageSafeDataset

  provenance_policy: ProvenancePolicy
  export_policy: ExportPolicy

  mode: string
  # sequential | parallel

  parallel_options:
    max_workers: int | null
    # Number of concurrent artifact-building jobs in parallel mode.
    # null uses the system default.
    fail_fast: bool
    # If true, abort remaining jobs on the first failure.
    # If false (default), isolate failures and continue.
```

---

# 4. Provenance Sub-Schemas

These schemas define the typed content of each provenance category within
`ProvenanceBundle`. Each is populated by traversing the corresponding
fields in the `LeakageSafeDataset` hierarchy.

## 4.1 PipelineIngestionProvenance

Sourced from `AnnotatedStructureWithPlugins.canonical_structure_result.provenance`
(representative sample — one record per unique provider/source_uri pair).

```yaml
PipelineIngestionProvenance:
  providers: list[string]
  # e.g. ["pdbe"]
  source_uris: list[string | null]
  retrieved_at_range:
    earliest: string | null
    latest: string | null
  # ISO 8601 timestamps. Range across all ingested structures.
  from_cache: bool | null
  # null if mixed (some from cache, some fetched live).
```

---

## 4.2 PipelineCanonicalizationProvenance

Sourced from `AnnotatedStructureWithPlugins.canonical_structure_result`.

```yaml
PipelineCanonicalizationProvenance:
  policy_id: string
  policy_name: string
  policy_version: string
  canonicalized_at_range:
    earliest: string | null
    latest: string | null
  # ISO 8601 timestamps.
```

---

## 4.3 PipelineMetadataProvenance

Sourced from `AnnotatedStructureWithPlugins.metadata_annotations.provenance_metadata`.

```yaml
PipelineMetadataProvenance:
  metadata_policy_id: string
  metadata_policy_name: string
  metadata_policy_version: string
  sources_used: list[string]
  # e.g. ["pdbe", "sifts", "uniprot", "taxonomy"]
  retrieved_at_range:
    earliest: string | null
    latest: string | null
```

---

## 4.4 PipelineAnnotationProvenance

Sourced from `AnnotatedStructureWithPlugins.applied_plugins` and
`derived_annotations`.

```yaml
PipelineAnnotationProvenance:
  annotation_plugin_policy_id: string | null
  annotation_plugin_policy_version: string | null

  plugin_records: list[PluginProvenanceRecord]

PluginProvenanceRecord:
  plugin_id: string
  plugin_name: string
  plugin_version: string
  output_layer_type: string
  config: object | null
```

---

## 4.5 PipelineCurationProvenance

Sourced from `LeakageSafeDataset.source_dataset`.

```yaml
PipelineCurationProvenance:
  curation_policy_id: string
  curation_policy_name: string
  curation_policy_version: string
  total_input: int
  total_selected: int
  total_excluded: int
  total_duplicates_removed: int
  filters_applied: list[AppliedFilterRecord]
  # References AppliedFilterRecord from Component 04.
  exclusions: list[ExclusionRecord]
  # References ExclusionRecord from Component 04.
  deduplication_strategy: string | null
```

---

## 4.6 PipelineSplittingProvenance

Sourced from `LeakageSafeDataset` and `LeakageSafeDataset.similarity_clusters`.

```yaml
PipelineSplittingProvenance:
  leakage_policy_id: string
  leakage_policy_name: string
  leakage_policy_version: string
  similarity_engines: list[string]
  similarity_engine_versions: list[string | null]
  sequence_similarity_threshold: float | null
  structure_similarity_threshold: float | null
  clustering_strategy: string | null
  split_ratios_requested:
    train: float
    validation: float
    test: float
  split_ratios_achieved:
    train: float
    validation: float
    test: float
  leakage_detected: bool
  split_at: string | null
  # ISO 8601 timestamp.
```

---

## 4.7 SourceReleaseProvenanceRecord

```yaml
SourceReleaseProvenanceRecord:
  pdbe_release: string | null
  sifts_release: string | null
  uniprot_release: string | null
  other_sources: list[ExternalSourceRelease]

ExternalSourceRelease:
  source_name: string
  release_version: string | null
  release_date: string | null
  # ISO 8601 date.
```

---

## 4.8 PolicyProvenanceRecord

```yaml
PolicyProvenanceRecord:
  canonicalization_policy_id: string | null
  canonicalization_policy_version: string | null
  metadata_policy_id: string | null
  metadata_policy_version: string | null
  annotation_plugin_policy_id: string | null
  annotation_plugin_policy_version: string | null
  curation_policy_id: string | null
  curation_policy_version: string | null
  leakage_policy_id: string | null
  leakage_policy_version: string | null
  provenance_policy_id: string | null
  provenance_policy_version: string | null
```

---

# 5. Manifest Schema

The `PandoraManifest` is the primary machine-readable artifact shared
between researchers for reproducibility. It is a compact summary — it does
not embed full structure data.

```yaml
PandoraManifest:
  manifest_id: string
  manifest_format: string
  # yaml | json

  pandora_version: string
  generated_at: string
  # ISO 8601 timestamp.

  artifact_id: string

  dataset_summary:
    dataset_id: string
    dataset_version: string
    total_structures: int
    train_count: int
    validation_count: int
    test_count: int
    train_fraction_achieved: float
    validation_fraction_achieved: float
    test_fraction_achieved: float

  source_releases: SourceReleaseProvenanceRecord

  policies: PolicyProvenanceRecord

  pipeline_steps: list[string]
  # Ordered list of pipeline steps applied.
  # e.g. ["ingestion", "canonicalization", "metadata_integration",
  #        "annotation", "curation", "splitting"]

  checksums:
    artifact_checksum: string | null
    manifest_checksum: string | null
    split_checksum: string | null
    checksum_algorithm: string
    # Always "SHA-256" in V1.
```

---

# 6. Output Schemas

## 6.1 Pandora artifact

```yaml
PandoraArtifact:
  artifact_id: string
  artifact_name: string | null

  leakage_safe_dataset: LeakageSafeDataset
  # Embedded by value in V1. See Section 2 (Artifact embedding contract).

  provenance_bundle: ProvenanceBundle

  manifest: PandoraManifest
  # When export_policy.emit_manifest_yaml and/or emit_manifest_json are true,
  # this manifest is also serialised to disk in the requested format(s).
  # Both formats can be emitted simultaneously.

  checksums:
    artifact_checksum: string | null
    manifest_checksum: string | null
    split_checksum: string | null
    checksum_algorithm: string
    # "SHA-256" in V1.

  reproducibility_report: ReproducibilityReport | null
  # null when export_policy.emit_provenance_report: false.

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  provenance:
    generated_at: string
    # ISO 8601 timestamp. Always populated.
    pandora_version: string
    # Pandora software version. Always populated (non-nullable).
```

---

## 6.2 Provenance bundle

```yaml
ProvenanceBundle:
  pipeline_provenance:
    ingestion: PipelineIngestionProvenance | null
    # null when provenance_policy.record_execution_timestamps: false
    # and no ingestion metadata is available.

    canonicalization: PipelineCanonicalizationProvenance | null
    metadata_integration: PipelineMetadataProvenance | null
    annotation_plugins: PipelineAnnotationProvenance | null
    dataset_curation: PipelineCurationProvenance | null
    leakage_splitting: PipelineSplittingProvenance | null
    # Each field is null when the corresponding record_* flag is false
    # in ProvenancePolicy, or when data is unavailable.

  source_release_provenance: SourceReleaseProvenanceRecord
  # Empty lists / null fields when record_source_releases: false.

  policy_provenance: PolicyProvenanceRecord
  # All fields null when record_policy_versions: false.

  annotation_provenance: PipelineAnnotationProvenance | null
  # null when record_annotation_plugin_versions: false.

  curation_provenance: PipelineCurationProvenance | null
  # null when record_curation_history: false.

  splitting_provenance: PipelineSplittingProvenance | null
  # null when record_split_history: false.
```

---

## 6.3 Reproducibility report

```yaml
ReproducibilityReport:
  report_id: string
  artifact_id: string

  summary:
    pipeline_steps: int
    # Number of pipeline steps recorded in provenance.
    source_count: int
    # Number of unique source providers used.
    policy_count: int
    # Number of distinct policies applied (one per component).
    plugin_count: int
    # Number of annotation plugins applied.

  lineage: list[string]
  # Ordered list of human-readable pipeline step summaries.
  # Format: "{step_number} {step_name}: {brief description}"
  # e.g. ["01 ingestion: 150 structures from pdbe",
  #        "02 canonicalization: 150 structures (policy: default_v1)",
  #        "03 metadata: pdbe, sifts, uniprot retrieved",
  #        "04 curation: 120 selected, 30 excluded",
  #        "05 splitting: train=84, val=18, test=18 (leakage: clean)"]

  reproducibility_risks: list[Diagnostic]
  # Warnings about missing provenance fields or reproducibility weaknesses.
  # e.g. missing engine version, null source release, incomplete timestamps.
```

---

## 6.4 Batch output

```yaml
ReproducibilityBatchResult:
  mode: string
  # sequential | parallel

  summary:
    total: int
    success: int
    warning: int
    failed: int

  results:
    - artifact_id: string
      status: string
      # success | warning | failed

      artifact: PandoraArtifact | null
      # null when status == "failed".

      diagnostics:
        warnings: list[Diagnostic]
        errors: list[Diagnostic]
```

---

# 7. Public Functions

## 7.1 `build_provenance_bundle()`

### Responsibility

Aggregate provenance from the leakage-safe dataset and all upstream
pipeline stages by traversing the `LeakageSafeDataset` object hierarchy.

### Internal Workflow

```text
1. Extract pipeline provenance:
   aggregate_pipeline_provenance(leakage_safe_dataset, provenance_policy)
   Traversal path per stage:
     ingestion:     leakage_safe_dataset
                    .source_dataset.structures[*]
                    .canonical_structure_result.provenance
     canonicalization: .canonical_structure_result.applied_policy
                       .canonical_structure_result.provenance.canonicalized_at
     metadata:      .metadata_annotations.provenance_metadata.sources
     annotation:    .applied_plugins + .derived_annotations[*].plugin_id
     curation:      leakage_safe_dataset.source_dataset
                    (.applied_policy, .selection_summary, .deduplication_report)
     splitting:     leakage_safe_dataset
                    (.applied_policy, .partition_summary, .provenance)

2. Aggregate source releases (if record_source_releases: true):
   aggregate_source_releases(leakage_safe_dataset)

3. Aggregate policy versions (if record_policy_versions: true):
   aggregate_policy_versions(leakage_safe_dataset)

4. Aggregate plugin versions (if record_annotation_plugin_versions: true):
   aggregate_plugin_versions(leakage_safe_dataset)

5. Assemble and return ProvenanceBundle.
   Fields controlled by record_* flags are set to null when the
   corresponding flag is false.
```

### Input Schema

```yaml
build_provenance_bundle:
  leakage_safe_dataset: LeakageSafeDataset
  provenance_policy: ProvenancePolicy
```

### Output Schema

```yaml
build_provenance_bundle_result:
  provenance_bundle: ProvenanceBundle
```

---

## 7.2 `generate_manifest()`

### Responsibility

Generate a `PandoraManifest` describing the artifact, its dataset summary,
source releases, policy versions, and (optionally) checksums.

### Notes

When `export_policy.emit_manifest_yaml: true`, the manifest is serialised
to a `.yaml` file.

When `export_policy.emit_manifest_json: true`, the manifest is serialised
to a `.json` file.

Both flags can be `true` simultaneously — both files will be generated.
The `manifest_format` field in `PandoraManifest` records which formats
were requested (as a comma-separated string, e.g. `"yaml,json"`).

### Input Schema

```yaml
generate_manifest:
  leakage_safe_dataset: LeakageSafeDataset
  provenance_bundle: ProvenanceBundle
  export_policy: ExportPolicy
```

### Output Schema

```yaml
generate_manifest_result:
  manifest: PandoraManifest
```

---

## 7.3 `compute_checksums()`

### Responsibility

Compute SHA-256 integrity checksums for the artifact, manifest, and split
partition lists.

### Checksum contract

```yaml
checksum_algorithm: SHA-256

serialization:
  artifact_checksum:
    input: PandoraArtifact serialised to JSON with sorted keys and no
           whitespace, with the checksums.artifact_checksum field set to null
           before hashing (to avoid circular dependency).

  manifest_checksum:
    input: PandoraManifest serialised to JSON with sorted keys and no
           whitespace, with manifest.checksums.manifest_checksum set to null
           before hashing.

  split_checksum:
    input: JSON array of three sorted lists of entry_id strings:
           [sorted(train), sorted(validation), sorted(test)]
           serialised with sorted keys and no whitespace.

output_format: lowercase hex string (64 characters for SHA-256).
```

### Input Schema

```yaml
compute_checksums:
  manifest: PandoraManifest
  leakage_safe_dataset: LeakageSafeDataset
  export_policy: ExportPolicy
```

### Output Schema

```yaml
compute_checksums_result:
  artifact_checksum: string | null
  # null when export_policy.emit_checksum_bundle: false.

  manifest_checksum: string | null
  split_checksum: string | null
  checksum_algorithm: string
  # "SHA-256" in V1.
```

---

## 7.4 `build_pandora_artifact()`

### Responsibility

Construct the final `PandoraArtifact` from the leakage-safe dataset.
This is the main orchestrator for Component 06.

### Internal Workflow

```text
1. Build provenance bundle:
   build_provenance_bundle(leakage_safe_dataset, provenance_policy)
   → ProvenanceBundle

2. Generate manifest:
   generate_manifest(leakage_safe_dataset, provenance_bundle, export_policy)
   → PandoraManifest

3. Compute checksums (if record_checksums: true):
   compute_checksums(manifest, leakage_safe_dataset, export_policy)
   → checksums

4. Assemble PandoraArtifact with:
     leakage_safe_dataset, provenance_bundle, manifest, checksums,
     applied_policy, provenance (generated_at, pandora_version)

5. Export provenance report (if emit_provenance_report: true):
   export_provenance_report(artifact, export_policy)
   → ReproducibilityReport
   Attach to artifact.reproducibility_report.

6. Return PandoraArtifact.
```

### Input Schema

```yaml
build_pandora_artifact:
  input: ReproducibilityInput
```

### Output Schema

```yaml
build_pandora_artifact_result:
  artifact: PandoraArtifact
```

---

## 7.5 `export_provenance_report()`

### Responsibility

Generate a human-readable `ReproducibilityReport` summarising the pipeline
lineage and flagging any reproducibility risks.

### Input Schema

```yaml
export_provenance_report:
  artifact: PandoraArtifact
  export_policy: ExportPolicy
```

### Output Schema

```yaml
export_provenance_report_result:
  report: ReproducibilityReport
```

---

## 7.6 `build_pandora_artifact_many()`

### Responsibility

Build reproducibility artifacts for many leakage-safe datasets.

### Input Schema

```yaml
build_pandora_artifact_many:
  input: ReproducibilityBatchInput
```

### Output Schema

```yaml
build_pandora_artifact_many_result:
  result: ReproducibilityBatchResult
```

---

# 8. Internal Helper Functions

## 8.1 `aggregate_pipeline_provenance()`

### Responsibility

Traverse the `LeakageSafeDataset` hierarchy and collect per-stage
provenance records into the six `Pipeline*Provenance` sub-schemas.

### Input

```yaml
aggregate_pipeline_provenance:
  leakage_safe_dataset: LeakageSafeDataset
  provenance_policy: ProvenancePolicy
```

### Output

```yaml
aggregate_pipeline_provenance_result:
  ingestion: PipelineIngestionProvenance | null
  canonicalization: PipelineCanonicalizationProvenance | null
  metadata_integration: PipelineMetadataProvenance | null
  annotation_plugins: PipelineAnnotationProvenance | null
  dataset_curation: PipelineCurationProvenance | null
  leakage_splitting: PipelineSplittingProvenance | null
```

---

## 8.2 `aggregate_source_releases()`

### Responsibility

Collect archive and reference database release versions from the metadata
provenance records in the dataset.

### Input

```yaml
aggregate_source_releases:
  leakage_safe_dataset: LeakageSafeDataset
```

### Output

```yaml
aggregate_source_releases_result:
  source_release_provenance: SourceReleaseProvenanceRecord
```

---

## 8.3 `aggregate_policy_versions()`

### Responsibility

Collect all policy identifiers and versions applied across the pipeline.

### Input

```yaml
aggregate_policy_versions:
  leakage_safe_dataset: LeakageSafeDataset
```

### Output

```yaml
aggregate_policy_versions_result:
  policy_provenance: PolicyProvenanceRecord
```

---

## 8.4 `aggregate_plugin_versions()`

### Responsibility

Collect annotation plugin identities, versions, and configurations from
the annotated structures in the dataset.

### Input

```yaml
aggregate_plugin_versions:
  leakage_safe_dataset: LeakageSafeDataset
```

### Output

```yaml
aggregate_plugin_versions_result:
  annotation_provenance: PipelineAnnotationProvenance
```

---

## 8.5 `build_lineage_graph()`

### Responsibility

Construct the ordered `lineage` list for `ReproducibilityReport` from the
assembled `ProvenanceBundle`.

### Input

```yaml
build_lineage_graph:
  provenance_bundle: ProvenanceBundle
  leakage_safe_dataset: LeakageSafeDataset
```

### Output

```yaml
build_lineage_graph_result:
  lineage: list[string]
  # Ordered list of pipeline step summaries.
  # Format: "{step_number} {step_name}: {brief description}"
```

---

## 8.6 `summarize_reproducibility_risks()`

### Responsibility

Inspect the `ProvenanceBundle` for missing or null fields that weaken
reproducibility guarantees and emit a `list[Diagnostic]`.

### Input

```yaml
summarize_reproducibility_risks:
  provenance_bundle: ProvenanceBundle
  provenance_policy: ProvenancePolicy
```

### Output

```yaml
summarize_reproducibility_risks_result:
  reproducibility_risks: list[Diagnostic]
```

### V1 Risk Rules

```yaml
warning_rules:
  MISSING_ENGINE_VERSION:
    condition: "similarity_engine_versions contains a null entry."
    message: "Similarity engine version not recorded. Re-running may use a
              different version and produce different similarity scores."

  MISSING_SOURCE_RELEASE:
    condition: "pdbe_release, sifts_release, or uniprot_release is null."
    message: "Archive release version not recorded. Dataset identity cannot
              be fully verified against a specific archive snapshot."

  NULL_CANONICALIZATION_TIMESTAMP:
    condition: "canonicalization.canonicalized_at_range.earliest is null."
    message: "Canonicalization timestamps not recorded."

  MISSING_PANDORA_VERSION:
    condition: "PandoraArtifact.provenance.pandora_version is null."
    message: "Pandora version not recorded. This should never happen in V1."

  PARTIAL_METADATA_SOURCES:
    condition: "Not all requested metadata sources have status 'success'."
    message: "One or more metadata sources had partial or failed retrieval."
```

---

# 9. Provenance Fields Reference

This section maps each upstream component's provenance output to the
`ProvenanceBundle` field that records it.

## 9.1 Ingestion provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Provider | `MmCIFIngestionResult.provenance.provider` | `pipeline_provenance.ingestion.providers` |
| Source URI | `MmCIFIngestionResult.provenance.source_uri` | `pipeline_provenance.ingestion.source_uris` |
| Retrieval time | `MmCIFIngestionResult.provenance.retrieved_at` | `pipeline_provenance.ingestion.retrieved_at_range` |
| From cache | `MmCIFIngestionResult.provenance.from_cache` | `pipeline_provenance.ingestion.from_cache` |

---

## 9.2 Canonicalization provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Policy ID | `CanonicalStructureResult.applied_policy.policy_id` | `pipeline_provenance.canonicalization.policy_id` |
| Policy version | `CanonicalStructureResult.applied_policy.policy_version` | `pipeline_provenance.canonicalization.policy_version` |
| Timestamp | `CanonicalStructureResult.provenance.canonicalized_at` | `pipeline_provenance.canonicalization.canonicalized_at_range` |

---

## 9.3 Metadata provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Policy ID | `MetadataAnnotatedStructure.applied_metadata_policy.policy_id` | `pipeline_provenance.metadata_integration.metadata_policy_id` |
| Sources used | `AnnotatedStructureWithPlugins.provenance.metadata_sources` | `pipeline_provenance.metadata_integration.sources_used` |
| Retrieval time | `AnnotatedStructureWithPlugins.provenance.retrieved_at` | `pipeline_provenance.metadata_integration.retrieved_at_range` |
| Source versions | `ProvenanceMetadata.sources[*].source_version` | `source_release_provenance.*_release` |

---

## 9.4 Annotation provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Plugin IDs | `AnnotatedStructureWithPlugins.applied_plugins[*].plugin_id` | `annotation_provenance.plugin_records[*].plugin_id` |
| Plugin versions | `AnnotatedStructureWithPlugins.applied_plugins[*].plugin_version` | `annotation_provenance.plugin_records[*].plugin_version` |
| Plugin configs | `AnnotatedStructureWithPlugins.applied_plugins[*].config` | `annotation_provenance.plugin_records[*].config` |

---

## 9.5 Dataset curation provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Policy ID | `Dataset.applied_policy.policy_id` | `curation_provenance.curation_policy_id` |
| Filters applied | `Dataset.selection_summary.applied_filters` | `curation_provenance.filters_applied` |
| Exclusions | `Dataset.excluded_items` | `curation_provenance.exclusions` |
| Dedup strategy | `Dataset.deduplication_report.strategy` | `curation_provenance.deduplication_strategy` |

---

## 9.6 Leakage splitting provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Policy ID | `LeakageSafeDataset.applied_policy.policy_id` | `splitting_provenance.leakage_policy_id` |
| Engines | `LeakageSafeDataset.provenance.similarity_engines` | `splitting_provenance.similarity_engines` |
| Split time | `LeakageSafeDataset.provenance.split_at` | `splitting_provenance.split_at` |
| Achieved fractions | `LeakageSafeDataset.partition_summary.*_fraction_achieved` | `splitting_provenance.split_ratios_achieved` |
| Leakage detected | `LeakageSafeDataset.leakage_summary.leakage_detected` | `splitting_provenance.leakage_detected` |

---

# 10. Policy Schemas

## 10.1 ProvenancePolicy

```yaml
ProvenancePolicy:
  policy_id: string
  policy_name: string
  policy_version: string
  description: string

  record_software_versions: bool
  # If false, pandora_version is still recorded (non-nullable).
  # This flag controls whether dependency versions (e.g. gemmi, biopython)
  # are included.

  record_policy_versions: bool
  # If false, PolicyProvenanceRecord fields are null.

  record_source_releases: bool
  # If false, SourceReleaseProvenanceRecord fields are null.

  record_annotation_plugin_versions: bool
  # If false, PipelineAnnotationProvenance is null.

  record_curation_history: bool
  # If false, PipelineCurationProvenance is null.

  record_split_history: bool
  # If false, PipelineSplittingProvenance is null.

  record_checksums: bool
  # If false, all checksum fields are null.

  record_execution_timestamps: bool
  # If false, all *_at timestamp range fields are null.
  # generated_at in PandoraArtifact.provenance is always recorded.
```

---

## 10.2 ExportPolicy

```yaml
ExportPolicy:
  emit_manifest_yaml: bool
  # If true, serialise PandoraManifest to a .yaml file.

  emit_manifest_json: bool
  # If true, serialise PandoraManifest to a .json file.
  # emit_manifest_yaml and emit_manifest_json can both be true;
  # both files will be generated.

  emit_provenance_report: bool
  # If true, generate ReproducibilityReport and attach to the artifact.

  emit_checksum_bundle: bool
  # If true, compute and populate all checksum fields.
  # If false, all checksum fields are null.

  emit_lineage_graph: bool
  # If true, populate ReproducibilityReport.lineage with ordered
  # pipeline step summaries. Requires emit_provenance_report: true.
```

---

# 11. Non-Responsibilities

Component 06 is not responsible for:
  - ingestion
  - canonicalization
  - metadata_retrieval
  - annotation_computation
  - dataset_curation
  - leakage_safe_splitting
  - benchmark_creation
  - model_training

---

# 12. Component Definition

The Provenance & Reproducibility Layer aggregates pipeline provenance into reproducibility manifests, integrity checksums, lineage records, and exportable Pandora artifacts so that generated datasets can be audited and reproduced exactly.
