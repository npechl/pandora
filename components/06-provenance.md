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
* and exporting reproducibility bundles in either embedded or by-reference
  artifact mode.

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
  → PandoraDataset (Dataset | ChainDataset | InterfaceDataset | ResidueDataset)
  → LeakageSafeDataset
  → PandoraArtifact (embedded or by-reference)
```

---

# 2. Core Design Principles

## Provenance is cross-cutting

Every preceding component contributes provenance data.

Component 6 centralises that information into one reproducible artifact.
When the source dataset is in-memory, it traverses the
`LeakageSafeDataset` object hierarchy directly. When the source dataset is
materialized, it reads provenance metadata fields from the dataset without
loading atom-coordinate records.

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

## Two artifact storage modes

**Embedded mode** — default for small datasets (~≤10K structures).

`PandoraArtifact.leakage_safe_dataset` holds the full `LeakageSafeDataset`
as a nested Python object. The artifact is self-contained and shareable as
a single serialized file (e.g. JSON or pickle), but may become impractically
large for datasets with many structures.

**By-reference mode** — required for large datasets (>10K structures), or
when the source dataset was in materialized mode.

`PandoraArtifact.leakage_safe_dataset` is null. Instead, the artifact holds
an `ArtifactStoreRef` pointing to a directory on disk that contains:

```text
{store_root}/
  manifest.json       — PandoraManifest (always present)
  manifest.yaml       — Optional; emitted when emit_manifest_yaml: true
  provenance.json     — ProvenanceBundle serialised to JSON
  checksums.json      — ArtifactChecksums
  splits/
    train.parquet     — Item identifiers for the train split
    validation.parquet
    test.parquet
```

The `splits/` Parquet files contain only the partition item identifier
strings — not full atom-coordinate records. Structure data remains in the
C04 `DatasetStore` (referenced via `LeakageSafeDataset.source_dataset_ref`).

The artifact mode is set automatically based on whether
`LeakageSafeDataset.source_dataset` is non-null (embedded) or
`LeakageSafeDataset.source_dataset_ref` is non-null (by-reference). It can
also be overridden explicitly via `ExportPolicy.artifact_mode`.

---

## Provenance depth depends on granularity

For structure-level datasets, Component 06 traverses the full object
hierarchy to collect per-entry ingestion, canonicalization, metadata, and
annotation provenance.

For chain/interface/residue-level datasets, the upstream
`AnnotatedStructureWithPlugins` objects are not embedded in the dataset.
Only curation and splitting provenance are available. This is recorded as
an `UPSTREAM_PROVENANCE_NOT_EMBEDDED` warning in the reproducibility report.

---

# 3. ArtifactStoreRef Schema

```yaml
ArtifactStoreRef:
  artifact_id: string
  # Matches PandoraArtifact.artifact_id.

  store_root: string
  # Absolute path to the artifact store root directory.
  # Contains the directory layout described in Section 2 above.

  store_format: string
  # parquet
  # V1: split files are always Parquet.
```

---

# 4. Input Schemas

## 4.1 Reproducibility input

```yaml
ReproducibilityInput:
  leakage_safe_dataset: LeakageSafeDataset

  provenance_policy: ProvenancePolicy
  # See Section 10.1 for the full policy schema.

  export_policy: ExportPolicy
  # See Section 10.2 for the full policy schema.
```

---

## 4.2 Batch reproducibility input

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
    fail_fast: bool
```

---

# 5. Provenance Sub-Schemas

These schemas define the typed content of each provenance category within
`ProvenanceBundle`. Each is populated by traversing the corresponding
fields in the `LeakageSafeDataset` hierarchy.

## 5.1 PipelineIngestionProvenance

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

## 5.2 PipelineCanonicalizationProvenance

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

## 5.3 PipelineMetadataProvenance

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

## 5.4 PipelineAnnotationProvenance

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

## 5.5 PipelineCurationProvenance

Sourced from `LeakageSafeDataset.source_dataset` (in-memory) or
`LeakageSafeDataset.source_dataset_ref` metadata (materialized).

```yaml
PipelineCurationProvenance:
  curation_policy_id: string
  curation_policy_name: string
  curation_policy_version: string
  granularity: string
  total_input: int
  total_selected: int
  total_excluded: int
  total_duplicates_removed: int
  filters_applied: list[AppliedFilterRecord]
  exclusions: list[ExclusionRecord]
  deduplication_strategy: string | null
```

---

## 5.6 PipelineSplittingProvenance

Sourced from `LeakageSafeDataset`.

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

## 5.7 SourceReleaseProvenanceRecord

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

## 5.8 PolicyProvenanceRecord

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

# 6. Manifest Schema

The `PandoraManifest` is the primary machine-readable artifact shared
between researchers for reproducibility. It is a compact summary — it does
not embed full structure data in either artifact mode.

```yaml
PandoraManifest:
  manifest_id: string
  manifest_format: list[string]
  # ["json"] | ["yaml"] | ["json", "yaml"]
  # Records which formats were emitted.

  pandora_version: string
  generated_at: string
  # ISO 8601 timestamp.

  artifact_id: string
  artifact_mode: string
  # embedded | by_reference

  dataset_summary:
    dataset_id: string
    dataset_version: string
    granularity: string
    # Values: structure | chain | interface | residue
    total_items: int
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

  by_reference_layout: object | null
  # Populated in by_reference mode only.
  # Describes the file layout of the ArtifactStore:
  #   store_root: string
  #   files:
  #     manifest_json: string | null   (relative path)
  #     manifest_yaml: string | null   (relative path)
  #     provenance_json: string        (relative path)
  #     checksums_json: string         (relative path)
  #     train_parquet: string          (relative path)
  #     validation_parquet: string     (relative path)
  #     test_parquet: string           (relative path)

  checksums:
    artifact_checksum: string | null
    manifest_checksum: string | null
    split_checksum: string | null
    checksum_algorithm: string
    # Always "SHA-256" in V1.
```

---

# 7. Output Schemas

## 7.1 Pandora artifact

```yaml
PandoraArtifact:
  artifact_id: string
  artifact_name: string | null

  artifact_mode: string
  # embedded     — leakage_safe_dataset is populated; artifact_store_ref is null.
  # by_reference — leakage_safe_dataset is null; artifact_store_ref is populated.

  leakage_safe_dataset: LeakageSafeDataset | null
  # Non-null in embedded mode.
  # null in by_reference mode — data lives in artifact_store_ref.

  artifact_store_ref: ArtifactStoreRef | null
  # null in embedded mode.
  # Non-null in by_reference mode. Points to the store directory containing
  # manifest.json, provenance.json, checksums.json, and splits/ Parquet files.

  provenance_bundle: ProvenanceBundle

  manifest: PandoraManifest
  # Always populated in memory. Serialised to disk according to export_policy.

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
    # Always populated (non-nullable).
```

---

## 7.2 Provenance bundle

```yaml
ProvenanceBundle:
  pipeline_provenance:
    ingestion: PipelineIngestionProvenance | null
    canonicalization: PipelineCanonicalizationProvenance | null
    metadata_integration: PipelineMetadataProvenance | null
    annotation_plugins: PipelineAnnotationProvenance | null
    dataset_curation: PipelineCurationProvenance | null
    leakage_splitting: PipelineSplittingProvenance | null
    # Each field is null when the corresponding record_* flag is false
    # in ProvenancePolicy, or when data is unavailable.

  source_release_provenance: SourceReleaseProvenanceRecord
  policy_provenance: PolicyProvenanceRecord
  annotation_provenance: PipelineAnnotationProvenance | null
  curation_provenance: PipelineCurationProvenance | null
  splitting_provenance: PipelineSplittingProvenance | null
```

---

## 7.3 Reproducibility report

```yaml
ReproducibilityReport:
  report_id: string
  artifact_id: string

  summary:
    pipeline_steps: int
    source_count: int
    # Number of unique source providers used.
    policy_count: int
    # Number of distinct policies applied.
    plugin_count: int
    # Number of annotation plugins applied.

  lineage: list[string]
  # Ordered list of human-readable pipeline step summaries.
  # Format: "{step_number} {step_name}: {brief description}"
  # e.g. structure-level:
  #   ["01 ingestion: 150 structures from pdbe",
  #    "02 canonicalization: 150 structures (policy: default_v1)",
  #    "03 metadata: pdbe, sifts, uniprot retrieved",
  #    "04 curation (structure): 120 selected, 30 excluded",
  #    "05 splitting: train=84, val=18, test=18 (leakage: clean)"]
  #
  # chain-level:
  #   ["04 curation (chain): 847 chains selected, 23 excluded",
  #    "05 splitting: train=593, val=127, test=127 (leakage: clean)"]
  # Note: steps 01-03 will be null when source is a ChainDataset.

  reproducibility_risks: list[Diagnostic]
  # Warnings about missing provenance fields or reproducibility weaknesses.
```

---

## 7.4 Batch output

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

# 8. Public Functions

## 8.1 `build_provenance_bundle()`

### Responsibility

Aggregate provenance from the leakage-safe dataset and all upstream
pipeline stages.

In embedded mode, traverses the `LeakageSafeDataset` object hierarchy
directly.

In by-reference mode (when `source_dataset` is null and
`source_dataset_ref` is non-null), reads provenance metadata from the
dataset's schema fields without loading atom-coordinate records.

### Internal Workflow

```text
1. Extract pipeline provenance:
   aggregate_pipeline_provenance(leakage_safe_dataset, provenance_policy)

   When source_dataset.granularity == "structure" (in-memory only):
     ingestion:        source_dataset.structures[*]
                       .canonical_structure_result.provenance
     canonicalization: source_dataset.structures[*]
                       .canonical_structure_result.applied_policy
     metadata:         source_dataset.structures[*]
                       .metadata_annotations.provenance_metadata.sources
     annotation:       source_dataset.structures[*]
                       .applied_plugins + .derived_annotations[*]

   When source_dataset is null (materialized mode) OR
   when source_dataset.granularity in ("chain", "interface", "residue"):
     ingestion, canonicalization, metadata, annotation:
       null — recorded as null with UPSTREAM_PROVENANCE_NOT_EMBEDDED warning.

   For all cases:
     curation:   source_dataset (or source_dataset_ref metadata)
                 (.applied_policy, .selection_summary, .deduplication_report)
     splitting:  leakage_safe_dataset
                 (.applied_policy, .partition_summary, .provenance)

2. Aggregate source releases (if record_source_releases: true):
   aggregate_source_releases(leakage_safe_dataset)

3. Aggregate policy versions (if record_policy_versions: true):
   aggregate_policy_versions(leakage_safe_dataset)

4. Aggregate plugin versions (if record_annotation_plugin_versions: true):
   aggregate_plugin_versions(leakage_safe_dataset)

5. Assemble and return ProvenanceBundle.
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

## 8.2 `generate_manifest()`

### Responsibility

Generate a `PandoraManifest` describing the artifact, its dataset summary,
source releases, policy versions, and (optionally) checksums.

In by-reference mode, also populates `manifest.by_reference_layout` with
the expected file paths within the `ArtifactStore`.

When `export_policy.emit_manifest_yaml: true`, the manifest is serialised
to `{store_root}/manifest.yaml`.

When `export_policy.emit_manifest_json: true`, the manifest is serialised
to `{store_root}/manifest.json`.

Both flags can be `true` simultaneously.

### Input Schema

```yaml
generate_manifest:
  leakage_safe_dataset: LeakageSafeDataset
  provenance_bundle: ProvenanceBundle
  export_policy: ExportPolicy
  artifact_id: string
  artifact_mode: string
  # embedded | by_reference
  store_root: string | null
  # Required when artifact_mode == "by_reference". Path to the artifact
  # store root directory where files will be written.
```

### Output Schema

```yaml
generate_manifest_result:
  manifest: PandoraManifest
```

---

## 8.3 `compute_checksums()`

### Responsibility

Compute SHA-256 integrity checksums for the artifact, manifest, and split
partition lists.

### Checksum contract

```yaml
checksum_algorithm: SHA-256

serialization:
  artifact_checksum:
    input: >
      PandoraArtifact serialised to JSON with sorted keys and no
      whitespace, with checksums.artifact_checksum set to null before
      hashing (to avoid circular dependency).

  manifest_checksum:
    input: >
      PandoraManifest serialised to JSON with sorted keys and no
      whitespace, with manifest.checksums.manifest_checksum set to null
      before hashing.

  split_checksum:
    input: >
      JSON array of three sorted lists of item identifier strings:
      [sorted(train), sorted(validation), sorted(test)]
      serialised with sorted keys and no whitespace.
      In by_reference mode, the lists are read from the splits/ Parquet
      files before computing the checksum.

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

## 8.4 `build_pandora_artifact()`

### Responsibility

Construct the final `PandoraArtifact` from the leakage-safe dataset.
This is the main orchestrator for Component 06.

### Artifact mode selection

```text
If export_policy.artifact_mode is explicitly set:
    Use the specified mode.
Else:
    If leakage_safe_dataset.source_dataset is non-null:
        Use "embedded".
    Else (source_dataset is null, source_dataset_ref is non-null):
        Use "by_reference".
```

### Internal Workflow

```text
1. Determine artifact_mode from export_policy or dataset mode (see above).

2. Build provenance bundle:
   build_provenance_bundle(leakage_safe_dataset, provenance_policy)
   → ProvenanceBundle

3. If artifact_mode == "by_reference":
   a. Initialise artifact store directory at export_policy.store_root.
   b. Write split Parquet files:
      splits/train.parquet, splits/validation.parquet, splits/test.parquet
      (item identifier strings only — no atom coordinates).
   c. Write provenance.json.

4. Generate manifest:
   generate_manifest(leakage_safe_dataset, provenance_bundle,
                     export_policy, artifact_id, artifact_mode, store_root)
   → PandoraManifest

5. Compute checksums (if record_checksums: true):
   compute_checksums(manifest, leakage_safe_dataset, export_policy)
   → checksums

6. Write checksums.json to store_root (by_reference mode only).

7. Assemble PandoraArtifact:
   - embedded mode:
       leakage_safe_dataset = leakage_safe_dataset
       artifact_store_ref = null
   - by_reference mode:
       leakage_safe_dataset = null
       artifact_store_ref = ArtifactStoreRef(artifact_id, store_root)

8. Export provenance report (if emit_provenance_report: true):
   export_provenance_report(artifact, export_policy)
   → ReproducibilityReport
   Attach to artifact.reproducibility_report.

9. Return PandoraArtifact.
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

## 8.5 `export_provenance_report()`

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

## 8.6 `build_pandora_artifact_many()`

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

# 9. Internal Helper Functions

## 9.1 `aggregate_pipeline_provenance()`

Traverse the `LeakageSafeDataset` hierarchy and collect per-stage
provenance records.

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

## 9.2 `aggregate_source_releases()`

Collect archive and reference database release versions from metadata
provenance records.

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

## 9.3 `aggregate_policy_versions()`

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

## 9.4 `aggregate_plugin_versions()`

Collect annotation plugin identities and versions from the dataset.

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

## 9.5 `build_lineage_graph()`

Construct the ordered `lineage` list for `ReproducibilityReport`.

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
```

---

## 9.6 `summarize_reproducibility_risks()`

Inspect the `ProvenanceBundle` for missing or null fields that weaken
reproducibility guarantees.

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

  UPSTREAM_PROVENANCE_NOT_EMBEDDED:
    condition: "source_dataset.granularity != 'structure' OR
                leakage_safe_dataset.source_dataset is null."
    message: "Chain/interface/residue datasets and materialized datasets do
              not embed upstream AnnotatedStructureWithPlugins objects.
              Ingestion, canonicalization, metadata, and annotation provenance
              fields will be null. Only curation and splitting provenance are
              available."
```

---

# 10. Provenance Fields Reference

## 10.1 Ingestion provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Provider | `MmCIFIngestionResult.provenance.provider` | `pipeline_provenance.ingestion.providers` |
| Source URI | `MmCIFIngestionResult.provenance.source_uri` | `pipeline_provenance.ingestion.source_uris` |
| Retrieval time | `MmCIFIngestionResult.provenance.retrieved_at` | `pipeline_provenance.ingestion.retrieved_at_range` |
| From cache | `MmCIFIngestionResult.provenance.from_cache` | `pipeline_provenance.ingestion.from_cache` |

---

## 10.2 Canonicalization provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Policy ID | `CanonicalStructureResult.applied_policy.policy_id` | `pipeline_provenance.canonicalization.policy_id` |
| Policy version | `CanonicalStructureResult.applied_policy.policy_version` | `pipeline_provenance.canonicalization.policy_version` |
| Timestamp | `CanonicalStructureResult.provenance.canonicalized_at` | `pipeline_provenance.canonicalization.canonicalized_at_range` |

---

## 10.3 Metadata provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Policy ID | `MetadataAnnotatedStructure.applied_metadata_policy.policy_id` | `pipeline_provenance.metadata_integration.metadata_policy_id` |
| Sources used | `AnnotatedStructureWithPlugins.provenance.metadata_sources` | `pipeline_provenance.metadata_integration.sources_used` |
| Retrieval time | `AnnotatedStructureWithPlugins.provenance.retrieved_at` | `pipeline_provenance.metadata_integration.retrieved_at_range` |
| Source versions | `ProvenanceMetadata.sources[*].source_version` | `source_release_provenance.*_release` |

---

## 10.4 Annotation provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Plugin IDs | `AnnotatedStructureWithPlugins.applied_plugins[*].plugin_id` | `annotation_provenance.plugin_records[*].plugin_id` |
| Plugin versions | `AnnotatedStructureWithPlugins.applied_plugins[*].plugin_version` | `annotation_provenance.plugin_records[*].plugin_version` |
| Plugin configs | `AnnotatedStructureWithPlugins.applied_plugins[*].config` | `annotation_provenance.plugin_records[*].config` |

---

## 10.5 Dataset curation provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Policy ID | `Dataset.applied_policy.policy_id` | `curation_provenance.curation_policy_id` |
| Filters applied | `Dataset.selection_summary.applied_filters` | `curation_provenance.filters_applied` |
| Exclusions | `Dataset.excluded_items` | `curation_provenance.exclusions` |
| Dedup strategy | `Dataset.deduplication_report.strategy` | `curation_provenance.deduplication_strategy` |

---

## 10.6 Leakage splitting provenance

| Field | Source | ProvenanceBundle path |
|-------|--------|-----------------------|
| Policy ID | `LeakageSafeDataset.applied_policy.policy_id` | `splitting_provenance.leakage_policy_id` |
| Engines | `LeakageSafeDataset.provenance.similarity_engines` | `splitting_provenance.similarity_engines` |
| Split time | `LeakageSafeDataset.provenance.split_at` | `splitting_provenance.split_at` |
| Achieved fractions | `LeakageSafeDataset.partition_summary.*_fraction_achieved` | `splitting_provenance.split_ratios_achieved` |
| Leakage detected | `LeakageSafeDataset.leakage_summary.leakage_detected` | `splitting_provenance.leakage_detected` |

---

# 11. Policy Schemas

## 11.1 ProvenancePolicy

```yaml
ProvenancePolicy:
  policy_id: string
  policy_name: string
  policy_version: string
  description: string

  record_software_versions: bool
  # If false, pandora_version is still recorded (non-nullable).
  # This flag controls whether dependency versions (gemmi, etc.) are included.

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

## 11.2 ExportPolicy

```yaml
ExportPolicy:
  artifact_mode: string | null
  # embedded     — force embedded mode regardless of dataset mode.
  # by_reference — force by_reference mode regardless of dataset mode.
  # null         — auto-select based on dataset.source_dataset (see Section 8.4).

  store_root: string | null
  # Required when artifact_mode == "by_reference" or when auto-selection
  # results in by_reference mode.
  # Absolute path to the artifact store root directory.

  emit_manifest_yaml: bool
  # If true, serialise PandoraManifest to {store_root}/manifest.yaml.

  emit_manifest_json: bool
  # If true, serialise PandoraManifest to {store_root}/manifest.json.
  # Both emit_manifest_yaml and emit_manifest_json can be true.

  emit_provenance_report: bool
  # If true, generate ReproducibilityReport and attach to the artifact.

  emit_checksum_bundle: bool
  # If true, compute and populate all checksum fields.
  # If false, all checksum fields are null.

  emit_lineage_graph: bool
  # If true, populate ReproducibilityReport.lineage.
  # Requires emit_provenance_report: true.
```

---

# 12. Non-Responsibilities

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

# 13. Component Definition

The Provenance & Reproducibility Layer aggregates pipeline provenance into
reproducibility manifests, integrity checksums, lineage records, and
exportable Pandora artifacts. It supports two artifact storage modes:
embedded (for small datasets, where the full LeakageSafeDataset is held in
the artifact object) and by-reference (for large datasets, where the artifact
is a lightweight manifest and checksums pointing to split Parquet files in an
ArtifactStore directory on disk).
