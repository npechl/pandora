# Component 03 — Metadata Integration & Derived Annotation Layer

## Purpose

The Metadata Integration & Derived Annotation Layer enriches canonical
structures with:

* stable metadata,
* biological mappings,
* structural annotations,
* provenance-aware metadata,
* and optional plugin-style derived annotations.

This component does NOT modify canonical structures.

Instead, it attaches:

* metadata layers,
* and optional task-specific annotation layers.

This layer is:

* modular,
* policy-driven,
* provenance-aware,
* and extensible through user-defined plugins.

---

# 1. Architectural Role

```text
MmCIFIngestionResult
  → CanonicalStructureResult
  → MetadataAnnotatedStructure
  → AnnotatedStructureWithPlugins
```

---

# 2. Core Design Principles

## Metadata is attached, not embedded

The canonical structure remains unchanged after Component 02.

Metadata and derived annotations are attached as separate layers and never
written back into `CanonicalStructure`.

---

## Metadata retrieval is policy-driven

Users explicitly define:

* which metadata sources to retrieve,
* which metadata categories to include,
* retrieval behaviour,
* and provenance requirements.

---

## Derived annotations are plugin-based

Users may:

* use built-in annotation plugins,
* or define their own task-specific annotation plugins.

Examples:

* protein-protein interface annotations,
* contact maps,
* ligand-binding site annotations,
* pocket annotations,
* surface annotations,
* custom research annotations.

---

# 3. Input Schemas

## 3.1 Metadata integration input

```yaml
MetadataIntegrationInput:
  canonical_structure_result: CanonicalStructureResult
  # Must have status "success" or "warning".
  # "failed" results are rejected before retrieval begins.

  metadata_policy: MetadataIntegrationPolicy
  # See Section 9.1 for the full policy schema.
```

---

## 3.2 Annotation plugin input

```yaml
AnnotationPluginInput:
  structure:
    CanonicalStructureResult | MetadataAnnotatedStructure
  # If a CanonicalStructureResult is provided directly (without prior
  # metadata attachment), metadata_annotations in the resulting
  # AnnotatedStructureWithPlugins will contain null category objects.
  # Prefer MetadataAnnotatedStructure when metadata is available.

  plugin_policy: AnnotationPluginPolicy
  # See Section 9.2 for the full policy schema.

  plugins: list[Plugin]
  # Plugins to apply. See Section 4.3 for the Plugin schema.
```

---

## 3.3 Batch integration input

```yaml
MetadataIntegrationBatchInput:
  structures:
    - CanonicalStructureResult

  metadata_policy: MetadataIntegrationPolicy

  plugin_policy: AnnotationPluginPolicy | null

  mode: string
  # sequential | parallel

  parallel_options:
    max_workers: int | null
    # Number of concurrent workers in parallel mode.
    # null uses the system default (typically CPU count).
    # Ignored in sequential mode.

    fail_fast: bool
    # If true, abort remaining entries on the first failure.
    # If false (default), isolate failures and continue.
```

---

# 4. Core Object Schemas

## 4.1 MetadataRecord

The intermediate representation produced by `retrieve_metadata()` and
consumed by `attach_metadata()`.

```yaml
MetadataRecord:
  entry_id: string

  archive_metadata: ArchiveMetadata | null
  # null when pdbe source is disabled or retrieval failed.

  biological_mappings: BiologicalMappings | null
  # null when sifts and uniprot sources are both disabled or failed.

  structural_annotations: StructuralAnnotations | null
  # null when structural annotation retrieval is disabled or failed.

  provenance_metadata: ProvenanceMetadata
  # Always populated; records what was retrieved and from where.

  retrieval_status:
    archive_metadata: string | null
    # success | partial | failed | disabled
    biological_mappings: string | null
    structural_annotations: string | null

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

---

## 4.2 AnnotationLayer

The typed output of a single annotation plugin application.

```yaml
AnnotationLayer:
  layer_name: string
  # Human-readable name for this annotation layer.
  # e.g. "protein_protein_interfaces", "residue_contact_map"

  layer_type: string
  # Machine-readable type identifier.
  # Built-in values:
  #   interface_annotations | contact_map | ligand_binding_sites |
  #   pocket_annotations | secondary_structure | surface_exposure |
  #   domain_annotations | custom
  # User-defined plugins should use a namespaced type:
  #   e.g. "org.mylab.my_annotation_type"

  schema_version: string | null
  # Version of the data schema for this layer_type.
  # Allows consumers to handle schema evolution across plugin versions.

  data: object
  # Annotation data. The structure of data is defined per layer_type.
  # Built-in layer_type schemas are documented in Section 5.
  # User-defined plugins are responsible for documenting their data schema.
```

---

## 4.3 Plugin

The formal definition of a registered annotation plugin.

```yaml
Plugin:
  plugin_id: string
  # Unique identifier for this plugin.
  # Built-in plugins use the prefix "pandora.builtin.".
  # User-defined plugins should use a namespaced identifier.

  plugin_name: string
  plugin_version: string

  plugin_type: string
  # pandora.builtin.interfaces       — protein-protein interface annotations
  # pandora.builtin.contacts         — residue contact maps
  # pandora.builtin.ligand_binding   — ligand-binding site annotations
  # pandora.builtin.pockets          — pocket annotations
  # pandora.builtin.secondary_structure — secondary structure assignments
  # pandora.builtin.surface          — surface exposure annotations
  # pandora.builtin.domains          — domain annotations
  # user_defined                     — user-supplied plugin

  input_type: string
  # canonical_structure          — accepts CanonicalStructureResult
  # metadata_annotated_structure — accepts MetadataAnnotatedStructure
  #                                (required when plugin reads metadata)

  output_layer_type: string
  # The layer_type value this plugin will produce in its AnnotationLayer.

  description: string | null

  config: object | null
  # Plugin-specific configuration. Schema is defined per plugin_type.
```

---

# 5. Metadata Category Schemas

These schemas define the content of the four metadata categories populated
by `retrieve_metadata()` in V1.

## 5.1 ArchiveMetadata

Populated from PDBe.

```yaml
ArchiveMetadata:
  entry_id: string

  title: string | null
  # _struct.title

  experimental_method: string | null
  # e.g. "X-RAY DIFFRACTION", "SOLUTION NMR", "ELECTRON MICROSCOPY"
  # _exptl.method

  resolution: float | null
  # In Angstroms. null for NMR and other non-diffraction methods.
  # _refine.ls_d_res_high

  r_factor: float | null
  r_free: float | null
  # _refine.ls_R_factor_obs / ls_R_factor_R_free

  deposition_date: string | null
  # ISO 8601 date. _pdbx_database_status.recvd_initial_deposition_date

  release_date: string | null
  # ISO 8601 date. _pdbx_database_status.status_code_sf

  revision_date: string | null
  # ISO 8601 date of the most recent revision.

  keywords: list[string]
  # _struct_keywords.pdbx_keywords split by comma.

  space_group: string | null
  # e.g. "P 21 21 21". _symmetry.space_group_name_H-M

  unit_cell:
    a: float | null
    b: float | null
    c: float | null
    alpha: float | null
    beta: float | null
    gamma: float | null
  # _cell.length_a/b/c, _cell.angle_alpha/beta/gamma
  # All fields null for NMR and EM structures.
```

---

## 5.2 BiologicalMappings

Populated from UniProt and SIFTS.

```yaml
BiologicalMappings:
  uniprot_mappings: list[UniProtMapping]
  sifts_mappings: list[SIFTSMapping]
  taxonomy: TaxonomyRecord | null

UniProtMapping:
  canonical_chain_id: string
  # References CanonicalStructure.Chain.chain_id

  uniprot_accession: string
  # e.g. "P00533"

  uniprot_id: string | null
  # e.g. "EGFR_HUMAN"

  segment_start: int
  segment_end: int
  # Canonical residue number range covered by this mapping.

  identity: float | null
  # Sequence identity [0.0, 1.0] between chain segment and UniProt sequence.

SIFTSMapping:
  canonical_chain_id: string

  uniprot_accession: string | null
  pfam_id: string | null
  cath_id: string | null
  scop_id: string | null
  interpro_id: string | null

  residue_level_mappings: list[SIFTSResidueLevelMapping] | null
  # Per-residue cross-reference data. null if not requested.

SIFTSResidueLevelMapping:
  canonical_chain_id: string
  canonical_seq_id: int
  uniprot_position: int | null
  uniprot_accession: string | null
  pfam_id: string | null

TaxonomyRecord:
  ncbi_taxon_id: int | null
  organism_scientific: string | null
  # e.g. "Homo sapiens"
  organism_common: string | null
  # e.g. "Human"
  lineage: list[string]
  # Ordered from root to species.
  # e.g. ["cellular organisms", "Eukaryota", ..., "Homo sapiens"]
```

---

## 5.3 StructuralAnnotations

Populated from trusted archive-provided annotations (PDBe).

```yaml
StructuralAnnotations:
  preferred_assembly_id: string | null
  # The assembly_id (canonical) annotated as the preferred biological unit.

  entity_annotations: list[EntityAnnotation]

  chain_annotations: list[ChainAnnotation]

EntityAnnotation:
  canonical_entity_id: string
  entity_type: string
  # polymer | non-polymer | water | branched
  description: string | null

ChainAnnotation:
  canonical_chain_id: string
  chain_length: int | null
  # Number of residues (polymer only; null for non-polymer chains).
  is_polymer: bool
  uniprot_accession: string | null
```

---

## 5.4 ProvenanceMetadata

Always populated; records what was retrieved, from where, and when.

```yaml
ProvenanceMetadata:
  sources: list[MetadataSourceRecord]

MetadataSourceRecord:
  source_name: string
  # pdbe | sifts | uniprot | taxonomy

  retrieval_status: string
  # success | partial | failed | disabled

  source_version: string | null
  # Release or API version used. null when not available.

  retrieved_at: string | null
  # ISO 8601 timestamp. null when source was disabled.

  base_url: string | null
  # The API or data endpoint used for retrieval.
```

---

# 6. Output Schemas

## 6.1 Metadata-annotated structure

```yaml
MetadataAnnotatedStructure:
  canonical_structure_result: CanonicalStructureResult

  metadata_annotations:
    archive_metadata: ArchiveMetadata | null
    biological_mappings: BiologicalMappings | null
    structural_annotations: StructuralAnnotations | null
    provenance_metadata: ProvenanceMetadata

  applied_metadata_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]

  provenance:
    metadata_sources: list[string]
    # Provider names used: e.g. ["pdbe", "sifts", "uniprot", "taxonomy"]
    retrieved_at: string | null
    # ISO 8601 timestamp of the most recent retrieval.
```

---

## 6.2 Annotation plugin result

```yaml
AnnotationPluginResult:
  plugin_id: string
  plugin_name: string
  plugin_version: string

  status: string
  # success | warning | failed

  annotation_layer: AnnotationLayer | null
  # null when status == "failed".

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]

  provenance:
    source_structure_id: string
    # entry_id of the structure this annotation was computed from.
    executed_at: string | null
    # ISO 8601 timestamp.
```

---

## 6.3 Fully annotated structure

```yaml
AnnotatedStructureWithPlugins:
  canonical_structure_result: CanonicalStructureResult

  metadata_annotations:
    archive_metadata: ArchiveMetadata | null
    biological_mappings: BiologicalMappings | null
    structural_annotations: StructuralAnnotations | null
    provenance_metadata: ProvenanceMetadata | null
  # All category fields are null when the input was a CanonicalStructureResult
  # passed directly without prior metadata attachment.

  derived_annotations: list[AnnotationLayer]
  # One AnnotationLayer per successfully applied plugin.
  # Empty list when no plugins were applied or all plugins failed.

  applied_metadata_policy:
    policy_id: string
    policy_name: string
    policy_version: string | null
  # policy_version is null when no metadata policy was applied
  # (i.e. plugins were applied directly to a CanonicalStructureResult).

  applied_plugins: list[Plugin]
  # The Plugin records for every plugin that was applied.

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]

  provenance:
    metadata_sources: list[string]
    # Provider names used: e.g. ["pdbe", "sifts", "uniprot", "taxonomy"]
    # Empty list when no metadata was retrieved.
    retrieved_at: string | null
    annotation_history: list[string]
    # Ordered list of plugin_ids applied, e.g. ["pandora.builtin.contacts",
    # "org.mylab.custom_annotation"]. Preserves application order.
```

---

## 6.4 Batch integration result

```yaml
MetadataAndAnnotationBatchResult:
  mode: string
  # sequential | parallel

  summary:
    total: int
    success: int
    warning: int
    failed: int

  results:
    - entry_id: string

      status: string
      # success | warning | failed

      annotated_structure: AnnotatedStructureWithPlugins | null
      # null when status == "failed".

      diagnostics:
        warnings: list[Diagnostic]
        errors: list[Diagnostic]
```

---

# 7. Public Functions

## 7.1 `retrieve_metadata()`

### Responsibility

Retrieve metadata from selected metadata providers according to the policy.

### Internal Workflow

```text
for each source enabled in metadata_policy.include_sources:
    retrieve_{source}_metadata()     [pdbe | sifts | uniprot | taxonomy]
    record retrieval status per source

if any source failed:
    if fail_on_missing: true
        return error — do not proceed to attach
    if allow_partial: true
        continue with available metadata; emit PARTIAL_METADATA warning
    if allow_partial: false and fail_on_missing: false
        continue with null for that category; emit MISSING_METADATA warning

merge retrieved records into MetadataRecord
return MetadataRecord with per-source retrieval_status
```

### Input Schema

```yaml
retrieve_metadata:
  canonical_structure_result: CanonicalStructureResult
  metadata_policy: MetadataIntegrationPolicy
```

### Output Schema

```yaml
retrieve_metadata_result:
  metadata_record: MetadataRecord

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

### Notes

This function is source-aware and policy-driven. It records per-source
retrieval status in `MetadataRecord.retrieval_status`.

---

## 7.2 `validate_metadata()`

### Responsibility

Validate the completeness and consistency of retrieved metadata records.

### Input Schema

```yaml
validate_metadata:
  metadata_record: MetadataRecord
  metadata_policy: MetadataIntegrationPolicy
```

### Output Schema

```yaml
validate_metadata_result:
  validation_status: string
  # valid | warning | invalid

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

### V1 Validation Rules

```yaml
error_rules:
  REQUIRED_SOURCE_MISSING:
    condition: "A source with fail_on_missing: true returned status 'failed'."
    result_status: invalid

warning_rules:
  PARTIAL_METADATA:
    condition: "One or more sources failed but allow_partial: true; metadata
                is incomplete."
    result_status: warning

  MISSING_METADATA:
    condition: "A source is disabled or failed; that category is null."
    result_status: warning

  STALE_METADATA:
    condition: "Retrieved metadata is older than a configurable threshold."
    result_status: warning
```

---

## 7.3 `attach_metadata()`

### Responsibility

Attach a validated `MetadataRecord` to a `CanonicalStructureResult`,
producing a `MetadataAnnotatedStructure`.

### Internal Workflow

```text
metadata_record = retrieve_metadata(canonical_structure_result, policy)
validate_metadata(metadata_record, policy)

if validation fails and fail_on_missing: true:
    raise error

return MetadataAnnotatedStructure(
    canonical_structure_result=canonical_structure_result,
    metadata_annotations={
        archive_metadata: metadata_record.archive_metadata,
        biological_mappings: metadata_record.biological_mappings,
        structural_annotations: metadata_record.structural_annotations,
        provenance_metadata: metadata_record.provenance_metadata,
    },
    applied_metadata_policy=policy,
    ...
)
```

This function does NOT alter the canonical structure itself.

### Input Schema

```yaml
attach_metadata:
  canonical_structure_result: CanonicalStructureResult
  metadata_policy: MetadataIntegrationPolicy
```

### Output Schema

```yaml
attach_metadata_result:
  metadata_annotated_structure: MetadataAnnotatedStructure
```

---

## 7.4 `register_annotation_plugin()`

### Responsibility

Register a user-defined annotation plugin into the plugin registry.

### Input Schema

```yaml
register_annotation_plugin:
  plugin: Plugin
```

### Output Schema

```yaml
register_annotation_plugin_result:
  status: string
  # registered        — Plugin registered successfully.
  # already_registered — A plugin with this plugin_id already exists.
  #                      Registration is skipped; existing plugin unchanged.
  # rejected          — Plugin definition is invalid (missing required fields,
  #                     unknown plugin_type, etc.). See diagnostics for details.

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

---

## 7.5 `apply_annotation_plugin()`

### Responsibility

Apply one annotation plugin to one structure and return a single
`AnnotationPluginResult`.

### Input Schema

```yaml
apply_annotation_plugin:
  structure:
    CanonicalStructureResult | MetadataAnnotatedStructure

  plugin: Plugin

  config: object | null
  # Runtime config overrides for this specific application.
  # Merged with plugin.config; runtime values take precedence.
```

### Output Schema

```yaml
apply_annotation_plugin_result:
  annotation_result: AnnotationPluginResult
```

---

## 7.6 `apply_annotation_plugins()`

### Responsibility

Apply multiple annotation plugins and produce an
`AnnotatedStructureWithPlugins` containing all resulting annotation layers.

### Internal Workflow

```text
for each plugin in plugins:
    result = apply_annotation_plugin(structure, plugin, config)

    if result.status == "failed":
        if plugin_policy.fail_on_plugin_error: true
            abort and return error
        else:
            record failure in diagnostics; continue

    if result.status in ("success", "warning"):
        if plugin_policy.allow_partial_annotation: true or result is success
            append result.annotation_layer to annotation_layers

assemble AnnotatedStructureWithPlugins(
    canonical_structure_result=...,
    metadata_annotations=...,
    derived_annotations=annotation_layers,
    applied_plugins=plugins,
    ...
)
```

In parallel mode, plugins are dispatched concurrently up to `max_workers`.
If `fail_fast` is true, remaining plugins are cancelled on the first failure.

### Input Schema

```yaml
apply_annotation_plugins:
  structure:
    CanonicalStructureResult | MetadataAnnotatedStructure

  plugins: list[Plugin]

  mode: string
  # sequential | parallel

  parallel_options:
    max_workers: int | null
    fail_fast: bool

  plugin_policy: AnnotationPluginPolicy
```

### Output Schema

```yaml
apply_annotation_plugins_result:
  annotated_structure: AnnotatedStructureWithPlugins
```

---

## 7.7 `validate_annotation_output()`

### Responsibility

Validate the output of an annotation plugin.

### Input Schema

```yaml
validate_annotation_output:
  annotation_result: AnnotationPluginResult
```

### Output Schema

```yaml
validate_annotation_output_result:
  validation_status: string
  # valid | warning | invalid

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

### V1 Validation Rules

```yaml
error_rules:
  PLUGIN_OUTPUT_NULL:
    condition: "annotation_layer is null but plugin status is 'success'."
    result_status: invalid

  UNKNOWN_LAYER_TYPE:
    condition: "annotation_layer.layer_type is not a known built-in type and
                is not namespaced (does not contain a '.')."
    result_status: invalid

warning_rules:
  MISSING_SCHEMA_VERSION:
    condition: "annotation_layer.schema_version is null."
    result_status: warning
```

---

## 7.8 `attach_metadata_many()`

### Responsibility

Run metadata integration and annotation plugins over many structures.

### Input Schema

```yaml
attach_metadata_many:
  structures:
    - CanonicalStructureResult

  metadata_policy: MetadataIntegrationPolicy

  plugin_policy: AnnotationPluginPolicy | null

  mode: string
  # sequential | parallel

  parallel_options:
    max_workers: int | null
    fail_fast: bool
```

### Output Schema

```yaml
attach_metadata_many_result:
  result: MetadataAndAnnotationBatchResult
```

---

# 8. Internal Helper Functions

## 8.1 `retrieve_pdbe_metadata()`

### Responsibility

Retrieve PDBe-specific archive metadata and structural annotations.

### Input

```yaml
retrieve_pdbe_metadata:
  entry_id: string
  include_categories:
    archive_metadata: bool
    structural_annotations: bool
  use_cache: bool
```

### Output

```yaml
retrieve_pdbe_metadata_result:
  archive_metadata: ArchiveMetadata | null
  structural_annotations: StructuralAnnotations | null
  retrieval_status: string
  # success | partial | failed
  diagnostics: list[Diagnostic]
```

---

## 8.2 `retrieve_uniprot_mappings()`

### Responsibility

Retrieve UniProt mappings for a canonical structure.

### Input

```yaml
retrieve_uniprot_mappings:
  entry_id: string
  canonical_chain_ids: list[string]
  use_cache: bool
```

### Output

```yaml
retrieve_uniprot_mappings_result:
  uniprot_mappings: list[UniProtMapping]
  retrieval_status: string
  diagnostics: list[Diagnostic]
```

---

## 8.3 `retrieve_sifts_mappings()`

### Responsibility

Retrieve SIFTS cross-reference mappings.

### Input

```yaml
retrieve_sifts_mappings:
  entry_id: string
  canonical_chain_ids: list[string]
  include_residue_level: bool
  use_cache: bool
```

### Output

```yaml
retrieve_sifts_mappings_result:
  sifts_mappings: list[SIFTSMapping]
  retrieval_status: string
  diagnostics: list[Diagnostic]
```

---

## 8.4 `retrieve_taxonomy_metadata()`

### Responsibility

Retrieve taxonomy metadata for the deposited organism.

### Input

```yaml
retrieve_taxonomy_metadata:
  entry_id: string
  use_cache: bool
```

### Output

```yaml
retrieve_taxonomy_metadata_result:
  taxonomy: TaxonomyRecord | null
  retrieval_status: string
  diagnostics: list[Diagnostic]
```

---

## 8.5 `merge_metadata_layers()`

### Responsibility

Merge per-source metadata records into a single `MetadataRecord`.

### Input

```yaml
merge_metadata_layers:
  entry_id: string
  archive_metadata: ArchiveMetadata | null
  biological_mappings_parts:
    uniprot_mappings: list[UniProtMapping]
    sifts_mappings: list[SIFTSMapping]
    taxonomy: TaxonomyRecord | null
  structural_annotations: StructuralAnnotations | null
  source_records: list[MetadataSourceRecord]
```

### Output

```yaml
merge_metadata_layers_result:
  metadata_record: MetadataRecord
```

---

## 8.6 `merge_annotation_layers()`

### Responsibility

Merge multiple `AnnotationPluginResult` outputs into a unified list of
`AnnotationLayer` records for inclusion in `AnnotatedStructureWithPlugins`.

### Input

```yaml
merge_annotation_layers:
  plugin_results: list[AnnotationPluginResult]
  allow_partial: bool
```

### Output

```yaml
merge_annotation_layers_result:
  annotation_layers: list[AnnotationLayer]
  diagnostics: list[Diagnostic]
```

---

# 9. Policy Schemas

## 9.1 Metadata integration policy

```yaml
MetadataIntegrationPolicy:
  policy_id: string
  policy_name: string
  policy_version: string
  description: string

  include_sources:
    pdbe: bool
    sifts: bool
    uniprot: bool
    taxonomy: bool

  include_categories:
    archive_metadata: bool
    biological_mappings: bool
    structural_annotations: bool
    provenance_metadata: bool
    # provenance_metadata is always populated regardless of this flag.

  retrieval_rules:
    eager: bool
    # If true, all enabled sources are fetched concurrently on first access.
    # If false, sources are fetched lazily when their category is first read.

    use_cache: bool
    fail_on_missing: bool
    # If true, any source failure causes the entire metadata retrieval to fail.
    # Takes precedence over allow_partial.

    allow_partial: bool
    # If true, partial results are returned when some sources fail.
    # Only applies when fail_on_missing: false.

  provenance_rules:
    record_source_versions: bool
    record_retrieval_time: bool
    record_mapping_history: bool
```

---

## 9.2 Annotation plugin policy

```yaml
AnnotationPluginPolicy:
  policy_id: string
  policy_name: string
  policy_version: string
  description: string

  execution_rules:
    mode: string
    # sequential | parallel

    parallel_options:
      max_workers: int | null
      fail_fast: bool

    fail_on_plugin_error: bool
    # If true, any plugin failure causes the entire annotation step to fail.

    allow_partial_annotation: bool
    # If true, successfully completed plugin layers are retained even when
    # other plugins fail.
    # Only applies when fail_on_plugin_error: false.

  provenance_rules:
    record_plugin_versions: bool
    record_plugin_configuration: bool
    record_execution_history: bool
```

---

# 10. Non-Responsibilities

Component 03 is not responsible for:
  - canonical structure modification
  - dataset filtering
  - leakage-safe splitting
  - benchmark generation
  - embeddings
  - graph generation
  - model training
  - automatic biological inference outside plugin scope

---

# 11. Component Definition

The Metadata Integration & Derived Annotation Layer attaches provenance-aware metadata to canonical structures and optionally applies user-defined plugin-style derived annotations for task-specific structural analyses.
