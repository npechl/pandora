# Component 04 — Dataset Construction & Curation Layer

## Purpose

The Dataset Construction & Curation Layer transforms collections of annotated
structures into reproducible, curated datasets suitable for downstream analyses.

This component is responsible for:

* dataset construction,
* structure selection and filtering,
* exact deduplication,
* inclusion and exclusion criteria,
* dataset provenance tracking,
* and reproducible dataset generation.

This component is where Pandora transitions from being structure-centric to
dataset-centric.

---

# 1. Architectural Role

```text
MmCIFIngestionResult
  → CanonicalStructureResult
  → MetadataAnnotatedStructure
  → AnnotatedStructureWithPlugins
  → Dataset
```

---

# 2. Core Design Principles

## Dataset-centric abstraction

Components 1–3 focus on individual structures.

Component 4 introduces the concept of a `Dataset` as a first-class object.

---

## Policy-driven curation

Dataset construction is fully reproducible.

Users explicitly define:

* selection criteria,
* filtering criteria,
* exact deduplication strategies,
* and provenance requirements.

---

## Composable filtering

Internally, dataset construction is implemented through reusable filtering
operations:

```text
apply_selection_rules()
apply_quality_filters()
apply_content_filters()
apply_organism_filters()
deduplicate_dataset()
```

Publicly, these are orchestrated through a single `build_dataset()` call
using a `DatasetCurationPolicy`.

---

## Content rules vs selection rules

**Selection rules, quality rules, and organism rules** determine *which
structures enter the dataset*. A structure failing any of these rules is
excluded entirely and recorded in `excluded_items`.

**Content rules** determine *what is kept within each retained structure*
(e.g. whether ligands or waters are included in the structure's content).
They do not exclude whole structures.

---

## Separation from leakage prevention

Component 4 determines:

```text
Which structures enter the dataset?
```

Component 5 determines:

```text
How is the dataset partitioned into train/validation/test splits?
```

These responsibilities are strictly separate. Similarity-based deduplication
(sequence identity, structure identity) requires pairwise similarity
computation and belongs to Component 5. Component 4 supports only exact
deduplication (`entry_id` and `exact_hash`).

---

# 3. Input Schemas

## 3.1 Dataset construction input

```yaml
DatasetConstructionInput:
  annotated_structures:
    - AnnotatedStructureWithPlugins
  # Structures with status "failed" from upstream components are silently
  # skipped and recorded in excluded_items with reason_code UPSTREAM_FAILURE.

  curation_policy: DatasetCurationPolicy
  # See Section 10 for the full policy schema.
```

---

## 3.2 Batch dataset construction input

```yaml
DatasetConstructionBatchInput:
  batches:
    - DatasetConstructionInput
  # Each element is an independent dataset construction job producing
  # one Dataset. This is NOT one large pool of structures — it is
  # N separate DatasetConstructionInputs producing N separate Datasets.

  mode: string
  # sequential | parallel

  parallel_options:
    max_workers: int | null
    # Number of concurrent dataset-construction jobs in parallel mode.
    # null uses the system default.
    # Ignored in sequential mode.

    fail_fast: bool
    # If true, abort remaining jobs on the first failure.
    # If false (default), isolate failures and continue.
```

---

# 4. Core Object Schemas

## 4.1 ExclusionRecord

Records a structure that was removed from the dataset and the reason for
its removal.

```yaml
ExclusionRecord:
  entry_id: string

  reason_code: string
  # UPSTREAM_FAILURE      — input had status "failed" from C01/C02/C03.
  # RESOLUTION_THRESHOLD  — resolution exceeded max_resolution.
  # NULL_RESOLUTION       — resolution is null and null_resolution_behavior
  #                         is "exclude".
  # CHAIN_TOO_SHORT       — all polymer chains are shorter than min_chain_length.
  # INCOMPLETE_CHAIN      — chain fails allow_incomplete_chains rule.
  # MISSING_RESIDUES      — missing residues not permitted by policy.
  # MISSING_ATOMS         — missing atoms not permitted by policy.
  # METHOD_EXCLUDED       — experimental method not in include_experimental_methods
  #                         or present in exclude_experimental_methods.
  # SOURCE_EXCLUDED       — source not in include_sources or in exclude_sources.
  # BIOMOLECULE_EXCLUDED  — biomolecule type not in include_biomolecules.
  # ORGANISM_EXCLUDED     — organism taxon in exclude_taxa, or not in include_taxa
  #                         when include_taxa is non-empty.
  # MISSING_TAXONOMY      — taxonomy metadata unavailable and organism filtering
  #                         is required.
  # DUPLICATE             — removed by deduplication.

  reason_message: string
  # Human-readable description including the offending value.
  # e.g. "Resolution 3.2 Å exceeds max_resolution threshold of 2.5 Å"

  applied_rule: string | null
  # The policy field that triggered this exclusion.
  # e.g. "quality_rules.max_resolution"
```

---

## 4.2 DeduplicationReport

```yaml
DeduplicationReport:
  enabled: bool

  strategy: string
  # entry_id   — exact match on entry_id string.
  # exact_hash — SHA-256 of the serialised CanonicalStructure content.

  duplicates_found: int

  removed_items: list[ExclusionRecord]
  # All items have reason_code == "DUPLICATE".
  # reason_message identifies the retained representative entry.
```

---

## 4.3 AppliedFilterRecord

Records a single filter operation applied during dataset construction,
including how many structures it removed.

```yaml
AppliedFilterRecord:
  filter_name: string
  # The policy field name. e.g. "max_resolution", "include_taxa"

  filter_category: string
  # quality | selection | content | organism | deduplication

  filter_value: object
  # The configured threshold or value, serialised for provenance.
  # e.g. { max_resolution: 2.5 } or { include_taxa: ["9606"] }

  structures_excluded: int
  # Number of structures removed by this filter alone.
```

---

## 4.4 SelectionSummary

```yaml
SelectionSummary:
  applied_filters: list[AppliedFilterRecord]
  # Ordered list of all filters that were active, including those that
  # excluded zero structures.
```

---

# 5. Output Schemas

## 5.1 Dataset object

```yaml
Dataset:
  dataset_id: string
  dataset_name: string
  dataset_version: string

  structures: list[AnnotatedStructureWithPlugins]
  # Final curated set after all filtering and deduplication.

  counts:
    total_input: int
    total_selected: int
    total_excluded: int
    total_duplicates_removed: int

  selection_summary: SelectionSummary

  excluded_items: list[ExclusionRecord]
  # Every structure that was removed, with reason codes.

  deduplication_report: DeduplicationReport

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]

  provenance:
    created_at: string | null
    # ISO 8601 timestamp.

    source_count: int
    # Number of unique entry_ids in the input.

    input_sources: list[string]
    # Provider names from input structures.
    # e.g. ["pdbe", "pdb"]
```

---

## 5.2 Batch dataset construction result

```yaml
DatasetConstructionBatchResult:
  mode: string
  # sequential | parallel

  summary:
    total: int
    success: int
    warning: int
    failed: int

  results:
    - dataset_id: string

      status: string
      # success | warning | failed

      dataset: Dataset | null
      # null when status == "failed".

      diagnostics:
        warnings: list[Diagnostic]
        errors: list[Diagnostic]
```

---

# 6. V1 Filter Definitions

This section defines the behaviour contract for each filter type applied
during dataset construction.

## 6.1 Resolution filter

```yaml
resolution_filter:
  policy_field: quality_rules.max_resolution
  checks: archive_metadata.resolution
  excludes_when:
    - resolution > max_resolution
    - resolution is null AND null_resolution_behavior == "exclude"
  null_resolution_behavior: string
  # exclude (default) — structures with null resolution are excluded.
  #                     Emits NULL_RESOLUTION_EXCLUDED warning.
  # include           — structures with null resolution are retained.
  notes: >
    NMR and some EM structures have null resolution. The default
    "exclude" behaviour is conservative. Set to "include" to retain
    NMR/EM structures when max_resolution is active.
```

---

## 6.2 Chain length filter

```yaml
chain_length_filter:
  policy_field: quality_rules.min_chain_length
  checks: minimum polymer chain length across all chains in the structure
  excludes_when:
    - all polymer chains have fewer than min_chain_length residues
  notes: >
    A structure is excluded only if ALL polymer chains are shorter than
    the threshold. Structures with at least one qualifying chain are retained.
```

---

## 6.3 Completeness filters

```yaml
completeness_filters:
  policy_fields:
    - quality_rules.allow_incomplete_chains
    - quality_rules.allow_missing_residues
    - quality_rules.allow_missing_atoms
  checks: diagnostics from C01 (parse) and C02 (canonicalization)
  excludes_when:
    - allow_incomplete_chains: false AND structure has INCOMPLETE_CHAIN diagnostic
    - allow_missing_residues: false AND structure has SEQUENCE_GAP diagnostic
    - allow_missing_atoms: false AND structure has MISSING_ATOMS diagnostic
```

---

## 6.4 Experimental method filter

```yaml
method_filter:
  policy_fields:
    - selection_rules.include_experimental_methods
    - selection_rules.exclude_experimental_methods
  checks: archive_metadata.experimental_method
  valid_method_values:
    - "X-RAY DIFFRACTION"
    - "SOLUTION NMR"
    - "SOLID-STATE NMR"
    - "ELECTRON MICROSCOPY"
    - "ELECTRON CRYSTALLOGRAPHY"
    - "FIBER DIFFRACTION"
    - "NEUTRON DIFFRACTION"
    - "SOLUTION SCATTERING"
  excludes_when:
    - include_experimental_methods is non-empty AND method not in list
    - method is in exclude_experimental_methods
  notes: >
    If both include and exclude lists are non-empty and a method appears
    in both, exclusion takes precedence.
```

---

## 6.5 Organism filter

```yaml
organism_filter:
  policy_fields:
    - organism_rules.include_taxa
    - organism_rules.exclude_taxa
  checks: biological_mappings.taxonomy.ncbi_taxon_id
  taxon_id_format: >
    NCBI taxon IDs as strings (e.g. "9606" for Homo sapiens,
    "10090" for Mus musculus). Scientific names are not supported
    as filter values in V1 to avoid ambiguity.
  excludes_when:
    - include_taxa is non-empty AND taxon_id not in include_taxa
      (and no ancestor taxon_id in include_taxa)
    - taxon_id in exclude_taxa (or any ancestor taxon_id in exclude_taxa)
    - taxonomy metadata is unavailable AND either filter list is non-empty
      (reason_code: MISSING_TAXONOMY)
  notes: >
    Taxon matching is ancestry-aware: specifying "9606" (Homo sapiens)
    also matches any entry whose taxonomy lineage includes 9606.
    Specify a higher-level taxon ID to match all organisms in that clade.
    If include_taxa is empty, all organisms are included (subject to
    exclude_taxa).
```

---

## 6.6 Biomolecule type filter

```yaml
biomolecule_filter:
  policy_field: selection_rules.include_biomolecules
  checks: entity types in canonical_structure.entities
  excludes_when:
    - include_biomolecules.proteins: false AND structure has only protein entities
    - include_biomolecules.rna: false AND structure has only RNA entities
    - include_biomolecules.dna: false AND structure has only DNA entities
    - include_biomolecules.complexes: false AND structure has entities of
      multiple types
  notes: >
    A structure is retained if it contains at least one entity type that
    is enabled in include_biomolecules. All flags false results in an
    empty dataset (a configuration error).
```

---

## 6.7 Content filters

```yaml
content_filters:
  policy_field: content_rules
  behaviour: >
    Content rules do NOT exclude structures. They control what is retained
    within each structure's content after the structure has been selected.
  actions:
    keep_ligands: false     → remove Ligand records where is_water: false,
                              is_ion: false
    keep_waters: false      → remove Ligand records where is_water: true
    keep_ions: false        → remove Ligand records where is_ion: true
    keep_nonpolymer_entities: false → remove all non-polymer entity content
```

---

# 7. Public Functions

## 7.1 `build_dataset()`

### Responsibility

Construct a curated dataset from a collection of annotated structures.
This is the main orchestrator for Component 04.

### Internal Workflow

```text
1. Reject upstream failures:
   Structures with status "failed" from C01/C02/C03 are moved to
   excluded_items with reason_code UPSTREAM_FAILURE.

2. Apply selection_rules:
   apply_selection_rules()
     — filter by source, biomolecule type, experimental method

3. Apply quality_rules:
   apply_quality_filters()
     — resolution, chain length, completeness (missing atoms/residues/chains)

4. Apply organism_rules:
   apply_organism_filters()
     — include/exclude by NCBI taxon ID

5. Apply content_rules:
   apply_content_filters()
     — strip ligands, waters, ions from retained structures per policy

6. Apply deduplication_rules (if enabled):
   deduplicate_dataset()
     — entry_id: remove entries with identical entry_id strings
     — exact_hash: remove entries with identical SHA-256 of CanonicalStructure

7. Validate:
   validate_dataset()

8. Assemble Dataset with counts, SelectionSummary, ExclusionRecords,
   DeduplicationReport, and provenance.
```

### Input Schema

```yaml
build_dataset:
  input: DatasetConstructionInput
```

### Output Schema

```yaml
build_dataset_result:
  dataset: Dataset
```

---

## 7.2 `filter_dataset()`

### Responsibility

Apply one or more explicit filter operations to a structure list or existing
dataset.

### Input Schema

```yaml
filter_dataset:
  structures: list[AnnotatedStructureWithPlugins]

  filters: list[FilterSpec]
```

```yaml
FilterSpec:
  filter_type: string
  # resolution | chain_length | completeness | experimental_method |
  # organism | biomolecule_type | source
  # See Section 6 for behaviour contracts per filter_type.

  parameters: object
  # Type-specific parameters matching the corresponding policy field.
  # e.g. { max_resolution: 2.5, null_resolution_behavior: "exclude" }
```

### Output Schema

```yaml
filter_dataset_result:
  retained: list[AnnotatedStructureWithPlugins]
  excluded: list[ExclusionRecord]

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

---

## 7.3 `deduplicate_dataset()`

### Responsibility

Remove exact duplicate structures according to the deduplication strategy.

Supported strategies are `entry_id` and `exact_hash` only.
Similarity-based deduplication (sequence identity, structure identity)
requires pairwise similarity computation and is handled by Component 05.

### Input Schema

```yaml
deduplicate_dataset:
  structures: list[AnnotatedStructureWithPlugins]

  deduplication_rules:
    enabled: bool

    strategy: string
    # entry_id   — exact match on entry_id string.
    #              The first occurrence is retained; subsequent duplicates
    #              are removed.
    # exact_hash — SHA-256 of the serialised CanonicalStructure content.
    #              Structurally identical entries are removed regardless of
    #              entry_id.
```

### Output Schema

```yaml
deduplicate_dataset_result:
  deduplicated_structures: list[AnnotatedStructureWithPlugins]
  deduplication_report: DeduplicationReport
```

---

## 7.4 `validate_dataset()`

### Responsibility

Validate dataset consistency, provenance completeness, and policy compliance.

### V1 Validation Rules

```yaml
error_rules:
  EMPTY_DATASET:
    condition: "Dataset.structures is empty after all filtering and deduplication."
    result_status: invalid

  ALL_BIOMOLECULES_EXCLUDED:
    condition: "All include_biomolecules flags are false in the policy."
    result_status: invalid

  PROVENANCE_INCOMPLETE:
    condition: "A provenance record_* rule is true but the corresponding
                field in Dataset.provenance is null."
    result_status: invalid

warning_rules:
  HIGH_EXCLUSION_RATE:
    condition: "More than 50% of total_input structures were excluded."
    result_status: warning

  AGGRESSIVE_DEDUPLICATION:
    condition: "More than 30% of structures were removed by deduplication."
    result_status: warning

  EMPTY_SOURCE_CONTRIBUTION:
    condition: "An entry in include_sources contributed zero structures
                to the final dataset."
    result_status: warning

  NULL_RESOLUTION_EXCLUDED:
    condition: "max_resolution is set and structures with null resolution
                were excluded (NMR or EM structures may be affected)."
    result_status: warning

  MISSING_TAXONOMY_SKIPPED:
    condition: "Organism filtering was active but some structures had no
                taxonomy metadata and were excluded with MISSING_TAXONOMY."
    result_status: warning
```

### Status determination

```yaml
status_rules:
  failed:  Any error_rule fires.
  warning: No error_rules fire, but one or more warning_rules fire.
  valid:   No rules fire.
```

### Input Schema

```yaml
validate_dataset:
  dataset: Dataset
  policy: DatasetCurationPolicy
```

### Output Schema

```yaml
validate_dataset_result:
  validation_status: string
  # valid | warning | invalid

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

---

## 7.5 `build_dataset_many()`

### Responsibility

Construct multiple independent datasets from multiple
`DatasetConstructionInput` jobs.

Each element in `batches` is an independent job producing one `Dataset`.
This is not a single large pool of structures — it is N inputs producing
N datasets.

### Input Schema

```yaml
build_dataset_many:
  batches: list[DatasetConstructionInput]

  mode: string
  # sequential | parallel

  parallel_options:
    max_workers: int | null
    fail_fast: bool
```

### Output Schema

```yaml
build_dataset_many_result:
  result: DatasetConstructionBatchResult
```

---

# 8. Internal Helper Functions

## 8.1 `apply_selection_rules()`

### Responsibility

Filter structures by source, biomolecule type, and experimental method.

### Input

```yaml
apply_selection_rules:
  structures: list[AnnotatedStructureWithPlugins]
  selection_rules:
    include_sources: list[string]
    exclude_sources: list[string]
    include_biomolecules:
      proteins: bool
      rna: bool
      dna: bool
      complexes: bool
    include_experimental_methods: list[string]
    exclude_experimental_methods: list[string]
```

### Output

```yaml
apply_selection_rules_result:
  retained: list[AnnotatedStructureWithPlugins]
  excluded: list[ExclusionRecord]
```

---

## 8.2 `apply_quality_filters()`

### Responsibility

Filter structures by resolution, chain length, and completeness.

### Input

```yaml
apply_quality_filters:
  structures: list[AnnotatedStructureWithPlugins]
  quality_rules:
    max_resolution: float | null
    null_resolution_behavior: string
    min_chain_length: int | null
    allow_incomplete_chains: bool
    allow_missing_residues: bool
    allow_missing_atoms: bool
```

### Output

```yaml
apply_quality_filters_result:
  retained: list[AnnotatedStructureWithPlugins]
  excluded: list[ExclusionRecord]
```

---

## 8.3 `apply_content_filters()`

### Responsibility

Strip ligand, water, ion, and non-polymer content from retained structures
according to content rules. Does not exclude whole structures.

### Input

```yaml
apply_content_filters:
  structures: list[AnnotatedStructureWithPlugins]
  content_rules:
    keep_ligands: bool
    keep_waters: bool
    keep_ions: bool
    keep_nonpolymer_entities: bool
```

### Output

```yaml
apply_content_filters_result:
  structures: list[AnnotatedStructureWithPlugins]
  # Same list length as input; content of individual structures may differ.
```

---

## 8.4 `apply_organism_filters()`

### Responsibility

Filter structures by NCBI taxon ID using inclusion and exclusion lists.

### Input

```yaml
apply_organism_filters:
  structures: list[AnnotatedStructureWithPlugins]
  organism_rules:
    include_taxa: list[string]
    # NCBI taxon IDs as strings. Empty list = include all.
    exclude_taxa: list[string]
    # NCBI taxon IDs as strings. Empty list = exclude none.
```

### Output

```yaml
apply_organism_filters_result:
  retained: list[AnnotatedStructureWithPlugins]
  excluded: list[ExclusionRecord]
```

---

## 8.5 `summarize_exclusions()`

### Responsibility

Aggregate all `ExclusionRecord` lists from each filter step into a unified
list and compute per-filter exclusion counts for `SelectionSummary`.

### Input

```yaml
summarize_exclusions:
  exclusion_batches: list[list[ExclusionRecord]]
  applied_filter_records: list[AppliedFilterRecord]
```

### Output

```yaml
summarize_exclusions_result:
  all_exclusions: list[ExclusionRecord]
  selection_summary: SelectionSummary
```

---

## 8.6 `build_dataset_provenance()`

### Responsibility

Construct the provenance record for the `Dataset` from input structures
and applied policy.

### Input

```yaml
build_dataset_provenance:
  input_structures: list[AnnotatedStructureWithPlugins]
  policy: DatasetCurationPolicy
```

### Output

```yaml
build_dataset_provenance_result:
  provenance:
    created_at: string
    # ISO 8601 timestamp.
    source_count: int
    input_sources: list[string]
```

---

# 9. Dataset Unit

## V1 Dataset Unit

```yaml
dataset_unit:
  AnnotatedStructureWithPlugins
```

The initial dataset abstraction is structure-based.

Future downstream components may generate:

* chain-level datasets,
* interface-level datasets,
* residue-level datasets,
* task-specific datasets.

---

# 10. Policy Schema

## 10.1 DatasetCurationPolicy

```yaml
DatasetCurationPolicy:
  policy_id: string
  policy_name: string
  policy_version: string
  description: string

  selection_rules:
    include_sources: list[string]
    # e.g. ["pdbe", "pdb"]. Empty list = include all sources.

    exclude_sources: list[string]
    # e.g. []. Empty list = exclude no sources.

    include_biomolecules:
      proteins: bool
      rna: bool
      dna: bool
      complexes: bool

    include_experimental_methods: list[string]
    # Valid values defined in Section 6.4.
    # Empty list = include all methods.

    exclude_experimental_methods: list[string]
    # Empty list = exclude no methods.

  quality_rules:
    max_resolution: float | null
    # null = no resolution filter applied.

    null_resolution_behavior: string
    # exclude (default) | include
    # Behaviour when max_resolution is set and a structure has null resolution.

    min_chain_length: int | null
    # null = no chain length filter applied.

    allow_incomplete_chains: bool
    allow_missing_residues: bool
    allow_missing_atoms: bool

  content_rules:
    keep_ligands: bool
    keep_waters: bool
    keep_ions: bool
    keep_nonpolymer_entities: bool

  organism_rules:
    include_taxa: list[string]
    # NCBI taxon IDs as strings. Empty list = include all organisms.
    exclude_taxa: list[string]
    # NCBI taxon IDs as strings. Empty list = exclude no organisms.

  deduplication_rules:
    enabled: bool

    strategy: string
    # entry_id   — remove entries with identical entry_id.
    # exact_hash — remove entries with identical CanonicalStructure content
    #              (SHA-256 hash).
    #
    # Note: similarity-based deduplication (sequence_identity,
    # structure_identity) is NOT supported here. It requires pairwise
    # similarity computation and is handled by Component 05.

  provenance_rules:
    record_filters: bool
    record_exclusions: bool
    record_dataset_version: bool
    record_input_sources: bool
```

---

# 11. Non-Responsibilities

Component 04 is not responsible for:
  - train_test_splitting
  - leakage_safe_partitioning
  - similarity_clustering
  - similarity_based_deduplication   # belongs to Component 05
  - benchmark_generation
  - model_training
  - graph_generation
  - embeddings
  - task_specific_labels

---

# 12. Component Definition

The Dataset Construction & Curation Layer transforms collections of annotated structures into reproducible curated datasets using explicit curation policies and provenance-aware filtering operations.
