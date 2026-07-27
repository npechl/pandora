# Component 04 — Dataset Construction & Curation Layer

## Purpose

The Dataset Construction & Curation Layer transforms collections of annotated
structures into reproducible, curated datasets suitable for downstream
analyses.

This component is responsible for:

* dataset construction at structure, chain, interface, and residue granularity,
* structure selection and filtering,
* exact deduplication,
* inclusion and exclusion criteria,
* dataset provenance tracking,
* materialization of large-scale datasets to disk,
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
  → Dataset                  [structure-level]
      → ChainDataset         [chain-level]
          → ResidueDataset   [residue-level]
      → InterfaceDataset     [interface-level]
          → ResidueDataset   [residue-level] (via chain_record_1/2)
```

Extraction is always a downstream step from the structure-level `Dataset`.
Chain and interface datasets are independent siblings; residue datasets
can be extracted from either chain or structure granularity.

---

# 2. Core Design Principles

## Dataset-centric abstraction

Components 1–3 focus on individual structures.

Component 4 introduces the concept of a `Dataset` as a first-class object
and provides four granularity levels: structure, chain, interface, and
residue.

---

## Multi-level granularity

Different ML tasks require different units of analysis:

* **Structure-level** — fold classification, global structure quality,
  structure-based drug design at entry level.
* **Chain-level** — sequence models, secondary structure prediction,
  homology detection, protein language model pre-training.
* **Interface-level** — protein-protein interaction prediction, binding
  affinity models, antibody-antigen studies, docking benchmarks.
* **Residue-level** — binding site prediction, mutation effect prediction,
  contact map prediction, solvent accessibility tasks.

Extraction always starts from a curated structure-level `Dataset`. Each
level is a typed, self-contained record that embeds the data needed for
downstream ML without requiring access to the parent structure.

---

## Two execution modes

**In-memory mode** — default for datasets up to ~10K structures.

Annotated structures are passed as a list. Filtered records accumulate in
memory. The resulting `Dataset` has `mode == "in_memory"` and `structures`
is fully populated.

**Materialized (streaming) mode** — required for datasets above ~10K
structures.

Annotated structures are supplied as an iterable/generator and processed
one at a time. Per-structure selection and quality filters are applied on
the fly. Records that pass are written immediately to a `DatasetStore`
(Parquet files on disk). The resulting `Dataset` has `mode == "materialized"`,
`structures` is empty, and `store` holds a `DatasetStoreRef` pointing to
the written data.

Deduplication in streaming mode is performed as a secondary pass over the
stored record IDs: exact-hash deduplication streams the stored IDs,
computes the deduplication set, then removes duplicates from the store.
No full structure objects are loaded during this pass.

All extraction functions (`extract_chain_records`, `extract_interface_records`,
`extract_residue_records`) support both modes.

---

## Policy-driven curation

Dataset construction is fully reproducible.

Users explicitly define:

* selection criteria,
* filtering criteria,
* exact deduplication strategies,
* extraction granularity and rules,
* execution mode,
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
Which structures (or chains, interfaces, residues) enter the dataset?
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
    list[AnnotatedStructureWithPlugins] | Iterable[AnnotatedStructureWithPlugins]
  # In in-memory mode: pass a list.
  # In streaming mode: pass any iterable or generator. Structures are
  # consumed one at a time; do not re-iterate after build_dataset() returns.
  # Structures with status "failed" from upstream components are silently
  # skipped and recorded in excluded_items with reason_code UPSTREAM_FAILURE.

  curation_policy: DatasetCurationPolicy
  # See Section 11 for the full policy schema.

  store: DatasetStoreRef | null
  # null = in-memory mode (default).
  # non-null = streaming/materialized mode. The store must be initialised
  # before calling build_dataset(). Use init_dataset_store() to create one.
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
    fail_fast: bool
```

---

# 4. Core Object Schemas

## 4.1 DatasetStoreRef

A typed reference to a materialized dataset stored on disk.

```yaml
DatasetStoreRef:
  store_id: string
  # Unique identifier for this store. Generated by init_dataset_store().

  store_format: string
  # parquet
  # V1 supports Parquet only.

  store_path: string
  # Absolute path to the store root directory.
  # Contains one Parquet file per record type:
  #   structures.parquet   — for Dataset
  #   chains.parquet       — for ChainDataset
  #   interfaces.parquet   — for InterfaceDataset
  #   residues.parquet     — for ResidueDataset
  # Multiple granularity files may coexist in the same store directory
  # when extraction was performed in streaming mode.

  granularity: string
  # structure | chain | interface | residue
  # Indicates which record type is present at this store.

  item_count: int
  # Total number of records written to the store.

  created_at: string | null
  # ISO 8601 timestamp.
```

---

## 4.2 ExclusionRecord

Records a unit (structure, chain, interface, or residue) that was removed
from the dataset and the reason for its removal.

```yaml
ExclusionRecord:
  unit_id: string
  # entry_id for structures; chain_id for chains;
  # interface_id for interfaces; residue_id for residues.

  granularity: string
  # structure | chain | interface | residue

  reason_code: string
  # ── Structure-level ────────────────────────────────────────────────
  # UPSTREAM_FAILURE      — input had status "failed" from C01/C02/C03.
  # RESOLUTION_THRESHOLD  — resolution exceeded max_resolution.
  # NULL_RESOLUTION       — resolution is null and null_resolution_behavior
  #                         is "exclude".
  # CHAIN_TOO_SHORT       — all polymer chains are shorter than min_chain_length.
  # INCOMPLETE_CHAIN      — chain fails allow_incomplete_chains rule.
  # MISSING_RESIDUES      — missing residues not permitted by policy.
  # MISSING_ATOMS         — missing atoms not permitted by policy.
  # METHOD_EXCLUDED       — experimental method excluded by policy.
  # SOURCE_EXCLUDED       — source excluded by policy.
  # BIOMOLECULE_EXCLUDED  — biomolecule type excluded by policy.
  # ORGANISM_EXCLUDED     — organism taxon excluded by policy.
  # MISSING_TAXONOMY      — taxonomy metadata unavailable.
  # DUPLICATE             — removed by deduplication.
  # ── Chain-level ────────────────────────────────────────────────────
  # CHAIN_TYPE_EXCLUDED   — chain type not enabled in include_chain_types.
  # CHAIN_LENGTH_EXCLUDED — chain shorter than chain_extraction_rules.min_chain_length.
  # ── Interface-level ────────────────────────────────────────────────
  # MISSING_INTERFACE_ANNOTATIONS — interface annotations not found in
  #                         derived_annotations. Run pandora.builtin.interfaces
  #                         plugin in C03 before extracting interfaces.
  # INTERFACE_TYPE_EXCLUDED       — interface type not enabled in interface_types.
  # INTERFACE_BELOW_MIN_CONTACTS  — fewer contact residues than min_contact_residues.
  # INTERFACE_BELOW_MIN_AREA      — interface area below min_interface_area.
  # ── Residue-level ──────────────────────────────────────────────────
  # RESIDUE_TYPE_EXCLUDED  — residue type not enabled in include_residue_types.
  # INCOMPLETE_BACKBONE    — residue excluded due to require_full_backbone: true.

  reason_message: string
  applied_rule: string | null
```

---

## 4.3 DeduplicationReport

```yaml
DeduplicationReport:
  enabled: bool
  strategy: string
  # entry_id | exact_hash
  duplicates_found: int
  removed_items: list[ExclusionRecord]
```

---

## 4.4 AppliedFilterRecord

```yaml
AppliedFilterRecord:
  filter_name: string
  filter_category: string
  # quality | selection | content | organism | deduplication | extraction
  filter_value: object
  structures_excluded: int
```

---

## 4.5 SelectionSummary

```yaml
SelectionSummary:
  applied_filters: list[AppliedFilterRecord]
```

---

## 4.6 ChainRecord

A self-contained chain-level record. One record per polymer chain per
structure. Embeds all data needed for downstream ML without access to
the parent structure.

```yaml
ChainRecord:
  chain_id: string
  # Canonical chain_id from the parent CanonicalStructure.

  entry_id: string
  # Parent PDB entry_id.

  chain_type: string
  # protein | rna | dna

  entity_id: string
  # Canonical entity_id of the polymer entity this chain belongs to.

  sequence: string
  # One-letter canonical amino acid or nucleotide sequence.
  # Sourced from Entity.sequence in the parent CanonicalStructure.

  chain_length: int
  # Number of residues in the sequence.

  residues: list[Residue]
  # All Residue records for this chain from the parent CanonicalStructure.

  metadata:
    archive_metadata: ArchiveMetadata | null
    # Inherited from the parent structure (resolution, method, dates, etc.).

    uniprot_mappings: list[UniProtMapping]
    # Chain-specific UniProt mappings, filtered from the parent
    # BiologicalMappings by this chain_id.

    sifts_mappings: list[SIFTSMapping]
    # Chain-specific SIFTS mappings.

    taxonomy: TaxonomyRecord | null
    # Organism of the source entity for this chain.

  derived_annotations: list[AnnotationLayer]
  # Annotation layers from the parent structure that apply to this chain.
  # Layers are filtered by chain_id where applicable.

  parent_entry_id: string
  # Explicit reference to the parent Dataset entry.

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string
```

---

## 4.7 InterfaceRecord

A self-contained interface-level record representing a physical contact
between two polymer chains. One record per chain pair per structure.

```yaml
InterfaceRecord:
  interface_id: string
  # Format: "{entry_id}_{chain_id_1}_{chain_id_2}"
  # chain_id_1 < chain_id_2 lexicographically.

  entry_id: string

  chain_id_1: string
  chain_id_2: string

  interface_type: string
  # protein_protein | protein_rna | protein_dna | protein_ligand

  chain_record_1: ChainRecord
  chain_record_2: ChainRecord
  # Full ChainRecord for each partner.

  interface_residues_chain_1: list[string]
  # Residue_ids (canonical) of chain_1 residues in contact with chain_2.

  interface_residues_chain_2: list[string]
  # Residue_ids (canonical) of chain_2 residues in contact with chain_1.

  contact_count: int
  # Number of inter-chain residue pairs in contact.

  interface_area: float | null
  # Buried surface area in Å². null if not computed by the interface plugin.

  source_annotation_layer: AnnotationLayer | null
  # The AnnotationLayer from pandora.builtin.interfaces that produced
  # this interface record.

  metadata:
    archive_metadata: ArchiveMetadata | null
    # Inherited from the parent structure.

  parent_entry_id: string

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string
```

---

## 4.8 ResidueRecord

A self-contained residue-level record. One record per residue per chain.
Can be extracted from a structure-level `Dataset` or a `ChainDataset`.

```yaml
ResidueRecord:
  residue_id: string
  # Canonical residue_id from the parent structure.

  entry_id: string
  chain_id: string

  seq_id: int | null
  # Canonical sequence number. null for non-polymer residues.

  comp_id: string
  # Three-letter chemical component code (e.g. "ALA", "GLY").

  residue_type: string
  # amino_acid | nucleotide | non_standard

  atoms: list[Atom]
  # All atom records for this residue.

  neighboring_residues: list[NeighborReference] | null
  # References to residues within context_radius Å.
  # null when context_radius is not set in extraction policy.

  metadata:
    archive_metadata: ArchiveMetadata | null
    # Inherited from the parent structure.

  derived_annotations: list[AnnotationLayer]
  # Residue-level annotation layers (e.g. secondary structure,
  # solvent exposure) from the parent structure.

  parent_entry_id: string
  parent_chain_id: string

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

NeighborReference:
  residue_id: string
  chain_id: string
  distance: float
  # Distance in Å between Cα atoms (or equivalent for non-standard residues).
```

---

# 5. Output Schemas

## 5.1 Dataset (structure-level)

```yaml
Dataset:
  dataset_id: string
  dataset_name: string
  dataset_version: string
  granularity: string
  # Always "structure" for this type.

  mode: string
  # in_memory    — structures list is populated; store is null.
  # materialized — structures list is empty; store holds a DatasetStoreRef.

  structures: list[AnnotatedStructureWithPlugins]
  # Populated in in_memory mode.
  # Empty list in materialized mode — records live on disk in store.

  store: DatasetStoreRef | null
  # null in in_memory mode.
  # DatasetStoreRef in materialized mode pointing to structures.parquet.

  counts:
    total_input: int
    total_selected: int
    total_excluded: int
    total_duplicates_removed: int

  selection_summary: SelectionSummary
  excluded_items: list[ExclusionRecord]
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
    source_count: int
    input_sources: list[string]
```

---

## 5.2 ChainDataset

```yaml
ChainDataset:
  dataset_id: string
  dataset_name: string
  dataset_version: string
  granularity: string
  # Always "chain" for this type.

  mode: string
  # in_memory | materialized

  chains: list[ChainRecord]
  # Populated in in_memory mode. Empty in materialized mode.

  store: DatasetStoreRef | null
  # null in in_memory mode.
  # DatasetStoreRef in materialized mode pointing to chains.parquet.

  source_dataset_id: string
  # dataset_id of the structure-level Dataset this was extracted from.

  counts:
    total_structures_input: int
    total_chains_extracted: int
    total_chains_excluded: int

  excluded_items: list[ExclusionRecord]
  # Chains excluded by chain_extraction_rules.

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]

  provenance:
    created_at: string
    source_dataset_id: string
```

---

## 5.3 InterfaceDataset

```yaml
InterfaceDataset:
  dataset_id: string
  dataset_name: string
  dataset_version: string
  granularity: string
  # Always "interface" for this type.

  mode: string
  # in_memory | materialized

  interfaces: list[InterfaceRecord]
  # Populated in in_memory mode. Empty in materialized mode.

  store: DatasetStoreRef | null
  # null in in_memory mode.
  # DatasetStoreRef in materialized mode pointing to interfaces.parquet.

  source_dataset_id: string

  counts:
    total_structures_input: int
    total_interfaces_extracted: int
    total_interfaces_excluded: int

  excluded_items: list[ExclusionRecord]

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]

  provenance:
    created_at: string
    source_dataset_id: string
```

---

## 5.4 ResidueDataset

```yaml
ResidueDataset:
  dataset_id: string
  dataset_name: string
  dataset_version: string
  granularity: string
  # Always "residue" for this type.

  mode: string
  # in_memory | materialized

  residues: list[ResidueRecord]
  # Populated in in_memory mode. Empty in materialized mode.

  store: DatasetStoreRef | null
  # null in in_memory mode.
  # DatasetStoreRef in materialized mode pointing to residues.parquet.

  source_dataset_id: string
  source_granularity: string
  # "structure" or "chain" — which level this was extracted from.

  counts:
    total_source_units_input: int
    total_residues_extracted: int
    total_residues_excluded: int

  excluded_items: list[ExclusionRecord]

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]

  provenance:
    created_at: string
    source_dataset_id: string
```

---

## 5.5 Batch dataset construction result

```yaml
DatasetConstructionBatchResult:
  mode: string

  summary:
    total: int
    success: int
    warning: int
    failed: int

  results:
    - dataset_id: string
      status: string
      dataset: Dataset | null
      diagnostics:
        warnings: list[Diagnostic]
        errors: list[Diagnostic]
```

---

# 6. V1 Filter Definitions

This section defines the behaviour contract for each filter type applied
during structure-level dataset construction.

## 6.1 Resolution filter

```yaml
resolution_filter:
  policy_field: quality_rules.max_resolution
  checks: archive_metadata.resolution
  excludes_when:
    - resolution > max_resolution
    - resolution is null AND null_resolution_behavior == "exclude"
  null_resolution_behavior: string
  # exclude (default) | include
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
```

---

## 6.3 Completeness filters

```yaml
completeness_filters:
  policy_fields:
    - quality_rules.allow_incomplete_chains
    - quality_rules.allow_missing_residues
    - quality_rules.allow_missing_atoms
  checks: diagnostics from C01 (parse) and C02 (canonicalisation)
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
    If both lists are non-empty and a method appears in both,
    exclusion takes precedence.
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
    NCBI taxon IDs as strings (e.g. "9606" for Homo sapiens).
    Matching is ancestry-aware.
  excludes_when:
    - include_taxa non-empty AND taxon_id not in lineage
    - taxon_id in exclude_taxa lineage
    - taxonomy unavailable AND either filter list is non-empty
```

---

## 6.6 Biomolecule type filter

```yaml
biomolecule_filter:
  policy_field: selection_rules.include_biomolecules
  checks: entity types in canonical_structure.entities
  excludes_when:
    - no entity type in structure matches any enabled biomolecule flag
```

---

## 6.7 Content filters

```yaml
content_filters:
  policy_field: content_rules
  behaviour: >
    Content rules do NOT exclude structures. They strip content
    from retained structures after selection.
  actions:
    keep_ligands: false          → remove non-water, non-ion Ligand records
    keep_waters: false           → remove Ligand records where is_water: true
    keep_ions: false             → remove Ligand records where is_ion: true
    keep_nonpolymer_entities: false → remove all non-polymer entity content
```

---

# 7. V1 Extraction Definitions

This section defines the behaviour contract for each extraction type
applied after structure-level curation.

## 7.1 Chain extraction

```yaml
chain_extraction:
  input: Dataset (structure-level), either mode
  output: ChainDataset, same mode as input Dataset
  policy_field: extraction_rules.chain_extraction_rules

  for_each_structure_in_dataset:
    for_each_chain_in_canonical_structure:
      if chain.chain_type not in include_chain_types:
        exclude with CHAIN_TYPE_EXCLUDED
      if chain_length < min_chain_length:
        exclude with CHAIN_LENGTH_EXCLUDED
      else:
        create ChainRecord:
          - embed sequence and residues from CanonicalStructure
          - filter UniProtMapping and SIFTSMapping to this chain_id
          - inherit archive_metadata and taxonomy from parent structure
          - filter derived_annotations to layers applicable to this chain

  materialized_mode_behaviour: >
    In materialized mode, structures are streamed one at a time from
    store.store_path/structures.parquet. Extracted ChainRecord objects
    are written immediately to store.store_path/chains.parquet.
    No more than one structure is held in memory at a time.

  notes: >
    Each polymer chain produces one ChainRecord regardless of
    how many chains share the same entity. A homo-dimer with two
    chains A and B of the same sequence produces two ChainRecords.
```

---

## 7.2 Interface extraction

```yaml
interface_extraction:
  input: Dataset (structure-level), either mode
  output: InterfaceDataset, same mode as input Dataset
  policy_field: extraction_rules.interface_extraction_rules

  precondition: >
    Each structure in the Dataset must have an AnnotationLayer with
    layer_type == "interface_annotations" in derived_annotations.
    This requires that pandora.builtin.interfaces was applied in C03.
    If the annotation is absent, the structure is skipped with
    MISSING_INTERFACE_ANNOTATIONS.

  for_each_structure_in_dataset:
    read interface_annotations layer from derived_annotations
    for_each_interface_in_annotation_layer:
      determine interface_type from chain entity types
      if interface_type not in enabled interface_types:
        exclude with INTERFACE_TYPE_EXCLUDED
      if contact_count < min_contact_residues:
        exclude with INTERFACE_BELOW_MIN_CONTACTS
      if interface_area < min_interface_area (when set):
        exclude with INTERFACE_BELOW_MIN_AREA
      else:
        create InterfaceRecord:
          - build ChainRecord for each partner chain
          - embed contact residue lists and geometry
          - inherit archive_metadata from parent structure

  materialized_mode_behaviour: >
    In materialized mode, structures are streamed from store. Extracted
    InterfaceRecord objects are written immediately to interfaces.parquet.
```

---

## 7.3 Residue extraction

```yaml
residue_extraction:
  input: Dataset (structure-level) OR ChainDataset, either mode
  output: ResidueDataset, same mode as input
  policy_field: extraction_rules.residue_extraction_rules

  for_each_source_unit:
    for_each_residue_in_unit:
      determine residue_type from comp_id
      if residue_type not in include_residue_types:
        exclude with RESIDUE_TYPE_EXCLUDED
      if require_full_backbone: true AND backbone atoms incomplete:
        exclude with INCOMPLETE_BACKBONE
      else:
        create ResidueRecord:
          - embed atoms for this residue
          - if context_radius is set:
              find all residues with Cα within context_radius Å
              populate neighboring_residues as list[NeighborReference]
          - inherit archive_metadata from parent
          - filter derived_annotations to residue-level layers

  materialized_mode_behaviour: >
    In materialized mode, source units are streamed from store.
    Extracted ResidueRecord objects are written immediately to
    residues.parquet.

  notes: >
    Backbone completeness check requires N, CA, C, O atoms all present
    for amino acid residues.
```

---

# 8. Public Functions

## 8.1 `build_dataset()`

### Responsibility

Construct a curated structure-level dataset from annotated structures.
This is the main orchestrator and must always be run before extraction.

Supports both in-memory and materialized execution modes.

### Internal Workflow

```text
In-memory mode (store is null):
1. Collect all structures into a list.
2. apply_selection_rules()
3. apply_quality_filters()
4. apply_organism_filters()
5. apply_content_filters()
6. deduplicate_dataset()  (if enabled)
7. validate_dataset()
8. Assemble Dataset with mode="in_memory", structures=retained.

Materialized mode (store is a DatasetStoreRef):
1. Initialise writer for store.store_path/structures.parquet.
2. For each structure in the input iterable:
   a. apply per-structure selection, quality, organism, content filters.
   b. If structure passes, write immediately to store.
   c. Record exclusions.
3. If deduplication enabled:
   Stream record IDs from store. Compute deduplication set.
   Remove duplicates from store (rewrite without duplicates).
4. validate_dataset()
5. Assemble Dataset with mode="materialized", structures=[], store=ref.
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

## 8.2 `materialize_dataset()`

### Responsibility

Convert an existing in-memory `Dataset` (or sub-type) to materialized mode
by writing all records to a `DatasetStore` on disk. The returned dataset
has `mode="materialized"`, `structures=[]` (or `chains=[]`, etc.), and
`store` pointing to the written data.

Use this function when an in-memory dataset has grown too large for the
next pipeline stage (e.g. before calling C05 on a large dataset).

### Input Schema

```yaml
materialize_dataset:
  dataset: Dataset | ChainDataset | InterfaceDataset | ResidueDataset
  # The in-memory dataset to materialize.

  store: DatasetStoreRef
  # Must be initialised before calling. Use init_dataset_store().
  # The store's granularity must match the dataset's granularity.
```

### Output Schema

```yaml
materialize_dataset_result:
  dataset: Dataset | ChainDataset | InterfaceDataset | ResidueDataset
  # Same type as input, with mode="materialized" and records cleared.
```

---

## 8.3 `init_dataset_store()`

### Responsibility

Initialise a new `DatasetStore` directory at the given path and return a
`DatasetStoreRef`. The store directory is created if it does not exist.

### Input Schema

```yaml
init_dataset_store:
  store_path: string
  # Absolute path to the store root directory.
  granularity: string
  # structure | chain | interface | residue
  store_format: string
  # parquet (V1 only)
```

### Output Schema

```yaml
init_dataset_store_result:
  store: DatasetStoreRef
```

---

## 8.4 `extract_chain_records()`

### Responsibility

Extract individual polymer chains from a structure-level `Dataset`,
producing a `ChainDataset`. One `ChainRecord` is created per qualifying
polymer chain per structure.

The output mode matches the input dataset mode.

### Input Schema

```yaml
extract_chain_records:
  dataset: Dataset
  extraction_rules:
    include_chain_types:
      proteins: bool
      rna: bool
      dna: bool
    min_chain_length: int | null

  store: DatasetStoreRef | null
  # Required when dataset.mode == "materialized".
  # Ignored when dataset.mode == "in_memory".
```

### Output Schema

```yaml
extract_chain_records_result:
  chain_dataset: ChainDataset
```

---

## 8.5 `extract_interface_records()`

### Responsibility

Extract chain-chain interfaces from a structure-level `Dataset`, producing
an `InterfaceDataset`. Requires that `pandora.builtin.interfaces` was run
in C03 for every structure in the dataset.

The output mode matches the input dataset mode.

### Input Schema

```yaml
extract_interface_records:
  dataset: Dataset
  extraction_rules:
    interface_types:
      protein_protein: bool
      protein_rna: bool
      protein_dna: bool
      protein_ligand: bool
    min_contact_residues: int | null
    min_interface_area: float | null

  store: DatasetStoreRef | null
  # Required when dataset.mode == "materialized".
```

### Output Schema

```yaml
extract_interface_records_result:
  interface_dataset: InterfaceDataset
```

### Notes

Structures without interface annotations are recorded in
`excluded_items` with reason_code `MISSING_INTERFACE_ANNOTATIONS`
and skipped. The run does not fail unless all structures are missing
annotations.

---

## 8.6 `extract_residue_records()`

### Responsibility

Extract individual residues from a structure-level `Dataset` or a
`ChainDataset`, producing a `ResidueDataset`.

The output mode matches the input dataset mode.

### Input Schema

```yaml
extract_residue_records:
  source: Dataset | ChainDataset

  extraction_rules:
    include_residue_types:
      standard_amino_acids: bool
      non_standard: bool
      nucleotides: bool
    require_full_backbone: bool
    context_radius: float | null
    # Å. Populate neighboring_residues within this radius.
    # null = no context enrichment.

  store: DatasetStoreRef | null
  # Required when source.mode == "materialized".
```

### Output Schema

```yaml
extract_residue_records_result:
  residue_dataset: ResidueDataset
```

---

## 8.7 `filter_dataset()`

### Responsibility

Apply explicit filter operations to a structure list.

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
  parameters: object
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

## 8.8 `deduplicate_dataset()`

### Responsibility

Remove exact duplicate structures. Supports `entry_id` and `exact_hash`
strategies. Similarity-based deduplication belongs to Component 05.

### Input Schema

```yaml
deduplicate_dataset:
  structures: list[AnnotatedStructureWithPlugins]
  deduplication_rules:
    enabled: bool
    strategy: string
    # entry_id | exact_hash
```

### Output Schema

```yaml
deduplicate_dataset_result:
  deduplicated_structures: list[AnnotatedStructureWithPlugins]
  deduplication_report: DeduplicationReport
```

---

## 8.9 `validate_dataset()`

### Responsibility

Validate dataset consistency, provenance completeness, and policy
compliance.

### V1 Validation Rules

```yaml
error_rules:
  EMPTY_DATASET:
    condition: "Dataset has zero records after filtering and deduplication
                (whether in_memory or materialized)."
    result_status: invalid
  ALL_BIOMOLECULES_EXCLUDED:
    condition: "All include_biomolecules flags are false."
    result_status: invalid
  PROVENANCE_INCOMPLETE:
    condition: "A record_* flag is true but the corresponding field is null."
    result_status: invalid

warning_rules:
  HIGH_EXCLUSION_RATE:
    condition: "More than 50% of total_input structures were excluded."
    result_status: warning
  AGGRESSIVE_DEDUPLICATION:
    condition: "More than 30% of structures were removed by deduplication."
    result_status: warning
  EMPTY_SOURCE_CONTRIBUTION:
    condition: "An include_source contributed zero structures."
    result_status: warning
  NULL_RESOLUTION_EXCLUDED:
    condition: "max_resolution is set and structures with null resolution
                were excluded."
    result_status: warning
  MISSING_TAXONOMY_SKIPPED:
    condition: "Organism filtering was active but some structures had no
                taxonomy metadata."
    result_status: warning
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

## 8.10 `build_dataset_many()`

### Responsibility

Construct multiple independent datasets from multiple
`DatasetConstructionInput` jobs.

### Input Schema

```yaml
build_dataset_many:
  batches: list[DatasetConstructionInput]
  mode: string
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

# 9. Internal Helper Functions

## 9.1 `apply_selection_rules()`

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

## 9.2 `apply_quality_filters()`

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

## 9.3 `apply_content_filters()`

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
```

---

## 9.4 `apply_organism_filters()`

### Input

```yaml
apply_organism_filters:
  structures: list[AnnotatedStructureWithPlugins]
  organism_rules:
    include_taxa: list[string]
    exclude_taxa: list[string]
```

### Output

```yaml
apply_organism_filters_result:
  retained: list[AnnotatedStructureWithPlugins]
  excluded: list[ExclusionRecord]
```

---

## 9.5 `summarize_exclusions()`

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

## 9.6 `build_dataset_provenance()`

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
    source_count: int
    input_sources: list[string]
```

---

## 9.7 `extract_chain_from_structure()`

Internal helper used by `extract_chain_records()`. Extracts one
`ChainRecord` from one structure for a given chain_id.

### Input

```yaml
extract_chain_from_structure:
  structure: AnnotatedStructureWithPlugins
  chain_id: string
  policy: DatasetCurationPolicy
```

### Output

```yaml
extract_chain_from_structure_result:
  chain_record: ChainRecord | null
  # null if the chain is excluded by extraction_rules.
  exclusion: ExclusionRecord | null
```

---

## 9.8 `extract_interfaces_from_structure()`

Internal helper used by `extract_interface_records()`. Extracts all
qualifying `InterfaceRecord` objects from one structure.

### Input

```yaml
extract_interfaces_from_structure:
  structure: AnnotatedStructureWithPlugins
  extraction_rules: object
  # interface_extraction_rules sub-schema
```

### Output

```yaml
extract_interfaces_from_structure_result:
  interface_records: list[InterfaceRecord]
  excluded: list[ExclusionRecord]
```

---

## 9.9 `extract_residues_from_chain()`

Internal helper used by `extract_residue_records()`. Extracts all
qualifying `ResidueRecord` objects from one `ChainRecord`.

### Input

```yaml
extract_residues_from_chain:
  chain_record: ChainRecord
  extraction_rules: object
  # residue_extraction_rules sub-schema
```

### Output

```yaml
extract_residues_from_chain_result:
  residue_records: list[ResidueRecord]
  excluded: list[ExclusionRecord]
```

---

# 10. Dataset Units

| Granularity | Type | Unit | Primary use cases |
|-------------|------|------|-------------------|
| Structure | `Dataset` | `AnnotatedStructureWithPlugins` | Fold classification, global quality, structure-based drug design |
| Chain | `ChainDataset` | `ChainRecord` | Sequence models, secondary structure, homology detection, LM pre-training |
| Interface | `InterfaceDataset` | `InterfaceRecord` | PPI prediction, binding affinity, antibody-antigen, docking benchmarks |
| Residue | `ResidueDataset` | `ResidueRecord` | Binding site prediction, mutation effects, contact maps, solvent exposure |

Extraction hierarchy:

```text
Dataset
  └── extract_chain_records()    → ChainDataset
        └── extract_residue_records() → ResidueDataset
  └── extract_interface_records() → InterfaceDataset
        └── (residues accessible via chain_record_1/chain_record_2)
  └── extract_residue_records()  → ResidueDataset (direct from structures)
```

---

# 11. Policy Schema

## 11.1 DatasetCurationPolicy

```yaml
DatasetCurationPolicy:
  policy_id: string
  policy_name: string
  policy_version: string
  description: string

  execution_mode: string
  # in_memory    — default; accumulate structures in memory.
  # materialized — stream structures to disk; requires store in input.
  # If not specified, defaults to in_memory when no store is provided in
  # DatasetConstructionInput, and materialized when a store is provided.

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

  quality_rules:
    max_resolution: float | null
    null_resolution_behavior: string
    # exclude (default) | include
    min_chain_length: int | null
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
    exclude_taxa: list[string]

  deduplication_rules:
    enabled: bool
    strategy: string
    # entry_id | exact_hash

  extraction_rules:
    granularity: string
    # structure  — no extraction; Dataset is the final output (default).
    # chain      — run extract_chain_records() after build_dataset().
    # interface  — run extract_interface_records() after build_dataset().
    # residue    — run extract_residue_records() after build_dataset().
    #              Set source_granularity to "structure" or "chain".

    chain_extraction_rules:
      include_chain_types:
        proteins: bool
        rna: bool
        dna: bool
      min_chain_length: int | null
      # Overrides quality_rules.min_chain_length at chain level when set.

    interface_extraction_rules:
      interface_types:
        protein_protein: bool
        protein_rna: bool
        protein_dna: bool
        protein_ligand: bool
      min_contact_residues: int | null
      min_interface_area: float | null
      # Å². null = no minimum area filter.

    residue_extraction_rules:
      source_granularity: string
      # structure — extract directly from structures in Dataset.
      # chain     — extract from ChainDataset (requires chain extraction first).
      include_residue_types:
        standard_amino_acids: bool
        non_standard: bool
        nucleotides: bool
      require_full_backbone: bool
      # If true, exclude residues missing any backbone atom (N, CA, C, O).
      context_radius: float | null
      # Å. Populate ResidueRecord.neighboring_residues within this radius.

  provenance_rules:
    record_filters: bool
    record_exclusions: bool
    record_dataset_version: bool
    record_input_sources: bool
```

---

# 12. Non-Responsibilities

```yaml
not_responsible_for:
  - train_test_splitting
  - leakage_safe_partitioning
  - similarity_clustering
  - similarity_based_deduplication   # belongs to Component 05
  - benchmark_generation
  - model_training
  - graph_generation
  - embeddings
  - task_specific_labels
  - interface_geometry_computation   # computed by pandora.builtin.interfaces in C03
```

---

# 13. Component Definition

The Dataset Construction & Curation Layer transforms collections of
annotated structures into reproducible curated datasets using explicit
curation policies and provenance-aware filtering operations. It supports
four granularity levels — structure, chain, interface, and residue — each
producing a self-contained typed dataset in either in-memory or materialized
mode. Materialized mode enables processing of large-scale datasets (>10K
structures) by streaming records to disk one at a time, avoiding full
in-memory accumulation.
