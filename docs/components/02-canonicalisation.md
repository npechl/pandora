# Component 02 — Canonical Structure Object Layer

## Purpose

The Canonical Structure Object Layer converts parsed structural records into
standardized and provenance-preserving canonical structure representations.

This component is responsible for:

* structural normalization,
* identifier normalization,
* alternate conformation handling,
* missing data handling,
* assembly normalization,
* and canonical structural consistency.

This layer uses explicit canonicalisation policies so dataset preparation becomes:

* reproducible,
* transparent,
* and configurable.

The output of this component becomes the canonical structure representation used by all downstream Pandora components.

---

# 1. Input Schemas

## 1.1 canonicalisation input

```yaml
canonicalisationInput:
  ingestion_result: MmCIFIngestionResult
  # Must have status "success" or "warning".
  # Ingestion results with status "failed" are rejected with a
  # configuration error before canonicalisation begins.

  policy: canonicalisationPolicy
  # See Section 4 for the full policy schema and strategy definitions.
```

---

## 1.2 Batch canonicalisation input

```yaml
canonicalisationBatchInput:
  ingestion_results:
    - MmCIFIngestionResult

  policy: canonicalisationPolicy

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

# 2. Canonical Structure Object Schemas

These schemas define the typed records produced by `canonicalise_structure()`.
`CanonicalStructure` reuses the `Atom`, `Residue`, `Chain`, `Entity`,
`Assembly`, and `Ligand` types defined in Component 01, but all identifier
fields (`chain_id`, `seq_id`, `assembly_id`, `entity_id`) reflect
**canonical values** after policy application — not the original archive values.
The original values are preserved exclusively in `canonical_mappings`.

## 2.1 CanonicalStructure

```yaml
CanonicalStructure:
  atoms: list[Atom]
  # Atom.chain_id and Atom.residue_id reflect canonical identifiers.
  # Atom.altloc is null for all atoms when altloc resolution was applied.

  residues: list[Residue]
  # Residue.seq_id reflects the canonical sequence number.
  # Residue.auth_seq_id and Residue.insertion_code preserve the original
  # archive values for reference (not for use as canonical keys).

  chains: list[Chain]
  # Chain.chain_id is the canonical chain identifier.
  # Chain.auth_chain_id preserves the original author chain ID.

  entities: list[Entity]
  # Entity.entity_id reflects the canonical entity identifier.

  assemblies: list[Assembly]
  # Assembly.assembly_id reflects the canonical assembly identifier.

  ligands: list[Ligand]
  # Populated according to ligand_rules in the applied policy.
  # May be empty if ligand_rules.strategy == "annotate_only".
```

---

## 2.2 ChainIdMapping

Records the correspondence between canonical and original chain identifiers
for every chain in the structure.

```yaml
ChainIdMapping:
  items:
    - canonical_chain_id: string
      # The chain_id used in CanonicalStructure.

      original_chain_id: string
      # The original label_asym_id from the archive.

      original_auth_chain_id: string
      # The original auth_asym_id from the archive.
```

---

## 2.3 ResidueNumberMapping

Records the correspondence between canonical and original residue numbers
for every residue in every chain.

```yaml
ResidueNumberMapping:
  items:
    - canonical_chain_id: string

      canonical_seq_id: int
      # The seq_id used in CanonicalStructure.Residue.

      original_chain_id: string

      original_seq_id: int | null
      # The original label_seq_id. null for non-polymer residues.

      original_auth_seq_id: string
      # The original auth_seq_id (may be non-numeric).

      original_insertion_code: string | null
      # The original PDB insertion code. null when absent.
```

---

## 2.4 AssemblyMapping

Records the correspondence between canonical and original assembly identifiers.

```yaml
AssemblyMapping:
  items:
    - canonical_assembly_id: string
      original_assembly_id: string
```

---

## 2.5 EntityMapping

Records the correspondence between canonical and original entity identifiers.
A canonical entity may correspond to multiple original entities when
`entity_rules.strategy == "merge_equivalent_entities"`.

```yaml
EntityMapping:
  items:
    - canonical_entity_id: string

      original_entity_ids: list[string]
      # Single-element list in all strategies except "merge_equivalent_entities".
```

---

## 2.6 AltlocSelectionMapping

Records which alternate conformation was selected for each residue where
altloc resolution was applied.

```yaml
AltlocSelectionMapping:
  items:
    - canonical_chain_id: string

      residue_id: string
      # The canonical residue_id of the affected residue.

      selected_altloc: string
      # The altloc label that was retained (e.g. "A", "B").

      available_altlocs: list[string]
      # All altloc labels that were present before resolution.

      selection_reason: string
      # best_occupancy    — highest occupancy altloc was selected.
      # first_alphabetical — lowest alphabetical label was selected.
      # last_alphabetical  — highest alphabetical label was selected.
      # user_defined       — user-specified altloc label was selected.
```

---

# 3. Output Schemas

## 3.1 Canonical structure result

```yaml
CanonicalStructureResult:
  entry_id: string

  status: string
  # success  — canonicalised and validated without issues.
  # warning  — canonicalised; non-fatal issues detected.
  # failed   — canonicalisation failed; canonical_structure is null.

  canonical_structure: CanonicalStructure | null
  # null when status == "failed".

  canonical_mappings:
    chain_id_mapping: ChainIdMapping
    residue_number_mapping: ResidueNumberMapping
    assembly_mapping: AssemblyMapping
    entity_mapping: EntityMapping
    altloc_selection_mapping: AltlocSelectionMapping
  # canonical_mappings is null when status == "failed".

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]

  provenance:
    provider: string
    source_uri: string | null
    retrieved_at: string | null
    # ISO 8601 timestamp. Propagated from the ingestion result.
    canonicalised_at: string | null
    # ISO 8601 timestamp of when canonicalisation was applied.
```

---

## 3.2 Batch canonicalisation result

```yaml
canonicalisationBatchResult:
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

      canonical_structure: CanonicalStructure | null
      # null when status == "failed".

      canonical_mappings:
        chain_id_mapping: ChainIdMapping | null
        residue_number_mapping: ResidueNumberMapping | null
        assembly_mapping: AssemblyMapping | null
        entity_mapping: EntityMapping | null
        altloc_selection_mapping: AltlocSelectionMapping | null
      # All mapping fields are null when status == "failed".

      applied_policy:
        policy_id: string
        policy_name: string
        policy_version: string

      diagnostics:
        warnings: list[Diagnostic]
        errors: list[Diagnostic]

      provenance:
        retrieved_at: string | null
        canonicalised_at: string | null
```

---

# 4. Policy Schema

## 4.1 canonicalisationPolicy

```yaml
canonicalisationPolicy:
  policy_id: string
  policy_name: string
  policy_version: string
  description: string

  identifier_rules:
    chain_id:
      strategy: string
      # preserve           — Keep original label_asym_id as canonical chain_id.
      # remap              — Reassign chain IDs to sequential single-character
      #                      labels (A, B, C, ... Z, AA, AB, ...).
      # use_auth_chain_id  — Use auth_asym_id as the canonical chain_id.

    residue_numbering:
      strategy: string
      # preserve      — Keep original label_seq_id as canonical seq_id.
      # use_auth_seq  — Use auth_seq_id (integer part only, insertion codes
      #                 stripped) as the canonical seq_id.
      # renumber      — Renumber residues sequentially from 1 within each
      #                 chain, ignoring gaps in the original numbering.

      preserve_insertion_codes: bool
      # If true, insertion codes are retained in Residue.insertion_code.
      # If false, insertion codes are stripped (collapsed into seq_id).
      # Only relevant when strategy is "preserve" or "use_auth_seq".

    assembly_id:
      strategy: string
      # preserve      — Keep original _pdbx_struct_assembly.id values.
      # remap         — Reassign assembly IDs to sequential integers (1, 2, ...).
      # standardize   — Use the preferred biological assembly as assembly "1";
      #                 renumber remaining assemblies sequentially.

  missing_data_rules:
    missing_atoms:
      strategy: string
      # preserve             — Keep residues with missing atoms as-is.
      # annotate             — Flag residues with missing atoms in diagnostics
      #                        but retain them unchanged.
      # drop_partial_residue — Remove any residue where one or more non-hydrogen
      #                        atoms are absent.
      # impute               — Attempt to fill missing atom coordinates using
      #                        ideal geometry. Requires allow_imputation: true.

      allow_imputation: bool
      # Must be true for strategy "impute". Ignored otherwise.

      record_missingness: bool
      # If true, emit a MISSING_ATOMS diagnostic for every affected residue.

    missing_residues:
      strategy: string
      # preserve           — Keep sequence gaps as-is.
      # annotate           — Flag gaps in diagnostics but retain the chain.
      # drop_chain_segment — Remove chain segments that contain sequence gaps.
      # impute             — Attempt to model missing residues using sequence
      #                      and ideal geometry. Requires allow_imputation: true.

      record_gaps: bool
      # If true, emit a SEQUENCE_GAP diagnostic for every detected gap.

    incomplete_chains:
      strategy: string
      # preserve                    — Keep chains regardless of completeness.
      # exclude                     — Remove chains with missing terminal or
      #                               internal segments.
      # truncate_to_complete_regions — Retain only contiguous complete segments;
      #                               discard terminal incomplete regions.

  altloc_rules:
    strategy: string
    # preserve              — Retain all altloc variants. All atoms with an
    #                         altloc label are kept. AltlocSelectionMapping
    #                         will be empty.
    # select_best_occupancy — Select the altloc with the highest occupancy per
    #                         residue. Apply tie_breaker on equal occupancies.
    # select_first          — Select the altloc with the lowest alphabetical
    #                         label (A before B).
    # select_user_defined   — Select the altloc label specified in
    #                         user_defined_altloc.

    tie_breaker: string
    # Used when strategy == "select_best_occupancy" and two altlocs have
    # equal occupancy.
    #
    # alphabetical_first — Select the altloc with the lowest alphabetical label.
    # alphabetical_last  — Select the altloc with the highest alphabetical label.
    # lowest_b_factor    — Select the altloc with the lowest average B-factor.
    # highest_b_factor   — Select the altloc with the highest average B-factor.

    user_defined_altloc: string | null
    # The specific altloc label to select when strategy == "select_user_defined".
    # e.g. "A", "B". null for all other strategies.

    record_selection: bool
    # If true, populate AltlocSelectionMapping for every residue where
    # resolution was applied.

  assembly_rules:
    strategy: string
    # preserve_as_reported          — Use assemblies exactly as reported in
    #                                 the mmCIF file.
    # standardize_biological_assembly — Prefer the assembly annotated as the
    #                                 preferred biological unit per
    #                                 preferred_assembly_source. Other
    #                                 assemblies are retained but renumbered.
    # select_first_assembly         — Retain only the first listed assembly
    #                                 (assembly_id "1" in the source file).
    #                                 All others are discarded.

    preferred_assembly_source: string
    # Used when strategy == "standardize_biological_assembly".
    #
    # author — Use the assembly flagged by the depositing author as preferred.
    # pdbe   — Use the PDBe-annotated preferred biological assembly.
    # pdbx   — Use the pdbx_struct_assembly preferred flag from the mmCIF file.
    # first  — Fall back to the first listed assembly if no preferred flag exists.

    record_original_assembly_mapping: bool
    # If true, populate AssemblyMapping for every assembly.

  entity_rules:
    strategy: string
    # preserve                  — Keep original entity IDs unchanged.
    # standardize               — Renumber entity IDs sequentially (1, 2, ...).
    # merge_equivalent_entities — Merge entities with identical canonical
    #                             sequences or compositions into a single entity.
    #                             EntityMapping will reflect the merge.

    preserve_original_entity_ids: bool
    # If true, the original entity_id is retained in EntityMapping even when
    # strategy == "standardize".

  ligand_rules:
    strategy: string
    # preserve      — Include all non-polymer, non-water entities as Ligand records.
    # filter        — Apply keep_waters / keep_ions / keep_nonpolymer_ligands to
    #                 determine which ligands to include.
    # annotate_only — Ligands are recorded in diagnostics but excluded from
    #                 CanonicalStructure.ligands.

    keep_waters: bool
    keep_ions: bool
    keep_nonpolymer_ligands: bool
    # Active when strategy == "filter". Ignored otherwise.

  validation_rules:
    strictness: string
    # strict     — All warning_rules are treated as errors.
    # moderate   — Default behaviour: errors are errors, warnings are warnings.
    # permissive — Error rules are downgraded to warnings where possible.

    fail_on_unresolved_issues: bool
    # If true, any unresolved error causes status == "failed".

    warnings_as_errors: bool
    # If true, any warning diagnostic causes status == "failed".
    # Equivalent to strictness == "strict" at the result level.

  provenance_rules:
    record_original_mappings: bool
    # If true, all mapping schemas are populated.
    # If false, mapping objects are empty (reduces output size).

    record_transforms: bool
    # If true, each applied transformation is logged in provenance.

    record_policy_application: bool
    # If true, the full applied policy is embedded in the result.

    emit_canonicalisation_report: bool
    # If true, a summary of all applied transformations is included
    # in the diagnostics context.
```

---

# 5. Public Functions

## 5.1 `canonicalise_structure()`

### Responsibility

Convert a parsed structure into a canonical structure according to a
canonicalisation policy.

### Internal Workflow

```text
1. Reject if ingestion_result.status == "failed".

2. Apply identifier rules (in order):
   a. normalize_chain_ids()        — if chain_id.strategy != "preserve"
   b. normalize_residue_numbering() — if residue_numbering.strategy != "preserve"
   c. normalize_assemblies()       — if assembly_id.strategy != "preserve"

3. Apply missing data rules:
   d. handle_missing_data()        — for missing_atoms, missing_residues,
                                     and incomplete_chains per policy

4. Apply altloc rules:
   e. resolve_altlocs()            — if altloc_rules.strategy != "preserve"

5. Apply entity and ligand rules:
   f. normalize_entities()         — per entity_rules
   g. filter_ligands()             — per ligand_rules

6. Validate:
   h. validate_canonical_structure()

7. If validation status == "failed" and fail_on_unresolved_issues:
      return CanonicalStructureResult(status="failed",
                                      canonical_structure=null,
                                      canonical_mappings=null, ...)

8. Assemble and return CanonicalStructureResult.
```

Steps 2a–2c must be completed before steps 3–5, as missing data handling
and altloc resolution operate on already-normalized identifiers.

### Input Schema

```yaml
canonicalise_structure:
  input: canonicalisationInput
```

### Output Schema

```yaml
canonicalise_structure_result:
  result: CanonicalStructureResult
```

### Notes

This function:

* applies canonicalisation policies,
* standardizes identifiers,
* handles missing data,
* resolves alternate conformations,
* normalizes assemblies,
* records mappings back to the original archive representation.

This function does NOT:

* add metadata,
* generate graphs,
* create embeddings,
* infer missing biological structures unless explicitly configured via
  `allow_imputation: true`.

---

## 5.2 `validate_canonical_structure()`

### Responsibility

Validate canonical structures against Pandora canonicalisation consistency
rules.

### V1 Validation Rules

```yaml
error_rules:
  CANONICAL_CHAIN_ID_COLLISION:
    condition: "Two or more chains share the same canonical chain_id after
                identifier normalization."
    result_status: failed

  RESIDUE_NUMBER_COLLISION:
    condition: "Two or more residues in the same chain share the same
                canonical seq_id after renumbering."
    result_status: failed

  UNRESOLVED_ALTLOC:
    condition: "Altloc resolution was requested but one or more residue
                groups could not be resolved (e.g. missing occupancy data
                for all altlocs)."
    result_status: failed

  MAPPING_INCOMPLETE:
    condition: "A chain or residue in the canonical structure has no
                corresponding entry in canonical_mappings, and
                record_original_mappings is true."
    result_status: failed

  EMPTY_CANONICAL_STRUCTURE:
    condition: "The canonical structure contains zero chains after
                applying all rules."
    result_status: failed

warning_rules:
  IMPUTATION_APPLIED:
    condition: "Missing atom or residue imputation was applied to one
                or more residues."
    result_status: warning

  ASSEMBLY_SOURCE_FALLBACK:
    condition: "The preferred_assembly_source was unavailable; fell back
                to the first listed assembly."
    result_status: warning

  ENTITY_MERGE_SKIPPED:
    condition: "Entity merge was attempted but skipped for one or more
                entity pairs due to sequence or composition mismatch."
    result_status: warning

  LIGANDS_EXCLUDED:
    condition: "One or more ligands were excluded by ligand_rules.strategy
                == 'annotate_only' or filter settings."
    result_status: warning

  INCOMPLETE_CHAIN_TRUNCATED:
    condition: "One or more chains were truncated due to
                incomplete_chains.strategy == 'truncate_to_complete_regions'."
    result_status: warning
```

### Status determination

```yaml
status_rules:
  failed:  Any error_rule fires.
  warning: No error_rules fire, but one or more warning_rules fire.
  valid:   No rules fire.

  Override rules:
    - If validation_rules.warnings_as_errors: true, any warning is
      promoted to an error and result_status becomes "failed".
    - If validation_rules.strictness: "permissive", error rules that
      can be downgraded produce "warning" instead of "failed" where noted.
```

### Input Schema

```yaml
validate_canonical_structure:
  canonical_structure: CanonicalStructure
  canonical_mappings: object
  policy: canonicalisationPolicy
```

### Output Schema

```yaml
validate_canonical_structure_result:
  validation_status: string
  # valid | warning | invalid

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

### Notes

Validation ensures:

* structural consistency — no collisions in canonical identifiers,
* mapping consistency — every canonical identifier has a corresponding
  original-archive entry in canonical_mappings,
* policy consistency — applied transformations are consistent with the
  policy that was declared,
* canonical representation integrity — the structure is non-empty and
  navigable after all rules are applied.

---

## 5.3 `canonicalise_many_structures()`

### Responsibility

Run canonicalisation for a list of parsed structures.

Supports sequential or parallel execution.

### Internal Workflow

```text
for each ingestion_result:
    canonicalise_structure(ingestion_result, policy)

In parallel mode:
    Entries are dispatched to up to max_workers concurrent workers.
    If fail_fast is true:
        Remaining entries are cancelled on the first failure.
    If fail_fast is false (default):
        All entries are attempted; per-entry failures are isolated.
```

### Input Schema

```yaml
canonicalise_many_structures:
  input: canonicalisationBatchInput
```

### Output Schema

```yaml
canonicalise_many_structures_result:
  result: canonicalisationBatchResult
```

### Notes

This function:

* orchestrates batch canonicalisation,
* isolates failures (unless `fail_fast` is true),
* aggregates diagnostics across all entries,
* and generates batch summaries.

---

# 6. Internal Helper Functions

## 6.1 `normalize_chain_ids()`

### Responsibility

Normalize chain identifiers according to `identifier_rules.chain_id.strategy`
and populate `ChainIdMapping`.

### Input

```yaml
normalize_chain_ids:
  parsed_structure: ParsedStructure
  strategy: string
  # preserve | remap | use_auth_chain_id
```

### Output

```yaml
normalize_chain_ids_result:
  updated_structure: ParsedStructure
  # Chain.chain_id fields reflect canonical values.
  chain_id_mapping: ChainIdMapping
```

---

## 6.2 `normalize_residue_numbering()`

### Responsibility

Normalize residue numbers according to
`identifier_rules.residue_numbering.strategy` and populate
`ResidueNumberMapping`.

### Input

```yaml
normalize_residue_numbering:
  parsed_structure: ParsedStructure
  strategy: string
  # preserve | use_auth_seq | renumber
  preserve_insertion_codes: bool
```

### Output

```yaml
normalize_residue_numbering_result:
  updated_structure: ParsedStructure
  # Residue.seq_id fields reflect canonical values.
  residue_number_mapping: ResidueNumberMapping
```

---

## 6.3 `resolve_altlocs()`

### Responsibility

Apply alternate conformation handling according to `altloc_rules` and
populate `AltlocSelectionMapping`.

### Input

```yaml
resolve_altlocs:
  parsed_structure: ParsedStructure
  strategy: string
  # preserve | select_best_occupancy | select_first | select_user_defined
  tie_breaker: string
  user_defined_altloc: string | null
  record_selection: bool
```

### Output

```yaml
resolve_altlocs_result:
  updated_structure: ParsedStructure
  # Atom.altloc is null for all atoms when strategy != "preserve".
  altloc_selection_mapping: AltlocSelectionMapping
```

---

## 6.4 `normalize_assemblies()`

### Responsibility

Normalize biological assembly representations according to
`assembly_rules` and populate `AssemblyMapping`.

### Input

```yaml
normalize_assemblies:
  parsed_structure: ParsedStructure
  strategy: string
  # preserve_as_reported | standardize_biological_assembly | select_first_assembly
  preferred_assembly_source: string
  record_original_assembly_mapping: bool
```

### Output

```yaml
normalize_assemblies_result:
  updated_structure: ParsedStructure
  # Assembly.assembly_id fields reflect canonical values.
  assembly_mapping: AssemblyMapping
  diagnostics: list[Diagnostic]
  # Includes ASSEMBLY_SOURCE_FALLBACK if preferred source was unavailable.
```

---

## 6.5 `handle_missing_data()`

### Responsibility

Apply policy-driven handling of missing atoms, missing residues, and
incomplete chains according to `missing_data_rules`.

### Input

```yaml
handle_missing_data:
  parsed_structure: ParsedStructure
  missing_atoms_strategy: string
  missing_residues_strategy: string
  incomplete_chains_strategy: string
  allow_imputation: bool
  record_missingness: bool
  record_gaps: bool
```

### Output

```yaml
handle_missing_data_result:
  updated_structure: ParsedStructure
  diagnostics: list[Diagnostic]
  # Includes MISSING_ATOMS, SEQUENCE_GAP, and IMPUTATION_APPLIED entries
  # according to record_missingness, record_gaps, and imputation settings.
```

---

# 7. Internal Architecture

## Single-structure workflow

```text
MmCIFIngestionResult  [status: success | warning]
  → normalize_chain_ids()
  → normalize_residue_numbering()
  → normalize_assemblies()
  → handle_missing_data()
  → resolve_altlocs()
  → normalize_entities()
  → filter_ligands()
  → validate_canonical_structure()
  → CanonicalStructureResult
```

---

## Batch workflow

```text
list of MmCIFIngestionResult
  → canonicalise_many_structures()
  → repeated canonicalise_structure()   [sequential or parallel]
  → canonicalisationBatchResult
```

---

# 8. Non-Responsibilities

Component 02 is not responsible for:
  - metadata integration
  - UniProt mapping
  - SIFTS integration
  - dataset filtering
  - leakage-safe splitting
  - graph generation
  - embeddings
  - ML-ready exports

---

# 9. Component Definition
The Canonical Structure Object Layer converts parsed structural records into standardized, provenance-preserving canonical structure representations using explicit canonicalisation policies.

