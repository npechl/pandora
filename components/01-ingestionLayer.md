# Component 01 — mmCIF Ingestion Layer

## Purpose

The mmCIF Ingestion Layer is the foundational component of Pandora.

Its responsibility is to:

* retrieve raw mmCIF files,
* parse structural archive data,
* validate parsed records,
* expose the structural hierarchy,
* and support both single-entry and batch ingestion workflows.

This component focuses strictly on ingestion and parsing.

It is **not** responsible for:

* canonicalization,
* normalization,
* metadata enrichment,
* dataset filtering,
* leakage-safe splitting,
* graph generation,
* embeddings,
* or ML-ready exports.

---

# 1. Input Schemas

## 1.1 Single-file ingestion input

```yaml
MmCIFIngestionInput:
  entry_id: string
  # Always required — used for provenance and cache keying even when
  # source_uri or raw_content is provided directly.

  provider: string
  # pdbe | pdb | local | custom
  # See Section 4 (Provider Contract) for routing rules.

  source_uri: string | null
  # Required when provider is "local" or "custom".
  # When provided for "pdbe" or "pdb", takes precedence over the
  # provider's default URL template.

  raw_content: string | null
  # Pre-loaded plain-text mmCIF content.
  # If provided, fetch_mmCIF() is skipped entirely.
  # Must be plain-text (not compressed). Callers are responsible
  # for decompression before passing raw_content.

  fetch_options:
    allow_partial: bool
    # If true, accept partially downloaded files.
    use_cache: bool
    # If true, use cached content when available.
    decompress: bool
    # If true, decompress gzip/bzip2 content after fetching.
    # Ignored when raw_content is provided directly.
```

### Example

```yaml
entry_id: "1abc"
provider: "pdbe"
source_uri: null
raw_content: null

fetch_options:
  allow_partial: true
  use_cache: true
  decompress: true
```

---

## 1.2 Batch ingestion input

```yaml
MmCIFBatchInput:
  entries:
    - entry_id: string
      provider: string
      source_uri: string | null
      raw_content: string | null

  mode: string
  # sequential | parallel

  fetch_options:
    allow_partial: bool
    use_cache: bool
    decompress: bool

  parallel_options:
    max_workers: int | null
    # Number of concurrent workers in parallel mode.
    # null uses the system default (typically CPU count).
    # Ignored in sequential mode.

    fail_fast: bool
    # If true, abort remaining entries on the first failure.
    # If false (default), isolate failures and continue processing
    # all remaining entries regardless of errors.
```

### Example

```yaml
entries:
  - entry_id: "1abc"
    provider: "pdbe"
    source_uri: null
    raw_content: null

  - entry_id: "2xyz"
    provider: "pdbe"
    source_uri: null
    raw_content: null

mode: "parallel"

fetch_options:
  allow_partial: true
  use_cache: true
  decompress: true

parallel_options:
  max_workers: 8
  fail_fast: false
```

---

# 2. Core Structural Object Schemas

These schemas define the typed records produced by `parse_mmCIF()`. Every
downstream component operates on these types.

## 2.1 Atom

```yaml
Atom:
  atom_id: string
  # Unique identifier within the structure.
  # Source: _atom_site.id

  atom_name: string
  # Atom name within the residue (e.g. "CA", "N", "CB").
  # Source: _atom_site.label_atom_id

  element: string
  # Element symbol (e.g. "C", "N", "O", "ZN").
  # Source: _atom_site.type_symbol

  x: float
  y: float
  z: float
  # Orthogonal Cartesian coordinates in Angstroms.
  # Source: _atom_site.Cartn_x / Cartn_y / Cartn_z

  occupancy: float
  # Fractional occupancy [0.0, 1.0].
  # Source: _atom_site.occupancy

  b_factor: float
  # Isotropic B-factor (temperature factor).
  # Source: _atom_site.B_iso_or_equiv

  altloc: string | null
  # Alternate location indicator (e.g. "A", "B").
  # null when no alternate conformation exists.
  # Source: _atom_site.label_alt_id

  residue_id: string
  # Reference to the parent Residue.residue_id.

  chain_id: string
  # Reference to the parent Chain.chain_id (label_asym_id).

  model_num: int
  # Model number for NMR ensembles and multi-model depositions.
  # Single-model structures use model_num: 1.
  # Source: _atom_site.pdbx_PDB_model_num
```

---

## 2.2 Residue

```yaml
Residue:
  residue_id: string
  # Internal identifier formed as "{chain_id}:{seq_id}:{comp_id}".
  # For non-polymer residues: "{chain_id}:null:{comp_id}:{auth_seq_id}"

  comp_id: string
  # Three-letter chemical component code (e.g. "ALA", "GLY", "ATP").
  # Source: _atom_site.label_comp_id

  seq_id: int | null
  # Sequential residue index within the polymer chain.
  # null for non-polymer residues (ligands, ions, waters).
  # Source: _atom_site.label_seq_id

  auth_seq_id: string
  # Author-assigned sequence number, preserved as string to handle
  # non-numeric values.
  # Source: _atom_site.auth_seq_id

  insertion_code: string | null
  # PDB insertion code (e.g. "A", "B").
  # null when absent.
  # Source: _atom_site.pdbx_PDB_ins_code

  chain_id: string
  # Reference to the parent Chain.chain_id.

  atoms: list[Atom]
  # All Atom records belonging to this residue,
  # including all altloc variants.

  is_polymer: bool
  # True for standard amino acids and nucleotides.
  # False for ligands, ions, and waters.
```

---

## 2.3 Chain

```yaml
Chain:
  chain_id: string
  # Label asymmetric unit identifier.
  # Source: _struct_asym.id

  auth_chain_id: string
  # Author-assigned chain identifier.
  # Source: _atom_site.auth_asym_id

  entity_id: string
  # Reference to the parent Entity.entity_id.

  chain_type: string
  # polymer | non-polymer | water | branched
  # Derived from the parent entity type.

  residues: list[Residue]
  # All Residue records in this chain (polymer residues,
  # ligands, and waters).
```

---

## 2.4 Entity

```yaml
Entity:
  entity_id: string
  # Source: _entity.id

  entity_type: string
  # polymer | non-polymer | water | branched
  # Source: _entity.type

  description: string | null
  # Source: _entity.pdbx_description

  chain_ids: list[string]
  # All Chain.chain_id values assigned to this entity.
  # Source: _struct_asym.id filtered by _struct_asym.entity_id

  sequence: string | null
  # One-letter canonical sequence for polymer entities.
  # null for non-polymer and water entities.
  # Source: _entity_poly.pdbx_seq_one_letter_code_can
```

---

## 2.5 Assembly

```yaml
Assembly:
  assembly_id: string
  # Source: _pdbx_struct_assembly.id

  details: string | null
  # Source: _pdbx_struct_assembly.details

  assembly_gen: list[AssemblyGen]
  # One or more generation instructions for this assembly.

AssemblyGen:
  asym_id_list: list[string]
  # Chain IDs (label_asym_id) included in this assembly segment.
  # Source: _pdbx_struct_assembly_gen.asym_id_list

  oper_expression: string
  # Symmetry operator expression referencing _pdbx_struct_oper_list.
  # Source: _pdbx_struct_assembly_gen.oper_expression
```

---

## 2.6 Ligand

```yaml
Ligand:
  ligand_id: string
  # Internal identifier.

  chem_comp_id: string
  # Three-letter chemical component code (e.g. "ATP", "HEM", "ZN").

  chain_id: string
  # Chain containing this ligand.

  residue_id: string
  # Reference to the Residue record for this ligand.

  is_water: bool
  # True when chem_comp_id == "HOH".

  is_ion: bool
  # True for monoatomic ions (e.g. "ZN", "CA", "MG", "FE").
  # Determined by the CCD component type field.
```

---

## 2.7 Diagnostic

All `diagnostics.warnings` and `diagnostics.errors` lists across this
component contain `Diagnostic` items with this schema.

```yaml
Diagnostic:
  code: string
  # Machine-readable identifier.
  # See Section 6.3 for defined V1 codes.
  # e.g. "MISSING_ATOM_SITE", "DISCONTINUOUS_SEQID"

  severity: string
  # warning | error

  message: string
  # Human-readable description of the issue.

  entry_id: string | null
  # The affected entry, when applicable.

  field: string | null
  # The mmCIF category or field path, when applicable.
  # e.g. "_atom_site.label_seq_id"

  context: object | null
  # Optional structured context to aid debugging.
  # e.g. { chain_id: "A", offending_value: "-1" }
```

---

## 2.8 ParsedStructure

```yaml
ParsedStructure:
  atoms: list[Atom]
  residues: list[Residue]
  chains: list[Chain]
  entities: list[Entity]
  assemblies: list[Assembly]
  ligands: list[Ligand]
```

---

# 3. Output Schemas

## 3.1 Single-file ingestion result

```yaml
MmCIFIngestionResult:
  entry_id: string

  status: string
  # success  — parsed and validated without issues.
  # warning  — parsed; non-fatal validation issues detected.
  # failed   — parse or fetch failure; parsed_structure is null.

  parsed_structure: ParsedStructure | null
  # null when status == "failed".
  # Populated (possibly with warnings) when status is "success" or "warning".

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
    # Merged diagnostics from both parse_mmCIF() and validate_mmCIF().

  provenance:
    provider: string
    source_uri: string | null
    retrieved_at: string | null
    # ISO 8601 timestamp. null when raw_content was provided directly.
    from_cache: bool
    # True when the content was served from the local cache.
    # False when raw_content was provided directly.
```

---

## 3.2 Batch ingestion result

```yaml
MmCIFBatchResult:
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

      parsed_structure: ParsedStructure | null
      # null when status == "failed".

      diagnostics:
        warnings: list[Diagnostic]
        errors: list[Diagnostic]

      provenance:
        provider: string
        source_uri: string | null
        retrieved_at: string | null
        from_cache: bool
```

---

# 4. Provider Contract

## 4.1 Defined providers

```yaml
providers:
  pdbe:
    description: "PDBe archive (EMBL-EBI)"
    default_url_template: "https://www.ebi.ac.uk/pdbe/entry-files/download/{entry_id}_updated.cif"
    requires_source_uri: false

  pdb:
    description: "RCSB PDB archive"
    default_url_template: "https://files.rcsb.org/download/{entry_id}.cif"
    requires_source_uri: false

  local:
    description: "Local filesystem path"
    requires_source_uri: true
    source_uri_format: "Absolute or relative local file path"

  custom:
    description: "Arbitrary remote URI"
    requires_source_uri: true
    source_uri_format: "Any valid HTTP/HTTPS/FTP URI"
```

---

## 4.2 Routing rules

```yaml
routing_rules:
  precedence:
    - If raw_content is not null: skip fetch entirely.
    - Else if source_uri is not null: fetch from source_uri regardless of provider.
    - Else: fetch using the provider's default_url_template with entry_id substituted.

  validation:
    - provider must be one of: pdbe | pdb | local | custom.
      Unknown provider values raise a configuration error.

    - If provider is "local" or "custom" and source_uri is null,
      fetch_mmCIF() raises a configuration error.

    - entry_id is always required for provenance and cache keying,
      even when source_uri or raw_content is used.
```

---

# 5. Caching Contract

## 5.1 Cache key

```yaml
cache_key_rules:
  - When source_uri is null (provider-based retrieval):
      cache_key: "{provider}:{entry_id}"

  - When source_uri is provided:
      cache_key: "{source_uri}"
```

Cache keys are computed before fetching. Content is stored after decompression
(i.e. the cache always holds plain-text mmCIF).

---

## 5.2 Cache options

The following fields extend `fetch_options` when fine-grained cache control
is required:

```yaml
cache_options:
  use_cache: bool
  # If false, always fetch from source and bypass the cache.
  # If true, serve cached content when available and not stale.

  max_age_seconds: int | null
  # Maximum age (in seconds) for a cached entry to be considered fresh.
  # null means no expiry: cached entries are always served if present.

  stale_behavior: string
  # use_stale  — serve the stale cached entry silently (default when
  #              max_age_seconds is null).
  # warn       — serve the stale entry and emit a STALE_CACHE diagnostic
  #              warning in the result.
  # fail       — raise an error if the cached entry exceeds max_age_seconds.
```

---

## 5.3 Cache storage contract

* The cache stores decompressed plain-text mmCIF only.
* Decompression occurs before writing to cache.
* Cache entries are keyed by the rules in Section 5.1.
* Cache invalidation in V1 is time-based only (max_age_seconds).
  Archive-release-aware invalidation is reserved for a future version.

---

# 6. Public Functions

## 6.1 `fetch_mmCIF()`

### Responsibility

Retrieve raw mmCIF content from PDBe, PDB, local filesystem, or custom URI,
following the routing and caching rules in Sections 4 and 5.

### Input Schema

```yaml
fetch_mmCIF:
  entry_id: string
  provider: string
  # pdbe | pdb | local | custom

  source_uri: string | null

  fetch_options:
    use_cache: bool
    allow_partial: bool
    decompress: bool
    max_age_seconds: int | null
    stale_behavior: string
```

### Output Schema

```yaml
fetch_mmCIF_result:
  raw_content: string
  # Plain-text mmCIF content (decompressed if decompress: true).

  provenance:
    provider: string
    source_uri: string | null
    retrieved_at: string
    # ISO 8601 timestamp of retrieval or cache read.
    from_cache: bool
    # True if the result was served from cache.
```

### Notes

This function only retrieves raw data.

It does NOT:

* parse structures,
* validate records,
* normalize data,
* or modify structural content.

---

# 6.2 `parse_mmCIF()`

### Responsibility

Parse raw mmCIF text into a `ParsedStructure` record.

### Parsed mmCIF categories (V1)

```yaml
required_categories:
  _atom_site:               # Atomic coordinates and per-atom identifiers.
  _entity:                  # Entity type and description.
  _struct_asym:             # Chain-to-entity assignment.

conditional_categories:
  _entity_poly:             # Polymer sequences. Present for polymer entities.
  _pdbx_struct_assembly:    # Biological assembly definitions.
  _pdbx_struct_assembly_gen: # Assembly generation instructions.
  _pdbx_entity_nonpoly:     # Non-polymer entity details. Present when applicable.

optional_categories:
  _struct_conn:             # Covalent bonds, disulfides, metal coordination.
  _pdbx_poly_seq_scheme:    # Author-to-label residue number mapping.
```

Categories not listed above are ignored in V1.

### Input Schema

```yaml
parse_mmCIF:
  raw_content: string
```

### Output Schema

```yaml
parse_mmCIF_result:
  parsed_structure: ParsedStructure | null
  # null when parse_status == "failed".

  parse_status: string
  # success | warning | failed

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

### Notes

This function:

* parses required and available mmCIF categories,
* builds the structural hierarchy defined in Section 2,
* exposes a navigable structural record.

This function does NOT:

* canonicalize structures,
* rename chains,
* normalize residue numbering,
* resolve alternate conformations,
* repair missing residues.

---

# 6.3 `validate_mmCIF()`

### Responsibility

Validate parsed structural records for consistency and completeness.

### V1 Validation Rules

```yaml
error_rules:
  MISSING_ATOM_SITE:
    condition: "No _atom_site records are present."
    result_status: failed

  MISSING_ENTITY:
    condition: "No _entity records are present."
    result_status: failed

  MISSING_STRUCT_ASYM:
    condition: "No _struct_asym records are present."
    result_status: failed

  ZERO_ATOMS:
    condition: "Total atom count across all chains is zero."
    result_status: failed

  UNRESOLVED_ENTITY_REF:
    condition: "A chain's entity_id references an entity not present in _entity."
    result_status: failed

  UNRESOLVED_ASSEMBLY_CHAIN_REF:
    condition: "An assembly generation record references a chain_id not present
                in _struct_asym."
    result_status: failed

warning_rules:
  EMPTY_CHAIN:
    condition: "A chain has zero residue records."
    result_status: warning

  MISSING_ASSEMBLY:
    condition: "No _pdbx_struct_assembly records are present."
    result_status: warning

  MISSING_ENTITY_POLY:
    condition: "A polymer entity is present but has no _entity_poly record."
    result_status: warning

  ALTLOC_WITHOUT_OCCUPANCY:
    condition: "Alternate conformation atoms are present but occupancy is
                unset (0.0) for one or more altloc groups."
    result_status: warning

  DISCONTINUOUS_SEQID:
    condition: "label_seq_id values are non-monotonic within a polymer chain."
    result_status: warning

  STALE_CACHE:
    condition: "Cached content exceeds max_age_seconds and stale_behavior
                is 'warn'."
    result_status: warning
```

### Status determination

```yaml
status_rules:
  failed:  Any error_rule fires, or the upstream parse_status is "failed".
  warning: No error_rules fire, but one or more warning_rules fire.
  valid:   No rules fire.
```

### Input Schema

```yaml
validate_mmCIF:
  parsed_structure: ParsedStructure
```

### Output Schema

```yaml
validate_mmCIF_result:
  validation_status: string
  # valid | warning | invalid

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

### Notes

Validation reports problems but does NOT repair them.

---

# 6.4 `ingest_mmCIF()`

### Responsibility

Run the complete ingestion workflow for a single entry.

### Internal Workflow

```text
if input.raw_content is not null:
    raw = input.raw_content
    # fetch_mmCIF() is skipped.
    # provenance.retrieved_at will be null.
    # provenance.from_cache will be false.
else:
    raw = fetch_mmCIF(entry_id, provider, source_uri, fetch_options)

parse_result = parse_mmCIF(raw)

if parse_result.parse_status == "failed":
    return MmCIFIngestionResult(
        status="failed",
        parsed_structure=null,
        diagnostics=parse_result.diagnostics,
        provenance=...
    )

validation_result = validate_mmCIF(parse_result.parsed_structure)

final_status = resolve_status(parse_result.parse_status,
                              validation_result.validation_status):
    "failed"  if either is "failed"
    "warning" if either is "warning" and neither is "failed"
    "success" otherwise

return MmCIFIngestionResult(
    status=final_status,
    parsed_structure=parse_result.parsed_structure,  # null if failed
    diagnostics=merge(parse_result.diagnostics,
                      validation_result.diagnostics),
    provenance=...
)
```

### Input Schema

```yaml
ingest_mmCIF:
  entry_id: string
  provider: string
  source_uri: string | null
  raw_content: string | null

  fetch_options:
    allow_partial: bool
    use_cache: bool
    decompress: bool
    max_age_seconds: int | null
    stale_behavior: string
```

### Output Schema

```yaml
ingest_mmCIF_result:
  result: MmCIFIngestionResult
```

### Notes

* If `raw_content` is provided, `fetch_mmCIF()` is skipped entirely.
* Diagnostics from both `parse_mmCIF()` and `validate_mmCIF()` are merged.
* The final `status` reflects the worst outcome across both steps.
* `parsed_structure` is null when `status == "failed"`.

---

# 6.5 `ingest_list_mmCIF()`

### Responsibility

Run the ingestion workflow for a list of entries.

Supports:

* sequential execution,
* or parallel execution.

### Internal Workflow

```text
for each entry in entries:
    ingest_mmCIF(entry)

In parallel mode:
    Entries are dispatched to up to max_workers concurrent workers.
    If fail_fast is true:
        Remaining in-flight entries are cancelled on the first failure.
        Partial results up to the failure point are returned.
    If fail_fast is false (default):
        All entries are attempted regardless of individual failures.
        Per-entry failures are recorded in the result with status "failed".
```

### Input Schema

```yaml
ingest_list_mmCIF:
  entries: list[MmCIFIngestionInput]

  mode: string
  # sequential | parallel

  fetch_options:
    allow_partial: bool
    use_cache: bool
    decompress: bool
    max_age_seconds: int | null
    stale_behavior: string

  parallel_options:
    max_workers: int | null
    fail_fast: bool
```

### Output Schema

```yaml
ingest_list_mmCIF_result:
  result: MmCIFBatchResult
```

### Notes

This function:

* orchestrates batch ingestion,
* isolates per-entry failures (unless `fail_fast` is true),
* aggregates diagnostics across all entries,
* and generates batch summaries.

---

# 7. Internal Helper Functions

## 7.1 `extract_structure_hierarchy()`

### Responsibility

Internal helper used by `parse_mmCIF()`.

Builds structural relationships from the flat `_atom_site` table and
related mmCIF categories:

* atom → residue,
* residue → chain,
* chain → entity,
* entity → assembly.

### Notes

This is an internal parsing helper and is NOT exposed as a public component API.

---

# 8. Internal Architecture

## Single-entry workflow

```text
raw_content provided:
    raw_content (plain-text)
      → parse_mmCIF()
      → validate_mmCIF()
      → MmCIFIngestionResult

raw_content not provided:
    entry_id + provider + source_uri
      → fetch_mmCIF()          [applies provider routing + caching]
      → parse_mmCIF()
      → validate_mmCIF()
      → MmCIFIngestionResult
```

---

## Batch workflow

```text
list of MmCIFIngestionInput
  → ingest_list_mmCIF()
  → repeated ingest_mmCIF()   [sequential or parallel]
  → MmCIFBatchResult
```

---

# 9. Non-Responsibilities

Component #1 is not responsible for:
  - canonicalization
  - normalization
  - metadata integration
  - residue renumbering
  - chain-ID standardization
  - assembly correction
  - dataset filtering
  - leakage-safe splitting
  - graph generation
  - embeddings
  - ML-ready exports

---

# 10. Component Definition

The mmCIF Ingestion Layer converts raw PDBe/PDB mmCIF files into parsed structural records with diagnostics and provenance.
