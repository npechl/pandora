# Policies

## Canonicalisation

`canonicalise_structure(structure, policy)` never guesses. Every normalization decision, how chains are named, how altlocs are resolved, whether ligands are kept, comes from a `canonicalisationPolicy` (`pandora.schemas.canonicalisation.canonicalisationPolicy`). Unset rules fall back to "preserve as reported", so an empty policy is a no-op pass that only records provenance.

```python
from pandora.schemas.canonicalisation import canonicalisationPolicy

policy = canonicalisationPolicy(
    policy_id="my-policy",
    policy_name="My Policy",
    policy_version="1.0.0",
    description="",  # optional, free text
)
```

### `identifier_rules`

**`chain_id.strategy`** (default `preserve`)

| Value | Description |
|---|---|
| `preserve` | Keep the original `label_asym_id` as the canonical chain ID. |
| `remap` | Reassign chain IDs to sequential labels (`A`, `B`, ... `Z`, `AA`, `AB`, ...). |
| `use_auth_chain_id` | Use `auth_asym_id` as the canonical chain ID (falls back to `label_asym_id` if no author ID is set). |

**`residue_numbering.strategy`** (default `preserve`)

| Value | Description |
|---|---|
| `preserve` | Keep `label_seq_id` as the canonical `seq_id`. |
| `use_auth_seq` | Use the numeric part of `auth_seq_id` as the canonical `seq_id`. |
| `renumber` | Renumber residues sequentially from 1 per chain, in order of first appearance, ignoring gaps in the original numbering. |

`residue_numbering.preserve_insertion_codes` (default `True`) only applies
when `strategy` is `preserve` or `use_auth_seq`: `True` keeps insertion
codes, `False` strips them. It has **no effect under `renumber`** —
renumbering always drops insertion codes.

**`assembly_id.strategy`** (default `preserve`)

| Value | Description |
|---|---|
| `preserve` | Keep original assembly IDs. |
| `remap` | Reassign assembly IDs to sequential integers (`1`, `2`, ...) in existing order. |
| `standardize` | Intended to put the preferred biological assembly first, then renumber the rest. |

*`remap` and `standardize` currently behave identically — both just renumber sequentially in file order. The "preferred assembly goes first" behavior implied by `standardize` isn't implemented yet.*

### `missing_data_rules`

**`missing_atoms.strategy`** (default `annotate`) — a residue is "missing atoms" if any backbone `N`/`CA`/`C`/`O` atom is absent.

| Value | Description |
|---|---|
| `preserve` | Keep residues with missing backbone atoms unchanged, no diagnostic. |
| `annotate` | Keep residues unchanged; emit a `MISSING_ATOMS` diagnostic per affected residue (if `record_missingness`). |
| `drop_partial_residue` | Remove residues missing any backbone atom (also diagnosed like `annotate`). |
| `impute` | Intended to fill missing atom coordinates from ideal geometry, gated by `allow_imputation`. *Not implemented yet — currently identical to `preserve`.* |

**`missing_residues.strategy`** (default `annotate`) — a "gap" is a jump of more than 1 in a chain's `seq_id` sequence.

| Value | Description |
|---|---|
| `preserve` | Keep sequence gaps as-is; emits `SEQUENCE_GAP` diagnostics only if `record_gaps`. |
| `annotate` | Same effect as `preserve` — gaps are flagged, never modified. |
| `drop_chain_segment` | Remove every atom belonging to a chain that contains any gap. |
| `impute` | Intended to model missing residues from sequence + ideal geometry, gated by `allow_imputation`. *Not implemented yet — no `impute` branch exists; falls through with the same behavior as `preserve`/`annotate`.* |

**`incomplete_chains.strategy`** (default `preserve`) — "incomplete" means the chain's residue range has gaps between its first and last residue.

| Value | Description |
|---|---|
| `preserve` | Keep chains regardless of completeness. |
| `exclude` | Drop every asym unit/atom belonging to an incomplete chain. |
| `truncate_to_complete_regions` | Keep only the single longest contiguous run of residues per incomplete chain. |

Extra flags:

| Field | Description | Default |
|---|---|---|
| `missing_atoms.allow_imputation` | Precondition for `strategy: impute`. Currently unused since `impute` isn't implemented. | `False` |
| `missing_atoms.record_missingness` | Emit a `MISSING_ATOMS` diagnostic for every affected residue. | `True` |
| `missing_residues.record_gaps` | Emit a `SEQUENCE_GAP` diagnostic for every detected gap. | `True` |

### `altloc_rules`

**`strategy`** (default `select_best_occupancy`)

| Value | Description |
|---|---|
| `preserve` | Keep every altloc variant; no atoms are dropped. |
| `select_best_occupancy` | Per residue, pick the altloc group with the highest mean `occupancy`; ties broken by `tie_breaker`. |
| `select_first` | Pick the alphabetically-first altloc label. |
| `select_user_defined` | Pick `user_defined_altloc` if it's one of the residue's available altlocs, otherwise fall back to the alphabetically-first one. |

**`tie_breaker`** (default `alphabetical_first`) — only consulted when `strategy: select_best_occupancy` has a tie.

| Value | Description |
|---|---|
| `alphabetical_first` | Lowest altloc label. |
| `alphabetical_last` | Highest altloc label. |
| `lowest_b_factor` | Lowest mean `B_iso_or_equiv` across the tied groups. |
| `highest_b_factor` | Highest mean `B_iso_or_equiv` across the tied groups. |

`user_defined_altloc` (default `None`) — the altloc code to select under `strategy: select_user_defined`.

`record_selection` (default `True`) — populate `AltlocSelectionMapping` for every residue with more than one altloc, regardless of strategy. This flag alone controls altloc mapping; it is independent of `provenance_rules.record_original_mappings` (which governs the other four mapping types).

### `assembly_rules`

**`strategy`** (default `preserve_as_reported`)

| Value | Description |
|---|---|
| `preserve_as_reported` | Keep assemblies as reported (aside from ID renumbering under `identifier_rules.assembly_id`). |
| `standardize_biological_assembly` | Intended to prefer the assembly annotated as biological, per `preferred_assembly_source`. *Not implemented yet — currently identical to `preserve_as_reported`.* |
| `select_first_assembly` | Keep only the first listed assembly, discard the rest. |

**`preferred_assembly_source`** (default `author`) — *currently unused: no strategy reads this field yet.*

| Value | Description |
|---|---|
| `author` | Use the assembly flagged by the depositing author as preferred. |
| `pdbe` | Use the PDBe-annotated preferred biological assembly. |
| `pdbx` | Use the `pdbx_struct_assembly` preferred flag from the mmCIF file. |
| `first` | Fall back to the first listed assembly if no preferred flag exists. |

`record_original_assembly_mapping` (default `True`) — populate `AssemblyMapping` for every assembly.

### `entity_rules`

**`strategy`** (default `preserve`)

| Value | Description |
|---|---|
| `preserve` | Keep entities and their original IDs. |
| `standardize` | Renumber entity IDs sequentially (`1`, `2`, ...). |
| `merge_equivalent_entities` | Merge entities that share an identical canonical polymer sequence (`pdbx_seq_one_letter_code_can`) into one entity ID. Non-polymer entities (no canonical sequence) are never merged with each other, even if their composition matches — only sequence-based merging for polymers is implemented. |

`preserve_original_entity_ids` (default `True`) — under `strategy: standardize`, controls whether `EntityMapping` records the pre-renumber ID (`True`) or just the new ID again (`False`).

### `ligand_rules`

**`strategy`** (default `preserve`)

| Value | Description |
|---|---|
| `preserve` | Include every non-polymer entity as-is, **including waters** — see warning below. |
| `filter` | Apply `keep_waters` / `keep_ions` / `keep_nonpolymer_ligands` to decide what to drop. |
| `annotate_only` | Record ligands in diagnostics but exclude them from the canonical structure. |

`keep_waters` / `keep_ions` / `keep_nonpolymer_ligands` (default `True`) — only consulted under `strategy: filter`.

*Ion detection is a keyword heuristic, not a chemical-component lookup: an entity counts as an ion only if its `pdbx_description` contains one of a fixed set of keywords (`ION`, `ZINC`, `CALCIUM`, `SODIUM`, `CHLORIDE`, ...). An ion with a description that doesn't match falls into `keep_nonpolymer_ligands` instead of `keep_ions`.*

!!! warning "`preserve` keeps waters too"
    `strategy: preserve` returns every atom unchanged, including waters —
    `keep_waters`/`keep_ions`/`keep_nonpolymer_ligands` are only consulted
    under `strategy: filter`. If you want waters excluded, use
    `strategy: filter` with `keep_waters: false`, not `preserve`.

### `validation_rules`

**`strictness`** (default `moderate`)

| Value | Description |
|---|---|
| `strict` | Promote every warning diagnostic to an error. |
| `moderate` | Default behavior — errors stay errors, warnings stay warnings. |
| `permissive` | Intended to downgrade error-level rules to warnings. *Not implemented yet — accepted but has no effect.* |

`fail_on_unresolved_issues` (default `True`) — if `True`, any unresolved error causes `status == "failed"`.
`warnings_as_errors` (default `False`) — if `True`, any warning also causes `status == "failed"`.

### `provenance_rules`

| Field | Description | Default |
|---|---|---|
| `record_original_mappings` | Populate the chain/residue/assembly/entity mapping objects; when `False`, those four come back empty (smaller output). **Does not affect `AltlocSelectionMapping`** — that's gated only by `altloc_rules.record_selection`, independently of this flag. | `True` |
| `record_transforms` | Intended to gate whether applied transforms are logged. *Currently unused — transforms are always recorded.* | `True` |
| `record_policy_application` | Intended to embed the full applied policy in the result. *Currently unused — the policy is never embedded.* | `True` |
| `emit_canonicalisation_report` | When `True`, include a warning/error count summary in `canonicalisationProvenance.report`. | `False` |

!!! warning "Naming drift from the design spec"
    The design spec in `policies/canonicalisation.yaml` names the report
    flag `emit_canonicalization_report` (American spelling); the
    implemented schema field is `emit_canonicalisation_report` (British,
    matching the rest of the codebase). Use the schema's spelling.

### Example: three policies, one structure

```python
from pandora.schemas.canonicalisation import (
    canonicalisationPolicy,
    IdentifierRules,
    ChainIdRules,
    ResidueNumberingRules,
    AltlocRules,
    AssemblyRules,
    EntityRules,
    LigandRules,
)

# Preserve everything as reported (including waters — see ligand_rules above).
preserve_policy = canonicalisationPolicy(
    policy_id="preserve",
    policy_name="Preserve As Reported",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        chain_id=ChainIdRules(strategy="preserve"),
        residue_numbering=ResidueNumberingRules(strategy="preserve"),
    ),
    altloc_rules=AltlocRules(strategy="select_best_occupancy"),
    ligand_rules=LigandRules(strategy="preserve"),
)

# Remap chain IDs and renumber residues sequentially.
remap_policy = canonicalisationPolicy(
    policy_id="remap",
    policy_name="Remap Chains And Renumber",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        chain_id=ChainIdRules(strategy="remap"),
        residue_numbering=ResidueNumberingRules(strategy="renumber"),
    ),
    altloc_rules=AltlocRules(
        strategy="select_best_occupancy", tie_breaker="lowest_b_factor"
    ),
    entity_rules=EntityRules(strategy="merge_equivalent_entities"),
    ligand_rules=LigandRules(
        strategy="filter", keep_waters=False, keep_ions=False
    ),
)

# Keep author-assigned identifiers, strict validation.
author_policy = canonicalisationPolicy(
    policy_id="author",
    policy_name="Author Identifiers",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        chain_id=ChainIdRules(strategy="use_auth_chain_id"),
        residue_numbering=ResidueNumberingRules(strategy="use_auth_seq"),
    ),
    altloc_rules=AltlocRules(strategy="select_first"),
    assembly_rules=AssemblyRules(
        strategy="standardize_biological_assembly",
        preferred_assembly_source="author",
    ),
    ligand_rules=LigandRules(strategy="annotate_only"),
)
```

(This is the same set of policies used in `examples/overview.py`, minus
the no-op `keep_waters=False` under the preserve policy — see the
warning above.)

## Provenance

Every rule that deviates from its "preserve" default is recorded as a transform label (e.g. `"chain_id:remap"`) in the returned `canonicalisationProvenance.transforms`, alongside the policy id/name/ version and, if `provenance_rules.emit_canonicalisation_report` is set, a report of warning/error counts. The original-to-canonical identifier mappings themselves come back separately as `CanonicalMappings`.
