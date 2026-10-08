# Quality and completeness: design

## Roadmap context

Sub-project 3 of the plan to let Pandora rebuild the four benchmark
specs in `examples/ppi/benchmarks/` (decisions in
`2026-10-07-contacts-design.md`, "Roadmap context"). The benchmarks are
equal peers; this spec covers the quality and completeness filters
they share.

| Need | Benchmarks | Before this spec |
|---|---|---|
| Cryo-EM resolution | ATOM3D, PPI v2.1, ProteinFlow, Pinder | Bug: only `_refine.ls_d_res_high` is read, so cryo-EM entries get `resolution=None` and any `max_resolution` excludes them as `NULL_RESOLUTION` |
| Resolution limit per method | PPI v2.1 | One `max_resolution` |
| Rfree, \|Rfree − Rwork\|, Rsym limits, X-ray only | PPI v2.1 (Pinder test set) | Rfree/Rwork extracted, no rules; Rsym/Rmerge not extracted |
| MSE→MET; exclude non-standard residues | PPI v2.1 | Not handled |
| Missing residues against SEQRES, tails vs middle | PPI v2.1, ProteinFlow | Only `label_seq_id` gap diagnostics, blind to missing tails |
| Maximum size | ProteinFlow (chain length), Spike 0 (atoms) | Minimums only |

Out of scope: the cryo-EM FSC ≤ 0.143 rule (it filters nothing,
PPI v2.1 limitation C3), interface gaps (pair layer, sub-project 5),
a size guard before parsing, and per-structure modelled-residue flags
(sub-project 4).

## Decisions

1. Chain-level rules **drop chains** and keep the rest of the entry.
   `curate_structure` returns the trimmed structure and one
   `ExclusionRecord` per dropped chain. An entry with no polymer chain
   left is excluded. This changes `curate_structure`'s return type.
2. Modified residues are mapped to their parent only when the parent is
   known from `_pdbx_struct_mod_residue`, or for MSE. Anything else stays
   unmapped and gets a diagnostic. No CCD lookup.
3. R-factor rules apply only to methods in `rfactor_methods` (default
   X-ray). A missing R-factor is **included** by default
   (`null_rfactor_behavior="include"`), because many X-ray entries lack
   Rsym.
4. `max_r_sym` reads `r_sym`, falling back to `r_merge`. They are the
   same statistic under two names.
5. Missing-residue fractions follow PPI v2.1: tails as a fraction of
   the SEQRES length; middle as a fraction of the SEQRES length minus
   the missing tails.

## Metadata: `QualityRecord`

`pandora/schemas/metadata.py`, `pandora/metadata/mmcif.py`.

- `resolution`: `_refine.ls_d_res_high`, falling back to
  `_em_3d_reconstruction.resolution`. Provenance lists
  `_em_3d_reconstruction` when the fallback was used.
- New `r_sym: float | None` from `_reflns.pdbx_Rsym_value`.
- New `r_merge: float | None` from `_reflns.pdbx_Rmerge_I_obs`.

## Canonicalisation: `modified_residue_rules`

New `ModifiedResidueRules` in `schemas/canonicalisation.py`, a
`modified_residue_rules` field on `canonicalisationPolicy`, and a new
module `pandora/canonicalisation/modified_residues.py`.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `strategy` | `"preserve" \| "map_to_parent"` | `"preserve"` | Whether modified residues are renamed to their parent |
| `comp_ids` | `list[str]` | `["MSE"]` | Which modified residues to map; empty means every one with a known parent |

Parents come from `_pdbx_struct_mod_residue` (`label_asym_id`,
`label_seq_id`, `label_comp_id` → `parent_comp_id`). MSE maps to MET
even when that category doesn't list it. Mapping a residue:

- sets `label_comp_id` (and `auth_comp_id`) to the parent;
- sets `group_PDB` to `"ATOM"`;
- for MSE only: renames atom `SE` to `SD` (`label_atom_id` and
  `auth_atom_id`) and sets `type_symbol` to `"S"`.

A residue named in `comp_ids` (or, with an empty list, any polymer
`HETATM` residue whose comp_id isn't standard) that has no known parent
gets a `MODIFIED_RESIDUE_UNMAPPED` warning. Each mapping is recorded in a new
`ModifiedResidueMapping` (`chain_id`, `seq_id`, `auth_seq_id`,
`original_comp_id`, `parent_comp_id`) on `CanonicalMappings`. The
transform label is `modified_residues:map_to_parent`.

The step runs first in `canonicalise_structure`, before chain-id
normalisation, so `label_asym_id` and `label_seq_id` still match
`_pdbx_struct_mod_residue`. Mapping items record the original chain id.
"Modified residue" means a polymer residue (`label_seq_id` set) whose
comp_id is not standard; the parser already marks polymer residues
`ATOM`, so `group_PDB` is set only for safety.

Known limit: canonicalisation diagnostics are not returned (see
CLAUDE.md, "Known issues"), so the unmapped warning surfaces only as a
count.

## Datasets: `chain_completeness`

New `ChainCompleteness` in `schemas/dataset.py`:

| Field | Meaning |
|---|---|
| `entry_id`, `chain_id` | The chain (`label_asym_id`) |
| `seqres_length` | Positions in `_entity_poly_seq` for the chain's entity |
| `missing_n_term` | Missing positions before the first present residue |
| `missing_c_term` | Missing positions after the last present residue |
| `missing_middle` | Missing positions between the first and last present residue |

New function in `pandora/datasets/completeness.py`:

```python
def chain_completeness(
    structure: Structure,
    definition: MissingResidueDefinition = "incomplete_backbone",
) -> tuple[list[ChainCompleteness], DiagnosticBundle]
```

- `MissingResidueDefinition = Literal["unobserved",
  "incomplete_backbone"]`. With `"incomplete_backbone"`, a residue with
  any backbone atom missing counts as missing. Backbone sets are the
  ones in `canonicalisation/missing_data.py` (N, CA, C, O for protein;
  the nucleic set for nucleic-acid entities), shared rather than copied.
- SEQRES positions are the distinct `num` values of `_entity_poly_seq`
  for the entity, so microheterogeneity counts once.
- Only polymer chains in `structure.asym_units` are reported.
- A chain whose entity has no `_entity_poly_seq` rows gets no record and
  a `NO_SEQRES` diagnostic. A chain with a `label_seq_id` outside
  `1..seqres_length`, or a residue whose comp_id differs from every
  `mon_id` at that position, gets no record and a `SEQRES_MISMATCH`
  diagnostic. A non-standard `mon_id` matches any comp_id, so a residue
  mapped by `modified_residue_rules` (MSE → MET) still lines up.
  Renumbering keeps ids within range, so the comp_id check is what
  catches it.
- A chain with no present residue has `missing_n_term = seqres_length`
  and zero elsewhere.

## Curation

`pandora/schemas/dataset.py`, `pandora/datasets/curation.py`.

### New `QualityRules` fields

| Field | Type | Default | Level |
|---|---|---|---|
| `max_resolution_by_method` | `dict[str, float]` | `{}` | entry |
| `rfactor_methods` | `list[str]` | `["X-RAY DIFFRACTION"]` | entry |
| `max_r_free` | `float \| None` | `None` | entry |
| `max_r_free_gap` | `float \| None` | `None` | entry |
| `max_r_sym` | `float \| None` | `None` | entry |
| `null_rfactor_behavior` | `"include" \| "exclude"` | `"include"` | entry |
| `exclude_nonstandard_residues` | `bool` | `False` | entry |
| `allowed_nonstandard_residues` | `list[str]` | `[]` | entry |
| `max_atoms` | `int \| None` | `None` | entry |
| `missing_residue_definition` | `MissingResidueDefinition` | `"incomplete_backbone"` | chain |
| `max_missing_tail_fraction` | `float \| None` | `None` | chain |
| `max_missing_middle_fraction` | `float \| None` | `None` | chain |
| `max_chain_length` | `int \| None` | `None` | chain |

Semantics:

- **Method matching** is case- and whitespace-insensitive throughout.
  `experimental_method` may list several methods joined by `"; "`.
- **Resolution limit**: the strictest of `max_resolution_by_method`
  values whose method the entry has; `max_resolution` if none match.
  `null_resolution_behavior` applies whenever a limit applies.
- **R-factor rules** run only if the entry has a method in
  `rfactor_methods`. For each set limit: a missing value follows
  `null_rfactor_behavior` (`NULL_RFACTOR`); otherwise `r_free >
  max_r_free` → `RFREE_THRESHOLD`, `|r_free − r_work| > max_r_free_gap`
  → `RFREE_GAP_THRESHOLD` (missing if either is missing), `(r_sym or
  r_merge) > max_r_sym` → `RSYM_THRESHOLD`.
- **Non-standard residues**: polymer residues present in the atoms
  (`label_seq_id` not None) whose comp_id is not one of the 20 standard
  amino acids, `UNK`, or the standard nucleotides (A, C, G, U, DA, DC,
  DG, DT, DU, N, DN), and not in `allowed_nonstandard_residues`, exclude
  the entry (`NONSTANDARD_RESIDUE`). Runs on the canonical structure,
  so mapped MSE is already MET.
- **`max_atoms`**: `len(structure.atoms) > max_atoms` →
  `TOO_MANY_ATOMS`.
- **Chain rules** drop a chain when its tail fraction exceeds
  `max_missing_tail_fraction` (`MISSING_TAILS`), its middle fraction
  exceeds `max_missing_middle_fraction` (`MISSING_MIDDLE`), or its
  observed residue count exceeds `max_chain_length` (`CHAIN_TOO_LONG`).
  If either missing-fraction limit is set and `chain_completeness`
  returned no record for a chain, the chain is dropped as `NO_SEQRES`.
  A middle fraction with a zero denominator is 0.

### `ExclusionRecord`

- New `chain_id: str | None = None`. `None` means the whole entry.
- New reason codes: `RFREE_THRESHOLD`, `RFREE_GAP_THRESHOLD`,
  `RSYM_THRESHOLD`, `NULL_RFACTOR`, `NONSTANDARD_RESIDUE`,
  `TOO_MANY_ATOMS`, `MISSING_TAILS`, `MISSING_MIDDLE`, `CHAIN_TOO_LONG`,
  `NO_SEQRES`, `NO_CHAINS_LEFT`.

### `curate_structure`

```python
def curate_structure(
    structure: Structure,
    metadata: MetadataRecord | None,
    policy: DatasetCurationPolicy,
) -> tuple[Structure | None, list[ExclusionRecord], CurationProvenance]
```

Order:

1. Entry rules: resolution, method, R-factors, organism, non-standard
   residues, `max_atoms`. The first failure excludes the entry; the
   list holds that one record.
2. Chain rules: drop failing chains (atoms and asym units), one record
   each.
3. On the trimmed structure: if a chain was dropped and no polymer chain
   is left → `NO_CHAINS_LEFT`; then `min_chain_length`,
   `min_polymer_chains`. An entry exclusion here is appended after the
   chain records.
4. Content rules, as today.

Returns `(structure, chain_records, provenance)` when the entry is kept
and `(None, records, provenance)` when it is excluded.

Callers to update: `pandora/cli/app.py`,
`pandora/provenance/reproduce.py`, `pandora/provenance/manifest.py`,
`examples/dataset_pipeline.py`, `examples/ppi_dataset_pipeline.py`,
`tests/test_curation.py`, and the docs that show `curate_structure`.
`docs/usage/datasets.md` gets a short migration note.

## Policies and docs

- `docs/reference/policies.md`: every new field, including the fraction
  definitions and the X-ray-only scoping.
- `docs/policies/dataset.yaml`, `docs/policies/canonicalisation.yaml`:
  new field types.
- `docs/usage/datasets.md`: `chain_completeness` and chain-level
  exclusions. `docs/usage/canonicalisation.md`: modified residues.
- Regenerate the ERD diagrams (schema relationships change).

## Testing

Fixtures from `datasets/dev/mmcif/` only.

| Test | Fixture |
|---|---|
| Cryo-EM resolution falls back to `_em_3d_reconstruction` | 22jy (2.2 Å) |
| `r_merge` extracted, `r_sym` None; `max_r_sym` uses `r_merge` | 1a08 |
| MSE → MET with SE → SD, mapping recorded | 1b6w |
| Mod residue with parent in `_pdbx_struct_mod_residue` mapped when listed | 1a08 (FTY → TYR) |
| Non-standard residue excludes; allowed list lets it through | 1a08 |
| Completeness counts tails and middle | 10mv (275: 21 / 18 / 5), 1aui |
| `incomplete_backbone` counts CA-only residues as missing | 1p58 (CA trace) |
| Renumbered chain → `SEQRES_MISMATCH`, no record | fixture after `residue_numbering="renumber"` |
| Per-method resolution, strictest wins | hand-built `MetadataRecord`s |
| R-factor rules skip cryo-EM; null behaviour both ways | hand-built `MetadataRecord`s |
| Chain drop leaves the rest; all dropped → `NO_CHAINS_LEFT` | fixture with two protein chains |
| `max_atoms`, `max_chain_length` | any fixture |
| No input mutated | each new function |
