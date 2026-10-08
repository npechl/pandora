# Canonicalisation

`pandora.canonicalisation` normalizes a parsed `Structure` according to a policy. The steps always run in this order: modified residues, chain IDs, assemblies, missing data, altlocs, residue numbering, entities, ligands, then validation. Residue numbering runs after the missing-data checks so `renumber` can't hide sequence gaps from them. Each section below changes one rule group on a real fixture from `datasets/dev/mmcif/` and shows the output it prints. The [policy reference](../reference/policies.md#canonicalisation) lists every field and value; [Architecture](../contributing/architecture.md#canonicalisation-policy) shows how the rule groups relate, and [Functions](../reference/functions.md#pandora.canonicalisation) has full signatures.

Every example on this page is a script in `docs/examples/canonicalisation/`, and the test suite checks that each one still prints what is shown here.

## Canonicalise with the default policy

An empty `canonicalisationPolicy` isn't a no-op: several rule groups default to something other than `"preserve"` (missing-atom/residue flagging, best-occupancy altloc selection), so even the default policy records transforms.

=== "`library`"

    ```python
    --8 < --"canonicalisation/default_policy.py"
    ```

    ```text
    --8<-- "canonicalisation/default_policy.out"
    ```

=== "`cli`"

    The same minimal policy, as YAML:

    ```yaml
    policy_id: p1
    policy_name: Default
    policy_version: 1.0.0
    ```

    ```bash
    pandora canonicalise --input-dir raw/ --policy policy.yaml --output-dir canonical/
    # canonicalised 1 structures -> canonical/
    ```

    `canonical/canonicalisation_provenance.json` records the same `transforms` list per entry as `provenance.transforms` above.

## Choose chain IDs

`identifier_rules.chain_id` decides what each chain is called. `remap` relabels chains `A`, `B`, `C`, ... in the order they appear. `use_auth_chain_id` uses the author's chain IDs; when several asym units share one author ID (a protein and its ligands usually do), the later ones get a suffix.

```python
--8 < --"canonicalisation/chain_ids.py"
```

```text
--8<-- "canonicalisation/chain_ids.out"
```

## Renumber residues

`identifier_rules.residue_numbering` decides each residue's `seq_id`. `preserve` keeps the SEQRES position, `use_auth_seq` uses the author's numbering, and `renumber` counts observed residues from 1, closing gaps. Every change is recorded in `mappings.residue_number_mapping`, so you can always get back to the original numbers.

```python
--8 < --"canonicalisation/residue_numbering.py"
```

```text
--8<-- "canonicalisation/residue_numbering.out"
```

!!! warning "Renumbering breaks the link to SEQRES"
    After `use_auth_seq` or `renumber`, `seq_id` no longer indexes `_entity_poly_seq`. [`chain_completeness()`](datasets.md#inspect-chain-completeness) then reports such chains as `SEQRES_MISMATCH` instead of measuring them, and curation's missing-residue rules drop them. Keep `preserve` if you use those rules.

## Pick alternate locations

`altloc_rules` decides which conformer survives when a residue has several. `select_best_occupancy` keeps the most occupied one and breaks ties with `tie_breaker`; `select_user_defined` keeps the altloc you name. Each choice is recorded in `mappings.altloc_selection_mapping`.

```python
--8 < --"canonicalisation/altlocs.py"
```

```text
--8<-- "canonicalisation/altlocs.out"
```

## Expand to the biological assembly

`assembly_rules.strategy: standardize_biological_assembly` picks one assembly (see `preferred_assembly_source`) and applies its symmetry operators, adding each copy as a new chain named `<chain>_<operator>`. `mappings.assembly_mapping` records which chain each copy came from.

```python
--8 < --"canonicalisation/assembly.py"
```

```text
--8<-- "canonicalisation/assembly.out"
```

!!! bug "Known issue: chains outside the assembly are kept"
    Chains E–H above belong to assembly 2, not assembly 1, but they stay in the expanded structure. Until this is fixed, remove them yourself before computing contacts. `test_assembly_expansion_keeps_only_assembly_chains` pins the bug.

## Map modified residues

`modified_residue_rules.strategy: map_to_parent` renames modified polymer residues to their standard parent, for example selenomethionine (`MSE`) to methionine (`MET`), whose `SE` atom becomes `SD`. By default only `MSE` is mapped; list other residues in `comp_ids`, or leave it empty to map every modified residue whose parent the file names in `_pdbx_struct_mod_residue`. A residue with no known parent is left unchanged with a `MODIFIED_RESIDUE_UNMAPPED` warning; Pandora doesn't guess.

```python
--8 < --"canonicalisation/modified_residues.py"
```

```text
--8<-- "canonicalisation/modified_residues.out"
```

## Handle missing atoms and residues

`missing_data_rules` has three parts. `missing_atoms` flags or drops residues missing a backbone atom. `missing_residues` flags or drops chains with gaps. `incomplete_chains` can also cut a gapped chain down to its longest unbroken stretch.

```python
--8 < --"canonicalisation/missing_data.py"
```

```text
--8<-- "canonicalisation/missing_data.out"
```

To measure how much of each chain is missing against SEQRES, without changing the structure, use [`chain_completeness()`](datasets.md#inspect-chain-completeness).

## Normalise entities

`entity_rules.strategy: standardize` renumbers entity IDs `1`, `2`, ...; `merge_equivalent_entities` gives polymer entities with the same canonical sequence one ID. PDB entries already number entities sequentially and give identical chains one entity, so on most files neither changes anything; they matter for files from other sources. See the [policy reference](../reference/policies.md#entity_rules).

## Strip waters, ions and ligands

`ligand_rules.strategy: filter` drops waters, ions and other ligands according to `keep_waters`, `keep_ions` and `keep_nonpolymer_ligands`.

```python
--8 < --"canonicalisation/ligands.py"
```

```text
--8<-- "canonicalisation/ligands.out"
```

### Filter ligands directly

`filter_ligands()`, used at the end of the script above, is what `ligand_rules` calls internally. It's exported so you can filter ligands without the rest of canonicalisation; [curation's content rules](datasets.md#curate-one-structure) reuse it the same way. No CLI subcommand wraps it.

## Trace what changed

`canonicalise_structure` returns three things. The structure is the result; `provenance.transforms` lists the rule groups that changed something, in the order they ran; and `mappings` holds one list per kind of change, so any canonical ID can be traced back to the deposited one. With `provenance_rules.emit_canonicalisation_report`, `provenance.report` also counts the warnings and errors the run raised.

```python
--8 < --"canonicalisation/trace_changes.py"
```

```text
--8<-- "canonicalisation/trace_changes.out"
```

## Load a policy from YAML

Policies are plain data, so you can keep them in YAML next to your dataset and load them with [`load_policy()`](ingestion.md#load-a-canonicalisation-policy). `datasets/canonicalisation.yaml` remaps chain IDs, renumbers residues, merges equivalent entities and drops waters and ions.

=== "`library`"

    ```python
    --8 < --"canonicalisation/from_yaml.py"
    ```

    ```text
    --8<-- "canonicalisation/from_yaml.out"
    ```

=== "`cli`"

    ```bash
    pandora canonicalise --input-dir raw/ --policy datasets/canonicalisation.yaml --output-dir canonical/
    # canonicalised 1 structures -> canonical/
    ```
