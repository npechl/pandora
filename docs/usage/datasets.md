# Datasets

`pandora.datasets` filters structures for inclusion (`curate_structure`, `deduplicate_structures`), measures chain completeness (`chain_completeness`), and reshapes a canonical `Structure` into flat, ML-friendlier records (`extract_*_records`, `entry_sequences`). See [Functions](../reference/functions.md#pandora.datasets) for full signatures.

The curation sections below each apply one kind of rule to real fixtures from `datasets/dev/mmcif/` and show the output. They are scripts in `docs/examples/curation/`, and the test suite checks that each one still prints what is shown here. The [curation policy reference](../reference/policies.md#curation) lists every field. The sections after curation use this setup:

```python
from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.metadata import collect_metadata
from pandora.schemas.canonicalisation import canonicalisationPolicy

policy = canonicalisationPolicy(
    policy_id="p", policy_name="p", policy_version="1.0.0"
)
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
canonical, _, _ = canonicalise_structure(structure, policy)
metadata = collect_metadata(canonical)
```

## Curate one structure

`curate_structure()` applies quality, organism, and content rules and returns `(curated_structure, exclusions, provenance)`. A kept entry has `exclusions` empty, or one record per chain the chain rules removed; an excluded entry has `curated_structure=None` and ends with an entry-level record (`chain_id=None`). Provenance is always populated, regardless of outcome. `metadata` can be `None`; rules that need missing data (such as resolution) then apply their configured default instead of being skipped.

=== "`library`"

    ```python
    --8 < --"curation/curate_one.py"
    ```

    ```text
    --8<-- "curation/curate_one.out"
    ```

=== "`cli`"

    The strict policy above, as YAML:

    ```yaml
    policy_id: c2
    policy_name: Strict
    policy_version: 1.0.0
    quality_rules:
      max_resolution: 1.0
    ```

    ```bash
    pandora curate --input-dir canonical/ --policy curation.yaml --output-dir curated/
    # curated: 4 retained, 1 excluded -> curated/
    ```

    Exclusion records, including removed chains, are written to `curated/curation_exclusions.json`. The CLI always runs `collect_metadata()` on each structure, so there's no way to pass `metadata=None` from it.

!!! note "Migrating from the single-exclusion return"
    Before chain-level rules, `curate_structure` returned `(structure, exclusion | None, provenance)`. It now returns a list. Replace `if curated is None: excluded.append(exclusion)` with `excluded.extend(exclusions)`, which also keeps the records of chains removed from entries that were kept.

Rules run in this order: entry rules (resolution, method, R-factors, organism, non-standard residues, `max_atoms`), then chain rules, then `min_chain_length` and `min_polymer_chains` on the chains that are left, then content rules.

## Set resolution limits per method

`max_resolution` applies one limit to every entry. `max_resolution_by_method` sets a limit per experimental method and overrides `max_resolution` for entries with that method; an entry with several methods gets the strictest. Cryo-EM resolution comes from `_em_3d_reconstruction`. An entry whose method has no limit and no `max_resolution` isn't checked at all.

```python
--8 < --"curation/resolution_by_method.py"
```

```text
--8<-- "curation/resolution_by_method.out"
```

## Check X-ray R-factors

`max_r_free`, `max_r_free_gap` (\|Rfree − Rwork\|) and `max_r_sym` check the entries whose method is in `rfactor_methods` (X-ray by default), so cryo-EM and NMR entries never fail them. `max_r_sym` uses Rsym, or Rmerge when Rsym isn't reported. A missing value keeps the entry unless `null_rfactor_behavior` is `exclude`.

```python
--8 < --"curation/rfactors.py"
```

```text
--8<-- "curation/rfactors.out"
```

## Reject non-standard residues

`exclude_nonstandard_residues` excludes an entry with any polymer residue that isn't one of the 20 amino acids, `UNK`, or a standard nucleotide. `allowed_nonstandard_residues` lets named residues through. Curation sees the canonical structure, so residues [mapped to their parent](canonicalisation.md#map-modified-residues) during canonicalisation, such as MSE, already count as standard.

```python
--8 < --"curation/nonstandard_residues.py"
```

```text
--8<-- "curation/nonstandard_residues.out"
```

## Drop chains with missing residues

`max_missing_tail_fraction` and `max_missing_middle_fraction` remove chains with too many missing residues, measured against SEQRES by [`chain_completeness()`](#inspect-chain-completeness). The tail fraction is the missing N- and C-terminal residues over the SEQRES length; the middle fraction is the missing residues between the first and last present ones over the SEQRES length minus the missing tails. `missing_residue_definition` decides whether a residue missing part of its backbone counts as missing (the default, `incomplete_backbone`) or only one with no atoms at all (`unobserved`).

Removing a chain keeps the rest of the entry. `min_polymer_chains` and `min_chain_length` then check what is left.

```python
--8 < --"curation/missing_residues.py"
```

```text
--8<-- "curation/missing_residues.out"
```

When a missing-fraction rule is set, a chain that can't be measured (no SEQRES, or a renumbered chain) is removed as `NO_SEQRES`.

## Cap the size

`max_atoms` excludes a whole entry above an atom count; `max_chain_length` removes single chains with more observed residues than the limit. `max_atoms` is checked after parsing, so it protects later steps such as contacts and export but not the parse itself.

```python
--8 < --"curation/size_limits.py"
```

```text
--8<-- "curation/size_limits.out"
```

## Filter by organism

`organism_rules` keeps (`include_taxa`) or drops (`exclude_taxa`) entries by NCBI taxonomy ID, read from the metadata. Without taxonomy metadata, an active organism filter excludes the entry as `MISSING_TAXONOMY`.

```python
--8 < --"curation/organism.py"
```

```text
--8<-- "curation/organism.out"
```

## Keep only the polymer

`content_rules` never excludes anything. It strips ligands, waters or ions from structures that pass, the same way [`filter_ligands()`](canonicalisation.md#filter-ligands-directly) does.

```python
--8 < --"curation/content_rules.py"
```

```text
--8<-- "curation/content_rules.out"
```

## Collect exclusions across a batch

Curating many entries, extend one list with every call's records. It then holds entry records (`chain_id` is `None`) for entries left out, and chain records for chains removed from entries that were kept. Write it out with `write_records()` or pass it to [`build_dataset_manifest()`](provenance.md#assemble-a-dataset-manifest) as `excluded`.

```python
--8 < --"curation/batch_exclusions.py"
```

```text
--8<-- "curation/batch_exclusions.out"
```

For whole benchmark policies built from these rules, see the [benchmark quality filters recipe](../recipes/benchmark-quality-filters.md).

## Inspect chain completeness

`chain_completeness()` counts each polymer chain's missing residues against SEQRES (`_entity_poly_seq`): at the N-terminus, at the C-terminus, and in the middle. The missing-residue rules use it; you can also call it yourself, for example to pick the chain with the fewest missing residues. A CA-only model shows why the definition matters: every residue is observed, but none has a full backbone.

```python
--8 < --"curation/completeness.py"
```

```text
--8<-- "curation/completeness.out"
```

Run it on a structure whose `label_seq_id`s still index SEQRES. After `residue_numbering: renumber` or `use_auth_seq` they don't, and the chain is reported as `SEQRES_MISMATCH` instead of measured.

## Deduplicate a batch

`deduplicate_structures()` drops every structure after the first with a
given `entry_id`:

=== "`library`"

    ```python
    from pandora.datasets import deduplicate_structures
    from pandora.schemas.dataset import DeduplicationRules

    duplicate, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
    duplicate_canonical, _, _ = canonicalise_structure(duplicate, policy)

    retained, removed, provenance = deduplicate_structures(
        [canonical, duplicate_canonical], DeduplicationRules(enabled=True)
    )
    print(
        len(retained),
        "retained,",
        len(removed),
        "removed:",
        [r.reason_code for r in removed],
    )
    # 1 retained, 1 removed: ['DUPLICATE']
    ```

    `DeduplicationRules(enabled=False)` (the default) is a no-op — every
    structure passes through, `removed` is always `[]`.

=== "`cli`"

    ```bash
    pandora dedup --input-dir curated/ --output-dir deduped/
    # dedup: 5 retained, 0 removed -> deduped/
    ```

    `--disable` runs `DeduplicationRules(enabled=False)` instead — every
    structure passes through, but dedup provenance is still recorded.

## Reshape into flat records

Each `extract_*_records()` function turns one canonical `Structure`
into a list of flat, independently-serializable records — see
[Export](export.md) for writing them to disk. No CLI subcommand wraps
this reshaping step.

```python
from pandora.datasets import (
    entry_sequences,
    extract_chain_records,
    extract_residue_records,
    extract_interface_records,
)

# One representative sequence per entry (its longest polymer chain) —
# what compute_sequence_similarity()/entry-keyed clustering expect.
sequences = entry_sequences({canonical.entry_id: canonical})
print(sequences)
# {'104M': 'MVLSEGEWQLVLHVWAKVEAD...'}  (157 residues)

chains = extract_chain_records(canonical)
print([(c.chain_id, c.residue_count) for c in chains])
# [('A', 153)]

residues = extract_residue_records(canonical)
print(len(residues), residues[0].chain_id, residues[0].comp_id)
# 153 A VAL

structure2, _, _ = mmcif_to_structure("datasets/dev/mmcif/1a3n.cif")
canonical2, _, _ = canonicalise_structure(structure2, policy)
interfaces = extract_interface_records(canonical2, distance_cutoff=4.0)
print(
    len(interfaces),
    interfaces[0].chain_id_1,
    interfaces[0].chain_id_2,
    interfaces[0].contact_count,
)
# 5 A B 35
```

`extract_residue_records()` carries each residue's full atom list (coordinates, B-factor, ...) rather than a bare count — residue-level ML use cases (contact maps, solvent exposure) don't need to go back to the structure. `extract_interface_records()` reshapes [`annotate_chain_interfaces()`](annotation.md#chain-chain-interfaces)'s output; it doesn't compute contacts itself. It takes the same `atom_set` and `polymer_types` arguments, and each `InterfaceRecord` carries the `atom_set` used and the `residue_pairs`.
