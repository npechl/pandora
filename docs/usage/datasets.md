# Datasets

`pandora.datasets` filters structures for inclusion (`curate_structure`,
`deduplicate_structures`) and reshapes a canonical `Structure` into flat,
ML-friendlier records (`extract_*_records`, `entry_sequences`). See
[Functions](../reference/functions.md#pandora.datasets) for full
signatures.

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

`curate_structure()` applies quality, organism, and content rules and
returns either `(curated_structure, None, provenance)` or
`(None, exclusion, provenance)` — provenance is always populated,
regardless of outcome.

=== "`library`"

    ```python
    from pandora.datasets import curate_structure
    from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules

    curation_policy = DatasetCurationPolicy(
        policy_id="c1", policy_name="Default", policy_version="1.0.0"
    )
    curated, exclusion, provenance = curate_structure(
        canonical, metadata, curation_policy
    )
    print(curated is not None, exclusion)
    # True None
    ```

    A structure that fails a rule comes back excluded, with a
    machine-readable reason:

    ```python
    strict_policy = DatasetCurationPolicy(
        policy_id="c2",
        policy_name="Strict",
        policy_version="1.0.0",
        quality_rules=QualityRules(max_resolution=1.0),
    )
    curated, exclusion, _ = curate_structure(canonical, metadata, strict_policy)
    print(exclusion.reason_code, "-", exclusion.message)
    # RESOLUTION_THRESHOLD - resolution 1.71 exceeds max_resolution 1.0
    ```

    `content_rules` (waters/ions/other ligands) never excludes — it
    filters the retained structure's atoms in place, the same way
    [`filter_ligands()`](canonicalisation.md#filter-ligands-directly) does:

    ```python
    from pandora.schemas.dataset import ContentRules

    content_policy = DatasetCurationPolicy(
        policy_id="c3",
        policy_name="NoLigands",
        policy_version="1.0.0",
        content_rules=ContentRules(
            keep_ligands=False, keep_waters=False, keep_ions=False
        ),
    )
    curated, _, _ = curate_structure(canonical, metadata, content_policy)
    print(len(canonical.atoms), "->", len(curated.atoms))
    # 1450 -> 1217
    ```

    `metadata` can be `None` — quality/organism checks that depend on
    missing data (e.g. null resolution) then apply their configured
    default instead of being skipped.

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

    Excluded structures' reasons are written to
    `curated/curation_exclusions.json`; `metadata` is always
    `collect_metadata()`'d from each structure internally, so there's no
    way to pass `metadata=None` from the CLI.

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

`extract_residue_records()` carries each residue's full atom list
(coordinates, B-factor, ...) rather than a bare count — residue-level ML
use cases (contact maps, solvent exposure) don't need to go back to the
structure. `extract_interface_records()` reshapes
[`annotate_chain_interfaces()`](annotation.md#chain-chain-interfaces)'s
output; it doesn't compute contacts itself.
