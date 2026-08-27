# Annotation

`pandora.annotations` computes derived structural summaries — nothing
here is source-backed like [Metadata](metadata.md); every function
recomputes its result from `Structure.atoms`/`entities` each call. See
[Functions](../reference/functions.md#pandora.annotations) for full
signatures.

All examples below use `1a3n` (deoxyhemoglobin: 2 alpha + 2 beta
chains, 4 heme groups) so that the contact-based annotations have
something to find.

```python
from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.schemas.canonicalisation import (
    canonicalisationPolicy,
    AltlocRules,
    EntityRules,
    LigandRules,
)

policy = canonicalisationPolicy(
    policy_id="p",
    policy_name="p",
    policy_version="1.0.0",
    altloc_rules=AltlocRules(strategy="select_best_occupancy"),
    entity_rules=EntityRules(strategy="merge_equivalent_entities"),
    ligand_rules=LigandRules(strategy="filter", keep_waters=False),
)
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1a3n.cif")
canonical, _, _ = canonicalise_structure(structure, policy)
```

## Per-entry counts

`annotate_structure_counts()` summarizes atom/residue/chain/entity/
assembly counts and a few breakdowns — a cheap sanity check to run on
every structure in a batch.

=== "`library`"

    ```python
    from pandora.annotations import annotate_structure_counts

    counts = annotate_structure_counts(canonical)
    print(
        counts.data["asym_unit_count"],
        "chains,",
        counts.data["residue_count"],
        "residues",
    )
    # 8 chains, 576 residues
    ```

=== "`cli`"

    ```bash
    pandora annotate --input-dir deduped/ --layers structure_counts --output-dir annotations/
    # annotated 5 entries -> annotations/
    ```

    Each entry's `structure_counts` layer lands in
    `annotations/annotations.json`, keyed by entry id.

## Chain-chain interfaces

`annotate_chain_interfaces()` finds every polymer chain pair with at
least one contact within `distance_cutoff` angstroms.

=== "`library`"

    ```python
    from pandora.annotations import annotate_chain_interfaces

    interfaces = annotate_chain_interfaces(canonical, distance_cutoff=4.0)
    print(len(interfaces.data["interfaces"]), "interfaces")
    # 5 interfaces
    for interface in interfaces.data["interfaces"][:2]:
        print(
            interface["chain_id_1"],
            interface["chain_id_2"],
            interface["contact_count"],
        )
    # A B 35
    # A C 6
    ```

=== "`cli`"

    ```bash
    pandora annotate --input-dir deduped/ --layers chain_interfaces --output-dir annotations/
    # annotated 5 entries -> annotations/
    ```

    `distance_cutoff` isn't exposed as a CLI flag — the CLI always uses
    `annotate_chain_interfaces()`'s default.

## Ligand contacts

`annotate_ligand_contacts()` finds polymer residues near each
non-polymer ligand (set `include_waters=True` to treat waters as
ligands too).

=== "`library`"

    ```python
    from pandora.annotations import annotate_ligand_contacts

    contacts = annotate_ligand_contacts(canonical, distance_cutoff=4.0)
    print(len(contacts.data["ligands"]), "ligands")
    # 4 ligands — one per heme group
    ```

=== "`cli`"

    ```bash
    pandora annotate --input-dir deduped/ --layers ligand_contacts --output-dir annotations/
    # annotated 5 entries -> annotations/
    ```

    `--layers` accepts several names at once (`--layers structure_counts
    ligand_contacts chain_interfaces`) to compute them all in one pass —
    it defaults to all three when omitted.

## Pairwise sequence identity

`annotate_pairwise_sequence_identity()` compares every polymer entity
in one structure against every polymer entity in another with a
simple ungapped, position-wise comparison (no alignment step — see
`examples/dataset_pipeline.py` for where this genuinely fools itself
on a shifted reading frame), and reports the best-scoring pair.

=== "`library`"

    ```python
    from pandora.annotations import annotate_pairwise_sequence_identity

    s2, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
    c2, _, _ = canonicalise_structure(s2, policy)

    identity = annotate_pairwise_sequence_identity(canonical, c2)
    print(identity.data["best_identity"])
    # 0.168 — 1a3n (hemoglobin) and 104m (myoglobin) are only distantly related
    ```

=== "`cli`"

    `--pairwise` computes this layer for every pair in the batch,
    alongside whichever `--layers` you pick:

    ```bash
    pandora annotate --input-dir deduped/ --layers structure_counts --pairwise --output-dir annotations/
    # annotated 5 entries -> annotations/
    ```

    Each pair's result is attached to both entries' entries in
    `annotations/annotations.json`.

For a real alignment-based identity between *sequences* (not entities
within structures) at dataset scale, see
[`compute_sequence_similarity()`](similarity.md#sequence-similarity-mmseqs2).
