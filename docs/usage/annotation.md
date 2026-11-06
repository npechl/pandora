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
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/a3/1a3n.cif")
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
    pandora annotate compute --input-dir deduped/ --layers structure_counts --output-dir annotations/
    # annotated 5 entries -> annotations/
    ```

    Each entry's `structure_counts` layer lands in
    `annotations/annotations.json`, keyed by entry id.

## Chain-chain interfaces

`annotate_chain_interfaces()` finds every polymer chain pair with at least one pair of selected atoms within `distance_cutoff` angstroms. For each pair it reports the contacting residues on each side and the residue pairs themselves (`residue_pairs`).

| Argument | Default | Meaning |
|---|---|---|
| `distance_cutoff` | `4.0` | Contact distance in angstroms |
| `atom_set` | `"heavy"` | Which atoms count; see below |
| `polymer_types` | `None` | Keep only chains of these entity polymer types, e.g. `["polypeptide(L)"]`. `None` keeps every polymer chain; `[]` keeps none |

| `atom_set` | Atoms | Used by |
|---|---|---|
| `all` | Every atom, hydrogens included | — |
| `heavy` | Every atom except H and D | ATOM3D (heavy-atom option), PRODIGY |
| `backbone` | Protein N, CA, C, O; nucleic-acid O5′, C5′, C4′, C3′, O3′ | Pinder's 10 Å pair gate |
| `ca` | CA atoms only | ATOM3D's 8 Å Cα default |

!!! warning "Changed after 0.5.6: hydrogens no longer count"
    Up to 0.5.6 every atom counted, including hydrogens, so structures that model H atoms got extra contacts. The default is now `atom_set="heavy"`. Layers record `atom_set` in `parameters`, so `reproduce_dataset()` rebuilds them exactly. Layers from older manifests record only `distance_cutoff` and are rebuilt with `heavy`; pass `atom_set="all"` yourself to match the old output.

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
    pandora annotate compute --input-dir deduped/ --layers chain_interfaces --output-dir annotations/
    # annotated 5 entries -> annotations/
    ```

    `distance_cutoff`, `atom_set` and `polymer_types` aren't exposed as CLI flags. The CLI always uses `annotate_chain_interfaces()`'s defaults.

### Pair gate, then interface (Pinder-style)

Pinder first keeps chain pairs with backbone atoms within 10 Å, then defines the interface with a tighter cutoff. That's two calls:

```python
gate = annotate_chain_interfaces(
    canonical, distance_cutoff=10.0, atom_set="backbone"
)
interfaces = annotate_chain_interfaces(canonical, distance_cutoff=4.0)

gated = {(i["chain_id_1"], i["chain_id_2"]) for i in gate.data["interfaces"]}
kept = [
    i
    for i in interfaces.data["interfaces"]
    if (i["chain_id_1"], i["chain_id_2"]) in gated
]
```

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
    pandora annotate compute --input-dir deduped/ --layers ligand_contacts --output-dir annotations/
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

    s2, _, _ = mmcif_to_structure("datasets/dev/mmcif/04/104m.cif")
    c2, _, _ = canonicalise_structure(s2, policy)

    identity = annotate_pairwise_sequence_identity(canonical, c2)
    print(identity.data["best_identity"])
    # 0.168 — 1a3n (hemoglobin) and 104m (myoglobin) are only distantly related
    ```

=== "`cli`"

    `--pairwise` computes this layer for every pair in the batch,
    alongside whichever `--layers` you pick:

    ```bash
    pandora annotate compute --input-dir deduped/ --layers structure_counts --pairwise --output-dir annotations/
    # annotated 5 entries -> annotations/
    ```

    Each pair's result is attached to both entries' entries in
    `annotations/annotations.json`.

For a real alignment-based identity between *sequences* (not entities
within structures) at dataset scale, see
[`compute_sequence_similarity()`](similarity.md#sequence-similarity-mmseqs2).

## Importing a layer from another tool

Pandora doesn't run tools such as COCaDA, FreeSASA or DSSP. You run the tool on the structures of a step folder (for example the output of `pandora curate`), convert its output into one of Pandora's layer templates, and import it. Pandora links every row to the residue or chain it names, reports the rows it can't link, and records where the data came from.

### Templates

| Template | One row per | Required value column |
|---|---|---|
| `residue_contacts` | residue pair (sides `_1`, `_2`) | none (`contact_type`, `distance` are common extras) |
| `residue_values` | residue | `value` (number) |
| `residue_labels` | residue | `label` (text) |
| `chain_values` | chain | `value` (number) |
| `chain_labels` | chain | `label` (text) |
| `interface_labels` | chain pair (sides `_1`, `_2`) | `label` (text) |
| `free_form` | anything | none |

Every row needs `entry_id` (any case). Each side names its residue or chain in either naming:

| Naming | Residue columns | Chain columns |
|---|---|---|
| Pandora's (mmCIF `label_*`) | `chain_id`, `seq_id` | `chain_id` |
| Author's (as in PDB files) | `auth_chain_id`, `auth_seq_id`, `ins_code` (optional) | `auth_chain_id` |

Pair templates add `_1` / `_2` to every key column (`chain_id_1`, `auth_seq_id_2`, ...). Other columns are kept as they are. A `free_form` row is checked against whichever key columns it has. A chain named by author ID resolves among polymer chains only, because ligands and waters usually share their polymer's author chain ID.

Ligands and waters have no `seq_id`, so name them in the author's naming. Each layer holds one `value` or `label` column; import several measurements as several layers.

The stored layer carries **both** namings, filled in from the structure, so later steps can rely on Pandora's keys.

=== "`library`"

    ```python
    from pandora.annotations import check_layer_columns, resolve_layer_rows
    from pandora.parsing import mmcif_to_structure

    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/a2/1a2c.cif")
    rows = [
        {
            "entry_id": "1a2c",
            "auth_chain_id_1": "L",
            "auth_seq_id_1": "1",
            "ins_code_1": "H",
            "auth_chain_id_2": "H",
            "auth_seq_id_2": "16",
            "contact_type": "hydrogen_bond",
        },
        {
            "entry_id": "1a2c",
            "auth_chain_id_1": "L",
            "auth_seq_id_1": "1",
            "ins_code_1": "H",
            "auth_chain_id_2": "H",
            "auth_seq_id_2": "147",
            "ins_code_2": "A",
            "contact_type": "hydrophobic",
        },
    ]
    check_layer_columns({c for row in rows for c in row}, "residue_contacts")
    matched, unmatched = resolve_layer_rows(
        rows, "residue_contacts", {structure.entry_id: structure}
    )
    row = matched[0]
    print(
        row["chain_id_1"], row["seq_id_1"], "-", row["chain_id_2"], row["seq_id_2"]
    )
    # A 1 - B 1
    print(unmatched[0]["reason"])
    # no_coordinates  (H 147A is in the sequence but has no atoms)
    ```

=== "`cli`"

    Describe the layer once:

    ```yaml
    # cocada.yaml
    name: cocada
    template: residue_contacts
    tool: COCaDA
    tool_version: "1.0"
    parameters:        # recorded as given
      cutoff: 5.0
    ```

    Then import the converted table (Parquet or JSON Lines) against the step folder the tool ran on:

    ```bash
    pandora annotate import --structures-dir curated/ \
      --table cocada_contacts.jsonl --description cocada.yaml \
      --source-files cocada_out/ --output-dir annotations/cocada
    # imported: 1 rows, 1 unmatched -> annotations/cocada
    jq -c '.counts' annotations/cocada/layer.json
    # {"rows_in":2,"rows_matched":1,"rows_unmatched":1,"unmatched:no_coordinates":1}
    ```

### The layer folder

| File | Content |
|---|---|
| `layer.<fmt>` | rows that matched, with both namings filled in and every other column kept |
| `unmatched.<fmt>` | rows that didn't match, with a `reason` column |
| `source/` | the tool's original files (only with `--source-files`) |
| `layer.json` | step summary: the description, Pandora version, counts, and `previous_summary` pointing to the step folder's `curation.json` or `canonicalisation.json` |

Unmatched rows don't stop the import. A table that lacks a required column altogether, or a description naming an unknown template, stops it before anything is written. `--overwrite` replaces an existing layer.

| Reason | Meaning |
|---|---|
| `entry_not_found` | the entry isn't in the step folder |
| `missing_key` | neither naming is complete for a side |
| `chain_not_found` | no such chain |
| `residue_not_found` | no such residue |
| `no_coordinates` | the residue is in the sequence but has no atoms (for author naming, only while the structure keeps `_pdbx_poly_seq_scheme`; otherwise `residue_not_found`) |
| `naming_conflict` | the row gives both namings and they point at different residues or chains |
| `ambiguous_author_id` | the author IDs point at more than one residue or polymer chain |
| `bad_value` | `seq_id` isn't an integer, `value` isn't a number, or `label` is empty |
| `duplicate_key` | a second row for the same key; kept for `residue_contacts` and `free_form`, where repeats are allowed |
