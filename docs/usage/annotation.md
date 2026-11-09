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

## Contacts are one layer

Pandora's contact annotations return a `residue_contacts` layer: a list of `ResidueContactRow`s, one per pair of residues in contact. It is the same table shape that [`pandora annotate import`](#importing-a-layer-from-another-tool) stores for another tool's contacts (COCaDA, for example). Everything that reads contacts ([`extract_interface_records()`](datasets.md#reshape-into-flat-records), [interface clustering](similarity.md)) takes any `residue_contacts` layer, whoever produced it.

Each row names both residues in both namings (`chain_id_1`, `seq_id_1`, `auth_chain_id_1`, `auth_seq_id_1`, `ins_code_1`, `comp_id_1`, and the same with `_2`), plus `contact_type` and `distance` when known. Ligand residues have no `seq_id`.

On the command line, `pandora annotate compute` computes one layer for every structure of a step folder and writes the same layer folder as `annotate import`: `layer.<fmt>` with the rows and a `layer.json` step summary whose `parameters` record the settings used.

## Chain-chain interfaces

`annotate_chain_interfaces()` finds every pair of residues on two different polymer chains with at least one pair of selected atoms within `distance_cutoff` angstroms. Side 1 is the chain that sorts first, and `distance` is the closest atom pair of the two residues.

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
    Up to 0.5.6 every atom counted, including hydrogens, so structures that model H atoms got extra contacts. The default is now `atom_set="heavy"`; pass `atom_set="all"` to match the old output.

[`extract_interface_records()`](datasets.md#reshape-into-flat-records) groups the rows into one record per chain pair.

=== "`library`"

    ```python
    from pandora.annotations import annotate_chain_interfaces
    from pandora.datasets import extract_interface_records

    rows = annotate_chain_interfaces(canonical, distance_cutoff=4.0)
    print(len(rows), "contacting residue pairs")
    r = rows[0]
    print(
        r.chain_id_1,
        r.seq_id_1,
        r.comp_id_1,
        "-",
        r.chain_id_2,
        r.seq_id_2,
        r.comp_id_2,
        r.distance,
    )

    records = extract_interface_records(rows, distance_cutoff=4.0, atom_set="heavy")
    print(len(records), "interfaces")
    for record in records[:2]:
        print(record.chain_id_1, record.chain_id_2, record.contact_count)
    ```

    ```text
    135 contacting residue pairs
    A 31 ARG - B 122 PHE 2.961
    5 interfaces
    A B 35
    A C 6
    ```

=== "`cli`"

    ```bash
    pandora annotate compute --layer chain_interfaces --structures-dir curated/ \
      --distance-cutoff 4.0 --atom-set heavy --output-dir annotations/interfaces
    jq -c '.counts, .policy.parameters' annotations/interfaces/layer.json
    ```

    ```text
    computed: 135 rows for 1 entries -> annotations/interfaces
    {"entries":1,"rows":135}
    {"distance_cutoff":4.0,"atom_set":"heavy","polymer_types":null}
    ```

    `--polymer-types` is repeatable (`--polymer-types "polypeptide(L)"`).

### Pair gate, then interface (Pinder-style)

Pinder first keeps chain pairs with backbone atoms within 10 Å, then defines the interface with a tighter cutoff. That's two layers:

```python
gate = annotate_chain_interfaces(
    canonical, distance_cutoff=10.0, atom_set="backbone"
)
gated = {
    (r.entry_id, r.chain_id_1, r.chain_id_2)
    for r in extract_interface_records(gate)
}
kept = [r for r in records if (r.entry_id, r.chain_id_1, r.chain_id_2) in gated]
```

## Ligand contacts

`annotate_ligand_contacts()` finds polymer residues near each ligand. A ligand is a residue of a non-polymer chain; set `include_waters=True` to count waters too. A modified residue inside a polymer chain (selenomethionine, sulfotyrosine, ...) is not a ligand, even when the file writes it as `HETATM`: it counts on the polymer side.

Side 1 is the polymer residue, side 2 the ligand, named by its author IDs (it has no `seq_id`), and `distance` the nearest atom pair. Because ligand rows have no `seq_id` on side 2, interface records skip them.

!!! warning "Changed after 0.6.0: modified residues are no longer ligands"
    Earlier versions counted every non-water `HETATM` residue as a ligand, including modified residues inside polymer chains.

=== "`library`"

    ```python
    from pandora.annotations import annotate_ligand_contacts

    rows = annotate_ligand_contacts(canonical, distance_cutoff=4.0)
    ligands = {(r.chain_id_2, r.auth_seq_id_2) for r in rows}
    print(len(rows), "contacts with", len(ligands), "ligands")
    r = rows[0]
    print(
        r.chain_id_1,
        r.seq_id_1,
        r.comp_id_1,
        "-",
        r.auth_chain_id_2,
        r.auth_seq_id_2,
        r.comp_id_2,
        r.distance,
    )
    ```

    ```text
    74 contacts with 4 ligands
    A 32 MET - A 142 HEM 3.828
    ```

=== "`cli`"

    ```bash
    pandora annotate compute --layer ligand_contacts --structures-dir curated/ \
      --distance-cutoff 4.0 --output-dir annotations/ligands
    ```

    ```text
    computed: 74 rows for 1 entries -> annotations/ligands
    ```

    `--include-waters` counts waters. A flag that belongs to the other layer (`--atom-set` here) is an error.

## Sequence similarity

Pairwise sequence identity is no longer an annotation: similarity is its own [layer](similarity.md#the-similarity-layer), computed with MMseqs2 (`pandora similarity search`) or imported from any other source (`pandora similarity import`).

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

## Using layers in a dataset

Pandora doesn't merge layers into the record tables. Every layer is its own table, keyed by the same columns as the records (`entry_id`, `chain_id`, `seq_id`; `chain_id_1`/`chain_id_2` for interfaces), so you join the layers you want with pandas, polars or DuckDB. Residue records also carry the author's keys (`auth_chain_id`, `auth_seq_id`, `ins_code`), so either naming works.

### Contact counts

`count_residue_contacts()` and `count_chain_contacts()` turn any `residue_contacts` layer, Pandora's own or imported, into value layers:

| Function | Layer | `value` |
|---|---|---|
| `count_residue_contacts()` | `residue_values` | distinct partner residues of the residue |
| `count_chain_contacts()` | `chain_values` | distinct residues of the chain in contact |

A pair given several times (for example with several contact types), in either order, counts once. Residues and chains without contacts get no row, so fill in 0 after the join.

### Joining

There is no join function; these are plain pandas joins. The example uses `1a2c` (thrombin light chain L and heavy chain H, with a hirudin peptide I and the inhibitor aeruginosin J), whose light chain has insertion codes.

=== "`library`"

    ```python
    import pandas as pd

    from pandora.annotations import (
        annotate_chain_interfaces,
        count_chain_contacts,
        count_residue_contacts,
        resolve_layer_rows,
    )
    from pandora.datasets import (
        extract_chain_records,
        extract_interface_records,
        extract_residue_records,
    )
    from pandora.parsing import mmcif_to_structure

    structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/a2/1a2c.cif")
    contacts = annotate_chain_interfaces(structure)


    def table(rows, exclude=()):
        return pd.DataFrame([r.model_dump(exclude=set(exclude)) for r in rows])


    # Residue counts onto residue records, on Pandora's keys
    residues = table(extract_residue_records(structure), exclude={"atoms"})
    counts = table(count_residue_contacts(contacts))
    keys = ["entry_id", "chain_id", "seq_id"]
    joined = residues.merge(counts[[*keys, "value"]], on=keys, how="left")
    joined = joined.fillna({"value": 0}).astype({"value": int})
    print(
        joined.loc[
            joined.value > 0,
            ["chain_id", "seq_id", "auth_seq_id", "ins_code", "value"],
        ].head(3)
    )

    # ... or on the author's keys
    auth = ["entry_id", "auth_chain_id", "auth_seq_id", "ins_code"]
    by_auth = residues.merge(counts[[*auth, "value"]], on=auth, how="left")
    print((by_auth.value.fillna(0) == joined.value).all())

    # Chain counts onto chain records
    chain_counts = table(count_chain_contacts(contacts))
    keys = ["entry_id", "chain_id"]
    chains = table(extract_chain_records(structure)).merge(
        chain_counts[[*keys, "value"]], on=keys, how="left"
    )
    chains = chains.fillna({"value": 0}).astype({"value": int})
    print(chains[["chain_id", "auth_chain_id", "residue_count", "value"]])

    # An entry total
    print(chain_counts.groupby("entry_id")["value"].sum())

    # Interface labels onto interface records, in either chain order
    labels, _ = resolve_layer_rows(
        [
            {
                "entry_id": "1a2c",
                "chain_id_1": "B",
                "chain_id_2": "A",
                "label": "light_heavy",
            }
        ],
        "interface_labels",
        {structure.entry_id: structure},
    )
    labels = pd.DataFrame(labels)
    interfaces = table(extract_interface_records(contacts))
    pair = ["chain_id_1", "chain_id_2"]
    for df in (interfaces, labels):
        df[pair] = pd.DataFrame(
            df[pair].apply(sorted, axis=1).tolist(), index=df.index
        )
    keys = ["entry_id", *pair]
    interfaces = interfaces.merge(labels[[*keys, "label"]], on=keys, how="left")
    print(interfaces[[*pair, "contact_count", "label"]])
    ```

    ```text
      chain_id  seq_id auth_seq_id ins_code  value
    0        A       1           1        H      4
    1        A       2           1        G      2
    2        A       3           1        F      2
    True
      chain_id auth_chain_id  residue_count  value
    0        A             L             36     27
    1        B             H            252     65
    2        C             I             10      8
    3        D             J              4      4
    entry_id
    1A2C    104
    Name: value, dtype: int64
      chain_id_1 chain_id_2  contact_count        label
    0          A          B             61  light_heavy
    1          B          C             21          NaN
    2          B          D             22          NaN
    ```

=== "`cli`"

    `pandora annotate count-contacts` writes a count layer folder (`layer.<fmt>`, `layer.json`) from a contact layer folder, one level per run:

    ```bash
    pandora annotate compute --layer chain_interfaces --structures-dir curated/ \
      --output-dir annotations/interfaces
    pandora annotate count-contacts --contacts annotations/interfaces \
      --level residue --output-dir annotations/residue_contacts
    pandora annotate count-contacts --contacts annotations/interfaces \
      --level chain --output-dir annotations/chain_contacts
    jq -c '.counts, .policy.template, .policy.parameters' annotations/chain_contacts/layer.json
    ```

    ```text
    computed: 120 rows for 1 entries -> annotations/interfaces
    counted: 104 residue rows from 120 contacts -> annotations/residue_contacts
    counted: 4 chain rows from 120 contacts -> annotations/chain_contacts
    {"rows_in":120,"rows_out":4}
    "chain_values"
    {"level":"chain","contacts_layer":"chain_interfaces"}
    ```

    Then join the layer onto `curated/chains.parquet` (written by `pandora curate`):

    ```python
    import pandas as pd

    keys = ["entry_id", "chain_id"]
    chains = pd.read_parquet("curated/chains.parquet")
    counts = pd.read_parquet("annotations/chain_contacts/layer.parquet")
    chains = chains.merge(counts[[*keys, "value"]], on=keys, how="left")
    chains = chains.fillna({"value": 0}).astype({"value": int})
    ```

    No command writes residue records; build them with `extract_residue_records()` as in the library tab.

Per-entry values from another tool (an R-free, a clashscore) import as a `free_form` layer, which needs only `entry_id`, and join on `entry_id`.
