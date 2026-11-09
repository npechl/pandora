# CLI

The `pandora` command wraps each component as a subcommand. Every
example below was actually run against `datasets/dev/mmcif/` and a
handful of small fixture entries; run `pandora <subcommand> -h` for the
full argument list, or see `pandora/cli/README.md` for the directory/
JSON conventions between stages.

Structure-transforming stages (`canonicalise`, `curate`, `dedup`) write mmCIF files in the PDB archive layout (`04/104m.cif`). Every stage that reads structures accepts a flat folder, that layout, or a step folder. `canonicalise` and `curate` write a step folder: the structures in `structures/<mid>/<entry>.cif` (the PDB archive layout), one row per entry in `entries.<fmt>`, `failures.jsonl` for inputs that raised (the run carries on), and a step summary (`canonicalisation.json` / `curation.json`) recording the policy, Pandora and gemmi versions, input and counts. Tables are Parquet by default; `--table-format jsonl` writes JSON Lines instead (no `export` extra needed, readable with `jq`). A table with no rows is not written. Both refuse an `--output-dir` that already holds step outputs; `--overwrite` replaces them (other files in the folder are left alone). Output files are named by the parsed entry ID, whatever the input file is called (`.cif.gz` inputs are read directly). A second input with the same entry ID is skipped and recorded in `failures.jsonl` as `DuplicateEntry`. For other folder layouts (e.g. the PDBe release tree, `entry/00/100d/clean_mmcif/100d_updated.cif.gz`), list the files with `pandora ingest` and pass its `entries.txt` to `canonicalise --entries`. `ingest` writes `entries.txt` (one absolute path per line) and an `ingestion.json` step summary, without opening any file; `canonicalise --entries` records that summary as its `previous_summary`. `fetch` also writes an `entries.txt` of the files it saved. Other provenance-producing stages write one JSON file per stage into their `--output`/`--output-dir`, which `manifest` reads back in.

## fetch

Download raw mmCIF files (with on-disk caching):

```sh
pandora fetch 104m 112m 118l 138l 1ayi --output-dir raw/
# fetched 5/5 entries -> raw/
```

`raw/ingestion_provenance.json` is written alongside the files — keep
it, `manifest` needs it later if you want `reproduce` to work.

## canonicalise

Parse and canonicalise every mmCIF file in a folder (`--input-dir`, flat or PDB layout) or listed in a text file (`--entries`, one path per line) against a policy YAML (see [Policies](../reference/policies.md)):

```sh
pandora canonicalise --input-dir raw/ --policy datasets/canonicalisation.yaml --output-dir canonical/ --table-format jsonl
# canonicalised 5 structures -> canonical/
```

`canonical/` then holds `structures/`, `entries.jsonl` (per entry: structure path, policy, which transforms ran), `chain_copies.jsonl` (symmetry-copy chains and their source, which mmCIF has no field for), `mappings.jsonl` (one line per entry: original-to-canonical chain, residue, assembly, entity, altloc and modified-residue mappings), `failures.jsonl` and `canonicalisation.json`. [Canonicalisation](canonicalisation.md) shows how to read them with `jq`.

## curate

Apply quality/organism/content rules from a `DatasetCurationPolicy` YAML to a `canonicalise` step folder (or any folder of mmCIF files):

```sh
pandora curate --input-dir canonical/ --policy curation.yaml --output-dir curated/ --table-format jsonl
# curated: 5 retained, 0 excluded -> curated/
```

`curated/` then holds `structures/` (kept entries), `entries.jsonl` (per kept entry: structure path, method, resolution, chain count, policy), `chains.jsonl` (per kept chain: sequence, completeness, method, resolution, `copy_of` for symmetry copies), `exclusions.jsonl` (per dropped entry or chain, with its reason), `failures.jsonl` and `curation.json`, whose `previous_summary` points to `canonical/canonicalisation.json`.

where `curation.yaml` is a minimal policy:

```yaml
policy_id: curation-default
policy_name: Default
policy_version: 1.0.0
```

## dedup

Remove structures sharing the same `entry_id`:

```sh
pandora dedup --input-dir curated/ --output-dir deduped/
# dedup: 5 retained, 0 removed -> deduped/
```

## similarity

`similarity search` runs an all-vs-all search (requires `mmseqs`/`foldseek` on `PATH`) over a step folder's chains and writes the hits as a [similarity layer](similarity.md#the-similarity-layer) folder: `layer.parquet`, `unmatched.parquet` and `layer.json`. MMseqs2 reads the chain sequences from `--structures-dir`'s `chains` table (written by `curate`), leaving symmetry copies out; Foldseek searches the per-chain files in `--input-dir`, named `{entry}_{chain}.cif`. Similarity layers are Parquet and need the `export` extra.

```sh
pandora similarity search --engine mmseqs2 --structures-dir curated/ --max-seqs 300 --output-dir similarity/seq
# searched: 9 hits, 0 unmatched -> similarity/seq
```

`--max-seqs` (and `--exhaustive-search` for Foldseek) control how many hits each query keeps; see the [`max_seqs` warning](similarity.md#sequence-similarity-mmseqs2).

`similarity import` brings in an edge list computed any other way (Parquet, or TSV with a header row), described by a YAML with `template: similarity` and `items: entry | chain | interface`. Entry and chain items are checked against `--structures-dir`'s tables, interface items against a contact layer given with `--contacts`. See [Importing an edge list](similarity.md#importing-an-edge-list).

```sh
pandora similarity import --table r_graph_edges.tsv --description r_graph.yaml --structures-dir curated/ --output-dir similarity/r_graph
# imported: 2 rows, 1 unmatched -> similarity/r_graph
```

Both refuse an `--output-dir` that already holds a layer, unless `--overwrite` is given.

## cluster

Connected-component clustering of a similarity layer, filtered by a [`HitFilter`](similarity.md#filtering-hits-hitfilter). The items come from `--structures-dir` (entries or chains) or `--contacts` (interfaces), so items with no hits are single-item clusters:

```sh
pandora cluster --similarity similarity/seq --structures-dir curated/ --min-score 0.9 --output clusters.json
# 3 clusters from 4 edges -> clusters.json
```

`--hit-filter filter.yaml` loads a full `HitFilter`; `--min-score`/`--min-coverage` override its fields. One of the three is required; there is no default threshold. `--items entry` on a chain layer clusters entries, linking two entries when any of their chains are similar ([Clustering entries from a chain layer](similarity.md#clustering-entries-from-a-chain-layer)). For `min_interface_coverage` on a chain layer, pass a contact layer folder with `--contacts` alongside `--structures-dir` (see [interface coverage](similarity.md#interface-restricted-coverage-ppi-pairs)). `--pairs pairs.json` also writes per-pair cluster keys to `cluster_pairs.json`.

## partition

Leakage-safe train/val/test split of clusters:

```sh
pandora partition --clusters clusters.json --pct-train 0.6 --pct-val 0.2 --pct-test 0.2 --output splits.json
# split sizes: {'train': 4, 'val': 1, 'test': 0} -> splits.json
```

`--labels <layer folder>` takes a `chain_labels` or `interface_labels` layer and places clusters so every label reaches every split where keeping clusters whole allows. `partition_provenance.json` then records `label_counts`, `missing_labels` and `labels_layer`; the command prints how many items were unlabelled and warns for each split missing a label (see [Every kind in every split](similarity.md#every-kind-in-every-split)).

## annotate compute

One of Pandora's contact layers (`chain_interfaces` or `ligand_contacts`) for every structure of a step folder, written as a layer folder (see [Contacts are one layer](annotation.md#contacts-are-one-layer)):

```sh
pandora annotate compute --layer chain_interfaces --structures-dir curated/ --distance-cutoff 4.0 --atom-set heavy --output-dir annotations/interfaces
# computed: 135 rows for 1 entries -> annotations/interfaces
```

It writes `layer.<fmt>` (`residue_contacts` rows) and `layer.json`, whose `parameters` record the settings used. Settings: `--distance-cutoff` (both layers), `--atom-set` and `--polymer-types` (repeatable; `chain_interfaces`), `--include-waters` (`ligand_contacts`). A flag that doesn't apply to the chosen `--layer` is an error. `--overwrite` and `--table-format` work as on `annotate import`.

## annotate import

Import a table produced by an external tool, converted into one of Pandora's layer templates, and link it to the structures the tool ran on (see [Importing a layer from another tool](annotation.md#importing-a-layer-from-another-tool)):

```sh
pandora annotate import --structures-dir curated/ --table cocada_contacts.jsonl --description cocada.yaml --output-dir annotations/cocada
# imported: 1 rows, 1 unmatched -> annotations/cocada
```

It writes `layer.<fmt>`, `unmatched.<fmt>`, `layer.json` and, with `--source-files`, a copy of the tool's files in `source/`.

## annotate count-contacts

Count a `residue_contacts` layer's distinct partners per residue (`--level residue`, a `residue_values` layer) or residues in contact per chain (`--level chain`, a `chain_values` layer), see [Contact counts](annotation.md#contact-counts):

```sh
pandora annotate count-contacts --contacts annotations/interfaces --level chain --output-dir annotations/chain_contacts
# counted: 4 chain rows from 120 contacts -> annotations/chain_contacts
```

It writes `layer.<fmt>` and `layer.json`, whose `previous_summary` is the contact layer's `layer.json`. A folder that isn't a `residue_contacts` layer is an error. `--overwrite` and `--table-format` work as on `annotate import`.

## annotate filter-contacts

Drop the interfaces an `interface_labels` layer marks with a `--drop` label (repeatable) from a `residue_contacts` layer, writing a new contact layer; see [Dropping interfaces by label](annotation.md#dropping-interfaces-by-label):

```sh
pandora annotate filter-contacts --contacts annotations/interfaces --labels annotations/hbond_verdict --drop low_hbond --output-dir annotations/interfaces_strong
# filtered: 1 interfaces dropped, 2 kept, 0 unlabelled -> annotations/interfaces_strong
```

It writes `layer.<fmt>` and `layer.json`, whose `previous_summary` is the contact layer's `layer.json` and whose `counts` give rows in and out, interfaces in, dropped and unlabelled, and labelled interfaces not in the contact layer. Interfaces without a label are kept. It is an error, before anything is written, when a folder isn't a layer of the right template, `--output-dir` is an input folder, a `--drop` label is on no interface, or the two layers come from different structures folders (a layer whose structures folder can't be traced gives a warning). `--overwrite` and `--table-format` work as on `annotate import`.

## manifest

Assemble a `DatasetManifest` from every prior stage's output:

```sh
pandora manifest \
  --dataset-id dev-cli-demo \
  --dataset-name "CLI demo dataset" \
  --dataset-version 1.0.0 \
  --structures-dir deduped/ \
  --canonicalisation-policy datasets/canonicalisation.yaml \
  --curation-policy curation.yaml \
  --ingestion-provenance raw/ingestion_provenance.json \
  --canonicalisation-dir canonical/ \
  --curation-dir curated/ \
  --cluster-provenance cluster_provenance.json \
  --partition-provenance partition_provenance.json \
  --splits splits.json \
  --output manifest.json
# manifest for 5 structures -> manifest.json
```

Every `--*-provenance`/`--*-policy`/`--*-dir` flag is optional — only pass what the earlier stages actually produced. `--canonicalisation-dir` and `--curation-dir` take step folders and fill each structure's canonicalisation and curation provenance from their `entries` tables; `--exclusions` takes `exclusions` tables (`.jsonl`, `.json` or `.parquet`). `--ingestion-provenance` is the
one flag worth not skipping: without it, `reproduce` (below) can't
re-fetch that structure later.

## export

Convert one mmCIF file to mmCIF or JSON:

```sh
pandora export --input deduped/04/104m.cif --output 104m.json
# exported -> 104m.json
```

The output format is inferred from `--output`'s suffix.

## reproduce

Replay a `DatasetManifest` from scratch — re-fetches every structure,
then re-runs canonicalisation/curation/dedup/similarity/clustering/
partition/annotation exactly as the manifest recorded them:

```sh
pandora reproduce --manifest manifest.json --output-dir reproduced/
# reproduced 5 structures -> reproduced/
```

Requires every structure's manifest entry to carry `ingestion` provenance (see `manifest` above), and — if the manifest recorded clustering — `ClusteringProvenance.search` and `.hit_filter` to be set (a precomputed hit file is re-used only if it still exists; interface coverage filters can't be reproduced automatically). It's a best-effort re-run, not a guaranteed byte-identical rebuild; diff `reproduced/reproduced_manifest.json` against the input to see what changed.
