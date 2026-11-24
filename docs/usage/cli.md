# CLI

The `pandora` command groups its subcommands by the kind of data they write: `entries`, `structures`, `annotations`, `similarity`, `clusters`, `splits` and `pairs`. Every example below was run against `datasets/dev/mmcif/` and a handful of small fixture entries; run `pandora <group> <command> -h` for the full argument list, or see `pandora/cli/README.md` for the input and output conventions.

```text
pandora entries      fetch | scan
pandora structures   canonicalise | curate | drop-contained | convert
pandora annotations  contacts | ligand-contacts | contact-counts | interface-contacts | dimer-length | interface-gaps | filter | import
pandora similarity   mmseqs | foldseek | interface-overlap | import
pandora clusters     build | representatives
pandora splits       assign
pandora pairs        positives | negatives | assemble
```

## Conventions {: #conventions }

Every command reads as `pandora <group> <command> MAIN_INPUT... OUT [options]`: what the command processes comes first, then where the result goes, and every other input is a named flag (`--structures`, `--contacts`, `--labels`, `--clusters`, `--similarity`, ...).

`OUT` is always a folder. It holds the data under a fixed name plus a `summary.json` (the step name, policy, Pandora and gemmi versions, input and counts). The one exception is `structures convert`, which writes exactly the file you name. A command refuses an `OUT` that already holds its outputs; `--overwrite` replaces them (other files in the folder are left alone). `entries fetch` is the exception: it adds to the folder it is given, which is also its download cache.

Every input accepts the folder a Pandora command wrote or the bare file in it (`clusters.json`, `layer.parquet`, a `.parquet` table, ...). Nothing follows `previous_summary` to find a missing input: when a command needs the structures behind a layer, you pass `--structures`.

Structures come in two kinds. The entry points (`structures canonicalise`, `structures curate`, `annotations contacts`, `annotations ligand-contacts` and the `--structures` of `annotations interface-gaps`) accept a folder of `.cif`/`.cif.gz` files (flat or in the PDB archive layout), an `entries.txt` list or a single file. Every other structures input needs the structures folder written by `structures curate`, with its `entries` and `chains` tables.

Tables are Parquet by default; `--table-format jsonl` writes JSON Lines instead (no `export` extra needed, readable with `jq`). A table with no rows is not written. Warnings go to stderr, and a bad input exits with status 1 and a one-line `<group> <command>: ...` message.

## entries fetch {: #entries-fetch }

Download raw mmCIF files (with on-disk caching):

```sh
pandora entries fetch 104m 112m 118l 138l 1ayi raw/
# fetched 5/5 entries -> raw/
```

IDs can come from the command line, from `--ids FILE` (one per line, `#` comments), or both. `raw/` then holds the mmCIF files, `entries.txt`, `ingestion_provenance.json` and `summary.json`. Fetching into a folder that already holds a fetch downloads only the missing files and merges the lists.

## entries scan {: #entries-scan }

List mmCIF files you already have on disk, without downloading, copying or opening them:

```sh
pandora entries scan datasets/dev/mmcif scan/ --path-template '{mid}/{id}.cif' --ids ids.txt
# scanned 5 entries -> scan/
```

`scan/entries.txt` has one absolute path per line; pass it to `structures canonicalise` or `structures curate`. `--glob` lists files by pattern instead of a template, and `--ids` works with both.

## structures canonicalise {: #structures-canonicalise }

Parse and canonicalise every mmCIF file of a folder (flat or PDB layout), an `entries.txt` list or a single file against a policy YAML (see [Policies](../reference/policies.md)):

```sh
pandora structures canonicalise raw/ canonical/ --policy datasets/canonicalisation.yaml --table-format jsonl
# canonicalised 5 structures -> canonical/
```

`canonical/` then holds `structures/<mid>/<entry>.cif` (the PDB archive layout), `entries.jsonl` (per entry: structure path, policy, which transforms ran), `chain_copies.jsonl` (symmetry-copy chains and their source, which mmCIF has no field for), `mappings.jsonl` (one line per entry: original-to-canonical chain, residue, assembly, entity, altloc and modified-residue mappings), `failures.jsonl` (inputs that raised; the run carries on) and `summary.json`. Output files are named by the parsed entry ID, whatever the input file is called, and a second input with the same entry ID is skipped and recorded in `failures.jsonl` as `DuplicateEntry`. When the canonicalisation policy sets [`max_atoms`](../reference/policies.md#max_atoms), the command estimates each entry's atom count before parsing it and records an entry above the limit as `TooLarge` (counted in `entries_too_large`). [Canonicalisation](canonicalisation.md) shows how to read the outputs with `jq`.

## structures curate {: #structures-curate }

Apply quality/organism/content rules from a `DatasetCurationPolicy` YAML to a `structures canonicalise` folder (or any structures input):

```sh
pandora structures curate canonical/ curated/ --policy curation.yaml --table-format jsonl
# curated: 5 retained, 0 excluded -> curated/
```

`curated/` then holds `structures/` (kept entries), `entries.jsonl` (per kept entry: structure path, method, resolution, chain count, policy), `chains.jsonl` (per kept chain: sequence, completeness, method, resolution, `copy_of` for symmetry copies), `exclusions.jsonl` (per dropped entry or chain, with its reason), `failures.jsonl` and `summary.json`. It reads `chain_copies` from a canonicalise folder when present. Repeated entry IDs are dropped here, as in `structures canonicalise`.

where `curation.yaml` is a minimal policy:

```yaml
policy_id: curation-default
policy_name: Default
policy_version: 1.0.0
```

## structures drop-contained {: #structures-drop-contained }

Remove entries contained in another, by chain counts per cluster of `--clusters` (the `clusters.json` of `clusters build` over the same folder's chains, or the folder holding it). It needs a curate folder with tables and writes a structures folder with a `contained` table and a `summary.json` ([Contained biounits](datasets.md#contained-biounits)):

```sh
pandora structures drop-contained work/curation work/uncontained --clusters clusters/seq90
# 29 retained, 2 contained -> work/uncontained
```

## structures convert {: #structures-convert }

Convert one mmCIF file to mmCIF or JSON. It takes exactly one input file and writes exactly the file you name, with no folder or summary:

```sh
pandora structures convert curated/structures/04/104m.cif 104m.json
# converted -> 104m.json
```

The output format is inferred from the output's suffix.

## annotations contacts and ligand-contacts {: #annotations-contacts }

One of Pandora's contact layers for every structure of a structures input, written as a layer folder (see [Contacts are one layer](annotation.md#contacts-are-one-layer)):

```sh
pandora annotations contacts curated/ annotations/interfaces --distance-cutoff 4.0 --atom-set heavy
# computed: 135 rows for 1 entries -> annotations/interfaces
```

It writes `layer.<fmt>` (`residue_contacts` rows) and `summary.json`, whose `parameters` record the settings used. `contacts` takes `--distance-cutoff`, `--atom-set` and `--polymer-types` (repeatable); `ligand-contacts` takes `--distance-cutoff` and `--include-waters`.

## annotations import {: #annotations-import }

Import a table produced by an external tool, converted into one of Pandora's layer templates, and link it to the structures the tool ran on (see [Importing a layer from another tool](annotation.md#importing-a-layer-from-another-tool)):

```sh
pandora annotations import cocada_contacts.jsonl annotations/cocada --structures curated/ --description cocada.yaml
# imported: 1 rows, 1 unmatched -> annotations/cocada
```

It writes `layer.<fmt>`, `unmatched.<fmt>`, `summary.json` and, with `--source-files`, a copy of the tool's files in `source/`.

## annotations contact-counts {: #annotations-contact-counts }

Count a `residue_contacts` layer's distinct partners per residue (`--level residue`, a `residue_values` layer) or residues in contact per chain (`--level chain`, a `chain_values` layer), see [Contact counts](annotation.md#contact-counts):

```sh
pandora annotations contact-counts annotations/interfaces annotations/chain_contacts --level chain
# counted: 4 chain rows from 120 contacts -> annotations/chain_contacts
```

It writes `layer.<fmt>` and `summary.json`. A folder that isn't a `residue_contacts` layer is an error.

## annotations interface-contacts, dimer-length and interface-gaps {: #annotations-interface-labels }

Each writes an `interface_labels` layer with one row per chain pair; see [Labelling interfaces by contacts and size](annotation.md#labelling-interfaces-by-contacts-and-size) and [Labelling interfaces with a hole](annotation.md#labelling-interfaces-with-a-hole).

```sh
pandora annotations interface-contacts annotations/interfaces annotations/few_contacts --min-residues-per-chain 5
pandora annotations dimer-length annotations/interfaces annotations/length --structures curated/ --max-total-residues 1000
pandora annotations interface-gaps annotations/interfaces annotations/gaps --structures curated/
```

`interface-contacts` labels `few_contacts` when either chain has fewer residues in contact than the threshold, else `ok`. `dimer-length` labels `too_large` when the two chains' sequence lengths add up to more than the limit (it reads `seqres_length` from the `chains` table of `--structures`), else `ok`. `interface-gaps` labels `interface_gap`, else `unknown`, else `ok`.

## annotations filter {: #annotations-filter }

Drop the interfaces that `interface_labels` layers mark with a `--drop` label from a `residue_contacts` layer, writing a new contact layer; see [Dropping interfaces by label](annotation.md#dropping-interfaces-by-label):

```sh
pandora annotations filter annotations/interfaces annotations/interfaces_strong --labels annotations/hbond_verdict --drop low_hbond
# filtered: 1 interfaces dropped, 2 kept, 0 unlabelled -> annotations/interfaces_strong
```

`--labels` and `--drop` are repeatable: an interface is dropped when any of the layers gives it a `--drop` label, so one pass replaces several. It writes `layer.<fmt>` and `summary.json`, whose `counts` give rows in and out, interfaces in, dropped and unlabelled, and labelled interfaces not in the contact layer. Interfaces without a label are kept. It is an error, before anything is written, when a folder isn't a layer of the right template, `OUT` is an input folder, or a `--drop` label is on no interface of any layer. When the `summary.json` files record different input structures, the command warns.

## similarity mmseqs and foldseek {: #similarity-search }

All-vs-all search (requires `mmseqs`/`foldseek` on `PATH`) over a curate folder's chains, written as a [similarity layer](similarity.md#the-similarity-layer) folder: `layer.parquet`, `unmatched.parquet` and `summary.json`. MMseqs2 reads the chain sequences from the `chains` table, leaving symmetry copies out; Foldseek searches per-chain files that Pandora cuts from the structures itself, or the `{entry}_{chain}.cif` files of `--chain-files`. `--chains` limits the search to a list of chain IDs or to the representatives of a `clusters representatives` folder. Similarity layers are Parquet and need the `export` extra.

```sh
pandora similarity mmseqs curated/ similarity/seq --max-seqs 300
# searched: 9 hits, 0 unmatched -> similarity/seq
```

`--max-seqs` (and `--exhaustive-search` for Foldseek) control how many hits each query keeps; see the [`max_seqs` warning](similarity.md#sequence-similarity-mmseqs2).

## similarity import {: #similarity-import }

Bring in an edge list computed any other way (Parquet, or TSV with a header row), described by a YAML with `template: similarity` and `items: entry | chain | interface`. Entry and chain items are checked against `--structures`, interface items against a contact layer given with `--contacts`. See [Importing an edge list](similarity.md#importing-an-edge-list).

```sh
pandora similarity import r_graph_edges.tsv similarity/r_graph --description r_graph.yaml --structures curated/
# imported: 2 rows, 1 unmatched -> similarity/r_graph
```

## similarity interface-overlap {: #similarity-interface-overlap }

Score every pair of a contact layer's interfaces whose chains share `--chain-clusters` (a `clusters.json` from `clusters build`) and are at least `--min-identity` (default 0.9) identical, by per-side Jaccard of their interface residues; `score` is the smaller side. Every compared pair is written. See [Comparing interfaces](similarity.md#comparing-interfaces).

```sh
pandora similarity interface-overlap annotations/interfaces similarity/iface_jaccard --chain-clusters clusters/seq90 --structures curated/
# scored: 6 interface pairs (4 zero, 0 below identity) from 5 interfaces -> similarity/iface_jaccard
```

## clusters build {: #clusters-build }

Connected-component clustering of one or more similarity layers, filtered by a [`HitFilter`](similarity.md#filtering-hits-hitfilter). The items come from `--structures` (entries or chains) or `--contacts` (interfaces), so items with no hits are single-item clusters:

```sh
pandora clusters build similarity/seq clusters/seq90 --structures curated/ --min-score 0.9
# 3 clusters from 4 edges -> clusters/seq90
```

`--hit-filter filter.yaml` loads a full `HitFilter`; `--min-score`/`--min-coverage` override its fields. One of the three is required; there is no default threshold. A similarity layer's item level comes from its `summary.json`; a bare table or a folder without one needs `--similarity-items {entry,chain,interface}`. `--items entry` on a chain layer clusters entries, linking two entries when any of their chains are similar ([Clustering entries from a chain layer](similarity.md#clustering-entries-from-a-chain-layer)). For `min_interface_coverage` on a chain layer, pass a contact layer folder with `--contacts` alongside `--structures` (see [interface coverage](similarity.md#interface-restricted-coverage-ppi-pairs)). The output folder holds `clusters.json` and `summary.json`.

## clusters representatives {: #clusters-representatives }

One representative per cluster, ranked by a `RepresentativePolicy` (default: fewest missing residues, then most hits inside the cluster, then lowest id). Clusters are limited to the chains of `--structures`; symmetry copies are never picked. The `n_hits` key needs `--similarity` and a hit filter ([Representatives](similarity.md#representatives)). With `--contacts` it picks interfaces instead of chains:

```sh
pandora clusters representatives clusters/seq90 clusters/seq90_reps --structures work/uncontained --similarity layers/seq --hit-filter policies/seq90.yaml
# 31 representatives from 51 chains -> clusters/seq90_reps
```

The output folder holds `representatives.json` and `summary.json`.

## splits assign {: #splits-assign }

Leakage-safe train/val/test split of clusters:

```sh
pandora splits assign clusters/seq90 splits/seq90 --pct-train 0.6 --pct-val 0.2 --pct-test 0.2
# split sizes: {'train': 4, 'val': 1, 'test': 0} -> splits/seq90
```

Clusters always go to a split whole. `--labels <layer folder>` takes a `chain_labels` or `interface_labels` layer and places clusters so every label reaches every split where keeping clusters whole allows; `--stratify <layer folder>` splits each label's clusters separately, so every split keeps the label mix (the two flags are mutually exclusive). `summary.json` then records `label_counts`, `missing_labels` and the labels layer; the command prints how many items were unlabelled and warns for each split missing a label (see [Every kind in every split](similarity.md#every-kind-in-every-split)). The output folder holds `splits.json` and `summary.json`.

## pairs {: #pairs }

`pairs positives`, `pairs negatives` and `pairs assemble` build the positive and negative dimer pairs and join them with the splits; see [Rebalancing homomers and heteromers](similarity.md#rebalancing-homomers-and-heteromers), [Sampling negatives](similarity.md#sampling-negatives) and [Splits and the final table](similarity.md#splits-and-the-final-table) for the full commands.

```sh
pandora pairs positives similarity/seq similarity/iface work/positives --contacts annotations/interfaces --structures curated/ --chain-clusters clusters/seq90 --hit-filter policies/dimer.yaml
pandora pairs negatives work/positives work/negatives --representatives clusters/seq90_reps --contacts annotations/any_contacts --seed 0
pandora splits assign work/positives splits/dimers --stratify work/positives/kinds
pandora pairs assemble work/positives work/dataset --negatives work/negatives --splits splits/dimers --representatives clusters/seq90_reps
```

`pairs positives` builds the dimer groups itself, as the connected components of the dimer graph of its similarity inputs, and writes them as `clusters.json` next to the `positives` table, so `splits assign` can take the positives folder directly. These groups match `clusters build --items interface` run on the same similarity inputs without `--structures`: symmetry copies are not linked to their source chain.
