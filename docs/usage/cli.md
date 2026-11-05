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

All-vs-all sequence or structure similarity (requires `mmseqs`/`foldseek` on `PATH`). The tool's hits are kept in `--output`, and the `SimilaritySearch` record that `cluster` reads goes next to it as `<output>.search.json`:

```sh
pandora similarity --input-dir deduped/ --engine mmseqs2 --max-seqs 300 --output mmseqs.tsv
# hits -> mmseqs.tsv; search record -> mmseqs.tsv.search.json
```

`--max-seqs` (and `--exhaustive-search` for Foldseek) control how many hits each query keeps; see the [`max_seqs` warning](similarity.md#sequence-similarity-mmseqs2). `--precomputed-hits <tsv>` wraps a hit file you ran yourself instead of searching ([Precomputed searches](similarity.md#precomputed-searches)).

`--engine foldseek` works the same way. When chaining `similarity
--engine foldseek` into `cluster` below, keep to this exact
`similarity` → `cluster` order on the *same* directory — see
[Keep ids consistent between stages](similarity.md#keep-ids-consistent-between-stages).

## cluster

Connected-component clustering of a search's hits, filtered by a [`HitFilter`](similarity.md#filtering-hits-hitfilter):

```sh
pandora cluster --input-dir deduped/ --search mmseqs.tsv.search.json --min-score 0.9 --output clusters.json
# 3 clusters from 4 edges -> clusters.json
```

`--hit-filter filter.yaml` loads a full `HitFilter`; `--min-score`/`--min-coverage` override its fields. One of the three is required — there is no default threshold.

## partition

Leakage-safe train/val/test split of clusters:

```sh
pandora partition --clusters clusters.json --pct-train 0.6 --pct-val 0.2 --pct-test 0.2 --output splits.json
# split sizes: {'train': 4, 'val': 1, 'test': 0} -> splits.json
```

## annotate

Entry-level and/or pairwise annotation layers:

```sh
pandora annotate --input-dir deduped/ --layers structure_counts ligand_contacts --pairwise --output-dir annotations/
# annotated 5 entries -> annotations/
```

`--layers` defaults to every entry-level annotator; list a subset to
skip the rest. `--pairwise` additionally computes
`annotate_pairwise_sequence_identity()` for every pair in the batch.

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
  --annotations annotations/annotations.json \
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
