# CLI

The `pandora` command wraps the whole pipeline as subcommands. Every
example below was actually run against `datasets/dev/mmcif/` and a
handful of small fixture entries; run `pandora <subcommand> -h` for the
full argument list, or see `pandora/cli/README.md` for the directory/
JSON conventions between stages.

Structure-transforming stages (`canonicalise`, `curate`, `dedup`) read
and write directories of mmCIF files, one per entry. Provenance-
producing stages write one JSON file per stage into their
`--output`/`--output-dir`, which downstream stages — `manifest` in
particular — read back in.

## fetch

Download raw mmCIF files (with on-disk caching):

```sh
pandora fetch 104m 112m 118l 138l 1ayi --output-dir raw/
# fetched 5/5 entries -> raw/
```

`raw/ingestion_provenance.json` is written alongside the files — keep
it, `manifest` needs it later if you want `reproduce` to work.

## canonicalise

Parse and canonicalise every `*.cif` in a directory against a policy
YAML (see [Policies](../reference/policies.md)):

```sh
pandora canonicalise --input-dir raw/ --policy datasets/canonicalisation.yaml --output-dir canonical/
# canonicalised 5 structures -> canonical/
```

## curate

Apply quality/organism/content rules from a `DatasetCurationPolicy`
YAML:

```sh
pandora curate --input-dir canonical/ --policy curation.yaml --output-dir curated/
# curated: 5 retained, 0 excluded -> curated/
```

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

All-vs-all sequence or structure similarity (requires `mmseqs`/
`foldseek` on `PATH`):

```sh
pandora similarity --input-dir deduped/ --engine mmseqs2 --output relationships.json
# computed 2 relationships -> relationships.json
```

`--engine foldseek` works the same way. When chaining `similarity
--engine foldseek` into `cluster` below, keep to this exact
`similarity` → `cluster` order on the *same* directory — see
[Keep ids consistent between stages](similarity.md#keep-ids-consistent-between-stages).

## cluster

Connected-component clustering of a relationship network:

```sh
pandora cluster --input-dir deduped/ --relationships relationships.json --threshold 0.9 --output clusters.json
# 3 clusters at threshold=0.9 -> clusters.json
```

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
  --canonicalisation-provenance canonical/canonicalisation_provenance.json \
  --cluster-provenance cluster_provenance.json \
  --partition-provenance partition_provenance.json \
  --splits splits.json \
  --annotations annotations/annotations.json \
  --output manifest.json
# manifest for 5 structures -> manifest.json
```

Every `--*-provenance`/`--*-policy` flag is optional — only pass what
the earlier stages actually produced. `--ingestion-provenance` is the
one flag worth not skipping: without it, `reproduce` (below) can't
re-fetch that structure later.

## export

Convert one mmCIF file to mmCIF or JSON:

```sh
pandora export --input deduped/104m.cif --output 104m.json
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

Requires every structure's manifest entry to carry `ingestion`
provenance (see `manifest` above), and — if the manifest recorded
clustering — `ClusteringProvenance.similarity_method` to be set. It's
a best-effort re-run, not a guaranteed byte-identical rebuild; diff
`reproduced/reproduced_manifest.json` against the input to see what
changed.
