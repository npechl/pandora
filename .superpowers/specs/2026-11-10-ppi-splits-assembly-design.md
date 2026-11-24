# Splits and final assembly for PPI datasets

Date: 2026-11-10. Related: #136 (`rebalance`), #137 (`negatives`),
PPI v2.1 benchmark doc (`examples/ppi/benchmarks/ppi-v2.1-benchmark.md`,
"Splits", "Output", limitations B1, E3, F).

## Why

The PPI v2.1 build now stops at two step folders, `work/positives` and
`work/negatives`. The dataset is meant for **training a model (user)**, so
it needs a leakage-safe train/val/test split and one final table. v2.1
itself specifies neither (B1, open decisions). Its only guidance, which
this design follows **(user)**: assign whole Step 6 dimer groups to
splits, and put each negative in its positive's split.

Decisions marked **(user)** came from the user; **(my call)** were
proposed by Claude and accepted, or left to Claude.

## Goals

1. `rebalance` also writes what `partition` needs: the dimer groups
   limited to the positives, and a homomer/heteromer label layer.
2. `partition --stratify`: whole groups to splits, with every split
   keeping roughly the overall homomer:heteromer ratio. **(user: "preserve
   balance")**
3. `pandora assemble`: one final table of positive and negative pairs with
   their split, plus a measured C-chain leakage count. **(user)**
4. `examples/ppi/v2.1/src/08_splits.sh` runs 2 and 3.

## Non-goals

- Split-aware negative sampling (drawing C per split so no C resembles a
  chain of another split). Deferred: build it only if the leakage count
  on the real build is high. **(user)**
- Labels beyond homomer/heteromer (DSSP, interface area, pH: not built).
- Export formats beyond the step-folder tables.

## Anti-goals

- No orchestrator; three separate commands as today.
- Never split a dimer group across splits.
- No silent drops: a positive in no split, or a negative whose positive is
  unknown, is an error naming it.

## 1. `rebalance` writes two more outputs

In its `--output-dir` (alongside `positives`, `dropped`,
`rebalance.json`), always written, no new flags **(my call: no
output-selection flags, per project convention)**:

- `groups.json`: the `--clusters` groups limited to the kept positives
  (`restrict_clusters`), groups with no positive left out. Same
  `SimilarityCluster` list format as `cluster`'s output.
- `kinds/`: an `interface_labels` layer folder that `partition --labels`
  reads: `layer.<fmt>` with one row per positive (`entry_id`,
  `chain_id_1`, `chain_id_2`, `label` = `homomer`/`heteromer`) and
  `layer.json` whose `policy.template` is `interface_labels` (build it the
  way `annotate label-interfaces` does). **(my call: a layer nested in the
  step folder, since `rebalance` derives it; the alternative, a separate
  `layers/` folder, needs an extra output flag)**

`--overwrite` must also clear `groups.json` and `kinds/`.

## 2. `partition --stratify`

```sh
pandora partition --clusters work/positives/groups.json \
  --labels work/positives/kinds --stratify --output splits.json
```

- Needs `--labels` and whole clusters (not with
  `--no-keep-similar-items`); exit otherwise.
- Every cluster must carry a single label; exit naming the first cluster
  with two. (True for `rebalance` output: a homomer-only group keeps one
  homomer, a mixed group keeps only heteromers.)
- Algorithm: run today's whole-cluster placement once **per label**, on
  that label's clusters only, with the same `pct_*` targets computed from
  that label's item count; concatenate the per-label splits. Labels in
  sorted order, so the result is deterministic.
- Provenance: `stratify: bool` on `PartitionProvenance`;
  `label_counts`/`missing_labels` already report the achieved mix. A label
  with fewer clusters than active splits leaves a split without it; that
  is the existing `missing_labels` warning, not an error.
- Library: `partition_dataset(..., stratify=False)`.

## 3. `pandora assemble`

```sh
pandora assemble --positives work/positives --negatives work/negatives \
  --splits splits.json \
  --representatives clusters/seq90/representatives.json \
  --similarity layers/tm --hit-filter policies/tm.yaml \
  --output-dir work/dataset
```

### Output: `work/dataset/`

- `pairs.<fmt>`, one row per pair, positives then negatives, each sorted:
  - `chain_a`, `chain_b`: chain item ids (positive: its two chains;
    negative: A, C).
  - `label`: 1 positive, 0 negative.
  - `kind`: the positive's `homomer`/`heteromer` (negatives inherit their
    positive's).
  - `cluster`: the positive's dimer group (its `group`; negatives inherit).
  - `split`: `train`/`val`/`test`; a negative takes its positive's. **(user)**
  - `positive`: the dimer id (for a positive, itself).
  - `c_leak_positive`: negatives only, true when C is structurally
    similar to a positive's chain in another split (below); null for
    positives.
  - `c_leak_negative`: negatives only, true when C is structurally
    similar to another split's negative C; null for positives.
- `assemble.json`: step summary. Settings: the input paths and the hit
  filter. Counts: pairs per split and label, `negatives_c_leak_positive`,
  `negatives_c_leak_negative`, and
  `positives_unpaired` (positives with no negative). `previous_summary` =
  `negatives.json`.

### The leakage count **(user: measure, don't prevent)**

Map every chain to its representative (`--representatives` members;
TM edge ends outside it kept as-is, like `negatives`). For each split,
collect two sets: the representatives of its positives' chains, and the
representatives of its negatives' C chains. A negative (A, C) in split s
has `c_leak_positive = true` when C's representative is TM-linked (after
the hit filter) to a positive-chain representative of another split, and
`c_leak_negative = true` when it is TM-linked to another split's C. Two
separate flags, so neither kind of leak inflates the other. **(my call)**

With the `--representatives` file `negatives` used, C can never *equal* a
positive's representative or another C (the pool excludes them and C is
drawn once), so only TM links count. If `--representatives` or the TM
layer path differs from the one recorded in `negatives.json`, print a
warning (not an error): the counts are then not comparable. **(my call)**

The `cluster` column comes from `PositiveDimer.group` (the original
group's first member, possibly a dropped dimer), not from `groups.json`,
whose first member can differ.

### Errors

- A positive in no split of `--splits`, or in two: error naming it.
- A split item that is not a positive: error naming it.
- A negative whose `positive` is not in the positives table: error.
- A chain that must be mapped and is in no `members`: error (as
  `negatives`).

## 4. Build script and docs

- `08_splits.sh`: the `partition` and `assemble` commands above, default
  ratios 0.6/0.2/0.2 **(my call; v2.1 gives none)**, then `jq .counts`.
- Usage docs: `partition --stratify` in the partitioning section; a
  "Final table" part in the "Positives and negatives (PPI)" section; real
  output, library and CLI tabs. CLI README. Benchmark doc: "With
  Pandora" notes under Splits and Output, and limitation F gets the
  leakage-count note.

## Tests (local fixtures)

- `rebalance` writes `groups.json` (positives only) and a `kinds` layer
  that `partition --labels` accepts.
- `partition_dataset(stratify=True)`: two labels with 10 and 5
  single-item clusters at 0.6/0.2/0.2 give 6/2/2 and 3/1/1; a two-label
  cluster raises; `stratify` without labels raises; CLI flag checks.
- `assemble`: rows, splits inherited by negatives, counts;
  `c_leak_positive` / `c_leak_negative` each true for a hand-made TM edge
  across splits (to a positive chain / to another C) and false within a
  split; each error above.
