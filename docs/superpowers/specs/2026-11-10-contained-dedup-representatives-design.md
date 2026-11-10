# Containment dedup and cluster representatives (#138)

Date: 2026-11-10. Issue: #138. Related: #141 (interface representatives),
#139 (`similarity search --chains`), #108 (entry-id dedup).

## Why

PPI v2.1 needs two selection steps that only existed in the old part B
script, which is gone:

- **Step 1 filter 1:** drop biounits contained in another biounit.
- **Step 2.4:** pick one representative chain per 90% sequence cluster.

`examples/ppi/v2.1/src/02_sequence_clusters.sh` has a "NOT YET IN THE
CLI (#138)" block for both. Until they exist, every later step uses
every curated chain.

How the outputs are used, following the v2.1 benchmark doc
(`examples/ppi/benchmarks/ppi-v2.1-benchmark.md`):

- Step 4 (dimers) and everything after it read **every chain of the
  uncontained biounits**. Containment changes the dataset, so it writes a
  `work/` step folder.
- The representatives are a **selection**, not a new dataset. Step 3
  (Foldseek) may run on representatives only (through #139), and Step 8
  uses them as its negative candidate pool. So `representatives` writes a
  file under `clusters/`, not a step folder.

This follows the `work/` (changes the dataset) versus `layers/` and
`clusters/` (adds information) split.

## Goals

1. `pandora dedup --contained`: a `work/uncontained` step folder.
2. `pandora representatives`: `clusters/seq90/representatives.json`,
   one representative and its members per cluster.
3. Ranking keys and directions in policies, not hard-coded.
4. `representatives` is shaped so `--items interface` (#141) can be
   added later without changing the chain behaviour or the output format.

## Non-goals

- `--items interface` and its ranking keys (#141).
- Limiting `similarity search` to a chain list (#139).
- Making `cluster` spread a representative's scores to its cluster
  members (the Step 3 on representatives → Step 5 mapping). That's a
  later follow-up.
- Changing plain `pandora dedup` (entry-id dedup) or its output files.
- Step 2.1 ("within each biounit, keep one of each repeated chain") as
  its own command. Connected components plus representatives already
  cover it.
- Matching the old part B script. It no longer exists, so the tests below
  replace that check.

## Anti-goals

- No orchestrator or flag that runs both steps together. They stay two
  commands.
- No output-selection flags.
- Reading clusters from anything other than `pandora cluster`'s
  `clusters.json`. Only that file puts symmetry copies in their source
  chain's cluster.

## Inputs both commands share

- A step folder (`--input-dir` / `--structures-dir`) with the `entries`
  and `chains` tables that `pandora curate` writes. `chains` has
  `entry_id`, `chain_id`, `missing_n_term`, `missing_c_term`,
  `missing_middle`, `resolution`, `copy_of`; `entries` has `entry_id`,
  `resolution`.
- `--clusters`: a `clusters.json` from `pandora cluster`, a list of
  `SimilarityCluster` (`components`, `n_components`). Chain items are
  keys `"{ENTRY_ID_UPPER}_{chain_id}"` (`_chain_key` in
  `pandora/similarity/clustering.py`), symmetry copies included (for
  example `1ABC_A-2`). Map rows to cluster ids with the same key
  function. Don't parse the ids.

## 1. `pandora dedup --contained`

```sh
pandora dedup --contained --input-dir work/curation \
  --clusters clusters/seq90/clusters.json \
  [--policy policies/dedup.yaml] --output-dir work/uncontained
```

### Rule (decided)

- **It runs after `curate`** (decided "for now"). The v2.1 doc lists
  containment before the quality filters. Running it on `work/curation`
  means chains and biounits that `curate` dropped no longer count. An AA
  that lost a chain counts as A, and a biounit can't be contained in one
  that failed quality. Record this order in the benchmark doc's Step 1.
  Following the doc's order would need a second search and `cluster` on
  the canonicalised chains; revisit that later if needed.

- Turn each biounit (entry) into a **count per cluster**: how many of its
  chains, symmetry copies included, fall in each seq90 cluster.
- X is **contained in** Y (X ≠ Y) when, for every cluster, X's count is
  at most Y's count. Examples:
  - A ⊂ AA, so the monomer is dropped.
  - AA ⊄ A.
  - AB ⊂ AAB.
- This **counts chains**, which goes beyond the doc's literal wording
  ("every chain has a ≥ 90% identical chain in the other"). It was chosen
  so that a homodimer is never dropped in favour of a monomer of the same
  protein. Record it as a Pandora decision in the benchmark doc's Step 1.
- **Equal counts** (X ⊂ Y and Y ⊂ X): keep one, ranked by the dedup
  policy. The default is lowest `resolution` first, then lowest entry id
  as the final tie-break, which is always applied and isn't configurable.
  A missing resolution ranks last.
- Every entry that is strictly contained in some kept or removed entry is
  removed. Containment is transitive, so the result doesn't depend on
  order. Removing a contained X never makes a different entry
  uncontained. Strict containers are never removed for containment; only
  equal ties are.
- A chain of the step folder that isn't in `--clusters` counts as its own
  singleton cluster. Count these in the summary (`chains_unclustered`)
  and print a warning, because this usually means the clusters came from
  a different folder.

### Output: a step folder like `curate`'s

- `structures/` holds the kept entries' files, copied unchanged with no
  re-parsing (`copy_source` in `cli/_io.py`).
- `entries.<fmt>` and `chains.<fmt>` are the input rows for the kept
  entries.
- `contained.<fmt>` has one row per removed entry: `entry_id`,
  `contained_in` (one kept entry that contains it; the lowest id among
  the candidates), `equal` (bool, true when the removal was a tie).
- `failures.jsonl` follows the convention (here it's normally empty).
- `dedup.json` is the step summary from `build_step_summary("dedup",
  ...)`, with counts `entries_in`, `entries_kept`, `entries_contained`,
  `chains_kept`, `chains_unclustered`, and `previous_summary` = the
  input's `curation.json`.
- It honours `--table-format`, `--overwrite`, and `prepare_step_dir`.

Plain `dedup` (without `--contained`) is unchanged. `--clusters` is
required with `--contained` and an error without it.

### Library

- In `pandora/datasets/curation.py` (or a new `datasets/containment.py`):
  `find_contained(cluster_counts: dict[str, Counter[str]], rank_values:
  dict[str, tuple], ...) -> list[ContainmentRecord]` or similar. It's a
  plain function on plain data, with no `Structure` objects. Exported in
  `pandora.datasets.__all__`. The signature is the builder's call.
- Expect up to around 10⁵ entries. Compare only against entries that
  share a cluster, not all pairs: index entries by cluster, and check
  only candidates that contain X's rarest cluster.

### Schemas (`pandora/schemas/dataset.py`)

- `DeduplicationRules` gets `contained_rank: list[RankKey] =
  [RankKey(key="resolution", order="min")]`. Keys: any numeric column of
  the `entries` table; an unknown key is an error naming it. The existing
  `enabled` field is unchanged and **ignored by `--contained`** (plain
  `dedup` doesn't read it either; it uses `--disable`). `--contained`
  always removes; `--contained --disable` writes the `contained` rows but
  keeps every entry, as plain `--disable` does.
- `ContainmentRecord(entry_id, contained_in, equal)`.
- `--policy` is an optional `DeduplicationRules` YAML; without it the
  defaults apply.

## 2. `pandora representatives`

```sh
pandora representatives --clusters clusters/seq90/clusters.json \
  --structures-dir work/uncontained \
  [--similarity layers/seq --hit-filter policies/seq90.yaml] \
  [--policy policies/representatives.yaml] [--items chain] \
  --output clusters/seq90/representatives.json
```

### Behaviour (decided)

- `--items` accepts only `chain` for now (the default). Any other value
  is an argparse error, until #141 adds `interface`.
- **Limit clusters to the step folder:** drop members whose chain isn't
  in `--structures-dir`'s `chains` table, for example because
  `dedup --contained` removed its biounit. Drop clusters left empty.
- **Who can be picked:** non-copy chains only (`copy_of` is null). Copies
  stay in `members`. If a cluster has only copies (it shouldn't), fall
  back to every member.
- **Ranking:** the policy's ordered keys, then the lowest item id as an
  always-applied, non-configurable final tie-break. Keys available for
  chains:
  - `missing_residues`: `missing_n_term + missing_c_term +
    missing_middle`, an **absolute count** (decided; the full length is in
    `seqres_length`, so a fraction can be derived later). The counts come
    from `curate`, before any rebuilding of missing residues, as the doc
    requires.
  - `n_hits`: the number of distinct other chains in the **same cluster,
    after limiting it to the step folder**, linked to this chain by a hit
    that passes `--hit-filter`. Chains whose biounit was removed don't
    count. This needs
    `--similarity` and `--hit-filter` (or `--min-score`/`--min-coverage`,
    as in `cluster`); it's an error if the policy uses `n_hits` without
    them.
  - Any numeric column of the `chains` table (for example `resolution`,
    `residue_count`). An unknown key is an error that names the key.
  - A null value ranks last for both `min` and `max`.
- Default policy, which is the v2.1 rule: `[missing_residues: min,
  n_hits: max]`. Because the default uses `n_hits`, the v2.1 call passes
  `--similarity`.

### Output

`representatives.json` is a JSON list, one object per cluster, sorted by
representative id:

```json
[{"representative": "1A3N_B", "members": ["1A3N_B", "1A3N_D", "1ABW_B", "1ABW_C"], "n_members": 4}]
```

Each `RepresentativeCluster` holds the members after limiting to the step
folder, copies included. A provenance JSON next to it
(`representatives_provenance.json`, following the
`partition_provenance.json` pattern) records the policy, the inputs, the
`items` kind and the counts.

The output must not change when #141 adds interfaces: the same model
with interface ids.

### Library and schemas

- `pandora/similarity/representatives.py`: `pick_representatives(clusters:
  list[SimilarityCluster], values: dict[str, dict[str, float | None]],
  policy: RepresentativePolicy, eligible: set[str] | None = None) ->
  list[RepresentativeCluster]`. It works on plain values, so #141 only has
  to supply interface values. A helper computes `n_hits` from edges and
  clusters. Both are exported. Exact signatures are the builder's call.
- `schemas/similarity.py`: `RepresentativePolicy(items: SimilarityItems
  = "chain", rank: list[RankKey])` and
  `RepresentativeCluster(representative, members, n_members)`.
- `RankKey(key: str, order: Literal["min", "max"])` lives in
  `schemas/common.py` and is shared with `DeduplicationRules`.

## Docs, policies, example

- `docs/policies/`: add a YAML listing field types for
  `RepresentativePolicy`, and extend the dedup rules. They must match the
  models.
- `docs/reference/policies.md`: describe `contained_rank`, every
  `RepresentativePolicy` field, the rank keys, and the null and final
  tie-break rules.
- Usage docs, next to the `cluster` docs in `docs/usage/similarity.md`
  and wherever `dedup` is documented: full code inline, library and CLI
  tabs, with the outputs shown and tested.
- `pandora/cli/README.md`: add `representatives` to the command list,
  plus one paragraph on both commands.
- `examples/ppi/v2.1/src/02_sequence_clusters.sh`: replace the "NOT YET"
  block with the two commands, and add
  `dataset/policies/representatives.yaml`. Steps 03 to 05 move from
  `work/curation` to `work/uncontained` (Step 3 still uses every chain
  until #139).
- Benchmark doc, Step 1: record the counting decision. Step 2: record
  absolute counts and the lowest-id tie-break.
- Regenerate the ERD
  (`uv run --extra docs python docs/scripts/generate_erd.py`), and build
  the docs site.

## Tests (local fixtures only, plain pytest)

Use hand-built `clusters.json` and tables, or a small step folder built
from `datasets/dev/mmcif/` with `curate`.

`tests/test_containment.py` (or extend `test_curation.py`):

- A and AA: A removed, AA kept.
- AB and AAB: AB removed.
- Equal sets: the better resolution kept; at equal resolution, the lower
  id kept; a missing resolution loses.
- An unclustered chain is a singleton, and is counted.
- CLI: the step folder layout, `contained` rows, `dedup.json` counts, and
  `--clusters` required with `--contained`.

`tests/test_representatives.py`:

- Fewest missing residues wins.
- A tie broken by `n_hits`, then by lowest id.
- A copy is never picked.
- Members outside the step folder are dropped; a cluster emptied that way
  disappears.
- A `max` order works; an unknown key is an error.
- `n_hits` without `--similarity` is an error.
- `--items interface` is rejected for now.
- Inputs are never mutated.

## Left to the builder

- Exact function names and signatures, and which module containment goes
  in.
- Whether `contained` and representatives' provenance use existing
  provenance models or new small ones (in `schemas/`).
- The data structure for the containment candidate search, within the
  "don't compare all pairs" constraint.

## Calls made without explicit sign-off

- The `RankKey` model is shared between the two policies, and the lowest
  id is always the final, non-configurable tie-break.
- `dedup --contained` writes a full step folder (tables plus
  `dedup.json`), while plain `dedup` keeps its older output files.
- `n_hits` counts only partners inside the same cluster, after the
  cluster is limited to the step folder.
- `--contained` ignores `DeduplicationRules.enabled`; `--disable` means
  "record but keep".
- `contained_rank` keys are numeric `entries` columns only (no per-entry
  aggregation of chain values).
- Copies are never picked as representatives.
- A null rank value always ranks last.
- An unclustered chain counts as a singleton, with a warning.
