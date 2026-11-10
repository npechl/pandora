# Positives, negatives and chain-limited search (#136, #137, #139)

Date: 2026-11-10. Issues: #136 (Step 7 rebalancing), #137 (Step 8
negative sampling), #139 (`similarity search --chains`). Related: #135
(dimer clusters), #138 (chain representatives), #141 (interface
representatives), #83 (no negative-sampling primitive).

## Why

These are the last three PPI v2.1 steps without a `pandora` command
(`examples/ppi/benchmarks/ppi-v2.1-benchmark.md`, Steps 3, 7 and 8).
`examples/ppi/v2.1/src/03_structure_similarity.sh` is a placeholder
waiting for #139, and there is no script for Steps 7 and 8. Once these
land, the v2.1 build runs end to end with `pandora` commands only.

Decisions marked **(user)** were made by the user in the brainstorm.
Decisions marked **(my call)** were proposed by Claude and accepted
without objection, or left to Claude.

## Goals

1. `pandora similarity search --chains <file>` limits the search to a
   list of chains; Foldseek no longer needs `--input-dir`.
2. `pandora rebalance` writes the positive set (Step 7) as a step folder.
3. `pandora negatives` writes the negative pairs (Step 8) as a step
   folder, reproducible from a seed.
4. `03_structure_similarity.sh` runs for real; new scripts
   `06_positives.sh` and `07_negatives.sh` run Steps 7 and 8.

## Non-goals

- A homomer/heteromer label layer for `partition --labels`. The
  positives table carries a `kind` column; a layer can come when a split
  needs it.
- More than one negative per positive, or both chains of a positive
  getting a negative.
- Step 9 (final assembly / export) and splits.
- Any change to `representatives` or to the dimer groups' contents.

## Anti-goals

- No orchestrator: three independent commands, each reading files the
  earlier steps wrote.
- No hidden randomness: the only random choices are in `negatives`, and
  they all come from one seeded generator recorded in the summary.
- No silent fallbacks: an id that can't be mapped (a chain missing from
  the clusters, a `--chains` id missing from the step folder) is an
  error naming it.

## Shared conventions

- Item ids: chain `E_C`, interface `E_C1_C2` (chains sorted). Build and
  parse them with `item_id` / `parse_item_id`, never by splitting
  strings. Symmetry copies look like `1AYI_A-2`.
- Step folders follow `pandora/cli/README.md`: tables in `--table-format`
  (default Parquet, `jsonl` allowed), refuse a non-empty output folder
  unless `--overwrite`, a step summary JSON built with
  `build_step_summary(...)` whose `previous_summary` points at the input
  step/layer summary.
- Logic in `pandora/similarity/` (plain functions, no I/O), models in
  `pandora/schemas/similarity.py`, CLI wiring in
  `pandora/cli/similarity.py`. The step-folder writing pattern
  (`build_step_summary`, `--table-format`) is in `pandora/cli/annotate.py`;
  `representatives` writes a single JSON instead. `_interface_values`
  hard-codes `"representatives"` in its messages; pass the command name
  in when reusing it.
- Tests in `tests/` with local fixtures. Fake Foldseek:
  `FAKE_FOLDSEEK` in `tests/test_structure_similarity.py` and
  `FAKE_FOLDSEEK_CLI` in `tests/test_cli.py`.

## 1. `similarity search --chains` (#139)

```sh
pandora similarity search --engine foldseek \
  --structures-dir work/uncontained \
  --chains clusters/seq90/representatives.txt \
  --exhaustive-search --output-dir layers/tm
```

- `--chains <file>`: chain item ids, one per line; blank lines and
  `#` comments skipped. Works for both engines. **(user)**
- Every id must be a chain of `--structures-dir`'s `chains` table;
  otherwise exit naming the first missing id. **(my call)**
- MMseqs2: filter `_chain_sequences(...)` to the listed ids. Copies are
  still left out (a copy in the list is ignored, not an error).
- Foldseek **without `--input-dir`**: Pandora writes one single-chain
  mmCIF per chain (all chains of the step folder, or the `--chains`
  list) into a temporary folder, and passes
  `compute_structure_similarity` an `{item_id: path}` mapping so hits
  come back named by item id. The folder is deleted afterwards.
  **(user)** Builder's choice how to cut the chain (e.g. read the step
  folder structure, keep one chain, `structure_to_mmcif`); only
  coordinates matter to Foldseek. Symmetry copies are left out, as for
  MMseqs2.
- Foldseek **with `--input-dir`** keeps today's behaviour. `--chains`
  together with `--input-dir` is an error (we can't map file names to
  chain ids reliably). **(my call)**
- The layer's items stay the step folder's chains, so `unmatched` and
  counts behave as today. `layer.json` records the `--chains` path and
  the number of chains searched.

Build script: `03_structure_similarity.sh` writes
`clusters/seq90/representatives.txt` with
`jq -r '.[].representative'` and runs the command above.

## 2. `pandora rebalance` (#136)

```sh
pandora rebalance --clusters clusters/dimer_groups/clusters.json \
  --contacts layers/contacts \
  --chain-clusters clusters/seq90/clusters.json \
  --chain-representatives clusters/seq90/representatives.json \
  --similarity layers/iface_sim --hit-filter policies/same_proteins.yaml \
  --similarity layers/tm --hit-filter policies/tm.yaml \
  --policy policies/dimer_representatives.yaml \
  --output-dir work/positives
```

### Inputs

- `--clusters`: the Step 6 dimer groups (`SimilarityCluster` list, interface
  items). Every dimer in a group is a candidate. **(user)**
- `--contacts`: the contact layer the groups were built from; gives the
  interface items and, through the structures folder it traces to, the
  rank values (reuse `_interface_values`).
- `--chain-clusters`: the 90% chain clusters (`pandora cluster` output).
- `--similarity` / `--hit-filter` / `--min-score` / `--min-coverage`:
  the same layers and filters `cluster` took in Step 5, with the same
  per-layer rules (`_hit_filters`).
- `--policy`: a `RepresentativePolicy` YAML with `items: interface`, used
  only to break ties. Default: the interface default (most interface
  residues, best resolution, fewest missing residues, lowest id).
  **(my call, user did not object)**

### Homomer / heteromer

A dimer is a **homomer** when both chains are in the same 90% chain
cluster, otherwise a **heteromer**. **(user)** A dimer chain missing
from `--chain-clusters` is an error naming it. Symmetry copies are in
their source chain's cluster, so `1AYI_A_A-2` is a homomer. Symmetry
copy dimers are normal candidates here, unlike in `representatives`.
**(my call)**

### The dimer graph

Degree needs the real dimer-to-dimer graph. **Don't reuse
`interface_edges`**: it is a connectivity-preserving shortcut (each
chain links only the first interface it appears in), so its degrees are
wrong. Build adjacency directly, within each group only:

- interface-level layer edge (e.g. `iface_sim`): the two dimers are
  adjacent;
- chain-level layer edge (c1, c2): every dimer containing c1 is
  adjacent to every dimer containing c2;
- two dimers sharing a chain are adjacent (a chain matches itself).

With `--chain-representatives` (section 4), every chain is first mapped
to its representative, **before all three rules**. So two dimers whose
chains share a 90% representative are adjacent, as TM ≈ 1 would make
them in v2.1. Symmetry copies are members of their source chain's
representative, so they need no separate copy handling. **(my call)**

Run with the same layers, filters and `--chain-representatives` as
Step 5, the groups are the connected components of this graph, so an
edge between two groups means the inputs don't match: exit with an
error naming both dimers. **(my call)** A dimer in `--clusters` missing
from `--contacts` is also an error. A pure function in
`pandora/similarity/`, e.g.
`dimer_adjacency(interface_ids, chain_edges, interface_edges)`.

### Selection, per group

- **Only homomers:** keep the homomer with the highest degree. **(v2.1)**
- **Any heteromer:** drop every homomer (reason `homomer_in_mixed_group`).
  Then, until no heteromer is left: compute each remaining heteromer's
  degree counting only remaining heteromers, keep the one with the
  lowest degree, remove it and its neighbours (reason
  `similar_to_kept`, with the kept dimer's id). Degrees are recounted
  every round. **(my call)**
- **Ties** (same degree): the `--policy` ranking; its last key is the
  lowest id, so the result is deterministic. **(my call)**
- Homomer-only groups drop the rest with reason `not_most_connected`.

### Outputs (`work/positives/`)

- `positives.<fmt>`: `dimer` (interface id), `entry_id`, `chain_id_1`,
  `chain_id_2`, `group` (the group's first component id), `kind`
  (`homomer` / `heteromer`), `degree` (at the time it was kept).
- `dropped.<fmt>`: `dimer`, `group`, `reason`, `kept` (the dimer that
  caused the drop, if any).
- `rebalance.json`: step summary with the policy, the layers and their
  hit filters, input paths, and counts (groups, dimers in, positives by
  kind, dropped by reason). `previous_summary` = the contact layer's
  `layer.json`.

## 3. `pandora negatives` (#137)

```sh
pandora negatives --positives work/positives \
  --representatives clusters/seq90/representatives.json \
  --similarity layers/tm --hit-filter policies/tm.yaml \
  --contacts layers/contacts_any --seed 0 --output-dir work/negatives
```

### Inputs

- `--positives`: the `rebalance` step folder.
- `--representatives`: the Step 2 representatives JSON. Its
  `representative`s are the pool; its `members` map every chain to its
  representative. **(user)** A chain that must be mapped (a positive's
  chain, a contact-layer chain) but appears in no `members` is an
  error naming it.
- `--similarity` / `--hit-filter` (and `--min-score` /
  `--min-coverage`): one **chain** layer of structure similarity; its
  edges after the filter define "TM ≥ 0.5". Edges are mapped to
  representatives; an edge naming a chain outside the representatives
  is mapped through `members` too.
- `--contacts`: the "any contact" `residue_contacts` layer. Two chains
  "bind" when they have any contact row; pairs are mapped to
  representatives (a homomer gives a self-binding representative).
- `--seed` (int, default 0).

### Algorithm

All lookups happen at representative level. **(user)** Let `rep(x)`
map a chain to its representative, `tm(S)` be the representatives
TM-linked to any of `S`, and `binders(S)` those binding any of `S`.

1. Pool = all representatives minus `rep(c)` for every chain `c` of
   every positive.
2. Process positives in sorted `dimer` order. **(my call)** For each:
   1. A = one of the two chains, chosen by the seeded generator.
      **(user, recommended by Claude)** `a = rep(A)`.
   2. set 1 = `{a} ∪ tm({a})`. `a` is always in set 1, even without a
      self-hit, so A's own partners (B among them) go through set 2.
   3. set 2 = `binders(set 1)`; set 3 = `tm(set 2)`.
   4. Candidates = pool − (set 1 ∪ set 2 ∪ set 3), sorted. The sets
      are excluded for this pick only, not removed from the pool.
      **(my call: v2.1 says only C is removed "for later negatives")**
   5. No candidates: record the positive as unpaired and go on.
      Otherwise C = a seeded choice from the sorted candidates; emit
      (A, C) and remove C from the pool.
3. One generator, `random.Random(seed)`, used in this order: A's choice,
   then C's, per positive.

Homomer positives (A = B) need nothing special. **(my call)**

### Outputs (`work/negatives/`)

- `negatives.<fmt>`: `chain_1` (A, the positive's actual chain),
  `chain_2` (C, a representative chain), `positive` (the dimer id).
- `negatives.json`: step summary with seed, hit filter, layer and input
  paths, counts (pool size at start, positives, negatives, unpaired) and
  the unpaired positives' ids. `previous_summary` = `rebalance.json`.

## 4. Step 5 when the TM layer covers representatives only

Decided: option A. **(user)**

After #139, `layers/tm` holds hits between 90% representatives only.
`05_dimer_clusters.sh` (Step 5) and `rebalance` read it as a chain layer
over **every** dimer, so a dimer whose chains aren't representatives
(for example `1A3N_C`, a member of `1A3N_A`'s cluster) gets no TM edges.
Its Step 6 group and its degree would be wrong. `negatives` is not
affected (it maps everything to representatives).

**Fix:** `cluster` and `rebalance` take `--chain-representatives
<representatives.json>` and map each chain to its representative
(through `members`) before using a chain layer. Chains with the same
representative then also count as similar, which is what TM ≥ 0.5 would
say anyway. One shared helper. A chain missing from every `members` is
an error naming it. Without the flag, behaviour is unchanged.

For `cluster`, `interface_edges` takes an optional `chain_map` and keys
every chain (of the interfaces and of the edges) by its representative.
Adding (member, representative) edges instead is not enough: a
representative that is in no dimer would be dropped by the rollup, and
its members would never link. `rebalance` uses the same map, so both
build the same graph. `--chain-representatives` needs `--items
interface`. `05_dimer_clusters.sh` passes
`--chain-representatives clusters/seq90/representatives.json`, and
`06_positives.sh` passes the same file.

The rejected option was keeping Step 3 on all chains in the build.

## Build scripts and docs

- `03_structure_similarity.sh`: real command (section 1).
- `06_positives.sh`, `07_negatives.sh`: one command each.
- `examples/ppi/v2.1/README.md` script table, and the benchmark doc's
  "With Pandora" notes for Steps 3, 7 and 8, recording the **(my call)**
  choices above.
- `pandora/cli/README.md` and `docs/usage/similarity.md`: the two new
  commands and `--chains`. Build the docs site.

## Tests

- `--chains`: MMseqs2 and fake-Foldseek runs search only the listed
  chains; an unknown id is an error; `--chains` + `--input-dir` is an
  error; Foldseek without `--input-dir` receives one file per chain.
- `dimer_adjacency`: shared chain, chain edge expansion, interface
  edge, mapping through representatives (copies included); a
  cross-group edge errors.
- `cluster --chain-representatives`: a TM edge between two
  representatives joins dimers built from their non-representative
  members.
- Rebalance: homomer-only group keeps the most connected; mixed group
  drops homomers; a path of heteromers A–B–C keeps A and C (recounting);
  a tie falls to the ranking; a chain missing from the clusters errors.
- Negatives: set 1/2/3 exclusion on a small hand-made graph; B is never
  A's negative; C is used once; pool exhaustion records unpaired; same
  seed gives the same output, a different seed may not.
