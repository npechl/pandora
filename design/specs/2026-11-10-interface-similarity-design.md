# Interface-vs-interface similarity (per-side Jaccard) — design

Issue #134. Priority 2 of 7 for building PPI v2.1 with Pandora
(`examples/ppi/benchmarks/ppi-v2.1-benchmark.md`, Step 4).

## Intent

PPI v2.1 Step 4 removes redundant dimers. Two dimers are redundant when
they are the same two proteins (each chain ≥ 90% sequence identity to
its counterpart) **and** they touch at the same spot. Pandora can find
interfaces today (a contact layer: every entry's chain pair with at
least one contact row) but cannot compare one interface with another.

This adds one thing: a command that scores every pair of interfaces
that are the same two proteins by how much their interface residues
overlap, and writes the scores as an interface-level `similarity`
layer. The existing `pandora cluster` then groups them.

The user does not yet know where the cut-off should sit (0.1 is
suggested, "pending a look at the distribution"). The first job of
this layer is to *show* the distribution, so it holds every compared
pair, including pairs scoring 0.

## Goals

- `pandora similarity interfaces` writes a `similarity` layer with
  `items: interface` from a `residue_contacts` layer and chain
  sequence clusters.
- Handles homomers: tries both chain orders and records which was used.
- Matches residues across non-identical (≥ 90%) chains through a
  pairwise sequence alignment.
- The output works unchanged with
  `pandora cluster --similarity ... --items interface --contacts ...`.

## Non-goals

- The representative set (one dimer per protein pair), the link from
  representatives to complete-set dimers, and choosing which dimer
  represents a group. Representative set: #135 (chain layer rolled up
  to interfaces). The rest is later work.
- Choosing the cut-off. The user picks it from the distribution and
  passes it to `cluster --min-score`.
- Other interface metrics (iDist, Foldseek interface mode). Those come
  in later through `similarity import`. No `--metric` flag until there
  is a second metric (the issue's sketch had one; dropped as YAGNI).
- Parallelisation, sharding, restart.

## Anti-goals

- A pair that is "the same spot" on only one side scoring as a match.
  One protein reusing one patch to bind its partner in two different
  places is two different interfaces.
- Residues silently mapped to the wrong position. A residue the
  alignment cannot place counts as not shared, never guessed.
- Writing only pairs above a threshold (it would hide the distribution
  the cut-off is chosen from).

## Decisions (made with the user)

1. **Input interfaces** come from a `residue_contacts` layer, normally
   the filtered `contacts_real` from `annotate filter-contacts`. An
   interface is an entry's unordered chain pair, as in
   `label_interfaces` (`pandora/annotations/counts.py`). Item id
   `1ABC_A_B` via `item_id("interface", ...)`
   (`pandora/similarity/layer.py`).
2. **Criterion 1 (same proteins)** comes from the user's Step 2 chain
   clusters at 90%: a `clusters.json` written by `pandora cluster`
   (list of `SimilarityCluster`, components are chain item ids
   `1ABC_A`). Clusters only decide which interfaces are *looked at*.
   The user regenerates Step 2 clusters with `pandora cluster` (identity
   and both-side coverage ≥ 0.9 via a `--hit-filter` YAML) instead of
   the example script's `chain_clusters.tsv`. All chains in the file
   take part, including those of biounits the example marked
   `entry_kept = no` (Step 4 uses all chains from Step 1).
3. **Pairwise identity gate (criterion 1 proper):** a chain order
   counts only if *each* chain is ≥ `--min-identity` (default 0.9)
   identical to its counterpart on the pairwise alignment. Single-
   linkage clusters can join chains below 90%; this check stops them.
4. **Residue matching:** identical sequences match by position
   (`seq_id`). Otherwise a pairwise alignment of the two chain
   sequences with gemmi (already a dependency;
   `gemmi.align_string_sequences`) maps `seq_id` positions.
5. **Score per side:** `jaccard_1` and `jaccard_2`, one per side.
6. **`score` = min(jaccard_1, jaccard_2).** Both sides must overlap.
   Both per-side columns are written so the mean or anything else can
   be looked at later.
7. **Homomers / chain order:** when both orders satisfy criterion 1,
   score both and keep the order with the higher `score`. Record it.
8. **Scope** stops at the score layer.

## Behaviour

### Interface residue sets

The interfaces are exactly those of `interface_items`
(`pandora/similarity/layer.py`): rows where both `seq_id`s are set and
the two chains differ. *(Builder-check call: `cluster --items interface
--contacts` builds its items this way and drops rows naming anything
else, so following `label_interfaces` instead, which also keeps `A_A`
and None-`seq_id` pairs, would break "works unchanged with
`cluster`".)* Other rows are ignored and counted.

For each interface `(entry, a, b)`: side `a` = distinct `seq_id` of
chain `a` residues in contact with chain `b`; side `b` likewise.

### Which pairs are compared

Let `cl(x)` be the cluster of chain `x`. Interfaces `I = (e, a, b)` and
`J = (f, c, d)`, `I ≠ J`, each unordered pair once:

- **direct** order fits when `cl(a) = cl(c)` and `cl(b) = cl(d)`;
- **swapped** order fits when `cl(a) = cl(d)` and `cl(b) = cl(c)`.

Score an order only if both of its chain pairings pass the identity
gate: identity = identical aligned columns / aligned (match or
mismatch) columns, from the same gemmi alignment used for mapping
(identical sequences: 1.0). Compare when at least one order fits and
passes; a pair where orders fit on clusters but none passes the gate
writes no row and is counted (`pairs_below_identity`). Pairs within the same entry are
compared too (two copies of a dimer in one assembly). Group interfaces
by their unordered cluster pair `{cl(a), cl(b)}` and compare all-vs-all
within each group. That is the only pairing loop.

A chain in no cluster has no cluster, so its interfaces are never
compared; they are counted in the summary.

### Mapping and Jaccard for one side

Comparing side `x` of `I` (residue set `Sx`, sequence `qx`) with side
`y` of `J` (`Sy`, `qy`):

- if `qx == qy`: `M(Sx) = Sx`;
- else align `qx` to `qy`; `M(Sx)` = the `qy` positions aligned to a
  residue of `Sx` (match or mismatch columns; gap columns drop out).
- `jaccard = |M(Sx) ∩ Sy| / (|Sy| + |Sx| − |M(Sx) ∩ Sy|)`. The
  denominator uses `|Sx|`, not `|M(Sx)|`, so residues the alignment
  could not place count against the overlap (anti-goal: never
  guessed). Empty union → 0.

Cache alignments per unordered `(sequence, sequence)` pair; the same
two entity sequences recur across many interfaces.

### Output rows

One row per compared pair, columns fitting the `similarity` template
(`SimilarityRow`, extra columns allowed):

| column | meaning |
|---|---|
| `item_id_1`, `item_id_2` | interface ids, `item_id_1 < item_id_2` |
| `score` | `min(jaccard_1, jaccard_2)` for the kept order |
| `jaccard_1` | first chain of `item_id_1` vs its matched chain in `item_id_2` |
| `jaccard_2` | second chain of `item_id_1` vs its matched chain |
| `chain_order` | `"direct"` or `"swapped"` |
| `identity_1`, `identity_2` | pairwise identity of each matched chain pair (kept order) |

Chain order within an id follows `label_interfaces`: smaller chain id
first. "direct" means item 1's first chain ↔ item 2's first chain.
Sorted by `(item_id_1, item_id_2)`.

### Sequences

Chain sequences come from the `chains` table of the structures folder
the contact layer was computed from, found the same way as
`annotate label-interfaces` (`_layer_origin` + `find_table(origin,
"chains")`, `DatasetChainRecord.sequence`). Position `i` (1-based) of
the canonical one-letter sequence is `seq_id = i`. An interface whose
chain has no sequence is skipped and counted.

*My call (builder check):* `seq_id = i` only holds when
canonicalisation kept `label_seq_id` (`residue_numbering: preserve`, as
the PPI v2.1 policy does; the repo default `renumber` breaks it), and
when `sequence` is the `_can` one-letter code (the non-`_can` fallback
has `(MSE)` groups). Guard instead of trusting it: every contact row
carries `comp_id_s`; if the one-letter code of `comp_id` (gemmi) does
not equal `sequence[seq_id - 1]`, or `seq_id` is out of range, the
interface is skipped and counted (`interfaces_sequence_mismatch`). A
renumbered dataset then fails loudly with large counts, never silently.
The usage docs state that the structures must keep `label_seq_id`.

gemmi's `align_string_sequences` takes lists of residue names: pass
`list(seq)` with default scoring (plain match/mismatch, fine at
≥ 90% identity).

## Library API

New module `pandora/similarity/interfaces.py`, exported from
`pandora/similarity/__init__.py`:

```python
def interface_similarity(
    contacts: Iterable[ResidueContactRow],
    chain_clusters: Iterable[SimilarityCluster],
    sequences: Mapping[tuple[str, str], str],
    min_identity: float = 0.9,
) -> list[LayerRow]: ...
```

`sequences` is keyed by `(entry_id.upper(), chain_id)`, like
`label_interfaces`' `lengths`. Returns rows as above. Pure function, no
I/O. No new schema model is needed (`SimilarityRow` with extra columns
covers the rows). Full Google docstring.

Whether the function also returns skip counts (or the CLI derives them)
is the builder's call; keep it simple.

## CLI

```sh
pandora similarity interfaces \
  --contacts layers/contacts_real \
  --chain-clusters clusters/seq90.json \
  --output-dir layers/iface_jaccard \
  [--min-identity 0.9] [--overwrite] [--table-format parquet|...]
```

Then, unchanged (note `--min-score` is inclusive, `>=`; the benchmark
says "> 0.1"):

```sh
pandora cluster --similarity layers/iface_jaccard --items interface \
  --contacts layers/contacts_real --min-score 0.1 \
  --output clusters/dimers_complete.json
```

Follows `annotate label-interfaces` for structure: `check_contact_layer`,
the origin lookup and its error messages, `prepare_output_dir`,
`write_table(..., "layer", ...)`, `layer.json` built with
`build_step_summary` from a `LayerDescription(name="iface_jaccard",
template="similarity", items="interface", tool="pandora", ...)` with
`parameters` recording the contact layer name, structures dir and the
chain clusters file path. Output must read back with the existing
similarity layer readers (`read_similarity_table`), so `cluster` takes
it as-is.

Summary counts in `layer.json` and one printed line: `interfaces`,
`interfaces_without_cluster`, `interfaces_without_sequence`,
`non_polymer_rows`, `interfaces_sequence_mismatch`, `pairs_compared`,
`pairs_below_identity`, `pairs_swapped`,
`pairs_zero_score`, `alignments`.

## Constraints

- No new dependency (gemmi alignment, stdlib otherwise).
- No new policy file: the only rule (`min` of sides) is fixed and both
  inputs to it are written out. The identity gate is a CLI flag with a
  default, recorded in `layer.json` `parameters`. Coverage is not
  re-checked (the Step 2 clusters already required it).
- Never mutate inputs.

## Known ceiling

All-vs-all inside each cluster-pair group is quadratic. A popular
protein pair (thousands of interfaces of the same complex) gives
millions of rows. Accepted for now; mark with a `ponytail:` comment.
Upgrade path: stream rows to the table in batches, or compare only
across entries' representatives.

## Tests (`tests/test_interface_similarity.py`, local only)

Small hand-built `ResidueContactRow` lists and sequences, plus one CLI
run on a fixture-derived contact layer:

- identical heterodimer interfaces → score 1.0, `direct`;
- disjoint contacts on one side only → score 0 (min, not mean);
- homodimer where only the swapped order lines up → `swapped`, score
  from the swapped order;
- non-identical sequences (an insertion) → residues mapped through the
  alignment;
- same cluster but < 90% pairwise identity → no row, counted;
- interfaces in different clusters → not compared; chain in no cluster
  → counted, not compared;
- CLI output clusters with `pandora cluster --items interface
  --min-score`.

## Docs

A usage page section with inline library + CLI tabs and tested output
(project convention), next to the other similarity commands; CLI
reference entry. Build the docs site.

## Left to the builder

- Alignment scoring parameters for gemmi (defaults are fine).
- How skip counts travel from the function to the CLI.
- Table writing details beyond "reads back with `read_similarity_table`".
