# Positives, Negatives and Chain-Limited Search Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `similarity search --chains`, `cluster --chain-representatives`, `pandora rebalance` and `pandora negatives`, so the PPI v2.1 build runs Steps 3–8 with `pandora` commands only (#139, #136, #137).

**Architecture:** Pure logic in `pandora/similarity/` (new `rebalance.py`, `negatives.py`; small additions to `representatives.py` and `clustering.py`), models in `pandora/schemas/similarity.py`, CLI wiring for the two new commands in a new `pandora/cli/sampling.py` that reuses the private helpers of `pandora/cli/similarity.py`. Step folders follow the `pandora/cli/annotate.py` pattern (`check_table_format`, `prepare_output_dir`, `write_table`, `build_step_summary`).

**Tech Stack:** Python 3, pydantic v2, argparse, pytest, stdlib `random`.

**Spec:** `docs/superpowers/specs/2026-11-10-positives-negatives-chain-search-design.md`

## Global Constraints

- `from __future__ import annotations`; type hints; line length 80; `ruff format`.
- Public functions: Google-style docstring (`Args:`, `Returns:`, `Raises:`) and an `__all__` entry in `pandora/similarity/__init__.py`.
- Never mutate inputs; no new dependencies; catch named exceptions only.
- Item ids built with `item_id(...)` / parsed with `parse_item_id(...)`, never by string splitting.
- Every unmappable id is an error naming it (no silent fallback).
- Tests use local fixtures only; Foldseek is faked (`FAKE_FOLDSEEK_CLI` in `tests/test_cli.py`).
- `examples/` is untracked: edit the scripts, don't commit them, add no pytest for them.

## Review Focus

- Empty inputs: a group with one dimer, a contacts layer with no rows, a positives table with no rows. Expect: one positive, empty outputs written, no crash.
- Duplicate or reversed contact rows (B–A and A–B) must give one binding pair, not two.
- `--chains` file with trailing whitespace, blank lines and `#` comments.
- A TM layer that names chains outside `--representatives` members (Step 3 run on all chains): must map through `members`, not error.
- Re-running with `--overwrite` replaces the step folder; without it, refuses.

Each line has its test in the owning task (Tasks 1, 4, 5, 6).

---

### Task 1: `similarity search --chains` and per-chain Foldseek inputs (#139)

**Files:**
- Modify: `pandora/cli/similarity.py` (`_cmd_search`, `_chain_sequences`, `register`)
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: `_read_chain_list(path: Path) -> list[str]` and `_chain_files(structures_dir: Path, chains: set[str] | None, out: Path) -> dict[str, Path]` (private, in `pandora/cli/similarity.py`).

- [ ] **Step 1: Write the failing tests** in `tests/test_cli.py`, next to `test_similarity_search_foldseek_matches_file_stems`:
  - `test_similarity_search_mmseqs_limits_to_chains`: `_similarity_step` folder, `--chains` file `"# reps\n104M_A \n\n"`; a fake `mmseqs` that records its input FASTA ids to a file; assert only `104M_A` was searched and `layer.json`'s counts include `chains_searched == 1`.
  - `test_similarity_search_chains_unknown_id_exits`: `--chains` with `9XXX_A` → `SystemExit` matching `"9XXX_A"`.
  - `test_similarity_search_chains_rejects_input_dir`: foldseek with `--chains` and `--input-dir` → `SystemExit` matching `"--chains"`.
  - `test_similarity_search_foldseek_writes_chain_files`: a real step folder (canonicalise + curate one fixture from `datasets/dev/mmcif/`, as other CLI tests do), foldseek without `--input-dir`; the fake binary lists `"$2"` into a file; assert one file per non-copy chain, named `<item_id>.cif`, each with atoms of one `label_asym_id` only.
- [ ] **Step 2: Run** `uv run pytest tests/test_cli.py -k "similarity_search" -v`. Expected: the four new tests FAIL.
- [ ] **Step 3: Implement.**
  - Flag: `--chains`, help `"file of chain item ids (one per line, # comments) to limit the search to"`.
  - `_read_chain_list`: strip lines, skip blank and `#`.
  - Check every id is in the step folder's chain items (`_item_set` result); exit `similarity search: --chains names <id>, not a chain of <dir>`.
  - `--chains` with `--input-dir` → exit. Remove the "foldseek needs --input-dir" exit.
  - MMseqs2: filter `_chain_sequences` to the list.
  - Foldseek without `--input-dir`: in a `tempfile.TemporaryDirectory()`, `_chain_files` reads `read_structures_dir(structures_dir)` (`# ponytail: parses every entry; limit to the listed entries if it gets slow`), and for each non-copy chain (`copy_of` empty in the chains table) in the list writes `structure.model_copy(update={"atoms": <atoms with label_asym_id == chain>, "asym_units": <that asym>, "connections": [], "raw": {}})` with `structure_to_mmcif` to `<out>/<item_id>.cif`; pass the dict to `compute_structure_similarity`.
  - Add `chains_searched` to the counts; record the `--chains` path as the summary `input` when given.
- [ ] **Step 4: Run** the same command. Expected: PASS, and `uv run pytest tests/test_cli.py -v` all pass.
- [ ] **Step 5: Update** `examples/ppi/v2.1/src/03_structure_similarity.sh`: `jq -r '.[].representative' $D/clusters/seq90/representatives.json > $D/clusters/seq90/representatives.txt`, then the command from spec section 1 (`--overwrite`).
- [ ] **Step 6: Commit** `feat(cli): similarity search --chains; foldseek writes per-chain files (#139)`.

### Task 2: chain map and `cluster --chain-representatives`

**Files:**
- Modify: `pandora/similarity/representatives.py`, `pandora/similarity/clustering.py`, `pandora/similarity/__init__.py`, `pandora/cli/similarity.py`
- Test: `tests/test_representatives.py`, `tests/test_clustering.py`, `tests/test_cli.py`

**Interfaces:**
- Produces:
  - `member_representatives(clusters: Iterable[RepresentativeCluster]) -> dict[str, str]` (chain item id → representative id; every member, the representative included).
  - `interface_edges(edges, interface_ids, chain_map: Mapping[str, str] | None = None)`: every chain (from the interfaces and from the edges) is keyed by `chain_map.get(key, key)`.
  - `_chain_map(command: str, path: Path, chains: Iterable[str]) -> dict[str, str]` in `pandora/cli/similarity.py`: loads the JSON (`load_json_list(RepresentativeCluster, ...)`), exits `"<command>: chain <id> is in no members list of <path>"` for the first missing chain.

- [ ] **Step 1: Write the failing tests.**
  - `test_member_representatives_maps_every_member`: `[RepresentativeCluster(representative="1A_A", members=["1A_A", "1A_A-2", "2B_C"], n_members=3)]` → `{"1A_A": "1A_A", "1A_A-2": "1A_A", "2B_C": "1A_A"}`.
  - `test_interface_edges_chain_map_links_members`: interfaces `["1AAA_A_B", "2BBB_C_D"]`, edges `[]`, map `{"1AAA_A": "9R_X", "2BBB_C": "9R_X"}` → yields a link between the two interfaces (shared representative, even though `9R_X` is in no interface).
  - `test_interface_edges_chain_map_maps_edge_ends`: edge `("9R_X", "8R_Y")`, map sends `1AAA_A → 9R_X` and `2BBB_D → 8R_Y` → the two interfaces are linked.
  - CLI `test_cluster_chain_representatives_joins_members`: reuse `_rollup_contacts`; a representatives JSON mapping `1AAA_A` → `9R_X` and `2BBB_D` → `8R_Y` (every other interface chain to itself), and a TM layer with the single edge `9R_X–8R_Y` → `[["1AAA_A_B", "2BBB_C_D"], ["3CCC_E_F"]]`.
  - `test_cluster_chain_representatives_missing_chain_exits` (match the chain id) and `test_cluster_chain_representatives_needs_interface_items`.
- [ ] **Step 2: Run** `uv run pytest tests/test_representatives.py tests/test_clustering.py tests/test_cli.py -k "member_representatives or chain_map or chain_representatives" -v`. Expected: FAIL.
- [ ] **Step 3: Implement** the three functions; in `_cmd_cluster` add flag `--chain-representatives` (help: `"representatives.json (pandora representatives): map each chain to its representative before rolling chain layers up to interfaces"`), exit unless `items_level == "interface"`, build the map over the chains of `item_ids` and pass it to `interface_edges`. Record the path in `cluster_provenance.json` only if `ClusteringProvenance` has a natural field; otherwise add `chain_representatives: str | None = None` to it.
- [ ] **Step 4: Run** the same command, then `uv run pytest`. Expected: PASS.
- [ ] **Step 5: Update** `05_dimer_clusters.sh`: add `--chain-representatives $D/clusters/seq90/representatives.json` to `pandora cluster`.
- [ ] **Step 6: Commit** `feat(cli): cluster --chain-representatives maps chains to their 90% representative`.

### Task 3: rebalance logic (#136)

**Files:**
- Create: `pandora/similarity/rebalance.py`
- Modify: `pandora/schemas/similarity.py`, `pandora/similarity/__init__.py`
- Test: `tests/test_rebalance.py`

**Interfaces:**
- Consumes: `parse_item_id`, `SimilarityCluster`.
- Produces:
  - Models: `PositiveDimer(dimer: str, entry_id: str, chain_id_1: str, chain_id_2: str, group: str, kind: Literal["homomer", "heteromer"], degree: int)`; `DroppedDimer(dimer: str, group: str, reason: Literal["homomer_in_mixed_group", "similar_to_kept", "not_most_connected"], kept: str | None = None)`.
  - `dimer_adjacency(interface_ids: Iterable[str], chain_edges: Iterable[Sequence[Any]], dimer_edges: Iterable[Sequence[Any]], chain_map: Mapping[str, str] | None = None) -> dict[str, set[str]]`: every interface id is a key; rules from spec section 2 "The dimer graph" (shared chain, chain-edge expansion, interface edges; ids mapped by `chain_map.get(k, k)` before all of them); no self-loops.
  - `rebalance(groups: Sequence[SimilarityCluster], adjacency: Mapping[str, set[str]], homomers: set[str], sort_key: Callable[[str], Any]) -> tuple[list[PositiveDimer], list[DroppedDimer]]`. `group` = the group's first component. Raises `ValueError` naming both dimers if an adjacency edge crosses groups.

- [ ] **Step 1: Write the failing tests** in `tests/test_rebalance.py`:
  - `test_adjacency_shared_chain`: `1A_A_B`, `1A_B_C` adjacent, `1A_D_E` not.
  - `test_adjacency_expands_chain_edge`: chain edge `(1A_A, 2B_C)` links all dimers holding `1A_A` to all holding `2B_C`.
  - `test_adjacency_interface_edge` and `test_adjacency_chain_map` (two dimers whose chains share a representative are adjacent; copies `1A_A-2` mapped to `1A_A`).
  - `test_homomer_only_group_keeps_most_connected`: star graph centre kept, others `not_most_connected`.
  - `test_mixed_group_drops_homomers`.
  - `test_heteromer_path_keeps_both_ends`: path `X–Y–Z` → keeps `X` and `Z` (`Y` dropped `similar_to_kept`, `kept="X"`), `degree` of `Z` is 0 (recounted).
  - `test_tie_uses_sort_key`: two isolated heteromers with equal degree; `sort_key` puts `B` first → `B` kept first (order of `positives` follows pick order).
  - `test_single_dimer_group`: one positive, no drops.
  - `test_cross_group_edge_raises`: `pytest.raises(ValueError, match="1A_A_B")`.
- [ ] **Step 2: Run** `uv run pytest tests/test_rebalance.py -v`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement.** Heteromer loop body (the one algorithmic choice):

```python
remaining = set(heteromers)
while remaining:
    degree = {d: len(adjacency[d] & remaining) for d in remaining}
    kept = min(remaining, key=lambda d: (degree[d], sort_key(d)))
    # record kept with degree[kept]; drop adjacency[kept] & remaining
    remaining -= adjacency[kept] | {kept}
```

  Homomer-only: `max` by degree within the group, ties by `sort_key` (smallest wins). `# ponytail: O(n²) per group; fine for group sizes seen, a heap if groups reach thousands.`
- [ ] **Step 4: Run** `uv run pytest tests/test_rebalance.py -v`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(similarity): dimer adjacency and homomer/heteromer rebalancing`.

### Task 4: `pandora rebalance` CLI

**Files:**
- Create: `pandora/cli/sampling.py` (with `register()`; add it to `pandora/cli/app.py` and `pandora/cli/__init__.py` like the other modules)
- Modify: `pandora/cli/similarity.py` (`_interface_values` takes `command: str`), `pandora/schemas/similarity.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: Task 2 `_chain_map`, `member_representatives`; Task 3 `dimer_adjacency`, `rebalance`, models; from `pandora/cli/similarity.py`: `_hit_filters`, `_layer_edges`, `_similarity_level`, `_interface_values`, `_contact_layer`; `rank_sort_key` from `pandora._util`.
- Produces: step folder `positives.<fmt>`, `dropped.<fmt>`, `rebalance.json` (step `"rebalance"`); model `RebalanceSettings(rank: RepresentativePolicy, layers: list[ClusteredLayer], clusters: str, chain_clusters: str, chain_representatives: str | None = None)` passed as the summary's policy.

- [ ] **Step 1: Write the failing tests** (build inputs with `_iface_rep_inputs` / `_rollup_contacts` helpers):
  - `test_rebalance_writes_positives_and_dropped`: one mixed group; check both tables' rows and `rebalance.json` counts (`groups`, `dimers_in`, `positives:homomer`, `positives:heteromer`, `dropped:<reason>`) and `previous_summary` = the contact layer's `layer.json`.
  - `test_rebalance_chain_missing_from_chain_clusters_exits` (match chain id).
  - `test_rebalance_dimer_missing_from_contacts_exits`.
  - `test_rebalance_refuses_existing_output` and `..._overwrite_replaces`.
  - `test_rebalance_inputs_mismatch_exits`: give a `--similarity` layer that links two groups → exit naming both dimers.
- [ ] **Step 2: Run** `uv run pytest tests/test_cli.py -k rebalance -v`. Expected: FAIL.
- [ ] **Step 3: Implement** `_cmd_rebalance` with flags exactly as spec section 2 plus `--min-score`, `--min-coverage`, `--table-format`, `--overwrite`. Homomers: map each chain item to its index in `--chain-clusters` (`load_json_list(SimilarityCluster, ...)`); missing chain → exit. Chain-level layers go to `chain_edges`, interface-level to `dimer_edges` (by `_similarity_level`). `sort_key = lambda d: rank_sort_key(d, values[d], policy.rank)`. `ValueError` from `rebalance` → `SystemExit("rebalance: ...")`.
- [ ] **Step 4: Run** the same command, then `uv run pytest`. Expected: PASS.
- [ ] **Step 5: Write** `examples/ppi/v2.1/src/06_positives.sh` with the spec section 2 command (same layers/filters as `05`).
- [ ] **Step 6: Commit** `feat(cli): rebalance picks positive dimers per group (#136)`.

### Task 5: negative sampling logic (#137)

**Files:**
- Create: `pandora/similarity/negatives.py`
- Modify: `pandora/schemas/similarity.py`, `pandora/similarity/__init__.py`
- Test: `tests/test_negatives.py`

**Interfaces:**
- Consumes: `PositiveDimer`, `item_id`.
- Produces:
  - Model `NegativePair(chain_1: str, chain_2: str, positive: str)`.
  - `sample_negatives(positives: Sequence[PositiveDimer], representative_of: Mapping[str, str], tm_edges: Iterable[Sequence[Any]], binding_pairs: Iterable[tuple[str, str]], seed: int) -> tuple[list[NegativePair], list[str]]`: negatives and unpaired positive ids. Raises `ValueError` naming a positive's chain or a binding-pair chain missing from `representative_of`. TM edge ends missing from `representative_of` are kept as-is (they can't be in the pool).

- [ ] **Step 1: Write the failing tests** in `tests/test_negatives.py` (hand-made reps `R1…R8`, identity map plus a few members):
  - `test_partner_never_negative`: positive (A, B); B's rep never chosen, whatever the seed (loop seeds 0–20).
  - `test_sets_one_to_three_excluded`: TM `R1–R2`, binding `R2–R3`, TM `R3–R4`; positive on `R1`'s chain → C ∉ {R1…R4}.
  - `test_c_used_once`: two positives, pool of exactly two candidates → two different C.
  - `test_pool_exhaustion_records_unpaired`.
  - `test_same_seed_same_output`; `test_reversed_and_duplicate_binding_rows_count_once` (result equals the deduplicated input's).
  - `test_tm_edge_outside_reps_mapped_through_members`.
  - `test_unmapped_positive_chain_raises`.
  - `test_empty_positives`: `([], [])`.
- [ ] **Step 2: Run** `uv run pytest tests/test_negatives.py -v`. Expected: FAIL.
- [ ] **Step 3: Implement** per spec section 3 "Algorithm": `rng = random.Random(seed)`; for positives sorted by `dimer`: `A = rng.choice([chain_1_id, chain_2_id])` (chain item ids), then `C = rng.choice(sorted(candidates))`; `chain_1` = A's chain item id, `chain_2` = C.
- [ ] **Step 4: Run** `uv run pytest tests/test_negatives.py -v`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(similarity): seeded negative sampling`.

### Task 6: `pandora negatives` CLI

**Files:**
- Modify: `pandora/cli/sampling.py`, `pandora/schemas/similarity.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: Task 4's step folder (`find_table(dir, "positives")`, read as `PositiveDimer`), Task 5 `sample_negatives`, `member_representatives`, `_hit_filters`, `_layer_edges`, `_contact_layer`.
- Produces: step folder `negatives.<fmt>`, `negatives.json` (step `"negatives"`, `previous_summary` = `rebalance.json`, counts `pool`, `positives`, `negatives`, `unpaired`); model `NegativeSamplingSettings(seed: int, hit_filter: HitFilter, similarity_layer: str, contacts: str, representatives: str, unpaired: list[str])` as the summary's policy.

- [ ] **Step 1: Write the failing tests:**
  - `test_negatives_writes_pairs_and_summary` (seed recorded, counts right).
  - `test_negatives_same_seed_reproduces` (two runs, identical tables).
  - `test_negatives_needs_chain_layer`: an interface layer as `--similarity` → exit.
  - `test_negatives_empty_contacts_layer` (layer with no rows → runs, all non-excluded reps eligible).
  - `test_negatives_refuses_existing_output`.
- [ ] **Step 2: Run** `uv run pytest tests/test_cli.py -k negatives -v`. Expected: FAIL.
- [ ] **Step 3: Implement** `_cmd_negatives` with flags from spec section 3 plus `--min-score`, `--min-coverage`, `--table-format`, `--overwrite`. Exactly one `--similarity`. Binding pairs from the contacts layer's `entry_id`, `chain_id_1`, `chain_id_2` → chain item ids. Missing positives table → exit naming `pandora rebalance`.
- [ ] **Step 4: Run** the same command, then `uv run pytest`. Expected: PASS.
- [ ] **Step 5: Write** `examples/ppi/v2.1/src/07_negatives.sh` with the spec section 3 command.
- [ ] **Step 6: Commit** `feat(cli): negatives samples one negative pair per positive (#137)`.

### Task 7: docs and the build

**Files:**
- Modify: `pandora/cli/README.md`, `docs/usage/similarity.md`, `docs/reference/policies.md` (only if it documents `RepresentativePolicy`/`HitFilter`: add `NegativeSamplingSettings.seed`), `examples/ppi/v2.1/README.md`, `examples/ppi/benchmarks/ppi-v2.1-benchmark.md`
- If `ClusteringProvenance` changed: `uv run --extra docs python docs/scripts/generate_erd.py`

- [ ] **Step 1:** CLI README: `--chains`, `--chain-representatives`, `rebalance`, `negatives`, in the existing paragraph style. Usage docs: full inline examples, library + CLI tabs (project convention).
- [ ] **Step 2:** Benchmark doc: "With Pandora (#139/#136/#137)" notes for Steps 3, 7 and 8, listing the spec's **(my call)** choices. v2.1 README: script table rows for 03, 06, 07; remove "waiting for #139".
- [ ] **Step 3: Run** `uv run zensical build --clean` (use the docs-build skill), `uv run ruff format --check .`, `uv run ruff check .`, `uv run pytest`. Expected: all clean.
- [ ] **Step 4: Run** the v2.1 scripts 01–07 on the dev fixtures with the real `mmseqs` (Foldseek if on `PATH`, else skip 03 and note it). Expected: `work/positives` and `work/negatives` written; report counts.
- [ ] **Step 5: Commit** `docs: rebalance, negatives, search --chains, cluster --chain-representatives`.
