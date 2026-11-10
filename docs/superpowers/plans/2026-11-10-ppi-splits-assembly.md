# PPI Splits and Assembly Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `rebalance` writes positive-only groups and a homomer/heteromer label layer; `partition --stratify` keeps the label mix per split; `pandora assemble` writes the final pairs table with splits and C-chain leak flags.

**Architecture:** Pure logic in `pandora/similarity/` (`partition.py` gains `stratify`; new `assembly.py`), models in `pandora/schemas/similarity.py`, CLI in `pandora/cli/sampling.py` (rebalance, assemble) and `pandora/cli/similarity.py` (partition).

**Tech Stack:** Python 3, pydantic v2, argparse, pytest.

**Spec:** `docs/superpowers/specs/2026-11-10-ppi-splits-assembly-design.md`

## Global Constraints

- `from __future__ import annotations`; type hints; line length 80; `ruff format`; ruff C901/PLR0915 limits (extract helpers).
- Public functions: Google-style docstring and `__all__` entry in `pandora/similarity/__init__.py`.
- Never mutate inputs; no new dependencies; item ids via `item_id` / `parse_item_id`.
- Errors name the offending id; no silent drops.
- Tests with local fixtures only; `examples/` untracked (edit scripts, don't commit).

## Review Focus

- A split with no negatives (all its positives unpaired): rows still written, counts 0.
- `--splits` from a plain `partition` run (not stratified) works too.
- `partition --stratify` with a label whose clusters are fewer than active splits: warning, no crash.
- `assemble` with an empty TM layer (no `layer.parquet`): all leak flags false.
- `--overwrite` on `rebalance` clears `groups.json` and `kinds/`.

Each has a test in the owning task (Tasks 1, 2, 3).

---

### Task 1: `rebalance` writes `groups.json` and `kinds/`

**Files:**
- Modify: `pandora/cli/sampling.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: in `--output-dir`: `groups.json` (`SimilarityCluster` list = `restrict_clusters(groups, positive ids)`), `kinds/layer.<fmt>` (rows `entry_id`, `chain_id_1`, `chain_id_2`, `label`) and `kinds/layer.json` (`build_step_summary("rebalance_kinds", LayerDescription(name="positive_kinds", template="interface_labels", tool="pandora", tool_version=__version__, parameters={}), ...)`, `previous_summary` = `rebalance.json`).

- [ ] **Step 1: Write failing tests** next to the other `test_rebalance_*`:
  - `test_rebalance_writes_groups_and_kinds`: `groups.json` components == the positives' dimers, grouped; `kinds/layer.jsonl` labels == positives' `kind`; then `main(["partition", "--clusters", groups, "--labels", kinds, "--output", ...])` succeeds.
  - `test_rebalance_overwrite_clears_groups_and_kinds`: after `--overwrite`, files rewritten (mtime or content), no stale files.
- [ ] **Step 2: Run** `uv run pytest tests/test_cli.py -k rebalance -q`. Expected: the two new tests FAIL.
- [ ] **Step 3: Implement** in `_cmd_rebalance` (extract a `_write_kinds(...)` helper). Add `"groups.json"` and `"kinds"` to `_REBALANCE_OUTPUTS`. The kinds table uses `--table-format`.
- [ ] **Step 4: Run** same command, then `uv run pytest -q`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(cli): rebalance writes positive-only groups and a homomer/heteromer layer`.

### Task 2: `partition --stratify`

**Files:**
- Modify: `pandora/similarity/partition.py`, `pandora/schemas/similarity.py` (`PartitionProvenance.stratify: bool = False`), `pandora/cli/similarity.py`
- Test: `tests/test_partition.py`, `tests/test_cli.py`

**Interfaces:**
- Produces: `partition_dataset(clusters, pct_train=0.6, pct_val=0.2, pct_test=0.2, keep_similar_items=True, labels=None, stratify=False)`. Raises `ValueError` when `stratify` without `labels`, with `keep_similar_items=False`, or when a cluster carries two labels (message names the cluster's first component).

- [ ] **Step 1: Write failing tests** in `tests/test_partition.py`:
  - `test_stratify_keeps_label_ratio`: 10 single-item clusters labelled `h`, 5 labelled `x` → per split `label_counts` train `{h: 6, x: 3}`, val `{h: 2, x: 1}`, test `{h: 2, x: 1}`; `prov.stratify is True`.
  - `test_stratify_needs_labels`, `test_stratify_rejects_mixed_cluster` (match the component id), `test_stratify_rejects_split_items`.
  - `test_stratify_rare_label_reports_missing`: 1 cluster labelled `x` → `missing_labels` has two splits lacking `x`, no exception.
  - CLI `test_partition_stratify_flag`: `--stratify` without `--labels` → `SystemExit` matching `--labels`.
- [ ] **Step 2: Run** `uv run pytest tests/test_partition.py tests/test_cli.py -k stratify -q`. Expected: FAIL.
- [ ] **Step 3: Implement.** With `stratify`: group clusters by their single label (sorted labels), run the existing whole-cluster placement per group with targets from that group's item count, merge the splits per name; label report as today. Flag `--stratify` (help: `"split each --labels label's clusters separately, keeping the label mix in every split"`); validation errors → `SystemExit("partition: ...")`.
- [ ] **Step 4: Run** same command, then `uv run pytest -q`. Expected: PASS.
- [ ] **Step 5: Commit** `feat(similarity): partition --stratify keeps each label's share per split`.

### Task 3: `pandora assemble`

**Files:**
- Create: `pandora/similarity/assembly.py`
- Modify: `pandora/schemas/similarity.py`, `pandora/similarity/__init__.py`, `pandora/cli/sampling.py`
- Test: `tests/test_assembly.py`, `tests/test_cli.py`

**Interfaces:**
- Consumes: `PositiveDimer`, `NegativePair`, `member_representatives`, `_positives` (sampling.py), `_layers`, `_layer_edges`.
- Produces:
  - Model `DatasetPair(chain_a: str, chain_b: str, label: int, kind: str, cluster: str, split: str, positive: str, c_leak_positive: bool | None = None, c_leak_negative: bool | None = None)`.
  - Model `AssemblySettings(positives: str, negatives: str, splits: str, representatives: str, similarity_layer: str, hit_filter: HitFilter)`.
  - `assemble_pairs(positives: Sequence[PositiveDimer], negatives: Sequence[NegativePair], splits: Mapping[str, Sequence[str]], representative_of: Mapping[str, str], tm_edges: Iterable[Sequence[Any]]) -> list[DatasetPair]`: positives then negatives, each sorted (by `positive`, then `chain_b`); `ValueError` for a positive in no split or two, a split item that's not a positive, a negative with an unknown positive, an unmapped chain.

- [ ] **Step 1: Write failing tests** in `tests/test_assembly.py` (positives `1P_A_B` in train, `2Q_A_B` in test; reps identity):
  - `test_negatives_inherit_split_kind_cluster`.
  - `test_c_leak_positive`: TM edge between train negative's C and a test positive's chain → `c_leak_positive is True` for that negative, `c_leak_negative is False`; positives have `None`.
  - `test_c_leak_negative`: TM edge between the train C and the test C → both negatives `c_leak_negative is True`.
  - `test_same_split_link_is_not_a_leak`.
  - `test_no_tm_edges_no_leaks`.
  - one `pytest.raises(ValueError, match=<id>)` per error in the spec.
- [ ] **Step 2: Run** `uv run pytest tests/test_assembly.py -q`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement** per spec section 3 "The leakage count": per split, sets of positive-chain reps and C reps; a negative leaks when `tm[rep(C)]` meets the other splits' set (TM adjacency undirected, ends mapped with `representative_of.get(x, x)`).
- [ ] **Step 4: Run** `uv run pytest tests/test_assembly.py -q`. Expected: PASS.
- [ ] **Step 5: Write failing CLI tests** reusing `_negatives_inputs` + a `partition` run on the rebalance `groups.json`: `test_assemble_writes_pairs_and_summary` (rows, `assemble.json` counts `pairs:<split>:<label>`, `negatives_c_leak_positive`, `negatives_c_leak_negative`, `positives_unpaired`; `previous_summary` = `negatives.json`), `test_assemble_warns_on_other_representatives` (capsys stderr has `warning`), `test_assemble_refuses_existing_output`.
- [ ] **Step 6: Implement** `_cmd_assemble` + parser (`--positives`, `--negatives`, `--splits`, `--representatives`, `_layer_flags` with one chain layer, `--output-dir`, `--table-format`, `--overwrite`; outputs `pairs.<fmt>`, `assemble.json`). Warning: compare `--representatives` and the layer's `layer.json` path with `negatives.json`'s `policy.representatives` / `policy.similarity_layer`; print `assemble: warning: ...` to stderr.
- [ ] **Step 7: Run** `uv run pytest -q`. Expected: PASS.
- [ ] **Step 8: Commit** `feat(cli): assemble writes the final pairs table with splits and leak flags`.

### Task 4: build script and docs

**Files:** `examples/ppi/v2.1/src/08_splits.sh` (untracked), `examples/ppi/v2.1/README.md`, `examples/ppi/benchmarks/ppi-v2.1-benchmark.md`, `docs/usage/similarity.md`, `pandora/cli/README.md`

- [ ] **Step 1:** `08_splits.sh`: `partition --clusters $D/work/positives/groups.json --labels $D/work/positives/kinds --stratify --output $D/splits.json`, then the spec's `assemble` command (`--table-format jsonl --overwrite`), then `jq -c .counts`.
- [ ] **Step 2: Run** scripts 06–08 on dev fixtures. Expected: `work/dataset/pairs.jsonl` written; report counts.
- [ ] **Step 3:** Docs with real output from a scratch copy (library + CLI tabs): `--stratify` under partitioning; "Final table" under "Positives and negatives (PPI)". CLI README sentences. Benchmark doc "With Pandora" under Splits and Output, leakage-count note under limitation F. v2.1 README rows.
- [ ] **Step 4: Run** `uv run zensical build --clean --strict`, `uv run ruff format --check .`, `uv run ruff check .`, `uv run pytest -q`; regenerate ERD (`uv run --extra docs python docs/scripts/generate_erd.py`). Expected: all clean.
- [ ] **Step 5: Commit** `docs: partition --stratify and assemble`.
