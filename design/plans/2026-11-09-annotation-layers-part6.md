# Annotation Layers Part 6 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `partition_dataset()` / `pandora partition` take a label per
item and place whole clusters so every label appears in every split
where leakage safety allows, reporting the labels that couldn't.

**Architecture:** A coverage pass in `partition_dataset()` places the
smallest clusters of the rarest labels first, then today's size pass
places the rest. `item_labels()` turns a `chain_labels` /
`interface_labels` layer into an item ID → label mapping using
`item_id()`. `pandora partition --labels <layer dir>` wires them up.

**Tech Stack:** Python 3, pydantic v2, argparse, pytest.

**Spec:** `design/specs/2026-11-06-annotation-layers-design.md`,
section "Part 6: every kind reaches every split".

## Global Constraints

- `from __future__ import annotations`; type hints; line length 80;
  `uv run ruff format`. Public functions: Google-style docstring,
  exported in the package `__init__.py` `__all__`.
- Never mutate inputs. No new dependency. Tests offline.
- Without `labels`, `partition_dataset()` output is identical to
  today's.
- Clusters are never split. Coverage never makes the partition fail:
  missing labels are reported, not raised.
- Items without a label are the label `"unlabelled"`; a real label
  `"unlabelled"` is the same kind.
- `labels` with `keep_similar_items=False` → `ValueError`.
- A split whose fraction is 0 is skipped by the coverage pass and never
  appears in `missing_labels`.
- No structures folder is read by `partition`; no SplitPolicy (#107).

## Review Focus

1. Clusters file at another level than the labels (chain clusters,
   interface labels) → runs, everything `unlabelled`, printout says so.
   Test in Task 3.
2. `--pct-val 0` with labels → val gets nothing from the coverage pass
   and isn't listed as missing anything. Test in Task 1.
3. Interface label written B–A, cluster item `1ABC_A_B` → matches.
   Test in Task 2.
4. Labels layer with no table (only `layer.json`) → every item
   `unlabelled`, no crash. Test in Task 3.
5. Extended entry IDs (`pdb_00001abc`) → `item_id()` handles them;
   label for `pdb_00001abc_A` matches. Test in Task 2.

---

### Task 1: Coverage pass in `partition_dataset`

**Files:**
- Modify: `pandora/similarity/partition.py`
- Modify: `pandora/schemas/similarity.py` (`PartitionProvenance`)
- Test: `tests/test_partition.py`

**Interfaces:**
- Produces:
  - `partition_dataset(clusters, pct_train=0.6, pct_val=0.2,
    pct_test=0.2, keep_similar_items=True, labels: Mapping[str, str] |
    None = None) -> tuple[dict[str, list[str]], PartitionProvenance]`
  - `PartitionProvenance` new fields, with docstring lines:
    `label_counts: dict[str, dict[str, int]] = {}` (split → label →
    items; only labels present in that split), `missing_labels:
    dict[str, list[str]] = {}` (split → sorted labels absent; splits
    missing none are left out), `labels_layer: str | None = None` (set
    by the CLI only). All empty/None without `labels`.

Test data (verified by simulation; today's partition leaves
`heteromer` out of train and val):

```python
KINDS = [
    [f"a{i}" for i in range(10)],
    [f"b{i}" for i in range(8)],
    [f"c{i}" for i in range(6)],
    ["h0", "h1"],
    ["h2"],
    ["h3"],
    ["s0"],
    ["s1"],
    ["s2"],
    ["s3"],
]
KIND_LABELS = {i: "homomer" for c in KINDS for i in c} | {
    h: "heteromer" for h in ("h0", "h1", "h2", "h3")
}
```

- [ ] **Step 1: Write the failing tests**

- `test_labels_reach_every_split`: plain partition of `KINDS` has a
  split with no heteromer item (assert it, so the test means
  something); with `labels=KIND_LABELS`, every split has both labels,
  every cluster is in exactly one split, all 32 items placed,
  `provenance.missing_labels == {}`, and
  `sum(provenance.label_counts[s]["heteromer"] for s in splits) == 4`.
- `test_rare_label_missing_is_reported`: one heteromer cluster only →
  no exception; `missing_labels` has exactly two splits, each
  `["heteromer"]`.
- `test_unlabelled_items_are_a_kind`: labels cover only some items →
  `"unlabelled"` appears in `label_counts`; an item labelled
  `"unlabelled"` counts the same.
- `test_no_labels_is_unchanged`: for `CLUSTERS` and `KINDS`, result
  with `labels=None` equals today's (compare against the hard-coded
  `[18, 6, 6]` and the existing tests still pass);
  `label_counts == {}`, `missing_labels == {}`.
- `test_zero_fraction_split_skipped`: `KINDS`, `pct_val=0`,
  `pct_train=0.8` → `split_sizes["val"] == 0`, `"val"` not in
  `missing_labels`, train and test both hold heteromer.
- `test_labels_need_keep_similar_items`: `pytest.raises(ValueError)`.
- `test_labels_input_not_mutated`: clusters and labels equal their
  deep copies after the call.

- [ ] **Step 2: Run them, expect failures**

Run: `uv run pytest tests/test_partition.py -v`

- [ ] **Step 3: Implement**

Coverage pass, before today's loop; today's loop then skips placed
clusters (track by index, not identity of equal lists):

```text
kinds[k]  = {labels.get(item, "unlabelled") for item in cluster k}
active    = splits whose fraction > 0
order     = all labels sorted by (number of clusters carrying it, name)
for label in order:
    for k in unplaced clusters carrying label,
             sorted by (n_components, components[0]):
        lacking = active splits whose covered set lacks label
        if not lacking: break
        split = max(lacking, key=target - current size)  # first wins ties
        place cluster k in split; covered[split] |= kinds[k]
```

Then `label_counts` from the final splits, `missing_labels` over
`active` splits. Docstring: describe `labels`, the coverage guarantee
("where leakage safety allows") and the new `Raises:` case.

- [ ] **Step 4: Run them, expect PASS** (whole file)

- [ ] **Step 5: Commit**

```bash
git add pandora/similarity/partition.py pandora/schemas/similarity.py tests/test_partition.py
git commit -m "feat(similarity): partition places every label in every split"
```

### Task 2: `item_labels`

**Files:**
- Modify: `pandora/similarity/layer.py`, `pandora/similarity/__init__.py`
- Test: `tests/test_similarity_layer.py`

**Interfaces:**
- Consumes: `item_id(level, entry_id, *chain_ids)` (same file).
- Produces: `item_labels(rows: Iterable[LayerRow | Mapping[str, Any]],
  template: str) -> dict[str, str]`. `chain_labels` rows →
  `item_id("chain", entry_id, chain_id)`; `interface_labels` →
  `item_id("interface", entry_id, chain_id_1, chain_id_2)`. Read rows
  with `dict(row)` as `filters.py` does. Raises `ValueError` for any
  other template, and for one item ID given two different labels
  (message names the ID and both labels). The same label twice is fine.

- [ ] **Step 1: Write the failing tests**

- `test_item_labels_chains`: rows `(1ABC, A, x)`, `(1ABC, B, y)` →
  `{"1ABC_A": "x", "1ABC_B": "y"}`.
- `test_item_labels_interface_order`: row `(1ABC, B, A, het)` →
  `{"1ABC_A_B": "het"}`.
- `test_item_labels_extended_entry`: `(pdb_00001abc, A, x)` →
  `{"pdb_00001abc_A": "x"}`.
- `test_item_labels_conflict_raises`: interface rows B–A `x` and A–B
  `y` → `ValueError` matching `"1ABC_A_B"`.
- `test_item_labels_wrong_template_raises`: `"residue_labels"` →
  `ValueError`.

- [ ] **Step 2: Run them, expect ImportError**

Run: `uv run pytest tests/test_similarity_layer.py -k item_labels -v`

- [ ] **Step 3: Implement and export `item_labels`**

- [ ] **Step 4: Run them, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add pandora/similarity tests/test_similarity_layer.py
git commit -m "feat(similarity): item_labels maps a label layer to item IDs"
```

### Task 3: `pandora partition --labels`

**Files:**
- Modify: `pandora/cli/similarity.py` (`_cmd_partition`, `register()`)
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `partition_dataset(..., labels=)`, `item_labels(rows,
  template)`, `find_table`, `read_records(table, LayerRow)`.

Behaviour:
- `--labels <dir>` optional. Not a layer folder (no `layer.json`) or
  `policy.template` not in `chain_labels` / `interface_labels` →
  `SystemExit` naming the template, before writing. `--labels` with
  `--no-keep-similar-items` → `SystemExit` before writing.
  `ValueError` from `item_labels` → `SystemExit`.
- No table in the layer folder → empty rows (all items unlabelled).
- Provenance written with `labels_layer=str(labels_dir)` via
  `model_copy`.
- Printout: today's line, then with `--labels`:
  `unlabelled items: <n> of <total>`, then one stderr line per split
  with missing labels: `partition: warning: <split> has no item
  labelled <a>, <b>`.

CLI tests write the labels layer folder by hand: `layer.json` =
`{"policy": {"name": "kinds", "template": "chain_labels"}}` plus
`layer.jsonl` (`_write_jsonl`); the clusters file is hand-written
JSON.

- [ ] **Step 1: Write the failing tests**

- `test_partition_labels_writes_counts`: chain clusters, chain labels
  → `partition_provenance.json` has `label_counts`, `missing_labels`,
  `labels_layer`; `splits.json` keys unchanged.
- `test_partition_labels_warns_missing(capsys)`: one heteromer cluster
  → exit 0, stderr has `has no item labelled heteromer` twice.
- `test_partition_labels_wrong_level(capsys)`: chain clusters,
  `interface_labels` layer → runs; stdout has `unlabelled items: N of
  N`.
- `test_partition_labels_empty_layer`: `layer.json` only → runs, all
  unlabelled.
- `test_partition_labels_rejects_template`: `residue_values` layer →
  `SystemExit`, no `splits.json`.
- `test_partition_labels_rejects_no_keep`: `SystemExit`, no
  `splits.json`.

- [ ] **Step 2: Run them, expect failures**

Run: `uv run pytest tests/test_cli.py -k partition_labels -v`

- [ ] **Step 3: Implement**

- [ ] **Step 4: Run them, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add pandora/cli/similarity.py tests/test_cli.py
git commit -m "feat(cli): partition --labels puts every label in every split"
```

### Task 4: Docs

**Files:**
- Modify: `docs/usage/similarity.md` (new `### Every kind in every
  split` at the end of `## Leakage-safe partitioning`)
- Modify: `docs/usage/cli.md`, `pandora/cli/README.md` (`partition
  --labels`)
- Modify: `docs/contributing/architecture.md` (`item_labels`: label
  rows in, item ID → label out; `partition_dataset` takes `labels`)

Page convention: full code inline, `library` / `cli` tabs, each
followed by a `text` block with the **real** output (run by hand on
the dev fixtures and paste). Show: a chain labels table (homomer /
heteromer) written in pandas, `annotate import` as `chain_labels`,
`partition --labels`, the provenance's `label_counts` /
`missing_labels` via `jq`. One sentence each on: clusters stay whole,
so coverage is best effort and missing labels are reported, not
refused; unlabelled items are their own kind; label symmetry-copy
chains too (they sit in their source chain's cluster), e.g. by
copying the source's label in pandas.

- [ ] **Step 1: Write the section and the reference edits**
- [ ] **Step 2: Run every new example by hand; paste real output**
- [ ] **Step 3: Regenerate diagrams:**
  `uv run --extra docs python docs/scripts/generate_erd.py`
- [ ] **Step 4: Build the docs** with the `docs-build` skill; no new
  warnings
- [ ] **Step 5: Final checks:** `uv run pytest`,
  `uv run ruff format --check .`, `uv run ruff check .`
- [ ] **Step 6: Commit**, then close #95 and comment on #107 that
  `labels` is the first field waiting for `SplitPolicy` (ask the user
  before touching GitHub).
