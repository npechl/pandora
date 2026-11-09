# Annotation Layers Part 5 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Drop labelled interfaces from a `residue_contacts` layer,
writing a new contact layer, so interface records and `cluster
--contacts` follow.

**Architecture:** Two plain functions in a new
`pandora/annotations/filters.py`: one returns the kept contact rows,
one returns the counts. A thin `pandora annotate filter-contacts`
command reads two layer folders, checks they come from the same
structures folder, and writes an ordinary layer folder, as
`count-contacts` does.

**Tech Stack:** Python 3, pydantic v2, argparse, pytest.

**Spec:** `design/specs/2026-11-06-annotation-layers-design.md`,
section "Part 5: interfaces are dropped by label".

## Global Constraints

- `from __future__ import annotations` at the top of new modules; type
  hints everywhere; line length 80; `uv run ruff format`.
- Public functions: Google-style docstring, exported in
  `pandora/annotations/__init__.py` `__all__`.
- Never mutate inputs. No new dependency. Tests offline, fixtures from
  `datasets/dev/mmcif/` only.
- An interface is `(entry_id, unordered {chain_id_1, chain_id_2})` on
  Pandora chain IDs. Entry ids compare case-insensitively (import
  matches them that way). No special case for intra-chain or ligand
  rows: (A, A) is a pair like any other.
- Unlabelled interfaces are **kept**. An interface is dropped if any of
  its label rows has a label in `drop`.
- A `drop` value present on no label row is a `ValueError` (library)
  and an exit before anything is written (CLI).
- Import is not changed.

## Review Focus

1. Contact layer with no rows (only `layer.json`, no table) → writes
   only `layer.json`; all counts 0; no crash. Test in Task 2.
2. Entry id case differs between the layers (`1A2C` labels, `1a2c`
   contacts) → still matches. Test in Task 1.
3. `--labels` is a layer of another template (e.g. `chain_labels`) →
   exits before writing, naming the template. Test in Task 2.
4. `--output-dir` equal to the `--labels` folder (not only the
   `--contacts` one) → exits before writing, input untouched. Test in
   Task 2.
5. Chained runs: the `--contacts` folder is itself a
   `filter-contacts` output → its origin is traced through it to the
   structures folder, so a second filter runs. Test in Task 2.

---

### Task 1: Filter functions

**Files:**
- Create: `pandora/annotations/filters.py`
- Modify: `pandora/annotations/__init__.py`
- Test: `tests/test_annotations.py` (reuse its `_contact()` helper)

**Interfaces:**
- Produces:
  - `drop_labelled_interfaces(contacts: Iterable[ResidueContactRow],
    labels: Iterable[LayerRow], drop: Collection[str]) ->
    list[ResidueContactRow]`: kept rows, input order. Raises
    `ValueError` naming the unknown `drop` values.
  - `count_interface_labels(contacts: Iterable[ResidueContactRow],
    labels: Iterable[LayerRow], drop: Collection[str]) -> dict[str,
    int]` with keys exactly `interfaces_in`, `interfaces_dropped`,
    `interfaces_unlabelled`, `labels_unmatched`. Does not raise for
    unknown `drop` values (the CLI calls the other function first).
  - internal `_interface(entry_id, chain_1, chain_2) -> tuple[str,
    frozenset[str]]` (upper-cased entry id); a label row's key comes
    from its `entry_id`, `chain_id_1`, `chain_id_2` attributes.

Definitions for the counts: `interfaces_in` = distinct interfaces in
`contacts`; `interfaces_dropped` = those of them with a `drop` label;
`interfaces_unlabelled` = those of them with no label row at all;
`labels_unmatched` = distinct labelled interfaces not in `contacts`.

- [ ] **Step 1: Write the failing tests**

Labels built as `LayerRow(entry_id=..., chain_id_1=..., chain_id_2=...,
label=...)`.

- `test_drop_labelled_interfaces_drops_both_orders`: contacts
  `[A1–B5, B6–A2, A1–C3]`, label (B, A) `low_hbond`, `drop={"low_hbond"}`
  → result is `[A1–C3]`.
- `test_drop_labelled_interfaces_keeps_order_and_unlabelled`: contacts
  `[A–C, A–B, B–C]`, label (A, B) `bad` → `[A–C, B–C]` in that order.
- `test_drop_labelled_interfaces_any_label_drops`: labels (A, B) `ok`
  and (B, A) `bad`, `drop={"bad"}` → no A–B rows.
- `test_drop_labelled_interfaces_self_pair_and_ligand`: contacts
  `[A1–A9, A1–J1]` (J a ligand chain, `seq_id_2=None`), labels (A, A)
  `x` and (A, J) `y`: `drop={"x"}` → `[A1–J1]`; `drop={"y"}` →
  `[A1–A9]`.
- `test_drop_labelled_interfaces_entry_case`: contacts with
  `entry_id="1abc"`, label (A, B) `bad` with `entry_id="1ABC"` → the
  A–B rows are dropped.
- `test_drop_labelled_interfaces_unknown_label_raises`:
  `pytest.raises(ValueError, match="typo")` for `drop={"bad", "typo"}`.
- `test_count_interface_labels`: contacts A–B, A–C, B–C; labels
  (A, B) `bad`, (A, C) `ok`, (C, D) `bad`; `drop={"bad"}` → `{
  "interfaces_in": 3, "interfaces_dropped": 1,
  "interfaces_unlabelled": 1, "labels_unmatched": 1}`.
- `test_filter_functions_do_not_mutate_input`: `model_dump()` of every
  input row is unchanged after both calls.

- [ ] **Step 2: Run them, expect ImportError**

Run: `uv run pytest tests/test_annotations.py -k "labelled or interface_labels or filter_functions" -v`

- [ ] **Step 3: Implement both functions in `filters.py`; export them**

Module docstring: `"""Filters: drop labelled interfaces from a
\`residue_contacts\` layer."""`. Build a `dict[key, set[label]]` from
`labels` once.

- [ ] **Step 4: Run them, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add pandora/annotations tests/test_annotations.py
git commit -m "feat(annotations): drop labelled interfaces from contact rows"
```

### Task 2: `pandora annotate filter-contacts`

**Files:**
- Modify: `pandora/cli/annotate.py` (new `_cmd_filter`, `_layer_origin`,
  `register()`)
- Test: `tests/test_cli.py` (reuse `_step_folder`, `_compute`; write
  label tables with `_write_jsonl`; import them with `main(["annotate",
  "import", ...])` against the same step folder)

**Interfaces:**
- Consumes: `drop_labelled_interfaces`, `count_interface_labels`
  (Task 1); `check_contact_layer`, `find_table`, `read_records`,
  `write_table`, `prepare_output_dir`, `build_step_summary`,
  `write_json`, `_LAYER_OUTPUTS` (existing).
- Produces: `_layer_origin(folder: Path) -> Path | None`.

Command: `pandora annotate filter-contacts --contacts <dir> --labels
<dir> --drop <label> [--drop ...] --output-dir <dir> [--overwrite]
[--table-format]`; `--drop` is `action="append"`, required.

Order inside `_cmd_filter`, everything before `prepare_output_dir`
exits with `SystemExit` and writes nothing:
1. `check_contact_layer("annotate filter-contacts", contacts)`.
2. `--labels` must have `layer.json` with `policy.template ==
   "interface_labels"`, else exit `"annotate filter-contacts: <dir> is
   a <template> layer, --labels needs interface_labels"`.
3. `--output-dir` resolving to either input folder → exit `"... <out>
   is the --contacts folder"` / `"... is the --labels folder"`.
4. Origins: if both `_layer_origin` values are not None and differ →
   exit `"annotate filter-contacts: --contacts comes from <a>,
   --labels from <b>; both layers must come from the same structures
   folder"`. If either is None → print a warning to stderr naming the
   layer, and carry on.
5. Read rows (`find_table` → `read_records(.., ResidueContactRow)` /
   `LayerRow`; no table → `[]`); call `drop_labelled_interfaces`,
   turning its `ValueError` into an exit.

`_layer_origin(folder)`: load `folder/layer.json`; while
`previous_summary` names a file called `layer.json`, follow it. Then
return the **parent folder** of that non-layer `previous_summary`
(resolved). If `previous_summary` is None: return the resolved
`input` when `step == "annotation_compute"`, else None.

Output folder, as `count-contacts`: `layer.<fmt>` (kept rows) and
`layer.json` = `build_step_summary("annotation_filter",
LayerDescription(name="filtered_contacts", template="residue_contacts",
tool="pandora", tool_version=__version__, parameters={
"contacts_layer": <contacts policy name>, "labels_layer": <labels
policy name>, "labels_dir": str(labels), "drop": sorted(drop)}),
input=contacts, output=out, started_at=..., counts={"rows_in",
"rows_out", **count_interface_labels(...)} sorted by key,
previous_summary=contacts / "layer.json")`.

Printed line: `filtered: <dropped> interfaces dropped, <kept> kept,
<unlabelled> unlabelled -> <out>` (kept = in − dropped).

- [ ] **Step 1: Write the failing tests**

Fixture setup helper `_filter_inputs(tmp_path, labels_rows)`: step
folder from `_step_folder` plus `curated/curation.json` = `"{}"`;
`_compute` a `chain_interfaces` layer into `annotations/interfaces`;
import `labels_rows` with a description of template
`interface_labels` into `annotations/verdict`; return the two folders.
Labels use the chain pairs the computed layer contains (read them from
its `layer.jsonl`).

- `test_annotate_filter_contacts_writes_layer`: label one interface
  `bad`, one `ok`; run with `--drop bad`. Written rows equal
  `drop_labelled_interfaces(...)` on the input rows; no row of the
  `bad` pair remains; summary `step == "annotation_filter"`,
  `policy.name == "filtered_contacts"`, `policy.template ==
  "residue_contacts"`, `policy.parameters["drop"] == ["bad"]`,
  `previous_summary == str(contacts / "layer.json")`, `counts` has the
  six keys with `interfaces_dropped == 1`.
- `test_annotate_filter_contacts_output_feeds_records`:
  `extract_interface_records` on the output rows has no record for
  the `bad` pair.
- `test_annotate_filter_contacts_chains`: a second run with
  `--contacts` = the first output and the same `--labels` succeeds
  (origin traced through) and drops nothing more.
- `test_annotate_filter_contacts_rejects_wrong_labels_template`: a
  `chain_labels` layer as `--labels` → `SystemExit` matching
  `"needs interface_labels"`; output folder absent.
- `test_annotate_filter_contacts_rejects_unknown_drop`: `--drop typo`
  → `SystemExit` matching `"typo"`; output folder absent.
- `test_annotate_filter_contacts_refuses_input_as_output`: parametrise
  over `--contacts` and `--labels` folders; `SystemExit`; that folder's
  table unchanged.
- `test_annotate_filter_contacts_rejects_different_origins`: import
  the labels against a copy of the step folder (`shutil.copytree` to
  `other/`, with its own `curation.json`) → `SystemExit` matching
  `"same structures folder"`; output folder absent.
- `test_annotate_filter_contacts_warns_without_origin`: import the
  labels against a step folder **without** `curation.json` → runs;
  `capsys` stderr contains `"warning"`.
- `test_annotate_filter_contacts_empty_contacts`: a hand-written
  contacts folder holding only `layer.json` (as
  `test_annotate_count_contacts_empty_layer`, plus `"step":
  "annotation_compute"`, `"input": <step folder>`) → output folder
  holds only `layer.json`; its `counts` has all six keys, each 0
  except `labels_unmatched` (= the number of labelled interfaces).

- [ ] **Step 2: Run them, expect failure** (unknown subcommand)

Run: `uv run pytest tests/test_cli.py -k filter_contacts -v`

- [ ] **Step 3: Implement `_layer_origin`, `_cmd_filter` and the
  parser entry** (help: `"Drop labelled interfaces from a contact
  layer."`)

- [ ] **Step 4: Run them, expect PASS**

- [ ] **Step 5: Commit**

```bash
git add pandora/cli/annotate.py tests/test_cli.py
git commit -m "feat(cli): annotate filter-contacts drops labelled interfaces"
```

### Task 3: Docs

**Files:**
- Modify: `docs/usage/annotation.md` (new `### Dropping interfaces by
  label` subsection at the end of `## Using layers in a dataset`)
- Modify: `docs/usage/cli.md`, `pandora/cli/README.md`
  (`annotate filter-contacts`)
- Modify: `docs/contributing/architecture.md` (both functions:
  `ResidueContactRow` + `LayerRow` in; kept rows / counts out)

Page convention: full code inline, `library` / `cli` tabs, each
followed by a `text` block with the **real** output (run by hand on
1a2c and paste). The section shows: a contact layer with hand-written
`contact_type` values (say so: the built-in layers leave it empty),
hydrogen bonds counted per interface in pandas with the chain pair
sorted, a labels table (`low_hbond` / `ok`), `annotate import` as
`interface_labels`, `annotate filter-contacts --drop low_hbond`, the
summary's `counts` via `jq`, and `cluster --contacts` /
`extract_interface_records()` reading the new layer. One sentence each
on: unlabelled interfaces are kept; one label per interface, so
several criteria are folded into one label or applied by chained
runs; both layers must come from the same structures folder.

- [ ] **Step 1: Write the section and the reference edits**
- [ ] **Step 2: Run every new example by hand; paste real output**
- [ ] **Step 3: Build the docs** with the `docs-build` skill; no new
  warnings
- [ ] **Step 4: Final checks:** `uv run pytest`,
  `uv run ruff format --check .`, `uv run ruff check .`
- [ ] **Step 5: Commit**

```bash
git add docs pandora/cli/README.md
git commit -m "docs: dropping interfaces by label, annotate filter-contacts"
```
