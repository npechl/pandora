# Annotation layers, part 1 (foundation) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `pandora annotate import` takes a user's table in one of
Pandora's layer templates, checks every row against a step folder's
structures, and writes a layer folder. `pandora annotate compute`
replaces today's `pandora annotate`.

**Architecture:** The models go in `pandora/schemas/annotation.py`. A
pure resolver in the new `pandora/annotations/layers.py` turns user rows
into matched rows (both namings filled in) and unmatched rows (with a
`reason`). It does no file I/O. The CLI in `pandora/cli/annotate.py`
becomes a group with `compute` and `import`, and reuses the existing
table, YAML and step-summary helpers.

**Tech Stack:** Python 3, pydantic v2, argparse, pytest. No new
dependencies.

**Spec:** `design/specs/2026-11-06-annotation-layers-design.md`, section
"Part 1: foundation".

## Global Constraints

- PEP 8, `ruff format`, line length 80. `from __future__ import
  annotations` at the top of new modules.
- Type hints everywhere. Public functions get a Google-style docstring
  (`Args:`, `Returns:`, `Raises:`). Internal helpers start with `_` and
  get a one-line docstring.
- Public functions are added to `pandora/annotations/__init__.py`
  `__all__`.
- Never mutate inputs. Build new dicts and lists for output rows.
- Catch named exception types only (Ruff BLE001, S110).
- No new dependency. Parquet stays behind the `export` extra, via the
  existing `read_records` / `write_records`.
- Tests are offline and use `datasets/dev/mmcif/` fixtures.
- The old form `pandora annotate --input-dir ...` stops working, with
  no transition period.
- Before calling the work done: `uv run pytest`, `uv run ruff check .`,
  `uv run ruff format --check .`, `uv run zensical build --clean`.

## Review Focus

1. Parquet tables turn integer columns that contain nulls into floats
   (`seq_id` 42.0, `auth_seq_id` 45.0). These must resolve as 42 and
   "45" (Task 2 test).
2. Tools often write lowercase entry IDs (`1a2c`). They must match the
   structure `1A2C`, and the stored row carries `1A2C` (Task 2 test).
3. Insertion codes written as `"?"`, `"."`, `""` or null all mean
   "none" (Task 2 test).
4. A chain-level row in author naming whose author chain is shared by
   protein, ligands and waters (1a2c chain `L`) must resolve to the
   polymer chain, not be reported as ambiguous (Task 3 test).
5. An empty table (0 rows) must still give a `layer.json` with zero
   counts and no tables, not a crash (Task 5 test).

## Fixture facts used in tests (`datasets/dev/mmcif/a2/1a2c.cif`)

Entry id `1A2C`. Label chain `A` = author chain `L` (polymer, entity 1).
Label chain `B` = author chain `H` (polymer, entity 2). Label chain `F`
holds waters with author chain `L`.

| label chain, seq_id | author chain, auth_seq_id, ins_code | residue |
|---|---|---|
| A, 1 | L, 1, H | THR |
| A, 2 | L, 1, G | PHE |
| A, 3 | L, 1, F | GLY |
| F, None | L, 471, none | HOH |
| B, 148 (no atoms) | H, 147, A | in `_pdbx_poly_seq_scheme` only |

---

### Task 1: Layer models and column check

**Files:**
- Modify: `pandora/schemas/annotation.py`
- Create: `pandora/annotations/layers.py`
- Modify: `pandora/annotations/__init__.py`
- Test: `tests/test_annotations.py`

**Interfaces:**
- Produces (in `pandora/schemas/annotation.py`):
  ```python
  LayerTemplate = Literal["residue_contacts", "residue_values",
      "residue_labels", "chain_values", "chain_labels",
      "interface_labels", "free_form"]
  UnmatchedReason = Literal["entry_not_found", "missing_key",
      "chain_not_found", "residue_not_found", "no_coordinates",
      "naming_conflict", "ambiguous_author_id", "bad_value",
      "duplicate_key"]
  class LayerDescription(BaseModel):  # model_config extra="forbid"
      name: str
      template: LayerTemplate
      tool: str
      tool_version: str | None = None
      parameters: dict[str, Any] = {}   # Field(default_factory=dict)
      description: str = ""
  class LayerRow(BaseModel):          # model_config extra="allow", no fields
  ```
  `LayerRow` exists so that arbitrary user tables go through the
  existing `read_records(path, LayerRow)` / `write_records`.
- Produces (in `pandora/annotations/layers.py`):
  `check_layer_columns(columns: Collection[str], template: LayerTemplate) -> None`
  raises `ValueError`.

Template rules (used here and in Tasks 2–3):

| Template | Sides (column suffix) | Level | Required value |
|---|---|---|---|
| `residue_contacts` | `_1`, `_2` | residue | none |
| `residue_values` | `""` | residue | `value` |
| `residue_labels` | `""` | residue | `label` |
| `chain_values` | `""` | chain | `value` |
| `chain_labels` | `""` | chain | `label` |
| `interface_labels` | `_1`, `_2` | chain | `label` |
| `free_form` | inferred from the columns present | inferred | none |

Naming columns per side `s`: Pandora naming = `chain_id{s}` (+
`seq_id{s}` at residue level). Author naming = `auth_chain_id{s}` (+
`auth_seq_id{s}`, optional `ins_code{s}` at residue level). For
`free_form`, side `""` exists if any unsuffixed naming column is present,
and sides `_1`/`_2` if any suffixed one is. A side is residue level if
its `seq_id`/`auth_seq_id` column is present, else chain level.

`check_layer_columns` raises `ValueError` when: `entry_id` is absent;
any side has neither naming's columns complete; or the required value
column is absent. The message names the missing column(s) and the
template.

- [ ] **Step 1: Write the failing tests**

```python
def test_layer_description_rejects_unknown_template_and_fields():
    with pytest.raises(ValidationError):
        LayerDescription(name="x", template="nope", tool="t")
    with pytest.raises(ValidationError):
        LayerDescription(name="x", template="free_form", tool="t", extra=1)

def test_layer_row_keeps_any_columns():
    assert LayerRow(a=1, b="x").model_dump() == {"a": 1, "b": "x"}

def test_check_layer_columns_accepts_either_naming():
    check_layer_columns({"entry_id", "chain_id", "seq_id", "value"}, "residue_values")
    check_layer_columns({"entry_id", "auth_chain_id", "auth_seq_id", "value"}, "residue_values")
    check_layer_columns({"entry_id", "chain_id_1", "seq_id_1",
                         "auth_chain_id_2", "auth_seq_id_2"}, "residue_contacts")
    check_layer_columns({"entry_id"}, "free_form")

@pytest.mark.parametrize("columns, template, missing", [
    ({"chain_id", "seq_id", "value"}, "residue_values", "entry_id"),
    ({"entry_id", "chain_id", "seq_id"}, "residue_values", "value"),
    ({"entry_id", "chain_id", "label"}, "residue_labels", "seq_id"),
    ({"entry_id", "chain_id_1", "seq_id_1"}, "residue_contacts", "_2"),
    ({"entry_id", "chain_id_1", "chain_id_2"}, "interface_labels", "label"),
])
def test_check_layer_columns_rejects_missing(columns, template, missing):
    with pytest.raises(ValueError, match=missing):
        check_layer_columns(columns, template)
```

- [ ] **Step 2: Run the tests and confirm they fail.** Run
  `uv run pytest tests/test_annotations.py -k "layer" -v`. Expected:
  ImportError / FAIL.
- [ ] **Step 3: Implement** the models and `check_layer_columns`. Keep
  the template rules in one private table in `layers.py` (e.g.
  `_TEMPLATES: dict[str, tuple[tuple[str, ...], str, str | None]]`),
  plus a private `_sides(template, columns)` that Task 2 reuses. Add
  `check_layer_columns` to `pandora/annotations/__init__.py` `__all__`.
- [ ] **Step 4: Run the tests and confirm they pass.** Same command.
  Expected: PASS.
- [ ] **Step 5: Regenerate the schema diagrams.** Run
  `uv run --extra docs python docs/scripts/generate_erd.py`.
- [ ] **Step 6: Commit.** Message: `feat(annotations): layer template models and column check`

### Task 2: Resolve one residue side (both namings)

**Files:**
- Modify: `pandora/annotations/layers.py`
- Modify: `pandora/annotations/__init__.py`
- Test: `tests/test_annotations.py`

**Interfaces:**
- Consumes: `LayerTemplate`, `_sides` from Task 1.
- Produces:
  ```python
  def resolve_layer_rows(
      rows: Sequence[Mapping[str, Any]],
      template: LayerTemplate,
      structures: Mapping[str, Structure],
  ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
      """(matched, unmatched)."""
  ```
  A matched row is the user's row (all columns kept) with `entry_id`
  set to `Structure.entry_id` and, per side `s`, `chain_id{s}`,
  `auth_chain_id{s}` and, at residue level, `seq_id{s}` (int or None
  for non-polymers), `auth_seq_id{s}` (str) and `ins_code{s}` (str or
  None) filled in from the structure. An unmatched row is the user's row
  unchanged plus `reason: UnmatchedReason`. Exactly one reason per row:
  the first check that fails, in the order the rules below (and in Task
  3) are listed.

This task covers residue level, single side (`residue_values` /
`residue_labels`, without the value checks, which are Task 3). Matching
rules:
- **Entry:** case-insensitive match on `entry_id`, else
  `entry_not_found`.
- **Which naming:** Pandora naming is used when `chain_id` and `seq_id`
  are both non-null. Author naming is used when `auth_chain_id` and
  `auth_seq_id` are both non-null. Neither → `missing_key`. Both →
  resolve each, and if they reach different residues → `naming_conflict`.
- **Normalising:** `seq_id` accepts int, or a float/str holding an
  integer (42.0, "42"). Anything else → `bad_value`. `auth_seq_id` is
  compared as a string, with integral floats written without the
  decimal (45.0 → "45"). `ins_code` None, `""`, `"?"` and `"."` all mean
  none.
- **Residue index** per entry, built once per call from
  `structure.atoms`. Pandora key = (`label_asym_id`, `label_seq_id`),
  polymer residues only. Author key = (`auth_asym_id`, `auth_seq_id`,
  ins code or ""), all residues including ligands and waters. An author
  key that maps to more than one residue → `ambiguous_author_id`.
- **Chain missing** (no such label asym id / author chain id in
  `asym_units`) → `chain_not_found`.
- **Residue missing:** if the key is in `structure.raw`'s
  `_pdbx_poly_seq_scheme` (Pandora: `asym_id` + `seq_id`; author:
  `pdb_strand_id` + `pdb_seq_num` + `pdb_ins_code`) → `no_coordinates`,
  else `residue_not_found`. If the category is absent (canonicalisation
  can drop it), every missing residue is `residue_not_found`.

- [ ] **Step 1: Write the failing tests** (fixture: `1a2c`, parsed once
  in a module-level helper; `S = {"1A2C": structure}`):

```python
def test_resolve_author_naming_fills_pandora_keys():
    m, u = resolve_layer_rows([{"entry_id": "1a2c", "auth_chain_id": "L",
        "auth_seq_id": "1", "ins_code": "G", "value": 0.5}], "residue_values", S)
    assert u == []
    assert m == [{"entry_id": "1A2C", "chain_id": "A", "seq_id": 2,
        "auth_chain_id": "L", "auth_seq_id": "1", "ins_code": "G", "value": 0.5}]

def test_resolve_pandora_naming_fills_author_keys():
    m, _ = resolve_layer_rows([{"entry_id": "1A2C", "chain_id": "A",
        "seq_id": 1, "value": 1.0}], "residue_values", S)
    assert (m[0]["auth_chain_id"], m[0]["auth_seq_id"], m[0]["ins_code"]) == ("L", "1", "H")

def test_resolve_parquet_floats_and_ins_code_spellings():
    rows = [{"entry_id": "1A2C", "chain_id": "A", "seq_id": 3.0, "value": 1},
            {"entry_id": "1A2C", "auth_chain_id": "L", "auth_seq_id": 471.0,
             "ins_code": "?", "value": 1}]
    m, u = resolve_layer_rows(rows, "residue_values", S)
    assert u == [] and m[0]["seq_id"] == 3
    assert (m[1]["chain_id"], m[1]["seq_id"]) == ("F", None)

@pytest.mark.parametrize("row, reason", [
    ({"entry_id": "9XYZ", "chain_id": "A", "seq_id": 1}, "entry_not_found"),
    ({"entry_id": "1A2C", "chain_id": "A"}, "missing_key"),
    ({"entry_id": "1A2C", "chain_id": "Z", "seq_id": 1}, "chain_not_found"),
    ({"entry_id": "1A2C", "chain_id": "A", "seq_id": 9999}, "residue_not_found"),
    ({"entry_id": "1A2C", "auth_chain_id": "H", "auth_seq_id": "147",
      "ins_code": "A"}, "no_coordinates"),
    ({"entry_id": "1A2C", "chain_id": "B", "seq_id": 148}, "no_coordinates"),
    ({"entry_id": "1A2C", "chain_id": "A", "seq_id": 1, "auth_chain_id": "L",
      "auth_seq_id": "1", "ins_code": "G"}, "naming_conflict"),
    ({"entry_id": "1A2C", "chain_id": "A", "seq_id": "x"}, "bad_value"),
])
def test_resolve_unmatched_reasons(row, reason):
    m, u = resolve_layer_rows([{**row, "label": "x"}], "residue_labels", S)
    assert m == [] and u == [{**row, "label": "x", "reason": reason}]

def test_resolve_ambiguous_author_id():
    # strip insertion codes, as preserve_insertion_codes=False does (#125)
    atoms = [a.model_copy(update={"pdbx_PDB_ins_code": None}) for a in S["1A2C"].atoms]
    s = S["1A2C"].model_copy(update={"atoms": atoms})
    _, u = resolve_layer_rows([{"entry_id": "1A2C", "auth_chain_id": "L",
        "auth_seq_id": "1", "label": "x"}], "residue_labels", {"1A2C": s})
    assert u[0]["reason"] == "ambiguous_author_id"

def test_resolve_does_not_mutate_rows():
    row = {"entry_id": "1a2c", "chain_id": "A", "seq_id": 1, "label": "x"}
    resolve_layer_rows([row], "residue_labels", S)
    assert row == {"entry_id": "1a2c", "chain_id": "A", "seq_id": 1, "label": "x"}
```

- [ ] **Step 2: Run the tests and confirm they fail.** Run
  `uv run pytest tests/test_annotations.py -k resolve -v`. Expected:
  FAIL.
- [ ] **Step 3: Implement `resolve_layer_rows`** with the rules above,
  for residue-level single-side templates. Structure the per-side
  resolution as a private function returning either the filled keys or
  a reason, so Task 3 can call it once per side. Add it to `__all__`.
- [ ] **Step 4: Run the tests and confirm they pass.** Same command.
  Expected: PASS.
- [ ] **Step 5: Commit.** Message: `feat(annotations): resolve layer rows against structures`

### Task 3: Pairs, chain level, values, duplicates, free-form

**Files:**
- Modify: `pandora/annotations/layers.py`
- Test: `tests/test_annotations.py`

**Interfaces:**
- Consumes: `resolve_layer_rows`, the per-side helper (Task 2), and
  `polymer_asym_ids` from `pandora/annotations/entry.py`.
- Produces: `resolve_layer_rows` handles all seven templates, with the
  same signature.

Rules:
- **Pair templates** resolve sides `_1` and `_2` independently. The
  row's reason is the first failure (side `_1` before `_2`).
- **Chain level:** Pandora naming = `chain_id{s}` must be an asym id in
  `asym_units`. Author naming = `auth_chain_id{s}` resolves among the
  **polymer** asym units (`polymer_asym_ids(structure)`) whose
  `auth_id` matches. None → `chain_not_found`, more than one →
  `ambiguous_author_id`.
- **Values:** `value` must be an int or a float, not bool and not NaN.
  `label` must be a non-empty str. Otherwise → `bad_value`.
- **Duplicates:** for every template except `residue_contacts` and
  `free_form`, a row whose resolved key (entry + every side's filled
  key columns) was already matched → `duplicate_key`. The first row is
  kept.
- **Free-form:** sides and levels are inferred per Task 1, but from
  **each row's own non-null columns**, not from the whole table. A row
  with only `entry_id` matches when the entry exists.

- [ ] **Step 1: Write the failing tests**

```python
def test_resolve_contacts_both_sides_and_repeats_kept():
    row = {"entry_id": "1A2C", "chain_id_1": "A", "seq_id_1": 1,
           "auth_chain_id_2": "L", "auth_seq_id_2": "1", "ins_code_2": "G",
           "contact_type": "hbond"}
    m, u = resolve_layer_rows([row, {**row, "contact_type": "vdw"}], "residue_contacts", S)
    assert u == [] and len(m) == 2
    assert (m[0]["auth_chain_id_1"], m[0]["chain_id_2"], m[0]["seq_id_2"]) == ("L", "A", 2)

def test_resolve_contacts_reports_first_failing_side():
    _, u = resolve_layer_rows([{"entry_id": "1A2C", "chain_id_1": "A",
        "seq_id_1": 1, "chain_id_2": "Z", "seq_id_2": 1}], "residue_contacts", S)
    assert u[0]["reason"] == "chain_not_found"

def test_resolve_chain_author_naming_picks_polymer():
    m, u = resolve_layer_rows([{"entry_id": "1A2C", "auth_chain_id": "L",
        "label": "light"}], "chain_labels", S)
    assert u == [] and m[0]["chain_id"] == "A"

def test_resolve_interface_labels():
    m, _ = resolve_layer_rows([{"entry_id": "1A2C", "chain_id_1": "A",
        "chain_id_2": "B", "label": "physiological"}], "interface_labels", S)
    assert (m[0]["auth_chain_id_1"], m[0]["auth_chain_id_2"]) == ("L", "H")

@pytest.mark.parametrize("value", [None, float("nan"), True, "high"])
def test_resolve_bad_value(value):
    _, u = resolve_layer_rows([{"entry_id": "1A2C", "chain_id": "A",
        "value": value}], "chain_values", S)
    assert u[0]["reason"] == "bad_value"

def test_resolve_duplicate_key_keeps_first():
    rows = [{"entry_id": "1A2C", "chain_id": "A", "seq_id": 1, "label": "a"},
            {"entry_id": "1A2C", "auth_chain_id": "L", "auth_seq_id": "1",
             "ins_code": "H", "label": "b"}]
    m, u = resolve_layer_rows(rows, "residue_labels", S)
    assert [r["label"] for r in m] == ["a"]
    assert u[0]["reason"] == "duplicate_key"

def test_resolve_free_form():
    rows = [{"entry_id": "1a2c", "note": "x"}, {"entry_id": "1a2c", "note": "y"},
            {"entry_id": "1A2C", "chain_id": "Z", "note": "z"}]
    m, u = resolve_layer_rows(rows, "free_form", S)
    assert [r["note"] for r in m] == ["x", "y"]
    assert u[0]["reason"] == "chain_not_found"
```

- [ ] **Step 2: Run the tests and confirm they fail.** Run
  `uv run pytest tests/test_annotations.py -k resolve -v`. Expected: the
  new tests FAIL, and Task 2's still pass.
- [ ] **Step 3: Implement** the rules above in `layers.py`.
- [ ] **Step 4: Run the tests and confirm they pass.** Run
  `uv run pytest tests/test_annotations.py -v`. Expected: PASS.
- [ ] **Step 5: Commit.** Message: `feat(annotations): pair, chain, value and duplicate rules for layers`

### Task 4: `pandora annotate` becomes a group; `compute` replaces the flat command

**Files:**
- Modify: `pandora/cli/annotate.py`
- Modify: `docs/usage/cli.md`, `docs/usage/annotation.md`, and any other
  hit of `grep -rn "pandora annotate" docs examples pandora README.md`
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: `register()` adds an `annotate` parser with required
  subcommands. `compute` has exactly today's arguments (`--input-dir`,
  `--layers`, `--pairwise`, `--output-dir`) and handler `_cmd_compute`
  (today's `_cmd_annotate`, renamed). Task 5 adds `import` to the same
  subparsers object.

- [ ] **Step 1: Write the failing tests**

```python
def test_annotate_compute_writes_annotations(tmp_path):
    src = tmp_path / "in"; src.mkdir()
    shutil.copy(FIXTURES_DIR / "04" / "104m.cif", src / "104m.cif")
    main(["annotate", "compute", "--input-dir", str(src),
          "--layers", "structure_counts", "--output-dir", str(tmp_path / "out")])
    data = json.loads((tmp_path / "out" / "annotations.json").read_text())
    assert list(data) == ["104M"]

def test_annotate_old_flat_form_is_gone(tmp_path):
    with pytest.raises(SystemExit):
        main(["annotate", "--input-dir", str(tmp_path), "--output-dir", str(tmp_path)])
```

- [ ] **Step 2: Run the tests and confirm they fail.** Run
  `uv run pytest tests/test_cli.py -k annotate -v`. Expected: the first
  test FAILs.
- [ ] **Step 3: Implement.** Use `add_subparsers(dest=...,
  required=True)` under `annotate`. Behaviour and output of `compute`
  are unchanged.
- [ ] **Step 4: Update the docs.** Every `pandora annotate --input-dir`
  becomes `pandora annotate compute --input-dir` (4 places in
  `docs/usage/annotation.md`, 1 in `docs/usage/cli.md`). In `cli.md`,
  rename the section to `annotate compute`.
- [ ] **Step 5: Run the tests and the docs build.** Run
  `uv run pytest tests/test_cli.py -v` and `uv run zensical build --clean`.
  Expected: PASS, and a build with no new warnings.
- [ ] **Step 6: Commit.** Message: `feat(cli)!: pandora annotate compute replaces flat annotate`

### Task 5: `pandora annotate import`

**Files:**
- Modify: `pandora/cli/_io.py`
- Modify: `pandora/cli/annotate.py`
- Modify: `pandora/cli/README.md`, `docs/usage/cli.md`,
  `docs/usage/annotation.md`
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `LayerDescription`, `LayerRow` (Task 1),
  `check_layer_columns`, `resolve_layer_rows` (Tasks 1–3), the `annotate`
  subparsers (Task 4). Also the existing `load_policy(path, model=)`,
  `read_records`, `write_table`, `check_table_format`,
  `read_structures_dir`, `build_step_summary`, `write_json` and
  `now_iso`.
- Produces, in `_io.py`:
  `prepare_output_dir(out: Path, outputs: Sequence[str], overwrite: bool) -> None`.
  `outputs` are names relative to `out`, files or folders. It refuses if
  any exists and `overwrite` is false (`SystemExit`, same message
  shape as today), otherwise deletes those that exist, then creates
  `out`. `prepare_step_dir` is rewritten to call it with its current
  list, with the same behaviour.

Command: `pandora annotate import --structures-dir DIR --table FILE
--description YAML [--source-files PATH] --output-dir DIR [--overwrite]
[--table-format {parquet,jsonl}]` (default `parquet`).

Behaviour, in order:
1. `started_at = now_iso()`, then `check_table_format`.
2. Load the description with `load_policy(path, model=LayerDescription)`.
   Read rows with `read_records(table, LayerRow)` and `model_dump()`
   each one. Run `check_layer_columns` over the union of row keys
   (skipped when there are 0 rows). A `ValueError` from any of these →
   `SystemExit` with its message, **before anything is written**.
3. `prepare_output_dir(out, ["layer.parquet", "layer.jsonl",
   "unmatched.parquet", "unmatched.jsonl", "source", "layer.json"],
   overwrite)`.
4. `read_structures_dir(structures_dir)`, then `resolve_layer_rows`.
5. `--source-files`: copy a file to `out/source/<name>`, or a folder's
   contents to `out/source/` (`shutil.copytree`).
6. Write `layer.<fmt>` and `unmatched.<fmt>` with `write_table` (rows as
   `LayerRow(**row)`). Empty tables are not written.
7. Write `layer.json` with `build_step_summary("annotation_import",
   description, input=table, output=out, started_at=..., counts=...,
   previous_summary=...)`. `counts` = `rows_in`, `rows_matched`,
   `rows_unmatched`, plus `unmatched:<reason>` per reason present,
   sorted. `previous_summary` = the first of `curation.json`,
   `canonicalisation.json` that exists in `--structures-dir`, else None.
8. Print `imported: {matched} rows, {unmatched} unmatched -> {out}`.

- [ ] **Step 1: Write the failing tests.** Helper `_layer_inputs(tmp_path)`:
  copy `1a2c.cif` to `tmp/curated/structures/a2/1a2c.cif`, write
  `tmp/curated/curation.json` as `{}`, write `tmp/cocada.yaml`
  (`name: t`, `template: residue_values`, `tool: T`,
  `tool_version: "1"`, `parameters: {cutoff: 5.0}`), write
  `tmp/raw/out.txt`, and write `tmp/t.jsonl` with these four rows:
  `{"entry_id":"1a2c","auth_chain_id":"L","auth_seq_id":"1","ins_code":"H","value":1.0}`,
  `{"entry_id":"1A2C","chain_id":"A","seq_id":2,"value":2.0}`,
  `{"entry_id":"1A2C","auth_chain_id":"H","auth_seq_id":"147","ins_code":"A","value":3.0}`,
  `{"entry_id":"9XYZ","chain_id":"A","seq_id":1,"value":4.0}`.

```python
def test_annotate_import_writes_layer_folder(tmp_path):
    args = _layer_inputs(tmp_path)  # list of CLI args incl. --source-files, jsonl
    main(["annotate", "import", *args])
    out = tmp_path / "annotations" / "t"
    matched = read_records(out / "layer.jsonl", LayerRow)
    assert [r.model_dump()["seq_id"] for r in matched] == [1, 2]
    unmatched = [r.model_dump()["reason"] for r in read_records(out / "unmatched.jsonl", LayerRow)]
    assert unmatched == ["no_coordinates", "entry_not_found"]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "annotation_import"
    assert summary["policy"]["template"] == "residue_values"
    assert summary["counts"] == {"rows_in": 4, "rows_matched": 2, "rows_unmatched": 2,
        "unmatched:entry_not_found": 1, "unmatched:no_coordinates": 1}
    assert summary["previous_summary"].endswith("curation.json")
    assert (out / "source" / "out.txt").exists()

def test_annotate_import_refuses_existing_layer_without_overwrite(tmp_path):
    args = _layer_inputs(tmp_path)
    main(["annotate", "import", *args])
    with pytest.raises(SystemExit):
        main(["annotate", "import", *args])
    main(["annotate", "import", *args, "--overwrite"])

def test_annotate_import_missing_value_column_writes_nothing(tmp_path):
    args = _layer_inputs(tmp_path)
    (tmp_path / "t.jsonl").write_text('{"entry_id":"1A2C","chain_id":"A","seq_id":1}\n')
    with pytest.raises(SystemExit, match="value"):
        main(["annotate", "import", *args])
    assert not (tmp_path / "annotations" / "t").exists()

def test_annotate_import_empty_table(tmp_path):
    args = _layer_inputs(tmp_path)
    (tmp_path / "t.jsonl").write_text("")
    main(["annotate", "import", *args])
    out = tmp_path / "annotations" / "t"
    assert json.loads((out / "layer.json").read_text())["counts"]["rows_in"] == 0
    assert not (out / "layer.jsonl").exists()
```

  Also add a test that `prepare_step_dir` still refuses and overwrites as
  before, if `tests/test_cli.py` does not already cover it (check with
  `grep -n "overwrite" tests/test_cli.py`).

- [ ] **Step 2: Run the tests and confirm they fail.** Run
  `uv run pytest tests/test_cli.py -k "annotate_import" -v`. Expected:
  FAIL.
- [ ] **Step 3: Implement** `prepare_output_dir`, refactor
  `prepare_step_dir` onto it, and add `_cmd_import` and its parser as
  described above.
- [ ] **Step 4: Run the full suite.** Run `uv run pytest`. Expected: PASS.
- [ ] **Step 5: Document it.**
  - `pandora/cli/README.md`: one paragraph on `annotate compute` /
    `annotate import`, covering the layer folder files and that
    unmatched rows don't stop the import.
  - `docs/usage/cli.md`: an `annotate import` section with the command
    and its printed line.
  - `docs/usage/annotation.md`: a section "Importing a layer from
    another tool" with the template table and the naming columns (from
    the spec), the description YAML, the command, the folder layout, and
    the unmatched reasons. Follow the usage-page convention: show the
    real command and paste the real output from running it on the
    fixture.
- [ ] **Step 6: Check.** Run `uv run ruff check .`,
  `uv run ruff format --check .` and `uv run zensical build --clean`.
  Expected: all clean.
- [ ] **Step 7: Commit.** Message: `feat(cli): pandora annotate import`
