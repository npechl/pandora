# Annotation Layers Part 4 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every layer joinable onto records by the user: residue
records gain author keys, contact layers can be turned into count
layers, and the docs show the joins.

**Architecture:** No join code. `ResidueRecord` gains `auth_chain_id`
and `ins_code`. Two plain functions in a new
`pandora/annotations/counts.py` turn `ResidueContactRow`s into
`residue_values` / `chain_values` rows. A thin `pandora annotate
count-contacts` command writes them as an ordinary layer folder. The
docs show the pandas joins.

**Tech Stack:** Python 3, pydantic v2, argparse, pytest, pandas (docs
examples only, `export` extra).

**Spec:** `design/specs/2026-11-06-annotation-layers-design.md`,
section "Part 4: layers are joined by the user".

## Global Constraints

- `from __future__ import annotations` at the top of new modules;
  type hints everywhere; line length 80; `uv run ruff format`.
- Public functions: Google-style docstring (`Args:`, `Returns:`,
  `Raises:`), exported in the package `__init__.py` `__all__`.
- Never mutate inputs.
- No new dependency. Tests offline, fixtures from
  `datasets/dev/mmcif/` only.
- No join command, no join function, no entry record, no
  `entry_values` template.
- Insertion code `None` for none (`""`, `"?"`, `"."` all mean none, as
  `_NO_INS_CODE` in `pandora/schemas/annotation.py`).
- Zero-contact residues and chains get **no row**.

## Review Focus

1. A contact layer folder with no rows (empty tables are not written,
   so only `layer.json` exists) → `count-contacts` writes no table,
   `counts == {"rows_in": 0, "rows_out": 0}`, no crash. Test in Task 3.
2. A layer read back from Parquet gives `seq_id` as float → read
   through `ResidueContactRow`, so counts key on ints. Covered by
   reading with `read_records(path, ResidueContactRow)` in Task 3.
3. A ligand side (`seq_id` None) must not merge with another ligand in
   the same chain → residue key includes `auth_seq_id` and
   `ins_code`. Test in Task 2.
4. A row whose two sides are the same residue → it counts as no
   partner and creates no row by itself. Test in Task 2.
5. Residues sharing `auth_seq_id` but with different insertion codes
   (1a2c chain L, 14A/14B/14C) → distinct records with distinct
   `ins_code`. Test in Task 1.

---

### Task 1: Residue records carry author keys

**Files:**
- Modify: `pandora/schemas/dataset.py` (`ResidueRecord`, ~line 64)
- Modify: `pandora/datasets/records.py` (`extract_residue_records`,
  ~lines 80-116)
- Test: `tests/test_curation.py` (holds the other record tests)

**Interfaces:**
- Produces: `ResidueRecord.auth_chain_id: str | None = None`,
  `ResidueRecord.ins_code: str | None = None`, both documented in the
  `Attributes:` block.

- [ ] **Step 1: Write the failing tests**

```python
def test_residue_records_carry_author_keys():
    structure = _load("1a2c")  # use the file's existing fixture loader
    records = extract_residue_records(structure)
    by_seq = {(r.chain_id, r.seq_id): r for r in records}
    assert by_seq[("A", 1)].auth_chain_id == "L"
    assert by_seq[("A", 1)].auth_seq_id == "1"
    assert by_seq[("A", 1)].ins_code == "H"
    codes = [by_seq[("A", s)].ins_code for s in (23, 24, 25)]
    assert codes == ["A", "B", "C"]
    assert all(by_seq[("A", s)].auth_seq_id == "14" for s in (23, 24, 25))
    assert any(r.ins_code is None for r in records)


def test_residue_record_without_author_keys_still_loads():
    record = ResidueRecord(
        entry_id="1abc", chain_id="A", seq_id=1, comp_id="ALA", atoms=[]
    )
    assert record.auth_chain_id is None and record.ins_code is None
```

- [ ] **Step 2: Run them, expect FAIL** (`AttributeError` /
  `ins_code` missing)

Run: `uv run pytest tests/test_curation.py -k "author_keys" -v`

- [ ] **Step 3: Implement**

Add the two fields after `auth_seq_id`. In `extract_residue_records`
keep the grouping key as is and fill both fields from the first atom
of each group: `auth_chain_id = atom.auth_asym_id`, `ins_code =
None if atom.pdbx_PDB_ins_code in _NO_INS_CODE else
atom.pdbx_PDB_ins_code`. Import `_NO_INS_CODE` from
`pandora.schemas.annotation` rather than redefining it.

- [ ] **Step 4: Run tests, expect PASS**, then `uv run pytest`
  (nothing else may break).

- [ ] **Step 5: Commit**

```bash
git add pandora/schemas/dataset.py pandora/datasets/records.py tests/test_curation.py
git commit -m "feat(datasets): residue records carry auth_chain_id and ins_code"
```

---

### Task 2: Count functions

**Files:**
- Create: `pandora/annotations/counts.py`
- Modify: `pandora/annotations/__init__.py` (import + `__all__`)
- Test: `tests/test_annotations.py`

**Interfaces:**
- Consumes: `ResidueContactRow` (`pandora/schemas/annotation.py`),
  `LayerRow`.
- Produces:
  - `count_residue_contacts(rows: Iterable[ResidueContactRow]) ->
    list[LayerRow]` — columns, in order: `entry_id`, `chain_id`,
    `seq_id`, `auth_chain_id`, `auth_seq_id`, `ins_code`, `comp_id`,
    `value` (int).
  - `count_chain_contacts(rows: Iterable[ResidueContactRow]) ->
    list[LayerRow]` — columns: `entry_id`, `chain_id`,
    `auth_chain_id`, `value` (int).
  - Both sorted by their key columns (`None` sorts first).

Definitions (from the spec):
- residue key = `(entry_id, chain_id, seq_id, auth_seq_id, ins_code)`
  of one side;
- residue `value` = number of **distinct other** residue keys it pairs
  with, over both sides; (A,B) and (B,A) and repeats count once; a
  side paired with itself adds nothing;
- chain `value` = number of distinct residue keys of that
  `(entry_id, chain_id)` that have at least one partner.

- [ ] **Step 1: Write the failing tests**

Build rows with a small helper `_contact(c1, s1, c2, s2, **kw)` that
fills both namings (`auth_chain_id_s = chain`, `auth_seq_id_s =
str(seq)`, `comp_id_s = "ALA"`, `entry_id = "1abc"`).

```python
def test_count_residue_contacts_counts_distinct_partners():
    rows = [
        _contact("A", 1, "B", 5, contact_type="hbond"),
        _contact("A", 1, "B", 5, contact_type="salt_bridge"),
        _contact("B", 5, "A", 1),
        _contact("A", 1, "B", 6),
    ]
    out = {(r.chain_id, r.seq_id): r.value for r in count_residue_contacts(rows)}
    assert out == {("A", 1): 2, ("B", 5): 1, ("B", 6): 1}


def test_count_residue_contacts_keeps_both_namings():
    row = count_residue_contacts([_contact("A", 1, "B", 5)])[0]
    assert row.model_dump() == {
        "entry_id": "1abc", "chain_id": "A", "seq_id": 1,
        "auth_chain_id": "A", "auth_seq_id": "1", "ins_code": None,
        "comp_id": "ALA", "value": 1,
    }


def test_count_residue_contacts_ligands_keyed_by_author_naming():
    rows = [
        _contact("A", 1, "C", None, auth_seq_id_2="301"),
        _contact("A", 2, "C", None, auth_seq_id_2="302"),
    ]
    ligands = [r for r in count_residue_contacts(rows) if r.chain_id == "C"]
    assert [(r.seq_id, r.auth_seq_id, r.value) for r in ligands] == [
        (None, "301", 1), (None, "302", 1)
    ]


def test_count_residue_contacts_self_pair_adds_nothing():
    assert count_residue_contacts([_contact("A", 1, "A", 1)]) == []


def test_count_chain_contacts_counts_distinct_residues():
    rows = [_contact("A", 1, "B", 5), _contact("A", 1, "B", 6),
            _contact("A", 2, "B", 5)]
    out = {r.chain_id: r.value for r in count_chain_contacts(rows)}
    assert out == {"A": 2, "B": 2}


def test_count_layers_pass_their_template_check():
    rows = [_contact("A", 1, "B", 5)]
    residue = count_residue_contacts(rows)[0].model_dump()
    chain = count_chain_contacts(rows)[0].model_dump()
    check_layer_columns(residue, "residue_values")  # raises on failure
    check_layer_columns(chain, "chain_values")


def test_count_functions_do_not_mutate_input():
    rows = [_contact("A", 1, "B", 5)]
    before = [r.model_copy(deep=True) for r in rows]
    count_residue_contacts(rows)
    count_chain_contacts(rows)
    assert rows == before
```

- [ ] **Step 2: Run, expect FAIL** (import error)

Run: `uv run pytest tests/test_annotations.py -k count -v`

- [ ] **Step 3: Implement both functions in
  `pandora/annotations/counts.py`** and export them. One pass building
  `partners: dict[residue_key, set[residue_key]]` plus a
  `names: dict[residue_key, (auth_chain_id, comp_id)]`; the chain
  count is derived from `partners` keys. No structures, no file I/O.

- [ ] **Step 4: Run, expect PASS**; `uv run pytest tests/test_annotations.py`.

- [ ] **Step 5: Commit**

```bash
git add pandora/annotations/counts.py pandora/annotations/__init__.py tests/test_annotations.py
git commit -m "feat(annotations): contact counts as residue_values / chain_values rows"
```

---

### Task 3: `pandora annotate count-contacts`

**Files:**
- Modify: `pandora/cli/_io.py` (add `check_contact_layer`)
- Modify: `pandora/cli/similarity.py` (`_contact_layer` calls it)
- Modify: `pandora/cli/annotate.py` (new subcommand)
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `count_residue_contacts`, `count_chain_contacts` (Task 2);
  `find_table`, `prepare_output_dir`, `write_table`,
  `check_table_format` (`_io.py`); `read_records`, `write_json`;
  `build_step_summary`; `LayerDescription`.
- Produces: `check_contact_layer(command: str, folder: Path) -> None`
  in `pandora/cli/_io.py`: the "not a layer folder" / "is a {template}
  layer, --contacts needs residue_contacts" `SystemExit`s, moved
  verbatim out of `similarity._contact_layer`, which then keeps only
  `check_contact_layer(...)` + its `find_table(...) or summary` return
  (behaviour unchanged).

Command:
`pandora annotate count-contacts --contacts DIR --level
{residue,chain} --output-dir DIR [--overwrite] [--table-format]`.

Writes `layer.<fmt>` (skipped when there are no rows, as other steps)
and `layer.json`:
- `build_step_summary("annotation_count", description, input=contacts,
  output=out, started_at=..., counts={"rows_in": n_in, "rows_out":
  n_out}, previous_summary=contacts / "layer.json")`
- `description = LayerDescription(name="contact_count",
  template="residue_values" | "chain_values", tool="pandora",
  tool_version=__version__, parameters={"level": level,
  "contacts_layer": <contacts layer.json policy name>})`
- prints `counted: {n_out} {level} rows from {n_in} contacts -> {out}`.

- [ ] **Step 1: Write the failing tests** (reuse `_step_folder` and
  `_compute` in `tests/test_cli.py`)

```python
def _count(contacts, out, level):
    main(["annotate", "count-contacts", "--contacts", str(contacts),
          "--level", level, "--output-dir", str(out),
          "--table-format", "jsonl"])


@pytest.mark.parametrize(
    "level, template, function",
    [("residue", "residue_values", count_residue_contacts),
     ("chain", "chain_values", count_chain_contacts)],
)
def test_annotate_count_contacts_writes_value_layer(tmp_path, level, template, function):
    contacts = tmp_path / "interfaces"
    _compute(_step_folder(tmp_path), contacts, "--layer", "chain_interfaces")
    out = tmp_path / f"counts_{level}"
    _count(contacts, out, level)

    rows = read_records(contacts / "layer.jsonl", ResidueContactRow)
    written = read_records(out / "layer.jsonl", LayerRow)
    assert [r.model_dump() for r in written] == [
        r.model_dump() for r in function(rows)
    ]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "annotation_count"
    assert summary["policy"]["template"] == template
    assert summary["policy"]["name"] == "contact_count"
    assert summary["policy"]["parameters"] == {
        "level": level, "contacts_layer": "chain_interfaces"}
    assert summary["counts"] == {"rows_in": len(rows), "rows_out": len(written)}
    assert summary["previous_summary"] == str(contacts / "layer.json")


def test_annotate_count_contacts_rejects_non_contact_layer(tmp_path):
    # a folder whose layer.json policy.template is "residue_values"
    # (e.g. the output of the test above) -> SystemExit matching
    # "needs residue_contacts"; the output dir is not created.


def test_annotate_count_contacts_empty_layer(tmp_path):
    # a contacts folder holding only layer.json (template
    # residue_contacts, no layer table) -> no layer table written,
    # counts == {"rows_in": 0, "rows_out": 0}.
```

- [ ] **Step 2: Run, expect FAIL** (`invalid choice: 'count-contacts'`)

Run: `uv run pytest tests/test_cli.py -k count_contacts -v`

- [ ] **Step 3: Implement** `check_contact_layer` in `_io.py`, switch
  `similarity._contact_layer` to it, then `_cmd_count` +
  parser in `annotate.py`. Read rows with
  `read_records(table, ResidueContactRow)` when `find_table(contacts,
  "layer")` returns a path, else `[]`. Check the layer before
  `prepare_output_dir`.

- [ ] **Step 4: Run, expect PASS**; then `uv run pytest` (the
  similarity `--contacts` tests must still pass).

- [ ] **Step 5: Commit**

```bash
git add pandora/cli/_io.py pandora/cli/similarity.py pandora/cli/annotate.py tests/test_cli.py
git commit -m "feat(cli): annotate count-contacts writes a contact count layer"
```

---

### Task 4: Docs and diagrams

**Files:**
- Modify: `docs/usage/annotation.md` (new `## Using layers in a
  dataset` section, after "Importing a layer from another tool")
- Modify: `docs/usage/cli.md`, `pandora/cli/README.md`
  (`annotate count-contacts`)
- Modify: `docs/usage/datasets.md` (`ResidueRecord` new fields, where
  residue records are described)
- Modify: `docs/contributing/architecture.md` (the two count
  functions: `ResidueContactRow` in, `LayerRow` out)
- Regenerate: ERD

Usage examples follow the page's convention: full code inline,
`=== "library"` / `=== "cli"` tabs where both exist, each followed by
a `text` block with the **real** output (run them by hand on the
1a2c fixture and paste). The section shows:
1. `count-contacts` at residue level, then residue values joined onto
   residue records with `pandas.merge` on `entry_id, chain_id,
   seq_id`, and the same join on `entry_id, auth_chain_id,
   auth_seq_id, ins_code`;
2. the chain count layer joined onto `chains.<fmt>` with
   `fillna(0)`;
3. interface labels onto interface records matching either chain
   order (a small imported `interface_labels` table; normalise both
   sides to a sorted chain pair before the merge);
4. an entry total: `groupby("entry_id")["value"].sum()` on the chain
   count layer;
5. one sentence: per-entry values from another tool import as a
   `free_form` layer and join on `entry_id`.

No library tab is needed for the joins themselves (they are pandas);
say so instead of an empty tab.

- [ ] **Step 1: Write the section and the reference edits.**
- [ ] **Step 2: Run every new example by hand; paste real output.**
- [ ] **Step 3: Regenerate the ERD:**
  `uv run --extra docs python docs/scripts/generate_erd.py`
- [ ] **Step 4: Build the docs** with the `docs-build` skill
  (`uv run zensical build --clean`); it must build with no new
  warnings, and the rendered HTML must contain the example code, not
  only the outputs.
- [ ] **Step 5: Final checks:** `uv run pytest`,
  `uv run ruff format --check .`, `uv run ruff check .`
- [ ] **Step 6: Commit**

```bash
git add docs pandora/cli/README.md
git commit -m "docs: joining layers onto records, annotate count-contacts"
```
