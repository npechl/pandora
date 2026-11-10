# Label Interfaces Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One command, `pandora annotate label-interfaces`, that writes an
`interface_labels` layer straight from a contact layer: per-pair contact
counts, dimer size, and a `too_large` / `few_contacts` / `ok` label
(issue #133).

**Architecture:** One plain function in `pandora/annotations/counts.py`,
`label_interfaces(contacts, lengths, ...)`, returning one labelled row
per interface. Chain lengths are not computed: `pandora curate` already
writes them as `seqres_length` in the step folder's `chains` table. A
thin CLI handler finds the step folder the contact layer traces to with
the existing `_layer_origin()`, reads its `chains` table, and writes an
ordinary layer folder, as `count-contacts` does. The output feeds
`filter-contacts --labels` and `partition --labels` unchanged.

**Tech Stack:** Python 3, pydantic v2, argparse, pytest.

**Spec:** issue #133, plus the design agreed in chat on 2026-11-10,
restated in full under Global Constraints.

## Global Constraints

- `from __future__ import annotations`; type hints; line length 80;
  `uv run ruff format`. Public functions: Google-style docstring,
  exported in `pandora/annotations/__init__.py` `__all__`. Never mutate
  inputs. No new dependency. Tests offline, fixtures from
  `datasets/dev/mmcif/` only.
- An interface is `(entry_id, unordered {chain_id_1, chain_id_2})` on
  Pandora chain ids; entry ids compare case-insensitively (as in
  `filters._interface`). Output rows put the smaller chain id in
  `chain_id_1`. No special case for (A, A) or ligand rows.
- `residues_1` = distinct residues of `chain_id_1` in a contact row with
  `chain_id_2`; `residues_2` the reverse. A residue is identified as in
  `counts._side` (`chain_id, seq_id, auth_seq_id, ins_code`); a pair
  given several times or in either order counts once.
- Chain length = `seqres_length` from the step folder's `chains` table
  (`DatasetChainRecord`, written by `pandora curate`): the full intended
  sequence, not modelled residues. A chain absent from the table, or
  with `seqres_length` null, has no length.
- `total_residues = length_1 + length_2`; `None` if either is unknown.
- Labels, one per interface, first match wins:
  1. `too_large` if `max_total_residues` is set and
     `total_residues > max_total_residues` (unknown total never matches);
  2. `few_contacts` if `min_residues_per_chain` is set and
     `min(residues_1, residues_2) < min_residues_per_chain`;
  3. `ok` otherwise.
- Both thresholds optional; omitted means that rule is off.
- Row columns, in this order: `entry_id, chain_id_1, chain_id_2,
  auth_chain_id_1, auth_chain_id_2, residues_1, residues_2, length_1,
  length_2, total_residues, label`. Rows sorted by
  `(entry_id, chain_id_1, chain_id_2)`.
- CLI flags: `--contacts` (dest `contacts_dir`, required),
  `--min-residues-per-chain INT`, `--max-total-residues INT`,
  `--output-dir` (required), `--overwrite`, `--table-format` (default
  parquet). No `--structures-dir`: the folder comes from `_layer_origin`.
- `layer.json`: step `annotation_label`; description `name="contact_labels"`,
  `template="interface_labels"`, `tool="pandora"`; parameters
  `{"contacts_layer", "structures_dir", "min_residues_per_chain",
  "max_total_residues"}` (thresholds `null` when off);
  `previous_summary = <contacts>/layer.json`; counts `rows_in`,
  `interfaces`, `ok`, `few_contacts`, `too_large`, `unknown_length`.

## Review Focus

1. Contact layer whose origin can't be traced (recorded `input` path not
   readable from here) → exits before writing, naming the problem.
   Test in Task 2.
2. A contact row's chain is missing from the `chains` table, or its
   `seqres_length` is null → `length_*`/`total_residues` are `None`,
   never `too_large`, counted in `unknown_length`. Tests in Tasks 1
   and 2. A step folder with no `chains` table at all → exits before
   writing, naming `pandora curate`. Test in Task 2.
3. Interface present only as B–A rows → one row with `chain_id_1 = A`,
   counts swapped to match. Test in Task 1.
4. Both rules match → `too_large`. Test in Task 1.
5. The output folder used as `filter-contacts --labels` → accepted
   (template and origin check pass) and drops `few_contacts`. Test in
   Task 2.

---

### Task 1: Library function

**Files:**
- Modify: `pandora/annotations/counts.py` (add `label_interfaces`; reuse
  `_side`)
- Modify: `pandora/annotations/__init__.py`
- Test: `tests/test_annotations.py` (reuse `_contact`)

**Interfaces:**
- Produces:
  - `label_interfaces(contacts: Iterable[ResidueContactRow],
    lengths: Mapping[tuple[str, str], int],
    min_residues_per_chain: int | None = None,
    max_total_residues: int | None = None) -> list[LayerRow]` —
    `lengths` keyed by `(entry_id.upper(), chain_id)`.

- [ ] **Step 1: Write the failing tests**

```python
def _lengths(**chains):
    return {("1ABC", c): n for c, n in chains.items()}


def test_label_interfaces_counts_each_side():
    rows = [
        _contact("B", 1, "A", 1),
        _contact("A", 2, "B", 1),
        _contact("A", 2, "B", 2),
        _contact("A", 1, "B", 1),  # repeat of the first, other order
    ]
    (row,) = label_interfaces(rows, _lengths(A=10, B=20))
    assert row.model_dump() == {
        "entry_id": "1abc",
        "chain_id_1": "A",
        "chain_id_2": "B",
        "auth_chain_id_1": "A",
        "auth_chain_id_2": "B",
        "residues_1": 2,
        "residues_2": 2,
        "length_1": 10,
        "length_2": 20,
        "total_residues": 30,
        "label": "ok",
    }


@pytest.mark.parametrize(
    "min_res, max_total, label",
    [
        (None, None, "ok"),
        (3, None, "few_contacts"),
        (2, None, "ok"),
        (None, 29, "too_large"),
        (None, 30, "ok"),
        (3, 29, "too_large"),  # both apply: too_large wins
    ],
)
def test_label_interfaces_rules(min_res, max_total, label):
    rows = [_contact("A", 1, "B", 1), _contact("A", 2, "B", 2)]
    (row,) = label_interfaces(rows, _lengths(A=10, B=20), min_res, max_total)
    assert row.label == label


def test_label_interfaces_unknown_length_is_never_too_large():
    rows = [_contact("A", 1, "B", 1)]
    (row,) = label_interfaces(rows, _lengths(A=10), max_total_residues=1)
    assert (row.length_2, row.total_residues, row.label) == (None, None, "ok")


def test_label_interfaces_entry_case_insensitive_and_sorted():
    rows = [_contact("C", 1, "B", 1), _contact("A", 1, "B", 1)]
    out = label_interfaces(rows, _lengths(A=1, B=1, C=1))
    assert [(r.chain_id_1, r.chain_id_2) for r in out] == [
        ("A", "B"),
        ("B", "C"),
    ]
    assert all(r.total_residues == 2 for r in out)


def test_label_interfaces_passes_template_check_and_keeps_input():
    rows = [_contact("A", 1, "B", 1)]
    before = [r.model_copy() for r in rows]
    out = label_interfaces(rows, _lengths(A=1, B=1))
    check_layer_columns(out[0].model_dump().keys(), "interface_labels")
    assert rows == before
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_annotations.py -k label_interfaces -v`
Expected: FAIL / ImportError on `label_interfaces`.

- [ ] **Step 3: Implement `label_interfaces`**

One pass over `contacts` building, per interface key
`(entry_id.upper(), frozenset)`, the set of residue keys per chain (via
`_side`), the first-seen `entry_id` and `auth_chain_id` per chain; then
emit rows per Global Constraints. Add it to `__all__`. The docstring
states the label order and the `lengths` key.

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_annotations.py -v`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add pandora/annotations tests/test_annotations.py
git commit -m "feat(annotations): label_interfaces labels interfaces by contacts and size"
```

---

### Task 2: `annotate label-interfaces` command

**Files:**
- Modify: `pandora/cli/annotate.py` (handler `_cmd_label`, parser entry
  after `count-contacts`)
- Modify: `pandora/cli/README.md` (one sentence in the annotate
  paragraph, after `count-contacts`)
- Test: `tests/test_cli.py` (reuse `_step_folder`, `_compute`,
  `_filter`, `_contact_pairs`, `_jsonl`)

**Interfaces:**
- Consumes: `label_interfaces` (Task 1); `_layer_origin`,
  `check_contact_layer`, `find_table`, `read_records`,
  `prepare_output_dir`, `write_table`, `build_step_summary` (existing);
  `DatasetChainRecord` (`pandora.schemas.dataset`).

- [ ] **Step 1: Write the failing tests**

`_step_folder` has no `chains` table, so the tests write a small one by
hand: every chain in the contact layer gets `seqres_length = 100`
except those listed in `skip`.

```python
def _label(contacts, out, *extra):
    main(
        ["annotate", "label-interfaces", "--contacts", str(contacts)]
        + ["--output-dir", str(out), "--table-format", "jsonl", *extra]
    )


def _label_inputs(tmp_path, skip=()):
    """1a2c step folder with a chains table, and its contact layer."""

    step = _step_folder(tmp_path)
    contacts = tmp_path / "interfaces"
    _compute(step, contacts, "--layer", "chain_interfaces")
    chains = sorted({c for pair in _contact_pairs(contacts) for c in pair})
    rows = [
        {"entry_id": "1A2C", "chain_id": c, "entity_id": "1"}
        | {"residue_count": 1, "atom_count": 1, "seqres_length": 100}
        for c in chains
        if c not in skip
    ]
    (step / "chains.jsonl").write_text(
        "".join(json.dumps(r) + "\n" for r in rows)
    )
    return step, contacts


def test_annotate_label_interfaces_writes_layer(tmp_path):
    step, contacts = _label_inputs(tmp_path)
    out = tmp_path / "labels"
    _label(contacts, out, "--min-residues-per-chain", "5")

    rows = read_records(contacts / "layer.jsonl", ResidueContactRow)
    lengths = {
        ("1A2C", c): 100
        for pair in _contact_pairs(contacts)
        for c in pair
    }
    written = read_records(out / "layer.jsonl", LayerRow)
    expected = label_interfaces(rows, lengths, 5, None)
    assert [r.model_dump() for r in written] == [
        r.model_dump() for r in expected
    ]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "annotation_label"
    assert summary["policy"]["template"] == "interface_labels"
    assert summary["policy"]["parameters"] == {
        "contacts_layer": "chain_interfaces",
        "structures_dir": str(step.resolve()),
        "min_residues_per_chain": 5,
        "max_total_residues": None,
    }
    assert summary["counts"]["interfaces"] == len(written)
    assert summary["counts"]["unknown_length"] == 0
    assert summary["previous_summary"] == str(contacts / "layer.json")


def test_annotate_label_interfaces_feeds_filter(tmp_path):
    _, contacts = _label_inputs(tmp_path)
    labels = tmp_path / "labels"
    _label(contacts, labels, "--min-residues-per-chain", "1000")
    out = tmp_path / "real"
    _filter(contacts, labels, out, "--drop", "few_contacts")
    assert read_records(out / "layer.jsonl", ResidueContactRow) == []


def test_annotate_label_interfaces_untraceable_origin_exits(tmp_path):
    _, contacts = _label_inputs(tmp_path)
    summary = json.loads((contacts / "layer.json").read_text())
    summary["input"] = str(tmp_path / "gone")
    (contacts / "layer.json").write_text(json.dumps(summary))
    with pytest.raises(SystemExit, match="structures folder"):
        _label(contacts, tmp_path / "labels")
    assert not (tmp_path / "labels").exists()


def test_annotate_label_interfaces_missing_chain_has_unknown_length(
    tmp_path,
):
    _, probe = _label_inputs(tmp_path / "probe")
    first = _contact_pairs(probe)[0][0]
    _, contacts = _label_inputs(tmp_path, skip={first})
    out = tmp_path / "labels"
    _label(contacts, out, "--max-total-residues", "150")
    written = read_records(out / "layer.jsonl", LayerRow)
    unknown = [r for r in written if first in (r.chain_id_1, r.chain_id_2)]
    assert unknown and all(r.total_residues is None for r in unknown)
    assert all(r.label == "ok" for r in unknown)
    assert all(r.label == "too_large" for r in written if r not in unknown)
    counts = json.loads((out / "layer.json").read_text())["counts"]
    assert counts["unknown_length"] == len(unknown)


def test_annotate_label_interfaces_without_chains_table_exits(tmp_path):
    step, contacts = _label_inputs(tmp_path)
    (step / "chains.jsonl").unlink()
    with pytest.raises(SystemExit, match="pandora curate"):
        _label(contacts, tmp_path / "labels")
    assert not (tmp_path / "labels").exists()
```

Also add the output-equals-input case to the existing parametrised
"output dir is the input" test if there is one; otherwise one test
asserting `SystemExit` matching `"--contacts folder"`.

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_cli.py -k label_interfaces -v`
Expected: FAIL (`invalid choice: 'label-interfaces'`).

- [ ] **Step 3: Implement `_cmd_label(args: argparse.Namespace) -> None`**

Order, mirroring `_cmd_count`: check table format and contact layer;
refuse `--output-dir` equal to `--contacts`; `origin =
_layer_origin(contacts)`, and if `None` raise `SystemExit("annotate
label-interfaces: cannot tell which structures folder <contacts> came
from")`; if `find_table(origin, "chains")` is `None`, raise
`SystemExit` saying `<origin>` has no chains table and that it is
written by `pandora curate`; read it with `DatasetChainRecord` and build
`lengths = {(entry_id.upper(), chain_id): seqres_length}` for rows with
a non-null `seqres_length`; read contact rows; call `label_interfaces`; `prepare_output_dir`, `write_table`, write
`layer.json` per Global Constraints (`structures_dir = str(origin)`).
`unknown_length` = rows with `total_residues is None`. Print
`labelled: N interfaces (X ok, Y few_contacts, Z too_large) -> <out>`.
Register the parser with the Global Constraints flags (`type=int` for
thresholds). Add the README sentence.

- [ ] **Step 4: Run tests**

Run: `uv run pytest tests/test_cli.py tests/test_annotations.py -v`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add pandora/cli tests/test_cli.py
git commit -m "feat(cli): annotate label-interfaces writes an interface_labels layer"
```

---

### Task 3: Usage docs

**Files:**
- Modify: `docs/usage/annotation.md` (new `### Labelling interfaces by
  contacts and size` after `### Dropping interfaces by label`)

- [ ] **Step 1: Write the section**

Two or three sentences: what the labels mean and their order, that
lengths are the `seqres_length` column of the `chains` table that
`pandora curate` wrote in the step folder the layer traces to, and that leaving out a threshold gives an all-`ok` layer of
counts for looking at the distribution. Then `library` and `cli` tabs
on `1a2c`, the style of `### Contact counts`. The CLI tab runs
`label-interfaces --min-residues-per-chain 4 --max-total-residues 5000`,
then `filter-contacts --drop few_contacts --drop too_large`, then `jq`
on the counts. Run both tabs by hand and paste the real output into the
`text` blocks.

- [ ] **Step 2: Build the docs**

Run: `uv run zensical build --clean` (follow the `docs-build` skill).
Expected: builds with no new warnings; grep the rendered HTML for
`label-interfaces` to confirm the code shows.

- [ ] **Step 3: Full check and commit**

Run: `uv run pytest && uv run ruff format --check . && uv run ruff check .`
Expected: all pass.

```bash
git add docs/usage/annotation.md
git commit -m "docs: label interfaces by contacts and dimer size"
```
