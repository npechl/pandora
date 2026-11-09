# Annotation layers, part 2 (contacts are one layer) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Pandora's chain interfaces and ligand contacts produce
`residue_contacts` rows, the same shape `annotate import` stores, and
interface records and interface clustering read any `residue_contacts`
rows. `annotate compute` writes a layer folder; `cluster --contacts`
reads one.

**Architecture:** One typed row model, `ResidueContactRow`, in
`pandora/schemas/annotation.py`. The two built-ins in
`pandora/annotations/entry.py` return `list[ResidueContactRow]`. The
consumers (`datasets/records.py`, `similarity/structure.py`) take an
iterable of those rows. The CLI reads and writes them through the
existing `read_records` / `write_table` helpers. `AnnotationLayer`
stays for pairwise identity and reproduce/manifest.

**Tech Stack:** Python 3, pydantic v2, argparse, pytest. No new
dependencies.

**Spec:** `design/specs/2026-11-06-annotation-layers-design.md`, section
"Part 2: contacts are one layer".

## Global Constraints

- PEP 8, `ruff format`, line length 80. `from __future__ import
  annotations` at the top of new modules.
- Type hints everywhere. Public functions get a Google-style docstring
  (`Args:`, `Returns:`, `Raises:`). Internal helpers start with `_` and
  get a one-line docstring.
- Public names in the package `__init__.py` `__all__`.
- Never mutate inputs.
- Catch named exception types only (Ruff BLE001, S110).
- No new dependency.
- Tests offline, `datasets/dev/mmcif/` fixtures or synthetic
  structures built with the helpers already in `tests/test_annotations.py`
  (`_atom`, `_structure`, `_two_chain_structure`, `_load`).
- No transition period for removed CLI flags or outputs (Pandora is 0.x).
- `examples/` is untracked: update the scripts there, but they get no
  tests and are not committed.
- Before calling the work done: `uv run pytest`, `uv run ruff check .`,
  `uv run ruff format --check .`, `uv run zensical build --clean`,
  `uv run --extra docs python docs/scripts/generate_erd.py`.

## Review Focus

1. A residue pair reported in both orders, (A:5, B:7) and (B:7, A:5),
   gives one interface and one residue pair (Task 2 test).
2. One residue pair on several rows (an imported layer with a hydrogen
   bond row and a hydrophobic row) counts once in `residue_pairs` and
   `contact_count` (Task 2 test).
3. A modified residue written as `HETATM` inside a polymer chain is
   neither a ligand nor dropped from the polymer side (Task 4 test).
4. A `layer.parquet` read back gives `seq_id` as float (42.0) and null
   `ins_code`; the row model must turn these into `42` and `None`
   (Task 1 test).
5. `cluster --contacts` on a layer folder whose `layer.json` template is
   not `residue_contacts` stops with a message, before clustering
   (Task 6 test).

## Fixture facts (`datasets/dev/mmcif/a2/1a2c.cif`, defaults: cutoff 4.0, heavy atoms)

Polymer chains (label → author): `A`→`L`, `B`→`H`, `C`→`I`, `D`→`J`.
Chain interfaces: `(A,B)` 74 residue pairs, `(B,C)` 21, `(B,D)` 25.
Ligand contacts: one ligand, `NA` in chain `E` (author `H` 626), with 2
contacting polymer residues. `TYS` (chain `C`) is parsed as `ATOM` in
this file, so the HETATM-in-polymer case needs a synthetic structure.

---

### Task 1: `ResidueContactRow`, `comp_id` on import, `distance` check

**Files:**
- Modify: `pandora/schemas/annotation.py`
- Modify: `pandora/annotations/layers.py` (`_Residue`, `_build_index`,
  `_resolve_residue`, `_bad_value`)
- Test: `tests/test_annotations.py`

**Interfaces:**
- Produces: `ResidueContactRow` (exported from `pandora.schemas.annotation`
  `__all__`, if the module has one, and re-exported where `LayerRow` is):

  ```python
  class ResidueContactRow(BaseModel):
      model_config = ConfigDict(extra="allow")
      entry_id: str
      chain_id_1: str
      seq_id_1: int | None = None
      auth_chain_id_1: str | None = None
      auth_seq_id_1: str | None = None
      ins_code_1: str | None = None
      comp_id_1: str | None = None
      chain_id_2: str
      seq_id_2: int | None = None
      auth_chain_id_2: str | None = None
      auth_seq_id_2: str | None = None
      ins_code_2: str | None = None
      comp_id_2: str | None = None
      contact_type: str | None = None
      distance: float | None = None
  ```

  A float `seq_id` that is a whole number (Parquet's 42.0) validates as
  `42` (a `field_validator(mode="before")`); `ins_code` `""`, `"?"`,
  `"."` become `None`.
- Produces: `resolve_layer_rows` fills `comp_id`, `comp_id_1`,
  `comp_id_2` on every matched residue side; for `residue_contacts`, a
  row whose `distance` is present but not a number (bool counts as not
  a number) is unmatched with reason `bad_value`.

- [ ] **Step 1: Write the failing tests**

```python
def test_residue_contact_row_reads_parquet_shapes():
    row = ResidueContactRow.model_validate(
        {
            "entry_id": "1A2C",
            "chain_id_1": "A",
            "seq_id_1": 42.0,
            "ins_code_1": "",
            "chain_id_2": "E",
            "seq_id_2": None,
            "ins_code_2": "?",
            "distance": "3.5",
            "tool_score": 7,
        }
    )
    assert row.seq_id_1 == 42
    assert row.ins_code_1 is None and row.ins_code_2 is None
    assert row.distance == 3.5
    assert row.model_dump()["tool_score"] == 7


def test_resolve_fills_comp_id():
    # 1a2c label A seq 1 is the first residue of author chain L.
    matched, _ = resolve_layer_rows(
        [
            {
                "entry_id": "1A2C",
                "chain_id_1": "A",
                "seq_id_1": 1,
                "chain_id_2": "B",
                "seq_id_2": 1,
            }
        ],
        "residue_contacts",
        {"1A2C": _load("1a2c")},
    )
    expected = _load("1a2c")  # read the comp ids from the atoms
    comp = {
        (a.label_asym_id, a.label_seq_id): a.label_comp_id
        for a in expected.atoms
    }
    assert matched[0]["comp_id_1"] == comp[("A", 1)]
    assert matched[0]["comp_id_2"] == comp[("B", 1)]


@pytest.mark.parametrize("distance", ["near", True])
def test_resolve_contacts_bad_distance(distance):
    _, unmatched = resolve_layer_rows(
        [
            {
                "entry_id": "1A2C",
                "chain_id_1": "A",
                "seq_id_1": 1,
                "chain_id_2": "B",
                "seq_id_2": 1,
                "distance": distance,
            }
        ],
        "residue_contacts",
        {"1A2C": _load("1a2c")},
    )
    assert unmatched[0]["reason"] == "bad_value"
```

Also extend `test_resolve_author_naming_fills_pandora_keys` /
`test_resolve_pandora_naming_fills_author_keys` expectations with the
`comp_id` key (they compare whole dicts).

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_annotations.py -k "residue_contact_row or comp_id or bad_distance or naming_fills" -v`
Expected: FAIL (`ResidueContactRow` not importable; no `comp_id` key).

- [ ] **Step 3: Implement**

`_Residue` gains `label_comp_id` as a sixth element; `_resolve_residue`
writes it as `comp_id{side}`. Numeric strings for `distance` are fine
(the row model coerces them); only non-numeric strings and bools are
`bad_value`.

- [ ] **Step 4: Run the test file**

Run: `uv run pytest tests/test_annotations.py tests/test_cli.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add pandora/schemas/annotation.py pandora/annotations/layers.py tests/test_annotations.py
git commit -m "feat(annotations): ResidueContactRow; import fills comp_id, checks distance"
```

---

### Task 2: Chain interfaces as rows; interface records from rows

**Files:**
- Modify: `pandora/annotations/entry.py` (`annotate_chain_interfaces`,
  `_chain_pair_contacts`)
- Modify: `pandora/datasets/records.py` (`extract_interface_records`)
- Modify: `pandora/schemas/dataset.py` (`InterfaceRecord`)
- Modify: `pandora/provenance/reproduce.py` (drop `chain_interfaces`
  from `ENTRY_ANNOTATION_DISPATCH`)
- Modify: `pandora/cli/annotate.py` (drop `chain_interfaces` from
  `_ENTRY_ANNOTATORS`; Task 5 rewrites this command)
- Test: `tests/test_annotations.py`

**Interfaces:**
- Consumes: `ResidueContactRow` (Task 1).
- Produces:

  ```python
  def annotate_chain_interfaces(
      structure: Structure,
      distance_cutoff: float = 4.0,
      atom_set: ContactAtomSet = "heavy",
      polymer_types: list[str] | None = None,
  ) -> list[ResidueContactRow]: ...


  def extract_interface_records(
      rows: Iterable[ResidueContactRow],
      distance_cutoff: float | None = None,
      atom_set: ContactAtomSet | None = None,
  ) -> list[InterfaceRecord]: ...
  ```

  `InterfaceRecord.distance_cutoff: float | None = None`;
  `InterfaceRecord.atom_set: ContactAtomSet | None = "all"`.

Row contents for `annotate_chain_interfaces`: one row per residue pair
across two kept polymer chains, side 1 the chain that sorts first;
`seq_id`, `auth_chain_id` (`auth_asym_id`), `auth_seq_id`, `ins_code`
(`pdbx_PDB_ins_code`, `None` for none), `comp_id` from the atoms;
`distance` = closest selected-atom distance for the pair, rounded to 3
decimals (as ligand contacts do); `contact_type` `None`. Rows sorted by
(chain_id_1, chain_id_2, seq_id_1, seq_id_2). Same `ValueError` /
`TypeError` as today.

`extract_interface_records` algorithm: keep rows where both `seq_id`s
are present and `chain_id_1 != chain_id_2`. Key each row by
`(entry_id, frozenset of chains)`; orient it so side 1 is the chain
that sorts first. Per key: `residue_pairs` = sorted set of
`("A:5", "B:7")` strings (`f"{chain_id}:{seq_id}"`),
`interface_residues_chain_1/2` = sorted sets of each side,
`contact_count` = their sizes summed. Records ordered by (entry_id,
chain_id_1, chain_id_2). Residue-id strings sort as strings, as today.

- [ ] **Step 1: Rewrite the chain-interface tests for rows**

Change the test helper `_layer_pairs` to read rows:

```python
def _layer_pairs(rows):
    """{(chain_1, chain_2): {residue pairs}} from chain-interface rows."""
    pairs = {}
    for r in rows:
        pairs.setdefault((r.chain_id_1, r.chain_id_2), set()).add(
            (f"{r.chain_id_1}:{r.seq_id_1}", f"{r.chain_id_2}:{r.seq_id_2}")
        )
    return pairs
```

Update every chain-interface test that reads `layer.data` or
`layer.parameters` to use rows and `_layer_pairs`. Delete
`test_chain_interfaces_rebuild_from_recorded_parameters` (rows carry no
parameters; Task 5 records them in `layer.json`). Rewrite
`test_interface_records_carry_atom_set_and_residue_pairs` as below and
add:

```python
def test_chain_interfaces_rows_on_fixture():
    rows = annotate_chain_interfaces(_load("1a2c"))
    counts = Counter((r.chain_id_1, r.chain_id_2) for r in rows)
    assert counts == {("A", "B"): 74, ("B", "C"): 21, ("B", "D"): 25}
    first = rows[0]
    assert first.auth_chain_id_1 == "L" and first.auth_chain_id_2 == "H"
    assert first.comp_id_1 and first.auth_seq_id_1
    assert all(0 < r.distance <= 4.0 for r in rows)


def test_chain_interfaces_distance_is_closest_atom_pair():
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1", x=0.0),
        _atom(id=2, label_asym_id="A", entity_id="1", x=1.0),
        _atom(id=3, label_asym_id="B", entity_id="2", x=3.5),
    ]
    rows = annotate_chain_interfaces(_two_chain_structure(atoms))
    assert [r.distance for r in rows] == [2.5]


def test_interface_records_from_rows_match_fixture():
    rows = annotate_chain_interfaces(_load("1a2c"))
    records = extract_interface_records(
        rows, distance_cutoff=4.0, atom_set="heavy"
    )
    assert [
        (r.chain_id_1, r.chain_id_2, len(r.residue_pairs)) for r in records
    ] == [("A", "B", 74), ("B", "C", 21), ("B", "D", 25)]
    assert records[0].distance_cutoff == 4.0
    assert records[0].atom_set == "heavy"
    assert records[0].contact_count == (
        len(records[0].interface_residues_chain_1)
        + len(records[0].interface_residues_chain_2)
    )


def test_interface_records_unordered_repeats_ligands_same_chain():
    base = {"entry_id": "X"}
    rows = [
        ResidueContactRow(
            **base,
            chain_id_1="A",
            seq_id_1=5,
            chain_id_2="B",
            seq_id_2=7,
            contact_type="hbond",
        ),
        ResidueContactRow(
            **base,
            chain_id_1="B",
            seq_id_1=7,
            chain_id_2="A",
            seq_id_2=5,
            contact_type="hydrophobic",
        ),
        ResidueContactRow(
            **base, chain_id_1="A", seq_id_1=5, chain_id_2="E", seq_id_2=None
        ),  # ligand
        ResidueContactRow(
            **base, chain_id_1="A", seq_id_1=5, chain_id_2="A", seq_id_2=9
        ),  # same chain
    ]
    (record,) = extract_interface_records(rows)
    assert (record.chain_id_1, record.chain_id_2) == ("A", "B")
    assert record.residue_pairs == [("A:5", "B:7")]
    assert record.contact_count == 2
    assert record.distance_cutoff is None
    assert record.atom_set is None
```

`test_interface_record_without_atom_set_loads_as_all` must keep passing
unchanged.

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_annotations.py -v`
Expected: FAIL (functions still return the old layer / take a structure).

- [ ] **Step 3: Implement both functions and the `InterfaceRecord` change**

`_chain_pair_contacts` returns `{(atom_a residue, atom_b residue):
closest distance}` instead of a set, keyed by the atoms' residue so the
row fields can be filled; keep the cell grid and box pruning as they
are. Update the comment on `InterfaceRecord` in `schemas/dataset.py`
(line ~119) to say it is built from `residue_contacts` rows.

- [ ] **Step 4: Run the suite**

Run: `uv run pytest -v`
Expected: PASS, except `tests/test_structure_similarity.py` and any
docs-snippet tests that still use the old shape — those are Task 3 and
Task 7. If any other test fails, fix it here.

- [ ] **Step 5: Commit**

```bash
git add pandora/annotations/entry.py pandora/datasets/records.py pandora/schemas/dataset.py pandora/provenance/reproduce.py pandora/cli/annotate.py tests/test_annotations.py
git commit -m "feat(annotations)!: chain interfaces as residue_contacts rows; interface records from any contact rows"
```

---

### Task 3: Interface clustering positions from rows

**Files:**
- Modify: `pandora/similarity/structure.py` (`interface_residues_from_annotation`)
- Test: `tests/test_structure_similarity.py`

**Interfaces:**
- Consumes: `ResidueContactRow` (Task 1). Reuse Task 2's
  qualifying-row rule (both `seq_id`s present, different chains) —
  factor it as `_interface_rows(rows)` in `datasets/records.py` if
  reusing it is cleaner; either is fine.
- Produces:

  ```python
  def interface_residues_from_annotation(
      structures: Mapping[str, Structure],
      rows: Iterable[ResidueContactRow],
  ) -> dict[str, set[int]]: ...
  ```

  Same return value as today: `{chain_item_id(entry, chain): {positions}}`.
  A row whose `entry_id` is not in `structures` raises `KeyError` naming
  the entry (as today's `structures[entry_id]` does).

- [ ] **Step 1: Rewrite the test for rows**

`test_interface_residues_from_annotation_maps_to_file_positions` builds
rows instead of an `AnnotationLayer`: `(A:6,B:2)` and `(A:7,B:2)`, plus
one ligand row `(A:6, E:None)` that must not add anything. Expected
result unchanged: `{"ENTRY1_A": {4, 5}, "ENTRY1_B": {2}}`.

- [ ] **Step 2: Run it to verify it fails**

Run: `uv run pytest tests/test_structure_similarity.py -v`
Expected: FAIL.

- [ ] **Step 3: Implement.** Drop the `AnnotationLayer` import if unused.

- [ ] **Step 4: Run**

Run: `uv run pytest tests/test_structure_similarity.py tests/test_annotations.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add pandora/similarity/structure.py tests/test_structure_similarity.py
git commit -m "feat(similarity)!: interface residues from residue_contacts rows"
```

---

### Task 4: Ligand contacts as rows; remove structure counts

**Files:**
- Modify: `pandora/annotations/entry.py` (`annotate_ligand_contacts`,
  `_is_ligand_atom`, `_residues_within_cutoff`; delete
  `annotate_structure_counts`)
- Modify: `pandora/annotations/__init__.py`, `pandora/metadata/__init__.py`
  (exports)
- Modify: `pandora/provenance/reproduce.py` (delete
  `ENTRY_ANNOTATION_DISPATCH` and its branch in `_reproduce_annotations`;
  only the pairwise branch remains)
- Modify: `pandora/cli/annotate.py` (drop the remaining entries of
  `_ENTRY_ANNOTATORS`; `compute` keeps only `--pairwise` until Task 5)
- Test: `tests/test_annotations.py`, `tests/test_provenance.py`,
  `tests/test_cli.py`

**Interfaces:**
- Consumes: `ResidueContactRow` (Task 1).
- Produces:

  ```python
  def annotate_ligand_contacts(
      structure: Structure,
      distance_cutoff: float = 4.0,
      include_waters: bool = False,
  ) -> list[ResidueContactRow]: ...
  ```

Rule: a ligand residue is any residue of a **non-polymer** chain
(`label_asym_id` not in `polymer_asym_ids(structure)`), waters
(`WATER_COMP_IDS`) only with `include_waters`. The polymer side is every
atom of a polymer chain, whatever its `group_PDB` (today it is
`group_PDB == "ATOM"`; a HETATM modified residue inside a chain must
still count as polymer). One row per (polymer residue, ligand residue):
side 1 polymer (both namings, `comp_id`), side 2 ligand (`chain_id_2` =
ligand `label_asym_id`, `seq_id_2` `None`, `auth_chain_id_2`,
`auth_seq_id_2`, `ins_code_2`, `comp_id_2`), `distance` = nearest,
rounded to 3 decimals. Sorted by (chain_id_2, auth_seq_id_2, chain_id_1,
seq_id_1).

- [ ] **Step 1: Rewrite / add the tests**

Rewrite `test_ligand_contacts_reports_nearest_polymer_residue` for rows
(one row: `seq_id_1 == 1`, `chain_id_2 == "B"`, `seq_id_2 is None`,
`auth_seq_id_2 == "101"`, `comp_id_2 == "ZN"`, `distance == 1.0`). Add:

```python
def test_ligand_contacts_fixture():
    rows = annotate_ligand_contacts(_load("1a2c"))
    assert {
        (r.chain_id_2, r.auth_chain_id_2, r.auth_seq_id_2, r.comp_id_2)
        for r in rows
    } == {("E", "H", "626", "NA")}
    assert len(rows) == 2


def test_ligand_contacts_hetatm_in_polymer_chain_is_not_a_ligand():
    # A: polymer with a HETATM modified residue (MSE) at seq 2.
    # B: second polymer chain. C: a real ligand near A:2.
    atoms = [
        _atom(id=1, label_asym_id="A", label_seq_id=1, entity_id="1", x=0.0),
        _atom(
            id=2,
            group_PDB="HETATM",
            label_asym_id="A",
            label_seq_id=2,
            label_comp_id="MSE",
            entity_id="1",
            x=10.0,
        ),
        _atom(id=3, label_asym_id="B", label_seq_id=1, entity_id="2", x=11.0),
        _atom(
            id=4,
            group_PDB="HETATM",
            label_asym_id="C",
            label_seq_id=None,
            label_comp_id="ZN",
            auth_seq_id="101",
            entity_id="3",
            x=9.0,
        ),
    ]
    entities = [
        EntityRecord(id="1", type="polymer"),
        EntityRecord(id="2", type="polymer"),
        EntityRecord(id="3", type="non-polymer"),
    ]
    asym_units = [
        AsymRecord(id=c, entity_id=e)
        for c, e in (("A", "1"), ("B", "2"), ("C", "3"))
    ]
    rows = annotate_ligand_contacts(_structure(atoms, entities, asym_units))
    assert {r.chain_id_2 for r in rows} == {"C"}  # MSE not a ligand
    assert ("A", 2) in {
        (r.chain_id_1, r.seq_id_1) for r in rows
    }  # MSE as polymer
    assert extract_interface_records(rows) == []
```

Delete `test_reproduce_dataset_regenerates_annotations` from
`tests/test_provenance.py` (it replays `structure_counts`, which no
longer exists; no entry-level replay remains). Change
`test_annotate_compute_writes_annotations` in `tests/test_cli.py` to use
`--pairwise` with two fixtures, or delete it — Task 5 replaces it.

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_annotations.py -k ligand -v`
Expected: FAIL.

- [ ] **Step 3: Implement, delete `annotate_structure_counts`**

`grep -rn structure_counts pandora tests` must return nothing afterwards.

- [ ] **Step 4: Run**

Run: `uv run pytest -v`
Expected: PASS (docs-snippet failures, if any, are Task 7).

- [ ] **Step 5: Commit**

```bash
git add -u pandora tests
git commit -m "feat(annotations)!: ligand contacts as residue_contacts rows; drop structure counts"
```

---

### Task 5: `pandora annotate compute` writes a layer folder

**Files:**
- Modify: `pandora/cli/annotate.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `annotate_chain_interfaces`, `annotate_ligand_contacts`
  (Tasks 2, 4); existing `prepare_output_dir`, `read_structures_dir`,
  `write_table`, `check_table_format`, `build_step_summary`,
  `write_json`, `LayerDescription`; `_LAYER_OUTPUTS` /
  `_PREVIOUS_SUMMARIES` already in this module.
- Produces: the command below. `--pairwise`, `--layers`, `--input-dir`
  and `annotations.json` are gone. `annotate_pairwise_sequence_identity`
  stays in the library but has no CLI.

```
pandora annotate compute --layer {chain_interfaces,ligand_contacts}
    --structures-dir DIR --output-dir DIR
    [--distance-cutoff FLOAT] [--atom-set {all,heavy,backbone,ca}]
    [--polymer-types TYPE]...  [--include-waters]
    [--overwrite] [--table-format {parquet,jsonl}]
```

- `--atom-set` / `--polymer-types` with `--layer ligand_contacts`, or
  `--include-waters` with `--layer chain_interfaces`: `SystemExit` with
  a message naming the flag, before anything is written. (Flags default
  to `None` so "not given" is detectable; the function defaults apply.)
- Writes `layer.<fmt>` (all entries' rows; not written when empty) and
  `layer.json`: `build_step_summary("annotation_compute",
  LayerDescription(name=args.layer, template="residue_contacts",
  tool="pandora", tool_version=<pandora version as other summaries get
  it>, parameters=<the settings actually used, defaults included>),
  input=structures_dir, output=out, started_at=..., counts={"entries":
  n, "rows": n}, previous_summary=<same lookup as import>)`.
- Refuses an existing layer without `--overwrite` (reuse
  `prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)`).
- Prints `computed: <rows> rows for <entries> entries -> <out>`.

- [ ] **Step 1: Write the failing tests** (replace
  `test_annotate_compute_writes_annotations`)

```python
def _step_folder(tmp_path):
    structures = tmp_path / "curated" / "structures" / "a2"
    structures.mkdir(parents=True)
    shutil.copy(FIXTURES_DIR / "a2" / "1a2c.cif", structures / "1a2c.cif")
    return tmp_path / "curated"


def test_annotate_compute_writes_contact_layer(tmp_path):
    out = tmp_path / "annotations" / "interfaces"
    main(
        [
            "annotate",
            "compute",
            "--layer",
            "chain_interfaces",
            "--structures-dir",
            str(_step_folder(tmp_path)),
            "--atom-set",
            "heavy",
            "--output-dir",
            str(out),
            "--table-format",
            "jsonl",
        ]
    )
    rows = read_records(out / "layer.jsonl", ResidueContactRow)
    assert len(rows) == 74 + 21 + 25
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "annotation_compute"
    assert summary["policy"]["template"] == "residue_contacts"
    assert summary["policy"]["tool"] == "pandora"
    assert summary["policy"]["parameters"] == {
        "distance_cutoff": 4.0,
        "atom_set": "heavy",
        "polymer_types": None,
    }
    assert summary["counts"] == {"entries": 1, "rows": 120}


def test_annotate_compute_rejects_flag_for_other_layer(tmp_path):
    out = tmp_path / "out"
    with pytest.raises(SystemExit, match="--include-waters"):
        main(
            [
                "annotate",
                "compute",
                "--layer",
                "chain_interfaces",
                "--include-waters",
                "--structures-dir",
                str(_step_folder(tmp_path)),
                "--output-dir",
                str(out),
            ]
        )
    assert not out.exists()


def test_annotate_compute_refuses_existing_layer(tmp_path):
    args = [
        "annotate",
        "compute",
        "--layer",
        "ligand_contacts",
        "--structures-dir",
        str(_step_folder(tmp_path)),
        "--output-dir",
        str(tmp_path / "out"),
        "--table-format",
        "jsonl",
    ]
    main(args)
    with pytest.raises(SystemExit):
        main(args)
    main([*args, "--overwrite"])


def test_compute_and_import_give_same_interface_records(tmp_path):
    folder = _step_folder(tmp_path)
    computed = tmp_path / "computed"
    main(
        [
            "annotate",
            "compute",
            "--layer",
            "chain_interfaces",
            "--structures-dir",
            str(folder),
            "--output-dir",
            str(computed),
            "--table-format",
            "jsonl",
        ]
    )
    # Re-import the same contacts as a tool would give them: author naming
    # only, sides swapped.
    table = tmp_path / "tool.jsonl"
    table.write_text(
        "".join(
            json.dumps(
                {
                    "entry_id": r.entry_id.lower(),
                    "auth_chain_id_1": r.auth_chain_id_2,
                    "auth_seq_id_1": r.auth_seq_id_2,
                    "ins_code_1": r.ins_code_2,
                    "auth_chain_id_2": r.auth_chain_id_1,
                    "auth_seq_id_2": r.auth_seq_id_1,
                    "ins_code_2": r.ins_code_1,
                }
            )
            + "\n"
            for r in read_records(computed / "layer.jsonl", ResidueContactRow)
        )
    )
    (tmp_path / "d.yaml").write_text(
        "name: tool\ntemplate: residue_contacts\ntool: T\n"
    )
    imported = tmp_path / "imported"
    main(
        [
            "annotate",
            "import",
            "--structures-dir",
            str(folder),
            "--table",
            str(table),
            "--description",
            str(tmp_path / "d.yaml"),
            "--output-dir",
            str(imported),
            "--table-format",
            "jsonl",
        ]
    )
    a = extract_interface_records(
        read_records(computed / "layer.jsonl", ResidueContactRow)
    )
    b = extract_interface_records(
        read_records(imported / "layer.jsonl", ResidueContactRow)
    )
    assert a == b and len(a) == 3
```

Keep `test_annotate_old_flat_form_is_gone`.

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_cli.py -k annotate -v`
Expected: FAIL (no `--layer`).

- [ ] **Step 3: Implement.** Writing rows: `write_table(rows, out,
  "layer", args.table_format)`; check whether `write_table` already skips
  empty tables (part 1 relies on it) and keep that.

- [ ] **Step 4: Run**

Run: `uv run pytest tests/test_cli.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add pandora/cli/annotate.py tests/test_cli.py
git commit -m "feat(cli)!: annotate compute writes one contact layer folder"
```

---

### Task 6: `pandora cluster --contacts --structures-dir`

**Files:**
- Modify: `pandora/cli/similarity.py` (cluster parser and handler)
- Modify: `pandora/cli/_io.py` (delete `load_interface_residues`; add
  a layer reader if the handler gets long)
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `interface_residues_from_annotation(structures, rows)`
  (Task 3), `ResidueContactRow`, `read_records`, `read_structures_dir`.
- Produces: `--contacts DIR` (a layer folder) and `--structures-dir
  DIR` (the step folder the layer refers to) replace
  `--interface-residues`. Behaviour:
  - `--contacts` without `--structures-dir` (or the reverse):
    `SystemExit` naming the missing flag.
  - Reads `<contacts>/layer.json`; its `policy.template` must be
    `residue_contacts`, else `SystemExit("cluster: <dir> is a <template>
    layer, --contacts needs residue_contacts")`. Missing `layer.json`:
    `SystemExit` too.
  - Reads `layer.parquet` or `layer.jsonl`, whichever exists; neither
    (an empty layer) means no interface residues, i.e. `{}`.
  - `--input-dir` is unchanged (the per-chain files Foldseek ran on).

- [ ] **Step 1: Write the failing tests**

```python
def test_cluster_reads_contact_layer(tmp_path, monkeypatch):
    folder = _step_folder(tmp_path)
    layer = tmp_path / "interfaces"
    main(
        [
            "annotate",
            "compute",
            "--layer",
            "chain_interfaces",
            "--structures-dir",
            str(folder),
            "--output-dir",
            str(layer),
            "--table-format",
            "jsonl",
        ]
    )
    captured = {}

    def fake_cluster(item_ids, search, hit_filter, interface_residues=None):
        captured["residues"] = interface_residues
        return [], ClusterProvenance(...)  # build the minimal valid one

    monkeypatch.setattr(
        "pandora.cli.similarity.cluster_similar_items", fake_cluster
    )
    input_dir, search = _cluster_inputs(tmp_path)
    main(
        [
            "cluster",
            "--input-dir",
            str(input_dir),
            "--search",
            str(search),
            "--min-score",
            "0.5",
            "--contacts",
            str(layer),
            "--structures-dir",
            str(folder),
            "--output",
            str(tmp_path / "c.json"),
        ]
    )
    assert set(captured["residues"]) == {"1A2C_A", "1A2C_B", "1A2C_C", "1A2C_D"}


def test_cluster_rejects_non_contact_layer(tmp_path):
    layer = tmp_path / "labels"
    layer.mkdir()
    (layer / "layer.json").write_text(
        json.dumps({"policy": {"template": "chain_labels"}})
    )
    input_dir, search = _cluster_inputs(tmp_path)
    with pytest.raises(SystemExit, match="residue_contacts"):
        main(
            [
                "cluster",
                "--input-dir",
                str(input_dir),
                "--search",
                str(search),
                "--min-score",
                "0.5",
                "--contacts",
                str(layer),
                "--structures-dir",
                str(tmp_path),
                "--output",
                str(tmp_path / "c.json"),
            ]
        )
```

Use whatever provenance object `cluster_similar_items` returns today for
the fake (check `pandora/similarity/clustering.py`); if building it is
awkward, let the real function run on an MMseqs2 search and instead
monkeypatch `interface_residues_from_annotation` to capture its rows.

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_cli.py -k cluster -v`
Expected: FAIL (`--contacts` unknown).

- [ ] **Step 3: Implement.**

- [ ] **Step 4: Run**

Run: `uv run pytest -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add pandora/cli/similarity.py pandora/cli/_io.py tests/test_cli.py
git commit -m "feat(cli)!: cluster --contacts reads a residue_contacts layer"
```

---

### Task 7: Docs, README, examples, ERD

**Files:**
- Modify: `docs/usage/annotation.md`, `docs/usage/datasets.md`,
  `docs/usage/similarity.md`, `docs/usage/cli.md`,
  `docs/usage/provenance.md`, `docs/getting-started/quickstart.md`,
  `docs/index.md`, `docs/recipes/ppi-01.md`,
  `docs/contributing/architecture.md`, `README.md`,
  `pandora/cli/README.md`
- Modify (untracked, not committed): `examples/multichain_ligand_complex.py`,
  `examples/ppi_dataset_pipeline.py`, `examples/overview.py`
- Regenerate: ERD via `docs/scripts/generate_erd.py`

Content changes:
- Every snippet reading `.data["interfaces"]` / `.data["ligands"]` /
  `.parameters` reads rows instead. `extract_interface_records(structure,
  4.0)` becomes `extract_interface_records(annotate_chain_interfaces(
  structure, 4.0), distance_cutoff=4.0)`.
- Every `annotate_structure_counts` mention removed; where it showed
  per-chain counts, point to chain records.
- `docs/usage/annotation.md`: one section saying built-in contacts and
  imported contacts are the same `residue_contacts` layer, with the
  `annotate compute` command (library + CLI tabs, per the docs
  convention); the ligand rule (non-polymer chains only; modified
  residues are not ligands); remove the admonition that cutoff/atom set
  aren't CLI flags.
- `docs/usage/similarity.md`: the interface-coverage recipe uses
  `annotate compute` + `cluster --contacts --structures-dir`.
- `docs/usage/cli.md`: `annotate compute` usage; remove
  `--annotations annotations/annotations.json` from the `manifest`
  example; `cluster` flags.
- `pandora/cli/README.md`: the `annotate compute` and `cluster`
  paragraphs.
- `docs/contributing/architecture.md`: function → model table rows.

- [ ] **Step 1: Edit the docs and README** (load the `docs-build` skill
  first for the project's doc conventions).
- [ ] **Step 2: Run the example scripts** — each must run without error:
  `uv run python examples/overview.py`, `.../multichain_ligand_complex.py`
  (skip `ppi_dataset_pipeline.py` if it needs network).
- [ ] **Step 3: Regenerate the ERD**

Run: `uv run --extra docs python docs/scripts/generate_erd.py`
Expected: diagram files updated (`InterfaceRecord` optional fields).

- [ ] **Step 4: Full check**

Run: `uv run pytest && uv run ruff check . && uv run ruff format --check . && uv run zensical build --clean`
Expected: all pass; `grep -rn "structure_counts\|annotations.json\|interface-residues\|data\[\"interfaces\"\]\|data\[\"ligands\"\]" docs README.md pandora` returns nothing (except
`manifest`/`reproduce` code paths that legitimately mention `AnnotationLayer`).

- [ ] **Step 5: Commit**

```bash
git add docs README.md pandora/cli/README.md
git commit -m "docs: contacts as residue_contacts layers; annotate compute; cluster --contacts"
```
