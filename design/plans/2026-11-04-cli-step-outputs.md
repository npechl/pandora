# CLI step outputs — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `pandora canonicalise` and `pandora curate` write the same step folders as the PPI scripts; the old JSON outputs go away; every other command reads the new layouts.

**Architecture:** Two small library additions (`read_records`, provenance-from-row), CLI file helpers that understand flat / mirror / step folders, rewritten `canonicalise` / `curate` / `manifest` handlers, then tests and docs.

**Tech Stack:** Python 3, pydantic v2, argparse, pandas/pyarrow (`export` extra), pytest, jq (docs examples).

**Spec:** `design/specs/2026-11-04-cli-step-outputs-design.md`

## Global Constraints

- CLAUDE.md coding standards (future import, type hints, 80 cols, Google docstrings, `__all__`, models in `schemas/`, no input mutation, no new dependencies, pandas imported lazily).
- Old files are not written anywhere: `canonicalisation_provenance.json`, `canonicalisation_mappings.json`, `curation_provenance.json`, `curation_exclusions.json`.
- Tables: `--table-format {parquet,jsonl}`, default `parquet`. A table with no rows is not written.
- Caught per entry: `ValueError, OSError, RuntimeError, KeyError, IndexError` → `failures.jsonl` row `{"path", "error_type", "message"}`.
- Printed lines stay as today so docs churn is small: `canonicalised {n} structures -> {dir}` and `curated: {kept} retained, {excluded} excluded -> {dir}`; append `, {k} failed` only when `k > 0`.
- Structure file stems: `canonicalise` / `curate` use the input file's stem as spelled; `write_structures_dir` (dedup, reproduce) keeps today's `entry_id.lower()`.
- **No commits.**

## Review Focus

1. **A failing entry** (fixture `datasets/dev/mmcif/2a/22ag.cif`, which fails canonicalisation validation) must not stop the run: the others are written, `failures.jsonl` has its row. Pinned in Task 4.
2. **Everything excluded** (strict policy): no `entries`/`chains` table files, summary counts say 0, no crash. Pinned in Task 5.
3. **`curate` on a plain flat folder** (no step folder, no `chain_copies`): works, `copy_of` is null, `previous_summary` is null. Pinned in Task 5.
4. **Parquet nulls** read back as `None`, not NaN. Pinned in Task 1.
5. **Downstream commands on a step folder**: `similarity`/`annotate`/`dedup`/`manifest` given `curate`'s output read its structures. Pinned in Task 3 (`read_structures_dir`) and Task 6 (`manifest`).

---

### Task 1: `read_records`

**Files:** Modify `pandora/export/records.py`, `pandora/export/__init__.py`. Test `tests/test_export.py`.

**Interfaces — Produces:** `read_records(path: str | Path, model: type[M]) -> list[M]` (M bound to `BaseModel`). Dispatch on suffix like `write_records`: `.jsonl` (one object per non-blank line), `.json` (array), `.parquet` (pandas, lazy import, same `RuntimeError` message style). Parquet: convert NaN/NA to `None` and numpy arrays (list columns) to lists before `model.model_validate`. Unsupported suffix → `ValueError`.

- [ ] **Step 1: Failing tests**

```python
@pytest.mark.parametrize("suffix", [".jsonl", ".json", ".parquet"])
def test_read_records_round_trips_with_nulls_and_lists(tmp_path, suffix):
    rows = [
        DatasetChainRecord(
            entry_id="1ABC",
            chain_id="A",
            entity_id="1",
            residue_count=1,
            atom_count=1,
            copy_of=None,
            resolution=None,
        ),
        DatasetChainRecord(
            entry_id="1ABC",
            chain_id="A_2",
            entity_id="1",
            residue_count=1,
            atom_count=1,
            copy_of="A",
            resolution=2.0,
        ),
    ]
    path = write_records(rows, tmp_path / f"t{suffix}")
    assert read_records(path, DatasetChainRecord) == rows


def test_read_records_list_column_parquet(tmp_path):
    row = ChainCopyRecord(
        entry_id="1", chain_id="A_2", source_chain_id="A", operator_id="2"
    )
    entry = CanonicalisationEntryRecord(
        entry_id="1",
        policy_id="p",
        policy_name="n",
        policy_version="1",
        canonicalised_at="t",
        transforms=["a:b"],
        dropped_raw_categories=[],
    )
    for r, m in ((row, ChainCopyRecord), (entry, CanonicalisationEntryRecord)):
        p = write_records([r], tmp_path / f"{m.__name__}.parquet")
        assert read_records(p, m) == [r]


def test_read_records_rejects_bad_suffix(tmp_path):
    with pytest.raises(ValueError):
        read_records(tmp_path / "x.csv", ChainCopyRecord)
```

- [ ] **Step 2: Run, expect FAIL** (ImportError): `uv run pytest tests/test_export.py -k read_records -q`
- [ ] **Step 3: Implement and export.**
- [ ] **Step 4: Run, expect PASS.**

---

### Task 2: Provenance from rows

**Files:** Modify `pandora/provenance/records.py`, `pandora/provenance/__init__.py`. Test `tests/test_provenance.py`.

**Interfaces — Produces:**
- `canonicalisation_provenance_from_record(record: CanonicalisationEntryRecord) -> canonicalisationProvenance`: `report` = `{"warnings": w, "errors": e}` when both are not None, else `{}`.
- `curation_provenance_from_record(record: CurationEntryRecord) -> CurationProvenance`.

- [ ] **Step 1: Failing test** — round trip: build `cprov`/`kprov` exactly as in `test_entry_records_write_to_parquet`, then

```python
assert (
    canonicalisation_provenance_from_record(
        canonicalisation_entry_record(canon, cprov)
    )
    == cprov
)
assert (
    curation_provenance_from_record(
        curation_entry_record(curated, metadata, kprov)
    )
    == kprov
)
reported = cprov.model_copy(update={"report": {"warnings": 2, "errors": 0}})
assert (
    canonicalisation_provenance_from_record(
        canonicalisation_entry_record(canon, reported)
    )
    == reported
)
```

- [ ] **Step 2: Run, expect FAIL.** `uv run pytest tests/test_provenance.py -k from_record -q`
- [ ] **Step 3: Implement and export.**
- [ ] **Step 4: Run, expect PASS.**

---

### Task 3: CLI file helpers understand step folders

**Files:** Modify `pandora/cli/_io.py`. Test `tests/test_cli.py`.

**Interfaces — Produces (internal to the CLI, `_io.py`):**
- `structure_files(folder: Path) -> list[Path]`: if `folder / "structures"` is a dir, `list_mirror_files(folder / "structures")`; else sorted union of `folder.glob("*.cif")` and `list_mirror_files(folder)`.
- `read_structures_dir(input_dir)` iterates `structure_files(input_dir)` (rest unchanged).
- `write_structures_dir(structures, output_dir)` writes `mirror_path(output_dir, entry_id.lower())`.
- `find_table(folder: Path, name: str) -> Path | None`: `folder/f"{name}.parquet"` or `.jsonl`, whichever exists (parquet first).
- `write_table(rows: list[BaseModel], folder: Path, name: str, table_format: str) -> None`: `write_records(rows, folder / f"{name}.{table_format}")` if `rows`, else nothing.

- [ ] **Step 1: Failing tests**

```python
def test_structure_files_reads_flat_mirror_and_step_folders(tmp_path):
    src = FIXTURES_DIR / "04" / "104m.cif"
    flat = tmp_path / "flat"
    flat.mkdir()
    shutil.copy(src, flat / "104m.cif")
    mirror = tmp_path / "mirror"
    (mirror / "04").mkdir(parents=True)
    shutil.copy(src, mirror / "04" / "104m.cif")
    step = tmp_path / "step"
    (step / "structures" / "04").mkdir(parents=True)
    shutil.copy(src, step / "structures" / "04" / "104m.cif")
    for folder in (flat, mirror, step):
        assert list(read_structures_dir(folder)) == ["104M"]


def test_write_structures_dir_uses_mirror_layout(tmp_path):
    s, _, _ = mmcif_to_structure(str(FIXTURES_DIR / "04" / "104m.cif"))
    write_structures_dir({"104M": s}, tmp_path)
    assert (tmp_path / "04" / "104m.cif").exists()
```

- [ ] **Step 2: Run, expect FAIL** (`[]` for mirror/step; flat file location). `uv run pytest tests/test_cli.py -k "structure_files or mirror_layout" -q`
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run, expect PASS.** Fix the existing chain test's `dedup` assertion to `deduped_dir / "04" / "104m.cif"` and its `export --input` path accordingly (Task 5 rewrites the rest of that test).

---

### Task 4: `pandora canonicalise`

**Files:** Rewrite `_cmd_canonicalise` and `register` in `pandora/cli/canonicalisation.py`. Test `tests/test_cli.py`.

**Interfaces — Consumes:** Task 3 helpers; `read_entry_paths`, `canonicalise_structure`, `structure_to_mmcif`, `mirror_path`, `canonicalisation_entry_record`, `chain_copy_records`, `build_step_summary`, `write_json`, `now_iso`.

**Arguments:** mutually exclusive required group `--input-dir DIR` | `--entries FILE`; `--policy` (required); `--output-dir` (required); `--table-format {parquet,jsonl}` default `parquet`.

**Behaviour:** inputs = `structure_files(DIR)` or `read_entry_paths(FILE)`. Clear nothing (the CLI never deletes; it overwrites files it writes). Per path: parse (None → `ValueError("could not parse …")`), canonicalise, write `mirror_path(out / "structures", path.stem)`, collect entry row (`structure_path` relative to `out`), copies, and one `mappings.jsonl` line `{"entry_id": ..., **mappings.model_dump(mode="json")}` written as the loop runs. Then tables `entries`, `chain_copies`; `canonicalisation.json` summary (`input` = DIR or FILE; counts `entries_in`, `entries_kept`, `entries_failed`, sorted).

- [ ] **Step 1: Failing tests**

```python
def test_canonicalise_writes_step_folder(tmp_path):
    out = tmp_path / "canon"
    raw = tmp_path / "raw"
    raw.mkdir()
    for p in ("04/104m.cif", "2a/22ag.cif"):
        shutil.copy(FIXTURES_DIR / p, raw / Path(p).name)
    main(
        [
            "canonicalise",
            "--input-dir",
            str(raw),
            "--policy",
            str(CANON_POLICY),
            "--output-dir",
            str(out),
            "--table-format",
            "jsonl",
        ]
    )
    assert (out / "structures" / "04" / "104m.cif").exists()
    entries = [
        json.loads(l) for l in (out / "entries.jsonl").read_text().splitlines()
    ]
    assert [e["entry_id"] for e in entries] == ["104M"]
    assert entries[0]["structure_path"] == "structures/04/104m.cif"
    mappings = [
        json.loads(l) for l in (out / "mappings.jsonl").read_text().splitlines()
    ]
    assert (
        mappings[0]["entry_id"] == "104M"
        and mappings[0]["chain_id_mapping"]["items"]
    )
    failures = [
        json.loads(l) for l in (out / "failures.jsonl").read_text().splitlines()
    ]
    assert [Path(f["path"]).name for f in failures] == ["22ag.cif"]
    summary = json.loads((out / "canonicalisation.json").read_text())
    assert summary["counts"] == {
        "entries_failed": 1,
        "entries_in": 2,
        "entries_kept": 1,
    }
    assert not (out / "canonicalisation_provenance.json").exists()
    assert not (out / "canonicalisation_mappings.json").exists()


def test_canonicalise_reads_entries_file_and_writes_parquet(tmp_path):
    entries = tmp_path / "entries.txt"
    entries.write_text(str(FIXTURES_DIR / "04" / "104m.cif") + "\n")
    out = tmp_path / "canon"
    main(
        [
            "canonicalise",
            "--entries",
            str(entries),
            "--policy",
            str(CANON_POLICY),
            "--output-dir",
            str(out),
        ]
    )
    assert (
        read_records(out / "entries.parquet", CanonicalisationEntryRecord)[
            0
        ].entry_id
        == "104M"
    )
```

(If `22ag` stops failing — its validation failure is a known open issue — swap in a file that does not parse, e.g. a `.cif` with the text `garbage`, and keep the assertions.)

- [ ] **Step 2: Run, expect FAIL.** `uv run pytest tests/test_cli.py -k canonicalise -q`
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run, expect PASS.**

---

### Task 5: `pandora curate`

**Files:** Rewrite `_cmd_curate` and its parser in `pandora/cli/datasets.py`. Test `tests/test_cli.py`.

**Interfaces — Consumes:** Task 1 `read_records`, Task 3 helpers; `collect_metadata`, `curate_structure`, `dataset_chain_records`, `curation_entry_record`, `build_step_summary`.

**Arguments:** `--input-dir`, `--policy`, `--output-dir`, `--table-format` as Task 4.

**Behaviour:** files = `structure_files(DIR)`; copies from `find_table(DIR, "chain_copies")` via `read_records(..., ChainCopyRecord)` → `{entry_id: {chain_id: source}}`. Per file: parse, metadata, curate; extend exclusions; if kept: write `mirror_path(out / "structures", path.stem)`, rows via `dataset_chain_records(curated, metadata, policy.quality_rules.missing_residue_definition, copy_of=..., structure_path=<relative>)`, entry row. Merge exclusions/rows into the shared lists only after the entry's `try` succeeds (review minor 4 of the last change). Tables `entries`, `chains`, `exclusions`; `curation.json` with `previous_summary` = `DIR / "canonicalisation.json"` if it exists else None; counts as in `examples/ppi/v2.1/curate.py` (sorted). Print line per Global Constraints; `excluded` = `entries_in − kept − failed`.

- [ ] **Step 1: Failing tests** — rewrite `test_canonicalise_curate_dedup_export_chain` to: canonicalise `104m` (jsonl) → curate the canonicalise folder → assert `curated/structures/04/104m.cif`, `curated/chains.jsonl` rows for 104M with `copy_of` key, `curation.json["previous_summary"] == str(canonical_dir / "canonicalisation.json")`, no `curation_provenance.json` / `curation_exclusions.json`; then `dedup --input-dir curated` → `deduped/04/104m.cif`; `export --input deduped/04/104m.cif`. Plus:

```python
def test_curate_flat_folder_and_everything_excluded(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copy(FIXTURES_DIR / "04" / "104m.cif", raw / "104m.cif")
    strict = tmp_path / "strict.yaml"
    strict.write_text(
        "policy_id: s\npolicy_name: s\npolicy_version: 1\n"
        "quality_rules:\n  max_resolution: 1.0\n"
    )
    out = tmp_path / "strict"
    main(
        [
            "curate",
            "--input-dir",
            str(raw),
            "--policy",
            str(strict),
            "--output-dir",
            str(out),
            "--table-format",
            "jsonl",
        ]
    )
    assert (
        not (out / "entries.jsonl").exists()
        and not (out / "chains.jsonl").exists()
    )
    ex = [
        json.loads(l)
        for l in (out / "exclusions.jsonl").read_text().splitlines()
    ]
    assert ex[0]["reason_code"] == "RESOLUTION_THRESHOLD"
    summary = json.loads((out / "curation.json").read_text())
    assert (
        summary["counts"]["entries_kept"] == 0
        and summary["previous_summary"] is None
    )


def test_curate_fills_copy_of_from_canonicalise_folder(tmp_path):
    # 1ayi's assembly 1 has a symmetry copy A_2 of chain A.
    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copy(FIXTURES_DIR / "ay" / "1ayi.cif", raw / "1ayi.cif")
    canon_policy = tmp_path / "canon.yaml"
    canon_policy.write_text(
        "policy_id: c\npolicy_name: c\npolicy_version: 1\n"
        "assembly_rules:\n  strategy: standardize_biological_assembly\n"
        "  preferred_assembly_source: first\n"
    )
    curation = tmp_path / "cur.yaml"
    curation.write_text("policy_id: k\npolicy_name: k\npolicy_version: 1\n")
    main(
        [
            "canonicalise",
            "--input-dir",
            str(raw),
            "--policy",
            str(canon_policy),
            "--output-dir",
            str(tmp_path / "c"),
            "--table-format",
            "jsonl",
        ]
    )
    main(
        [
            "curate",
            "--input-dir",
            str(tmp_path / "c"),
            "--policy",
            str(curation),
            "--output-dir",
            str(tmp_path / "k"),
            "--table-format",
            "jsonl",
        ]
    )
    rows = {
        r["chain_id"]: r
        for r in map(
            json.loads,
            (tmp_path / "k" / "chains.jsonl").read_text().splitlines(),
        )
    }
    assert rows["A_2"]["copy_of"] == "A" and rows["A"]["copy_of"] is None
```

- [ ] **Step 2: Run, expect FAIL.** `uv run pytest tests/test_cli.py -k "curate" -q`
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run, expect PASS.**

---

### Task 6: `pandora manifest`

**Files:** Modify `_cmd_manifest` and its parser in `pandora/cli/provenance.py`. Test `tests/test_cli.py`.

**Interfaces — Consumes:** Task 1, Task 2, Task 3 `find_table`.

**Arguments:** remove `--canonicalisation-provenance`; add `--canonicalisation-dir DIR` and `--curation-dir DIR`. `--exclusions` paths are read with `read_records(path, ExclusionRecord)` (any of `.json`, `.jsonl`, `.parquet`). `--structures-dir` uses `read_structures_dir` (step folders work via Task 3).

**Behaviour:** per-entry `canonicalisation` from `find_table(DIR, "entries")` rows via Task 2; same for `curation`; pass both to `build_provenance_bundle`.

- [ ] **Step 1: Failing test** — after canonicalise (jsonl) + curate (jsonl) of 104m:

```python
main(
    [
        "manifest",
        "--dataset-id",
        "d",
        "--dataset-name",
        "d",
        "--dataset-version",
        "1",
        "--structures-dir",
        str(curated),
        "--canonicalisation-dir",
        str(canonical),
        "--curation-dir",
        str(curated),
        "--exclusions",
        str(curated / "exclusions.jsonl"),
        "--output",
        str(tmp_path / "manifest.json"),
    ]
)
m = json.loads((tmp_path / "manifest.json").read_text())
bundle = m["structures"][0]
assert bundle["entry_id"] == "104M"
assert (
    bundle["canonicalisation"]["policy_id"]
    and bundle["curation"]["policy_id"] == "test"
)
```

(Use a curation policy that excludes a chain or nothing; if `exclusions.jsonl` is not written because nothing was excluded, drop `--exclusions` from this call and add a separate assertion that `--exclusions` reads an `exclusions.jsonl` from Task 5's strict run.)

- [ ] **Step 2: Run, expect FAIL.** `uv run pytest tests/test_cli.py -k manifest -q`
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run, expect PASS.** Then `uv run pytest -q` → all pass.

---

### Task 7: Docs, CLI README, PPI scripts

**Files:** `pandora/cli/README.md`; `docs/usage/cli.md`, `canonicalisation.md`, `datasets.md`, `provenance.md`, `ingestion.md`, `metadata.md`; `docs/recipes/benchmark-quality-filters.md`; `examples/ppi/v2.1/curate.py`, `b_sequence_clusters.py`; architecture table (`read_records`, the two `*_from_record` functions).

- [ ] **Step 1: Rewrite every CLI example** that runs `canonicalise` / `curate` / `manifest` / `dedup` / `reproduce` or reads an old file:
  - add `--table-format jsonl` to `canonicalise` and `curate` calls;
  - `jq '.["ID"].X' dir/canonicalisation_provenance.json` → `jq 'select(.entry_id=="ID") | .X' dir/entries.jsonl`;
  - `jq '.["ID"]…' dir/canonicalisation_mappings.json` → `jq 'select(.entry_id=="ID") | …' dir/mappings.jsonl`;
  - `jq '.[] | …' dir/curation_exclusions.json` → `jq '…' dir/exclusions.jsonl`;
  - `manifest --canonicalisation-provenance F` → `--canonicalisation-dir DIR`; mention `--curation-dir`;
  - prose describing the output files → the step-folder tables from the spec.
- [ ] **Step 2: Regenerate outputs.** Copy the removed docs runner from git (`git show 3afa338^:tests/test_docs_examples.py`) into the scratchpad, add the other pages to `PAGES`, run it with `uv run pytest <scratch>/test_docs_examples.py -q` from the repo root. Expected: every example passes. Where the printed output legitimately changed, paste the real output. Do not add the runner to `tests/` (the user removed it on purpose).
- [ ] **Step 3: PPI scripts:** replace `pd.read_parquet(...)` reads with `read_records` (`curate.py: read_copies`; part B rows); drop the NaN→None workaround. Rerun `curate.py` and part B; expected unchanged 31/53 and 47 chains.
- [ ] **Step 4: Build docs** (`docs-build` skill: strict build, duplicate-anchor and nav checks). `uv run ruff format --check . && uv run ruff check . && uv run pytest -q` → clean / all pass.
