# Modular canonicalisation and curation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give Pandora the pieces that let canonicalisation and curation
run as two standalone steps, each writing mmCIF structures plus its own
provenance, and rebuild PPI v2.1 part A as two thin scripts on top.

**Architecture:** New typed record models in `pandora/schemas/`, small
builder functions in the matching logic packages, all exported through
`__all__`. Writers are the existing `structure_to_mmcif`,
`export_chain_mmcif`, `write_records`, `write_json`, `write_fasta`. The
loop stays in the example scripts.

**Tech Stack:** Python 3, pydantic v2, gemmi, pandas/pyarrow (`export`
extra), pytest, Ruff, uv.

**Spec:** `design/specs/2026-11-04-modular-canonicalisation-curation-design.md`

## Global Constraints

- Start new modules with `from __future__ import annotations`; type hints everywhere; line length 80; `uv run ruff format .`.
- Public functions: Google-style docstring (`Args:`, `Returns:`, `Raises:`), exported in the package `__init__.py` and `__all__`. Helpers start with `_`.
- Models only in `pandora/schemas/`; no behaviour in them.
- Never mutate inputs; new objects via `model_copy(update=...)`.
- No new dependencies. pandas/pyarrow only through `write_records` (lazy import, as today).
- Tests use `datasets/dev/mmcif/` fixtures only; no network, no binaries. One test file per area.
- `examples/` gets no pytest (it is untracked).
- Entry ids are never re-cased. Paths use the id exactly as passed (scripts pass the input file's stem, e.g. `10kz`).
- **Do not commit.** The user asked for no commits for now; leave changes in the working tree.

## Deviations from the spec (decided while planning)

- **Spec item 2 (read a step's structures back) is dropped.** `list_mirror_files()` (Task 1) plus the existing `mmcif_to_structure`, which already returns `None` + diagnostics on failure, is one line in the script. A wrapper adds nothing.
- **Per-entry provenance rows are flat models, not `ProvenanceBundle`s.** Checked: `write_records` → Parquet fails on `canonicalisationProvenance.report == {}` ("Cannot write struct type 'report' with no child field"). The flat rows carry `report` as `warnings` / `errors` ints. Empty `list[str]` columns write fine.
- **Entry rows carry `structure_path`** (both steps), so the next step or user can find the file from the table.

## Review Focus

1. **Empty tables.** `write_records([])` raises `ValueError`. An input list with no symmetry copies, no exclusions, or nothing kept must not crash a script. Owned by Task 7 (scripts guard each write with `if rows:`). Task 4 pins that `chain_copy_records` returns `[]`.
2. **Chain without completeness** (no `_entity_poly_seq`) → its four completeness fields are `None`, not a `KeyError`. Pinned in Task 5.
3. **Metadata without quality** (`metadata.quality is None`) → `experimental_method` / `resolution` are `None`. Pinned in Task 5.
4. **Curating a re-read canonical mmCIF gives the same result as in memory** (the spec's "same result both ways"). Pinned in Task 5.
5. **Ids are not re-cased**, and chain stems with `_` (`1ayi_A_2`) still land under `ay/`. Pinned in Task 1.

---

### Task 1: Mirror-layout paths

**Files:**
- Modify: `pandora/ingestion/paths.py`, `pandora/ingestion/__init__.py`
- Test: `tests/test_ingestion.py`

**Interfaces:**
- Produces:
  - `mirror_path(root: str | Path, stem: str, suffix: str = ".cif") -> Path`: returns `Path(root) / stem[1:3] / f"{stem}{suffix}"`. No case change, no mkdir.
  - `list_mirror_files(root: str | Path, suffix: str = ".cif") -> list[Path]`: sorted `Path(root).glob(f"*/*{suffix}")`. Empty list if `root` does not exist.

- [ ] **Step 1: Write the failing tests**

```python
def test_mirror_path_keeps_case(tmp_path):
    assert mirror_path(tmp_path, "10kz") == tmp_path / "0k" / "10kz.cif"
    assert mirror_path(tmp_path, "10KZ") == tmp_path / "0K" / "10KZ.cif"
    assert mirror_path(tmp_path, "1ayi_A_2") == tmp_path / "ay" / "1ayi_A_2.cif"


def test_list_mirror_files_sorted_and_missing_root(tmp_path):
    for stem in ("1ayi", "10kz"):
        p = mirror_path(tmp_path, stem)
        p.parent.mkdir(parents=True)
        p.write_text("")
    assert list_mirror_files(tmp_path) == sorted(
        [tmp_path / "0k" / "10kz.cif", tmp_path / "ay" / "1ayi.cif"]
    )
    assert list_mirror_files(tmp_path / "nope") == []
```

- [ ] **Step 2: Run, expect FAIL** (`ImportError`): `uv run pytest tests/test_ingestion.py -k mirror -v`
- [ ] **Step 3: Implement both functions in `pandora/ingestion/paths.py`; export them.**
- [ ] **Step 4: Run, expect PASS:** `uv run pytest tests/test_ingestion.py -k mirror -v`

---

### Task 2: Step summary

**Files:**
- Modify: `pandora/schemas/provenance.py`, `pandora/provenance/manifest.py`, `pandora/provenance/__init__.py`
- Test: `tests/test_provenance.py`

**Interfaces:**
- Produces:
  - Model `StepSummary` (in `schemas/provenance.py`): `step: str`, `pandora_version: str`, `gemmi_version: str`, `policy: dict[str, Any]`, `input: str`, `previous_summary: str | None = None`, `output: str`, `started_at: str`, `finished_at: str`, `counts: dict[str, int] = {}`.
  - `build_step_summary(step: str, policy: BaseModel, *, input: str | Path, output: str | Path, started_at: str, counts: Mapping[str, int] = {}, previous_summary: str | Path | None = None) -> StepSummary`: fills `pandora_version` from `pandora.__version__`, `gemmi_version` from `gemmi.__version__`, `finished_at` from `pandora._util.now_iso()`, `policy` from `policy.model_dump(mode="json")`; paths stored as `str`.
  - Written with the existing `write_json(summary, path)`.

- [ ] **Step 1: Write the failing test**

```python
def test_build_step_summary_embeds_policy_and_versions(tmp_path):
    policy = DatasetCurationPolicy(
        policy_id="p", policy_name="n", policy_version="1"
    )
    summary = build_step_summary(
        "curation",
        policy,
        input=tmp_path / "in",
        output=tmp_path / "out",
        started_at="2026-11-04T00:00:00Z",
        counts={"entries_in": 2},
        previous_summary=tmp_path / "canonicalisation.json",
    )
    assert summary.step == "curation"
    assert summary.policy["policy_id"] == "p"
    assert summary.pandora_version == pandora.__version__
    assert summary.gemmi_version == gemmi.__version__
    assert summary.input == str(tmp_path / "in")
    assert summary.previous_summary == str(tmp_path / "canonicalisation.json")
    assert summary.counts == {"entries_in": 2}
    path = write_json(summary, tmp_path / "curation.json")
    assert StepSummary.model_validate_json(path.read_text()) == summary
```

(Construct `DatasetCurationPolicy` with whatever its required fields are; the assertion on `policy_id` is what matters.)

- [ ] **Step 2: Run, expect FAIL:** `uv run pytest tests/test_provenance.py -k step_summary -v`
- [ ] **Step 3: Add the model and `build_step_summary`; export both** (`StepSummary` is importable from `pandora.schemas.provenance`).
- [ ] **Step 4: Run, expect PASS.**

---

### Task 3: Per-entry provenance rows + curation slot in `ProvenanceBundle`

**Files:**
- Modify: `pandora/schemas/provenance.py`, `pandora/provenance/manifest.py`, `pandora/provenance/__init__.py`
- Create: `pandora/provenance/records.py`
- Test: `tests/test_provenance.py`

**Interfaces:**
- Consumes: `canonicalisationProvenance`, `CurationProvenance`, `MetadataRecord`, `Structure`.
- Produces:
  - Model `CanonicalisationEntryRecord`: `entry_id`, `structure_path: str | None = None`, `policy_id`, `policy_name`, `policy_version`, `canonicalised_at`, `transforms: list[str]`, `dropped_raw_categories: list[str]`, `warnings: int | None = None`, `errors: int | None = None`.
  - Model `CurationEntryRecord`: `entry_id`, `structure_path: str | None = None`, `experimental_method: str | None`, `resolution: float | None`, `chain_count: int`, `policy_id`, `policy_name`, `policy_version`, `curated_at`, `dropped_raw_categories: list[str]`.
  - `canonicalisation_entry_record(structure: Structure, provenance: canonicalisationProvenance, structure_path: str | Path | None = None) -> CanonicalisationEntryRecord`: `warnings`/`errors` from `provenance.report.get(...)`.
  - `curation_entry_record(structure: Structure, metadata: MetadataRecord | None, provenance: CurationProvenance, structure_path: str | Path | None = None) -> CurationEntryRecord`: `chain_count` = number of polymer `asym_units` (entity type `"polymer"`).
  - `ProvenanceBundle.curation: CurationProvenance | None = None`; `build_provenance_bundle(..., curation: CurationProvenance | None = None)`.

- [ ] **Step 1: Write the failing tests**

```python
def test_entry_records_write_to_parquet(tmp_path):
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / "ay" / "1ayi.cif"))
    canon, _, cprov = canonicalise_structure(
        structure, canonicalisationPolicy()
    )
    metadata = collect_metadata(canon)
    curated, _, kprov = curate_structure(
        canon, metadata, DatasetCurationPolicy()
    )
    c_row = canonicalisation_entry_record(
        canon, cprov, tmp_path / "ay/1ayi.cif"
    )
    k_row = curation_entry_record(curated, metadata, kprov)
    assert c_row.entry_id == canon.entry_id
    assert c_row.structure_path == str(tmp_path / "ay/1ayi.cif")
    assert c_row.policy_id == cprov.policy_id
    assert k_row.chain_count == len(extract_chain_records(curated))
    assert k_row.resolution == metadata.quality.resolution
    write_records(
        [c_row], tmp_path / "c.parquet"
    )  # report={} must not break it
    write_records([k_row], tmp_path / "k.parquet")


def test_provenance_bundle_carries_curation():
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / "ay" / "1ayi.cif"))
    _, _, kprov = curate_structure(
        structure, collect_metadata(structure), DatasetCurationPolicy()
    )
    assert build_provenance_bundle(structure, curation=kprov).curation == kprov
```

(If `DatasetCurationPolicy()` needs arguments or drops 1ayi, use the `_policy()` helper style from `tests/test_curation.py` with permissive settings; the entry must be kept.)

- [ ] **Step 2: Run, expect FAIL:** `uv run pytest tests/test_provenance.py -k "entry_records or carries_curation" -v`
- [ ] **Step 3: Add the models, `records.py` functions, the bundle field and parameter; export.**
- [ ] **Step 4: Run, expect PASS.**

---

### Task 4: Chain copy records

**Files:**
- Modify: `pandora/schemas/canonicalisation.py`, `pandora/canonicalisation/assemblies.py`, `pandora/canonicalisation/__init__.py`
- Test: `tests/test_canonicalisation.py`

**Interfaces:**
- Consumes: `CanonicalMappings.assembly_mapping.items[*].chain_copies` (`AssemblyChainCopy`: `canonical_chain_id`, `source_chain_id`, `operator_id`).
- Produces:
  - Model `ChainCopyRecord`: `entry_id: str`, `chain_id: str`, `source_chain_id: str`, `operator_id: str`.
  - `chain_copy_records(entry_id: str, mappings: CanonicalMappings) -> list[ChainCopyRecord]`: one per copy, in mapping order; `[]` when none.

- [ ] **Step 1: Write the failing tests**

```python
def test_chain_copy_records_from_expanded_assembly():
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / "ay" / "1ayi.cif"))
    policy = canonicalisationPolicy(
        assembly_rules={
            "strategy": "standardize_biological_assembly",
            "preferred_assembly_source": "first",
        }
    )
    canon, mappings, _ = canonicalise_structure(structure, policy)
    rows = chain_copy_records(canon.entry_id, mappings)
    assert [(r.chain_id, r.source_chain_id) for r in rows] == [("A_2", "A")]
    assert {r.chain_id for r in rows} <= {u.id for u in canon.asym_units}


def test_chain_copy_records_empty_without_copies():
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / "ay" / "1ayi.cif"))
    _, mappings, _ = canonicalise_structure(structure, canonicalisationPolicy())
    assert chain_copy_records(structure.entry_id, mappings) == []
```

(1ayi's assembly 1 yields one copy `A_2` under the PPI v2.1 canonicalisation policy, as seen in part A's output. If the exact list differs, assert against what `mappings` holds, but keep the `A_2 → A` pair.)

- [ ] **Step 2: Run, expect FAIL:** `uv run pytest tests/test_canonicalisation.py -k chain_copy_records -v`
- [ ] **Step 3: Implement and export.**
- [ ] **Step 4: Run, expect PASS.**

---

### Task 5: Dataset chain records

**Files:**
- Modify: `pandora/schemas/dataset.py`, `pandora/datasets/records.py`, `pandora/datasets/__init__.py`
- Test: `tests/test_curation.py`

**Interfaces:**
- Consumes: `extract_chain_records`, `chain_completeness`, `MetadataRecord.quality`.
- Produces:
  - Model `DatasetChainRecord(ChainRecord)`: adds `seqres_length: int | None = None`, `missing_n_term`, `missing_c_term`, `missing_middle` (same type), `experimental_method: str | None = None`, `resolution: float | None = None`, `structure_path: str | None = None`, `copy_of: str | None = None`.
  - `dataset_chain_records(structure: Structure, metadata: MetadataRecord | None, definition: MissingResidueDefinition = "incomplete_backbone", copy_of: Mapping[str, str] = {}, structure_path: str | Path | None = None) -> list[DatasetChainRecord]`: one per `extract_chain_records(structure)` row, in that order; `copy_of` maps chain_id → source chain id. Completeness diagnostics are dropped (as `extract_chain_records` reports none either); missing completeness gives `None`s.

- [ ] **Step 1: Write the failing tests**

```python
def test_dataset_chain_records_join_completeness_quality_copies(tmp_path):
    s = _load("1ayi")
    rows = dataset_chain_records(
        s,
        collect_metadata(s),
        copy_of={"B": "A"},
        structure_path=tmp_path / "x.cif",
    )
    by_id = {r.chain_id: r for r in rows}
    assert [r.chain_id for r in rows] == [
        r.chain_id for r in extract_chain_records(s)
    ]
    assert by_id["A"].seqres_length is not None
    assert by_id["A"].copy_of is None and by_id["B"].copy_of == "A"
    assert by_id["A"].structure_path == str(tmp_path / "x.cif")
    assert by_id["A"].resolution == collect_metadata(s).quality.resolution


def test_dataset_chain_records_without_seqres_or_quality():
    s = _load("1ayi")
    no_seqres = s.model_copy(
        update={
            "raw": {k: v for k, v in s.raw.items() if k != "_entity_poly_seq"}
        }
    )
    rows = dataset_chain_records(no_seqres, None)
    assert rows and all(
        r.seqres_length is None
        and r.missing_middle is None
        and r.resolution is None
        and r.experimental_method is None
        for r in rows
    )


def test_curating_reread_canonical_mmcif_matches_in_memory(tmp_path):
    s = _load("1ayi")
    canon, _, _ = canonicalise_structure(
        s,
        canonicalisationPolicy(
            assembly_rules={
                "strategy": "standardize_biological_assembly",
                "preferred_assembly_source": "first",
            }
        ),
    )
    reread, _, _ = mmcif_to_structure(
        str(structure_to_mmcif(canon, tmp_path / "1ayi.cif"))
    )
    policy = _policy()
    a = curate_structure(canon, collect_metadata(canon), policy)
    b = curate_structure(reread, collect_metadata(reread), policy)
    assert [(e.chain_id, e.reason_code) for e in a[1]] == [
        (e.chain_id, e.reason_code) for e in b[1]
    ]
    assert dataset_chain_records(
        a[0], collect_metadata(canon)
    ) == dataset_chain_records(b[0], collect_metadata(reread))
```

(1ayi has chains A and B; `copy_of={"B": "A"}` is a synthetic mapping, used only to check the join.)

- [ ] **Step 2: Run, expect FAIL:** `uv run pytest tests/test_curation.py -k "dataset_chain_records or reread" -v`
- [ ] **Step 3: Implement and export.**
- [ ] **Step 4: Run, expect PASS.** If the re-read test fails on a field that legitimately differs after mmCIF (e.g. `auth_chain_id` of copies), stop and report; don't loosen the assertion silently.

---

### Task 6: Docs, diagrams, full check

**Files:**
- Modify: `docs/contributing/architecture.md` (function → model table: rows for `mirror_path`, `list_mirror_files`, `build_step_summary`, `canonicalisation_entry_record`, `curation_entry_record`, `chain_copy_records`, `dataset_chain_records`; add `curation` to `build_provenance_bundle`'s inputs)
- Regenerate: ERD via `uv run --extra docs python docs/scripts/generate_erd.py`

- [ ] **Step 1: Update the table and regenerate the ERD.**
- [ ] **Step 2: Build the docs** with the `docs-build` skill (`uv run zensical build --clean`). Expected: builds, no new broken links.
- [ ] **Step 3: Full check.** `uv run pytest` → all pass (apart from existing strict xfails); `uv run ruff format --check .` and `uv run ruff check .` → clean.

---

### Task 7: PPI v2.1 scripts

No pytest (untracked examples). Verified by running on the 100 fixtures.

**Files:**
- Create: `examples/ppi/v2.1/canonicalise.py`, `examples/ppi/v2.1/curate.py`
- Delete: `examples/ppi/v2.1/a_prepare_entries.py`
- Modify: `examples/ppi/v2.1/b_sequence_clusters.py`, `examples/ppi/v2.1/c_structure_similarity.py`, `examples/ppi/v2.1/README.md`

**Interfaces:**
- Consumes: everything from Tasks 1–5, plus `read_entry_paths`, `load_policy`, `mmcif_to_structure`, `canonicalise_structure`, `collect_metadata`, `curate_structure`, `structure_to_mmcif`, `export_chain_mmcif`, `write_records`, `write_json`, `write_fasta`.

**`canonicalise.py`** (settings at top like today: `ENTRIES_FILE`, `OUT = HERE / "output/canonicalisation"`, `CANONICALISATION_POLICY`):
- `shutil.rmtree(OUT, ignore_errors=True)` first.
- Per path: `stem = path.stem`; parse (None → failure row), canonicalise, write `mirror_path(OUT / "structures", stem)`, collect a `CanonicalisationEntryRecord` (with that path) and `chain_copy_records`. Catch `(ValueError, OSError, RuntimeError, KeyError, IndexError)` per entry → `failures.jsonl` row (`path`, `error_type`, `message`), written while the loop runs.
- After: `entries.parquet`, `chain_copies.parquet` (each only `if rows:`), `canonicalisation.json` via `build_step_summary("canonicalisation", ...)` with counts `entries_in`, `entries_kept`, `entries_failed`.

**`curate.py`** (`IN = HERE / "output/canonicalisation"`, `OUT = HERE / "output/curation"`, `CURATION_POLICY`):
- `rmtree(OUT)`. Read `IN / "chain_copies.parquet"` if present → `{(entry_id, chain_id): source_chain_id}`.
- Per `list_mirror_files(IN / "structures")`: parse, `collect_metadata`, `curate_structure`. Write exclusions. If kept: write `mirror_path(OUT / "structures", stem)`, build `dataset_chain_records(curated, metadata, curation.quality_rules.missing_residue_definition, copy_of=..., structure_path=...)`, a `curation_entry_record`. For non-copy chains: `export_chain_mmcif(curated, chain_id, mirror_path(OUT / "chains", f"{stem}_{chain_id}"))` and a FASTA entry `f"{stem}_{chain_id}"`. Same failure handling as above.
- After: `entries.parquet`, `chains.parquet`, `exclusions.parquet` (each `if rows:`), `chains.fasta`, `curation.json` with `previous_summary=IN / "canonicalisation.json"` and counts `entries_in`, `entries_kept`, `entries_failed`, `chains_kept`, and `excluded_entry:<REASON>` / `excluded_chain:<REASON>`.

**Parts B and C:**
- B reads `output/curation/chains.parquet` (`pd.read_parquet(...).to_dict("records")`), skips rows with `copy_of` set, and builds chain ids as `f"{Path(r['structure_path']).stem}_{r['chain_id']}"` so they match the FASTA ids and file names. `resolution` and missing counts are unchanged.
- C builds chain file paths with `mirror_path(a_out / "chains", chain)` instead of `chain[1:3]`, with `a_out = output/curation`.

- [ ] **Step 1: Write `canonicalise.py`; run** `cd examples/ppi/v2.1 && uv run python canonicalise.py`. Expected: `output/canonicalisation/` has `structures/0k/10kz.cif` (lowercase), `entries.parquet`, `chain_copies.parquet` (includes `1AYI`/`A_2`), `canonicalisation.json`; `22AG` appears in `failures.jsonl` (known, spec side finding).
- [ ] **Step 2: Write `curate.py`; run** `uv run python curate.py`. Expected: 31 kept entries and 53 rows in `chains.parquet`, as in the old `output/a/counts.tsv` (22AG was excluded on resolution there, so its new canonicalisation failure doesn't change the kept numbers; only the `RESOLUTION_THRESHOLD` entry count drops by one). Report any other difference. `chains.fasta` has 53 − (rows with `copy_of`) records; no `_2`-style ids in it or in `chains/`.
- [ ] **Step 3: Update B and C; run B and C if MMseqs2/Foldseek are on `PATH`**, otherwise say they were not run.
- [ ] **Step 4: Delete `a_prepare_entries.py`; rewrite README part A** as two sections (Canonicalise, Curate) with the output tables from the spec. Leave the old `output/a/` folder alone and mention it to the user.
- [ ] **Step 5: `uv run ruff format examples/ppi/v2.1 && uv run ruff check examples/ppi/v2.1`** → clean.
