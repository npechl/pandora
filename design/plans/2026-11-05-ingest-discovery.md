# CLI `ingest` finds input files — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `pandora ingest` lists the input files of any archive layout into `entries.txt` + `ingestion.json`, without touching individual files. `canonicalise --entries` reads that list and links back to its summary.

**Architecture:** Two functions in `pandora/ingestion/paths.py` (`template_entry_paths`, `write_entry_paths`) and one model in `schemas/ingestion.py` (`IngestionDiscovery`). A rewritten `_cmd_ingest`, a `previous_summary` hook in `_cmd_canonicalise`, and `entries.txt` from `fetch`. The `--glob` option added earlier to `canonicalise`/`curate` is reverted.

**Tech Stack:** Python 3 stdlib (`pathlib`, `os`, `re`), pydantic v2, argparse, pytest.

**Spec:** `design/specs/2026-11-05-ingest-discovery-design.md`

## Global Constraints

- CLAUDE.md coding standards (future import, type hints, 80 cols, Google docstrings on public functions, `__all__`, models in `schemas/`, no new dependencies).
- `ingest` never opens, parses or stats a file below the ID segment, and never checks that a listed path exists.
- `ingest` writes exactly `entries.txt` (absolute paths, sorted, one per line, no header) and `ingestion.json` (`StepSummary`, `step="ingestion"`, `counts={"entries_found": n}`). No `failures.jsonl`, no `ingestion_provenance.json`.
- Placeholders: `{id}`, `{mid}` = `id[1:3]`. With `--ids`, IDs are lowercased; when the IDs come from listing, they keep the case found on disk.
- `manifest`, `reproduce`, `IngestionProvenance`, `ingest_local_mmcif` stay unchanged.
- Printed line: `ingested {n} entries -> {output_dir}`.
- **No commits** unless the user asks.

## Review Focus

1. **Stray files next to the target** (`entry/04/104m/other.cif`, `entry/04/104m/clean_mmcif/104m.cif`): the template must list only paths it builds. Pinned in Task 2.
2. **Files sitting at a pre-ID level** (`entry/README` next to the `<mid>` folders): skipped, not an error. Pinned in Task 2.
3. **A template matching nothing**: both files are written, count 0, exit 0. Pinned in Task 3.
4. **Relative `--input-dir`**: `entries.txt` holds absolute paths that `canonicalise --entries` resolves from another folder. Pinned in Task 3.
5. **`canonicalise --overwrite` on a folder that holds an `entries.txt`**: the list survives. Pinned in Task 3.

---

### Task 1: Revert `--glob` on `canonicalise` / `curate`

Keep the uncommitted entry-ID naming (`stem = entry_id.lower()`) and the `DuplicateEntry` check. Remove only the glob option.

**Files:**
- Modify: `pandora/cli/_io.py` (`structure_files` loses `pattern`; delete `GLOB_HELP`)
- Modify: `pandora/cli/canonicalisation.py`, `pandora/cli/datasets.py` (drop `--glob`, the `args.glob` calls and the `GLOB_HELP` import)
- Modify: `docs/usage/cli.md`, `pandora/cli/README.md` (drop the `--glob` sentence and keep the entry-ID naming sentence; mention `pandora ingest` for other layouts)
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: `structure_files(folder: Path) -> list[Path]`, unchanged from `main`.

- [ ] **Step 1:** In `tests/test_cli.py`, replace `test_step_glob_reads_nested_gzipped_layout` with `test_step_names_outputs_by_entry_id`: copy `04/104m.cif` gzipped to `raw/104m_updated.cif.gz`, write `entries.txt` with that one line, run `canonicalise --entries`, then assert that `canon/structures/04/104m.cif` exists and `entries.jsonl` has `entry_id == "104M"`.
- [ ] **Step 2:** Apply the code and doc changes above.
- [ ] **Step 3:** `uv run pytest tests/test_cli.py -q`. Expected: all pass, and `grep -rn "\-\-glob\|GLOB_HELP" pandora/cli` returns nothing.

### Task 2: Library: `template_entry_paths`, `write_entry_paths`, `IngestionDiscovery`

**Files:**
- Modify: `pandora/ingestion/paths.py`, `pandora/ingestion/__init__.py` (`__all__`)
- Modify: `pandora/schemas/ingestion.py`
- Modify: `docs/contributing/architecture.md` (two rows next to `read_entry_paths`)
- Test: `tests/test_ingestion.py`

**Interfaces:**
- Produces:
  - `template_entry_paths(root: str | Path, template: str, ids: list[str] | None = None) -> list[Path]`, sorted, each `Path(root) / template.format(id=..., mid=...)`. Raises `ValueError` if no segment contains `{id}` or the ID segment also contains `{mid}`.
  - `write_entry_paths(paths: Iterable[Path], path: str | Path) -> None`: one path per line with a trailing newline, parent folders created.
  - `IngestionDiscovery(BaseModel)` with `input_dir: str`, `glob: str | None = None`, `path_template: str | None = None`, `ids_file: str | None = None`, `source_label: str | None = None`, `release_date: str | None = None`. No behaviour.

Algorithm without `ids` (the part the signature doesn't fix):

```python
parts = template.split("/")
k = next(i for i, p in enumerate(parts) if "{id}" in p)  # else ValueError
dirs = [Path(root)]
for part in parts[:k]:  # pre-ID levels
    pattern = part.replace("{mid}", "*")
    dirs = [c for d in dirs for c in d.glob(pattern) if c.is_dir()]
id_re = re.compile(re.escape(parts[k]).replace(re.escape("{id}"), "(.+)"))
ids = [
    m.group(1)
    for d in dirs
    for name in os.listdir(d)
    if (m := id_re.fullmatch(name))
]
```

- [ ] **Step 1: Write the failing tests** in `tests/test_ingestion.py`:
  - `test_template_entry_paths_lists_only_down_to_the_id`: in `tmp_path/entry`, create `04/104m/clean_mmcif/104m_updated.cif.gz`, `0k/10kz/clean_mmcif/10kz_updated.cif.gz`, a stray `04/104m/other.cif`, a stray `04/104m/clean_mmcif/104m.cif`, and a file `README` in `entry/`. `template_entry_paths(root, "{mid}/{id}/clean_mmcif/{id}_updated.cif.gz") == [root/"04/104m/clean_mmcif/104m_updated.cif.gz", root/"0k/10kz/clean_mmcif/10kz_updated.cif.gz"]`.
  - `test_template_entry_paths_without_listing_when_ids_given`: `root = tmp_path / "missing"` (doesn't exist); `template_entry_paths(root, "{mid}/{id}/x/{id}.cif.gz", ["100D"]) == [root/"00/100d/x/100d.cif.gz"]`.
  - `test_template_entry_paths_on_plain_mirror`: on `FIXTURES_DIR`, `{mid}/{id}.cif` returns exactly `sorted(FIXTURES_DIR.glob("*/*.cif"))`.
  - `test_template_entry_paths_keeps_listed_case`: folder `AB/1ABC/` with template `{mid}/{id}` returns `[root/"AB/1ABC"]`.
  - `test_template_entry_paths_rejects_bad_templates`: `pytest.raises(ValueError)` for `"{mid}/x.cif"` and for `"{mid}{id}.cif"`.
  - `test_write_entry_paths_round_trips`: absolute paths written to `tmp_path/"sub"/"entries.txt"` come back unchanged from `read_entry_paths`.
- [ ] **Step 2:** `uv run pytest tests/test_ingestion.py -q`. Expected: FAIL (ImportError).
- [ ] **Step 3:** Implement the functions and the model.
- [ ] **Step 4:** `uv run pytest tests/test_ingestion.py -q`. Expected: PASS.
- [ ] **Step 5:** Regenerate the diagrams: `uv run --extra docs python docs/scripts/generate_erd.py`.

### Task 3: CLI `ingest`, and `canonicalise` linking to it

**Files:**
- Modify: `pandora/cli/ingestion.py` (`_cmd_ingest` and its parser)
- Modify: `pandora/cli/_io.py`: `prepare_step_dir(out, summary, overwrite, files: tuple[str, ...] = _STEP_FILES)`. `ingest` passes `files=("entries.txt",)`. Every other caller keeps the default, so `canonicalise --overwrite` never touches `entries.txt`.
- Modify: `pandora/cli/canonicalisation.py`: with `--entries FILE`, `previous_summary = FILE.parent / "ingestion.json"` if it exists, else `None`.
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `template_entry_paths`, `write_entry_paths`, `read_entry_paths`, `IngestionDiscovery`, `build_step_summary`, `write_json`, `prepare_step_dir`.
- Produces: `pandora ingest --input-dir D (--glob P | --path-template T) [--ids F] [--source-label L] [--release-date R] --output-dir O [--overwrite]`. `--glob` and `--path-template` are a required mutually exclusive group; `--ids` without `--path-template` causes `SystemExit` with a message. A `ValueError` from the template also becomes `SystemExit`. `--input-dir` is `.resolve()`d first; `--glob` is `sorted(root.glob(P))`.

- [ ] **Step 1: Write the failing tests** in `tests/test_cli.py`. Replace the old `test_ingest_local_directory`. Use one helper that builds the PDBe-like tree from Task 2 with real fixture bytes, gzipped.
  - `test_ingest_path_template_feeds_canonicalise`: `monkeypatch.chdir(tmp_path)` and pass `--input-dir entry` (relative), `--path-template '{mid}/{id}/clean_mmcif/{id}_updated.cif.gz'`, `--source-label lbl`, `--output-dir ing`. Assert that `ing/entries.txt` holds 2 absolute paths, and that `ingestion.json` has `step == "ingestion"`, `counts == {"entries_found": 2}` and `policy["path_template"]` / `policy["source_label"] == "lbl"`. Then `canonicalise --entries ing/entries.txt --output-dir canon` gives `canon/structures/04/104m.cif`, and `canonicalisation.json` has `previous_summary` ending with `ing/ingestion.json`.
  - `test_ingest_glob`: `--glob '*/*.cif'` on `FIXTURES_DIR` lists the same paths as `sorted(FIXTURES_DIR.glob("*/*.cif"))`.
  - `test_ingest_ids_file`: `--ids` with lines `# c`, empty, `104M` gives one path ending with `04/104m/clean_mmcif/104m_updated.cif.gz`.
  - `test_ingest_nothing_found`: a template that matches nothing gives exit 0, an empty `entries.txt` and `entries_found == 0`.
  - `test_ingest_rejects_ids_without_template`: `pytest.raises(SystemExit)`.
  - `test_canonicalise_overwrite_keeps_entries_txt`: run `ingest` and `canonicalise` with the same `--output-dir`, then `canonicalise --overwrite`. `entries.txt` still exists.
- [ ] **Step 2:** `uv run pytest tests/test_cli.py -q -k "ingest or overwrite_keeps"`. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** `uv run pytest -q`. Expected: all pass.

### Task 4: `fetch` writes `entries.txt`

**Files:**
- Modify: `pandora/cli/ingestion.py` (`_cmd_fetch`)
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `write_entry_paths`.
- Produces: `<output-dir>/entries.txt`, one line per fetched entry: `output_dir/<id>.cif.gz` if it exists, else `<id>.cif`. Absolute paths, sorted. `ingestion_provenance.json` is unchanged.

- [ ] **Step 1: Write the failing test** `test_fetch_writes_entries_txt`: monkeypatch `pandora.ingestion.mmcif._download` to return the bytes of `04/104m.cif`, then run `fetch 104m --provider pdb --no-cache --output-dir raw`. `read_entry_paths(raw/"entries.txt") == [(raw/"104m.cif").resolve()]`.
- [ ] **Step 2:** Run it. Expected: FAIL (no `entries.txt`).
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** `uv run pytest -q`. Expected: all pass.

### Task 5: Docs

**Files:**
- Modify: `docs/usage/ingestion.md`: rewrite "Ingest already-downloaded files" with a `library` tab (`template_entry_paths` + `write_entry_paths` on `datasets/dev/mmcif/` with `{mid}/{id}.cif`) and a `cli` tab (`pandora ingest ... --path-template '{mid}/{id}.cif'`, then `jq .counts`). Show the real outputs from running them. Add the PDBe release template as an example that isn't run, with one line on why the template is fast on NFS. In the `fetch` section, mention `entries.txt`.
- Modify: `docs/usage/cli.md`, `pandora/cli/README.md`: `ingest` writes `entries.txt` + `ingestion.json`, and `canonicalise --entries` links to it.

- [ ] **Step 1:** Run the two examples from the repo root into a scratch folder and paste their real output.
- [ ] **Step 2:** `uv run zensical build --clean`. Expected: `No issues found`. Grep the rendered `site/usage/ingestion/index.html` for `template_entry_paths` (the code itself, not only the output).
- [ ] **Step 3:** Final check: `uv run pytest -q`, `uv run ruff format --check .` and `uv run ruff check .` all pass.
