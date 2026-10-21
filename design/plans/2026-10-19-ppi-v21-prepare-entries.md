# PPI v2.1 part A: prepare entries — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A shardable script that turns PDB entry ids into filtered, protein-only assembly-1 chains (records, FASTA, per-chain and per-entry mmCIF) plus per-filter drop counts.

**Architecture:** One library fix (`preferred_assembly_source: first` picks assembly 1). Everything else is an example script, `examples/ppi/v2.1/a_prepare_entries.py`, that chains existing Pandora functions per entry, writes per-shard files, and merges shards. Thresholds live in two policy YAMLs next to it.

**Tech Stack:** Python 3.13, pydantic v2, gemmi, PyYAML, argparse, pytest, uv, Slurm (example job only).

**Spec:** `design/specs/2026-10-19-ppi-v21-prepare-entries-design.md`

## Global Constraints

- **Do not commit.** The user asked for no commits; skip every commit step and leave changes in the working tree.
- `from __future__ import annotations` at the top of every new module; type hints on all new code; line length 80; `uv run ruff format --check .` and `uv run ruff check .` pass.
- Catch named exceptions only (Ruff BLE001, S110).
- Never mutate inputs; derive new structures with `model_copy(update=...)`.
- Tests use `datasets/dev/mmcif/` only: no network, no MMseqs2/Foldseek.
- Every v2.1 threshold lives in a policy YAML, never in the script.
- No new dependencies.

## Review Focus

1. **Multi-line SEQRES.** `ChainRecord.sequence` comes from `pdbx_seq_one_letter_code_can`, which can contain newlines. A FASTA record with a newline in the sequence breaks MMseqs2. Strip all whitespace before writing. Pinned in Task 2 (`test_fasta_sequences_have_no_whitespace`).
2. **Entry ids in mixed case or with blank lines** in `--entries` (`1A3N`, a trailing newline). Expected: lower-cased, blanks skipped, duplicates collapsed. Pinned in Task 3 (`test_entry_list_normalised`).
3. **Assembly copies** (`A_2`, …) after expansion. Expected: each copy is its own chain in every output, with its own mmCIF file and FASTA record. Pinned in Task 2 via 1a3n's chain-id assertion (the test reads the actual ids; if 1a3n expands to copies, the assertion lists them).
4. **Kill mid-entry.** Expected: no partial rows. Rows go to disk only after the whole entry finished. Pinned in Task 3 (`test_restart_skips_done_and_failed`).
5. **An entry that fails to parse** (corrupt or truncated file). Expected: one `failures.jsonl` row, and the shard carries on. Pinned in Task 3 (`test_unparseable_entry_is_recorded`).

---

## File Structure

| Path | Responsibility |
|---|---|
| `pandora/canonicalisation/assemblies.py` (modify) | `first` returns the first listed assembly |
| `docs/reference/policies.md` (modify) | Describe the new `first` behaviour |
| `examples/ppi/v2.1/policies/canonicalisation.yaml` | Canonicalisation policy for part A |
| `examples/ppi/v2.1/policies/curation.yaml` | The v2.1 Step 1 thresholds |
| `examples/ppi/v2.1/a_prepare_entries.py` | Per-entry flow, shard runner, merge, CLI |
| `examples/ppi/v2.1/slurm/a_prepare_entries.sh` | Example Slurm array job |
| `examples/ppi/v2.1/README.md` | How to run locally and on Codon |
| `tests/test_canonicalisation.py` (modify) | Test for `first` |
| `tests/test_ppi_v21.py` | Tests for the script |

Fixture facts, measured with the Task 2 policies (2026-10-19):

| Fixture | Outcome |
|---|---|
| `1a3n` | kept; protein chains `A B C D` (entities 1, 2) |
| `12hh` | kept; chain `A` protein; `B`, `C` DNA → `NOT_PROTEIN` |
| `1a3e` | kept `A`, `B`; chain `C` → `MISSING_TAILS` |
| `1aci` | entry → `METHOD_EXCLUDED` (NMR) |
| `13qi` | entry → `RESOLUTION_THRESHOLD` (2.82 Å) |

---

### Task 1: `preferred_assembly_source: first` picks assembly 1

**Files:**
- Modify: `pandora/canonicalisation/assemblies.py:148-176` (`_select_biological_assembly`)
- Modify: `docs/reference/policies.md:130` (the `first` row)
- Test: `tests/test_canonicalisation.py`

**Interfaces:**
- Produces: with `preferred_assembly_source="first"`, `_select_biological_assembly` returns `assemblies[0]` whatever its flags. `author`, `pdbx` and `pdbe` behave as before.

- [ ] **Step 1: Write the failing test** in `tests/test_canonicalisation.py`, next to the existing `pdbe` test (~line 930), importing `_select_biological_assembly` as that file already imports private helpers:

```python
def test_first_assembly_source_ignores_author_flag():
    assemblies = [
        AssemblyRecord(id="1", software_determined=True),
        AssemblyRecord(id="2", author_determined=True),
    ]
    picked = _select_biological_assembly(
        assemblies, "first", DiagnosticBundle(), "test"
    )
    assert picked.id == "1"


def test_author_assembly_source_still_prefers_author_flag():
    assemblies = [
        AssemblyRecord(id="1", software_determined=True),
        AssemblyRecord(id="2", author_determined=True),
    ]
    picked = _select_biological_assembly(
        assemblies, "author", DiagnosticBundle(), "test"
    )
    assert picked.id == "2"
```

- [ ] **Step 2: Run it and check it fails**

Run: `uv run pytest tests/test_canonicalisation.py -k assembly_source -v`
Expected: `test_first_assembly_source_ignores_author_flag` FAILS (`'2' == '1'`), and the author test passes.

- [ ] **Step 3: Implement.** In `_select_biological_assembly`, return `assemblies[0]` straight after the empty check when `preferred_assembly_source == "first"`.

- [ ] **Step 4: Update the docs.** Replace the `first` row in `docs/reference/policies.md` with: "`first` | Use the first listed assembly (assembly 1 in PDB files), ignoring the author/software flags." Then run the `docs-build` skill (`uv run zensical build --clean`) and confirm the build succeeds.

- [ ] **Step 5: Run the area tests**

Run: `uv run pytest tests/test_canonicalisation.py -q`
Expected: all pass.

---

### Task 2: Policies and the per-entry flow

**Files:**
- Create: `examples/ppi/v2.1/policies/canonicalisation.yaml`, `examples/ppi/v2.1/policies/curation.yaml`
- Create: `examples/ppi/v2.1/a_prepare_entries.py` (per-entry part)
- Test: `tests/test_ppi_v21.py`

**Interfaces:**
- Consumes: Task 1's `first`.
- Produces (in `a_prepare_entries.py`):
  - `PROTEIN = "polypeptide(L)"`
  - `load_policies(canon_path: Path, curation_path: Path) -> tuple[canonicalisationPolicy, DatasetCurationPolicy]`. The canonicalisation policy comes from `pandora.ingestion.load_policy`; the curation policy is `DatasetCurationPolicy.model_validate(yaml.safe_load(...))`.
  - `@dataclass(frozen=True) class EntryResult` with fields `chains: list[dict]`, `exclusions: list[dict]`, `structure: Structure | None` (curated and protein-only, or `None` when nothing is left).
  - `prepare_entry(path: Path, canon: canonicalisationPolicy, curation: DatasetCurationPolicy) -> EntryResult`. Raises `ParseError` (defined in the script, subclass of `ValueError`) when `mmcif_to_structure` returns `None`.
  - `fasta_record(row: dict) -> str` → `">{entry_id}_{chain_id}\n{sequence}\n"`, with every whitespace character removed from the sequence.
  - `write_entry_structures(result: EntryResult, out: Path) -> None` writes `chains/<mid>/<entry>_<chain>.cif` with `export_chain_mmcif` and `biounits/<mid>/<entry>.cif`, where `mid = entry_id[1:3]`.

- [ ] **Step 1: Write the policy YAMLs.** `canonicalisation.yaml` sets:
  - `assembly_rules`: `strategy: standardize_biological_assembly`, `preferred_assembly_source: first`;
  - `identifier_rules`: chain id `preserve`, residue numbering `preserve`;
  - `entity_rules`: `strategy: preserve`;
  - `modified_residue_rules`: `strategy: map_to_parent`;
  - `ligand_rules`: `strategy: filter` with `keep_waters`, `keep_ions` and `keep_nonpolymer_ligands` all `false`;
  - `altloc_rules`: `strategy: select_best_occupancy`.

  `curation.yaml` sets `quality_rules`:
  - `max_resolution_by_method`: `{"X-RAY DIFFRACTION": 2.5, "ELECTRON MICROSCOPY": 2.0}`;
  - `include_experimental_methods`: `[X-RAY DIFFRACTION, ELECTRON MICROSCOPY]`;
  - `max_r_free_gap: 0.07`, `max_r_free: 0.35`, `max_r_sym: 0.10`;
  - `exclude_nonstandard_residues: true`;
  - `missing_residue_definition: incomplete_backbone`, `max_missing_tail_fraction: 0.30`, `max_missing_middle_fraction: 0.10`;
  - `max_atoms: 1000000`.

  Each value gets a comment naming its v2.1 filter number.

- [ ] **Step 2: Write the failing tests** in `tests/test_ppi_v21.py`. Load the script with `importlib.util.spec_from_file_location("a_prepare_entries", ROOT / "examples/ppi/v2.1/a_prepare_entries.py")`, cached in a module-level helper `_mod()`. Then add `FIX = ROOT / "datasets/dev/mmcif"` and a `policies` fixture that returns `load_policies(...)` on the two YAMLs.

```python
def test_kept_entry_has_protein_chains(policies):
    r = _mod().prepare_entry(FIX / "1a3n.cif", *policies)
    assert [c["chain_id"] for c in r.chains] == ["A", "B", "C", "D"]
    assert r.exclusions == []
    row = r.chains[0]
    for key in (
        "sequence",
        "seqres_length",
        "missing_n_term",
        "missing_c_term",
        "missing_middle",
        "method",
        "resolution",
    ):
        assert key in row


def test_non_protein_chains_dropped(policies):
    r = _mod().prepare_entry(FIX / "12hh.cif", *policies)
    assert [c["chain_id"] for c in r.chains] == ["A"]
    assert {(e["chain_id"], e["reason_code"]) for e in r.exclusions} == {
        ("B", "NOT_PROTEIN"),
        ("C", "NOT_PROTEIN"),
    }
    assert {a.label_asym_id for a in r.structure.atoms} == {"A"}


def test_chain_filter_reason(policies):
    r = _mod().prepare_entry(FIX / "1a3e.cif", *policies)
    assert [c["chain_id"] for c in r.chains] == ["A", "B"]
    assert ("C", "MISSING_TAILS") in {
        (e["chain_id"], e["reason_code"]) for e in r.exclusions
    }


@pytest.mark.parametrize(
    "entry,code",
    [("1aci", "METHOD_EXCLUDED"), ("13qi", "RESOLUTION_THRESHOLD")],
)
def test_entry_filter_reason(policies, entry, code):
    r = _mod().prepare_entry(FIX / f"{entry}.cif", *policies)
    assert r.chains == [] and r.structure is None
    assert [(e["chain_id"], e["reason_code"]) for e in r.exclusions] == [
        (None, code)
    ]


def test_fasta_sequences_have_no_whitespace():
    rec = _mod().fasta_record(
        {"entry_id": "1abc", "chain_id": "A", "sequence": "MK\nTA Y\n"}
    )
    assert rec == ">1abc_A\nMKTAY\n"


def test_structure_files_match_chains(policies, tmp_path):
    m = _mod()
    r = m.prepare_entry(FIX / "12hh.cif", *policies)
    m.write_entry_structures(r, tmp_path)
    assert (tmp_path / "chains/2h/12hh_A.cif").exists()
    s, _, _ = mmcif_to_structure(str(tmp_path / "biounits/2h/12hh.cif"))
    assert {a.label_asym_id for a in s.atoms} == {"A"}
```

- [ ] **Step 3: Run the tests and check they fail**

Run: `uv run pytest tests/test_ppi_v21.py -v`
Expected: FAIL, because the script doesn't exist yet.

- [ ] **Step 4: Implement the per-entry part of `a_prepare_entries.py`**
  - Module docstring: what part A does, with a pointer to the spec.
  - `prepare_entry` runs, in order:
    1. `mmcif_to_structure`;
    2. `canonicalise_structure`;
    3. `collect_metadata`;
    4. `curate_structure`, mapping each `ExclusionRecord` to a dict with `model_dump()`;
    5. protein ids: the asym units whose entity's `poly.type == PROTEIN`;
    6. one `NOT_PROTEIN` row (`entry_id`, `chain_id`, `reason_code`, `message`) per curated polymer chain that isn't protein;
    7. if no protein chain is left and curation kept the entry, one `NO_PROTEIN_CHAINS` row with `chain_id=None`;
    8. the protein-only structure: `model_copy` keeping only protein-chain atoms, connections whose two partners are both protein chains, and `raw={}`;
    9. chain rows: `extract_chain_records` merged with `chain_completeness` (by `chain_id`), plus `method` and `resolution` from the `MetadataRecord`'s quality record.

    A chain without a completeness record keeps `None` in the completeness fields.

- [ ] **Step 5: Run the tests and check they pass**

Run: `uv run pytest tests/test_ppi_v21.py -v`
Expected: PASS. If a fixture fact in the table above no longer holds, stop and report it; don't edit the test to fit.

---

### Task 3: Shard runner, restarts and the CLI

**Files:**
- Modify: `examples/ppi/v2.1/a_prepare_entries.py`
- Test: `tests/test_ppi_v21.py`

**Interfaces:**
- Consumes: Task 2's `load_policies`, `prepare_entry`, `fasta_record`, `write_entry_structures`, `ParseError`.
- Produces:
  - `read_entries(path: Path) -> list[str]`: lower-cased, blank lines skipped, de-duplicated, sorted.
  - `shard_of(entries: list[str], shard: int, n_shards: int) -> list[str]` → `entries[shard::n_shards]`.
  - `shard_dir(out: Path, shard: int, n_shards: int) -> Path` → `out / "shards" / f"{shard}-of-{n_shards}"`.
  - `locate(entry_id: str, mmcif_dir: Path) -> Path` returns `<id>.cif` or `<id>.cif.gz` if either exists. Otherwise it calls `fetch_mmcif(entry_id, "pdb", None, mmcif_dir)` and returns the file it wrote.
  - `run_shard(entries: list[str], mmcif_dir: Path, out: Path, shard: int, n_shards: int, canon_path: Path, curation_path: Path) -> None`. Structure files go under `out` (shared by all shards); row files go under `shard_dir(...)`.
  - Deviation from the spec, to keep things simple: local files are opened by path, without calling `ingest_local_mmcif`. The source directory is recorded once, in `manifest.json` args.
  - `main(argv: list[str] | None = None) -> int` with argparse options `--entries`, `--mmcif-dir`, `--out`, `--shard` (default 0), `--n-shards` (default 1), `--canonicalisation` and `--curation` (both defaulting to the YAMLs next to the script), and `--merge` (Task 4).

- [ ] **Step 1: Write the failing tests**

```python
ENTRIES = ["1a3n", "12hh", "1aci"]


def _run(tmp_path, entries=ENTRIES, shard=0, n=1, mmcif_dir=FIX):
    lst = tmp_path / "entries.txt"
    lst.write_text("\n".join(entries) + "\n")
    return _mod().main(
        [
            "--entries",
            str(lst),
            "--mmcif-dir",
            str(mmcif_dir),
            "--out",
            str(tmp_path / "out"),
            "--shard",
            str(shard),
            "--n-shards",
            str(n),
        ]
    )


def _jsonl(p):
    return [json.loads(l) for l in p.read_text().splitlines()]


def test_entry_list_normalised(tmp_path):
    p = tmp_path / "e.txt"
    p.write_text("1A3N\n\n12hh\n1a3n\n")
    assert _mod().read_entries(p) == ["12hh", "1a3n"]


def test_shard_outputs(tmp_path):
    assert _run(tmp_path) == 0
    d = tmp_path / "out/shards/0-of-1"
    chains = _jsonl(d / "chains.jsonl")
    ids = [f"{c['entry_id']}_{c['chain_id']}" for c in chains]
    headers = [
        l[1:]
        for l in (d / "chains.fasta").read_text().splitlines()
        if l.startswith(">")
    ]
    assert headers == ids
    assert sorted((d / "done.txt").read_text().split()) == sorted(ENTRIES)
    assert any(
        e["reason_code"] == "METHOD_EXCLUDED"
        for e in _jsonl(d / "exclusions.jsonl")
    )
    manifest = json.loads((d / "manifest.json").read_text())
    assert manifest["counts"]["entries_done"] == 3


def test_restart_skips_done_and_failed(tmp_path):
    _run(tmp_path)
    d = tmp_path / "out/shards/0-of-1"
    before = (d / "chains.jsonl").read_text()
    _run(tmp_path)
    assert (d / "chains.jsonl").read_text() == before


def test_unparseable_entry_is_recorded(tmp_path):
    bad = tmp_path / "mirror"
    bad.mkdir()
    (bad / "9zzz.cif").write_text("this is not mmcif\n")
    shutil.copy(FIX / "1a3n.cif", bad / "1a3n.cif")
    assert _run(tmp_path, ["9zzz", "1a3n"], mmcif_dir=bad) == 0
    d = tmp_path / "out/shards/0-of-1"
    [fail] = _jsonl(d / "failures.jsonl")
    assert fail["entry_id"] == "9zzz"
    assert {"stage", "error_type", "message"} <= fail.keys()
    assert "1a3n" in (d / "done.txt").read_text().split()
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `uv run pytest tests/test_ppi_v21.py -k "shard or restart or unparseable or normalised" -v`
Expected: FAIL (`AttributeError`, no `main` or `read_entries`).

- [ ] **Step 3: Implement**
  - `run_shard` reads the entry ids already in `done.txt` and `failures.jsonl` and skips them.
  - For each remaining entry, `locate` it and call `prepare_entry`, catching `ParseError`, `ValueError`, `OSError`, `RuntimeError` and `gemmi`'s `RuntimeError` (the same class).
    - On error, append one failures row (`entry_id`, `stage` as `"locate"` or `"prepare"`, `error_type`, `message`).
    - On success, call `write_entry_structures`, then append to `chains.jsonl`, `chains.fasta` and `exclusions.jsonl`, and append the id to `done.txt` last.
  - Open the files in append mode and flush after each entry.
  - At the end, write `manifest.json` (overwritten on every run) with:
    - args, plus both policies via `model_dump()`;
    - `pandora.__version__` and `gemmi.__version__`;
    - start and end timestamps;
    - `counts`: `entries_in_shard`, `entries_done`, `entries_failed`, `chains_kept`.
  - `main` returns 0.

- [ ] **Step 4: Run the tests and check they pass**

Run: `uv run pytest tests/test_ppi_v21.py -v`
Expected: PASS.

---

### Task 4: Merge and per-filter counts

**Files:**
- Modify: `examples/ppi/v2.1/a_prepare_entries.py`
- Test: `tests/test_ppi_v21.py`

**Interfaces:**
- Consumes: Task 3's shard layout.
- Produces: `merge(out: Path) -> None`, run by `main(["--merge", "--out", ...])`. It ignores the entry options.
  - It concatenates, in sorted shard-directory order, `chains.jsonl`, `chains.fasta`, `exclusions.jsonl` and `failures.jsonl` into `out/`.
  - It writes `out/counts.tsv` with the header `level\treason\tcount`, plus these rows:
    - `entry\tin`, `entry\tkept`, `entry\tfailed`, `chain\tkept`;
    - one `entry\t<reason>` per entry-level exclusion code (`chain_id` is null);
    - one `chain\t<reason>` per chain-level code;
    - all sorted by level, then reason.

  `entry in` is the sum of `entries_done` and `entries_failed` across the shard manifests. `entry kept` is the number of distinct `entry_id` values in the merged `chains.jsonl`.

- [ ] **Step 1: Write the failing tests**

```python
def _sorted_lines(p):
    return sorted(p.read_text().splitlines())


def test_two_shards_merge_like_one(tmp_path):
    one, two = tmp_path / "one", tmp_path / "two"
    one.mkdir()
    two.mkdir()
    _run(one)
    _run(two, shard=0, n=2)
    _run(two, shard=1, n=2)
    for d in (one, two):
        assert _mod().main(["--merge", "--out", str(d / "out")]) == 0
    for name in (
        "chains.jsonl",
        "chains.fasta",
        "exclusions.jsonl",
        "counts.tsv",
    ):
        assert _sorted_lines(one / "out" / name) == _sorted_lines(
            two / "out" / name
        )


def test_counts_table(tmp_path):
    _run(tmp_path)
    _mod().main(["--merge", "--out", str(tmp_path / "out")])
    rows = (tmp_path / "out/counts.tsv").read_text().splitlines()
    assert rows[0] == "level\treason\tcount"
    assert "entry\tin\t3" in rows
    assert "entry\tkept\t2" in rows
    assert "entry\tMETHOD_EXCLUDED\t1" in rows
    assert "chain\tNOT_PROTEIN\t2" in rows
    assert "chain\tkept\t5" in rows
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `uv run pytest tests/test_ppi_v21.py -k "merge or counts" -v`
Expected: FAIL.

- [ ] **Step 3: Implement `merge`**, counting with `collections.Counter`.

- [ ] **Step 4: Run the tests and check they pass**

Run: `uv run pytest tests/test_ppi_v21.py -v`
Expected: PASS.

---

### Task 5: Slurm example and README

**Files:**
- Create: `examples/ppi/v2.1/slurm/a_prepare_entries.sh`
- Create: `examples/ppi/v2.1/README.md`

- [ ] **Step 1: Write the Slurm script**
  - Header:
    - `#SBATCH --array=0-199`;
    - `--mem=16G` (with the comment `# spike: 6.5 GB to parse 3j3q (2.44M atoms)`);
    - `--cpus-per-task=1`, `--time=24:00:00`;
    - `--output=logs/a_%A_%a.out`.
  - `set -euo pipefail`.
  - Variables `ENTRIES`, `MMCIF_DIR` and `OUT`, each with a placeholder path.
  - The run line: `uv run python examples/ppi/v2.1/a_prepare_entries.py --entries "$ENTRIES" --mmcif-dir "$MMCIF_DIR" --out "$OUT" --shard "$SLURM_ARRAY_TASK_ID" --n-shards 200`.
  - A trailing comment: after all tasks finish, run `--merge --out "$OUT"`.

- [ ] **Step 2: Write `README.md`**
  - The part A–G table from the spec, with only A built.
  - A local run on the fixtures (the exact command from Task 3's test), plus `--merge`.
  - Codon:
    - point `--mmcif-dir` at a local PDB mirror if EBI has one, else pre-fetch and don't fetch from 200 tasks at once;
    - submit with `sbatch`;
    - check progress with `wc -l "$OUT"/shards/*/done.txt`;
    - rerun killed tasks by resubmitting the same array index;
    - run entries killed for memory again on a high-memory partition, with an `--entries` file listing just those entries.
  - The output file table from the spec.

- [ ] **Step 3: Check the script parses**

Run: `bash -n examples/ppi/v2.1/slurm/a_prepare_entries.sh`
Expected: no output, exit 0.

---

### Task 6: Whole-change verification

- [ ] **Step 1:** `uv run pytest`. Expected: all pass, including the existing suite.
- [ ] **Step 2:** `uv run ruff format --check . && uv run ruff check .`. Expected: clean.
- [ ] **Step 3:** Run the script by hand on all 100 fixtures with `--n-shards 4` (four runs), then `--merge`. Read `counts.tsv` and report it to the user as the first per-filter count. Don't commit.
