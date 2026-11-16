# Canonicalise Size Guard Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `pandora canonicalise` skips entries whose estimated atom count
(after assembly expansion) exceeds a policy limit, before building the
pydantic `Structure`, and records them as `TooLarge` failures.

**Architecture:** A new public function estimates atoms from a gemmi read
of the file, using the same assembly selection as canonicalisation. A
new top-level `max_atoms` field on `canonicalisationPolicy` turns the
check on. The CLI calls the estimate before `mmcif_to_structure`.

**Tech Stack:** Python 3, pydantic v2, gemmi, pytest, argparse.

**Spec:** `design/specs/2026-11-16-canonicalise-size-guard-design.md`

## Global Constraints

- Policy field: `max_atoms: int | None = None` on `canonicalisationPolicy`, top level; `None` = no check.
- Failure line `error_type` is exactly `"TooLarge"`; new count key `entries_too_large` (also counted in `entries_failed`).
- The estimate never builds pydantic structure records.
- v2.1 example policy sets `max_atoms: 1000000`.
- Tests use `datasets/dev/mmcif/` fixtures only.
- Module starts with `from __future__ import annotations`; Google-style docstring; exported via `__all__`.

## Review Focus

- Entry with no assemblies under `standardize_biological_assembly` → estimate is model-1 AU atoms (canonicalisation leaves it unexpanded). Test in Task 1.
- `.cif.gz` input → estimate works (gemmi reads gzip). Covered by gemmi; no extra test.
- Unreadable file in the estimate → gemmi raises `RuntimeError`/`ValueError`, caught by `ENTRY_ERRORS` as a normal failure line, not a crash. Test in Task 2.
- Entry exactly at the limit → kept (check is `>`). Test in Task 2.
- `max_atoms` unset → no gemmi pre-read at all (no extra cost for existing users). Test in Task 2.

---

### Task 1: Policy field and atom estimate

**Files:**
- Modify: `pandora/schemas/canonicalisation.py` (`canonicalisationPolicy`, ~line 514: add field + `Attributes:` docstring line)
- Create: `pandora/canonicalisation/size.py`
- Modify: `pandora/canonicalisation/__init__.py` (export)
- Modify: `docs/reference/policies.md` (Canonicalisation section, new `### max_atoms` before `identifier_rules`), `docs/policies/canonicalisation.yaml` (add `max_atoms: integer | null` with comment)
- Test: `tests/test_canonicalisation.py`

**Interfaces:**
- Produces: `estimate_atom_count(path: str | Path, policy: canonicalisationPolicy) -> int` in `pandora.canonicalisation`.

- [ ] **Step 1: Write the failing tests**

```python
FIXTURES = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"


def _policy(strategy, **extra):
    return canonicalisationPolicy(
        policy_id="t",
        policy_name="t",
        policy_version="1",
        assembly_rules={
            "strategy": strategy,
            "preferred_assembly_source": "first",
        },
        altloc_rules={"strategy": "preserve"},
        **extra,
    )


def test_estimate_atom_count_expands_assembly():
    # 1ar2: 857 AU atoms, assembly 1 = 2 operators over A, B
    path = FIXTURES / "ar" / "1ar2.cif"
    policy = _policy("standardize_biological_assembly")
    assert estimate_atom_count(path, policy) == 1714
    parsed, _, _ = mmcif_to_structure(str(path))
    canon, _, _ = canonicalise_structure(parsed, policy)
    assert estimate_atom_count(path, policy) == len(canon.atoms)


def test_estimate_atom_count_sixty_operators():
    path = FIXTURES / "au" / "1auy.cif"  # 4045 AU atoms × 60
    assert (
        estimate_atom_count(path, _policy("standardize_biological_assembly"))
        == 242700
    )


def test_estimate_atom_count_without_expansion_is_asymmetric_unit():
    path = FIXTURES / "ar" / "1ar2.cif"
    assert estimate_atom_count(path, _policy("preserve_as_reported")) == 857


def test_policy_max_atoms_defaults_to_none():
    assert _policy("preserve_as_reported").max_atoms is None
```

Plus one test for an entry with no assemblies under `standardize_biological_assembly`: pick any fixture where `gemmi.read_structure(p).assemblies` is empty (search `datasets/dev/mmcif`), assert the estimate equals its model-1 atom count. If none exists, drop this test and note it in the commit.

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_canonicalisation.py -k "estimate or max_atoms" -v`
Expected: FAIL (ImportError on `estimate_atom_count`).

- [ ] **Step 3: Add the field and implement `estimate_atom_count`**

Field: `max_atoms: int | None = None`, docstring "Skip an entry whose estimated atom count (model 1, after assembly expansion, before ligand/water/altloc filtering) is above this. None: no check."

Estimate in `size.py`: `st = gemmi.read_structure(str(path))`; atoms per subchain of `st[0]` (`ch.subchains()`, `sc.subchain_id()`, sum of `len(residue)`). If `policy.assembly_rules.strategy != "standardize_biological_assembly"` or the selected assembly is `None`, return the model-1 total. Otherwise sum over `asm.generators` of `len(gen.operators) * sum(atoms[s] for s in gen.subchains)`. Pick the assembly by calling the existing `_select_biological_assembly(list(st.assemblies), policy.assembly_rules.preferred_assembly_source, DiagnosticBundle(), entry_id)`; it only reads `author_determined` / `software_determined` and indexes, which gemmi assemblies provide. Widen its `assemblies` type hint to `Sequence[Any]`-style duck typing in a one-line comment; no behaviour change.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_canonicalisation.py -v`
Expected: all PASS.

- [ ] **Step 5: Docs**

`docs/reference/policies.md`: one table row/section for `max_atoms` (default `null`, checked by `pandora canonicalise` before parsing; upper bound including waters, ligands and all altlocs, so not the same number as curation's `max_atoms`; over the limit → `TooLarge` failure). `docs/policies/canonicalisation.yaml`: add the field. Run `uv run zensical build --clean`; expected: builds with no new warnings.

- [ ] **Step 6: Commit**

```bash
git add pandora/schemas/canonicalisation.py pandora/canonicalisation/size.py pandora/canonicalisation/__init__.py pandora/canonicalisation/assemblies.py docs/reference/policies.md docs/policies/canonicalisation.yaml tests/test_canonicalisation.py
git commit -m "feat(canonicalisation): max_atoms policy field and estimate_atom_count"
```

### Task 2: CLI guard and v2.1 policy

**Files:**
- Modify: `pandora/cli/canonicalisation.py:56-63` (inside the per-path `try`, before `mmcif_to_structure`)
- Modify: `examples/ppi/v2.1/dataset/policies/canonicalisation.yaml` (add `max_atoms: 1000000` with a comment naming 10ic)
- Modify: `docs/usage/cli.md` if it lists `failures.jsonl` error types (grep `DuplicateEntry`; add `TooLarge` next to it)
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `estimate_atom_count(path, policy) -> int` from Task 1.

- [ ] **Step 1: Write the failing tests**

Copy `ar/1ar2.cif` (estimate 1714) and `b6/1b6w.cif` (estimate 1064) into `tmp_path/raw/`, write a policy YAML with `assembly_rules.strategy: standardize_biological_assembly`, `preferred_assembly_source: first`, `max_atoms: <limit>`, run `main(["canonicalise", "--input-dir", raw, "--policy", yaml, "--output-dir", out, "--table-format", "jsonl"])`.

```python
def test_canonicalise_skips_entry_over_max_atoms(tmp_path):
    out = _run_size_guard(tmp_path, limit=1500)
    failures = _jsonl(out / "failures.jsonl")
    assert [f["error_type"] for f in failures] == ["TooLarge"]
    assert failures[0]["path"].endswith("1ar2.cif")
    assert "1714" in failures[0]["message"] and "1500" in failures[0]["message"]
    summary = json.loads((out / "canonicalisation.json").read_text())
    assert summary["counts"]["entries_too_large"] == 1
    assert summary["counts"]["entries_failed"] == 1
    assert summary["counts"]["entries_kept"] == 1

def test_canonicalise_keeps_entry_at_max_atoms(tmp_path):
    out = _run_size_guard(tmp_path, limit=1714)
    assert _jsonl(out / "failures.jsonl") == []

def test_canonicalise_without_max_atoms_skips_estimate(tmp_path, monkeypatch):
    # monkeypatch pandora.cli.canonicalisation.estimate_atom_count to raise
    # AssertionError; run with max_atoms unset; both entries kept.

def test_canonicalise_unreadable_file_is_a_failure_line(tmp_path):
    # raw/ also holds broken.cif with text "not a cif"; limit=1500;
    # run completes; failures has one line for broken.cif whose
    # error_type is not "TooLarge".
```

Check the exact key path of counts in `canonicalisation.json` (`build_step_summary`) and adjust `summary["counts"]` if it differs.

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_cli.py -k "max_atoms or unreadable" -v`
Expected: the over-limit test FAILS (no `TooLarge` line); others may pass already.

- [ ] **Step 3: Implement the guard**

Before `mmcif_to_structure`: if `policy.max_atoms is not None`, compute `n = estimate_atom_count(path, policy)`; if `n > policy.max_atoms`, increment `entries_failed` and `entries_too_large`, write `failure_line(path, "TooLarge", f"estimated {n} atoms > max_atoms {policy.max_atoms}")`, `continue`. Initialise `entries_too_large=0` in the `Counter` only when `max_atoms` is set, so existing summaries don't change.

- [ ] **Step 4: Run the full suite and lint**

Run: `uv run pytest && uv run ruff format --check . && uv run ruff check .`
Expected: all pass.

- [ ] **Step 5: v2.1 policy**

Add to `examples/ppi/v2.1/dataset/policies/canonicalisation.yaml`:

```yaml
# Size guard before parsing (10ic ran out of memory on codon). Estimate
# includes waters/ligands, so it is an upper bound on curation's count.
max_atoms: 1000000
```

Run `bash src/01_prepare_entries.sh` from `examples/ppi/v2.1`; expected: same `entries_kept` counts as before (no dev fixture is near 1M).

- [ ] **Step 6: Commit**

```bash
git add pandora/cli/canonicalisation.py tests/test_cli.py examples/ppi/v2.1/dataset/policies/canonicalisation.yaml docs/usage/cli.md
git commit -m "feat(cli): canonicalise skips entries over max_atoms before parsing"
```
