# Quality and Completeness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let curation express the quality and completeness filters the four benchmarks share: cryo-EM resolution, per-method resolution, R-factor limits, modified/non-standard residues, missing residues against SEQRES, and size limits.

**Architecture:** Metadata gains the missing quality values. Canonicalisation gains a `modified_residue_rules` step. `pandora/datasets` gains a pure `chain_completeness` function, and `curate_structure` gains entry rules and chain rules that drop chains and return one `ExclusionRecord` per dropped chain.

**Tech Stack:** Python 3.11+, pydantic v2, pytest, Ruff, uv, Zensical.

**Spec:** `design/specs/2026-10-08-quality-completeness-design.md`

## Global Constraints

- Every module starts with `from __future__ import annotations`; line length 80; `uv run ruff format .` and `uv run ruff check .` clean.
- New models in `pandora/schemas/`, logic in the matching package. Never mutate inputs; use `.model_copy(update=...)`.
- Public functions: type hints, Google-style docstring with `Args:`/`Returns:`, exported in the package `__all__`. Internal helpers start with `_` and get a one-line docstring.
- Tests use only `datasets/dev/mmcif/` fixtures. No network.
- No new dependencies.
- Catch named exceptions only (Ruff BLE001, S110).
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

- **Entry with no metadata and R-factor rules set:** no method is known, so R-factor rules are skipped (resolution rules still apply). Pinned in Task 5 (`test_rfactor_rules_skip_entries_without_metadata`).
- **A chain failing two chain rules:** exactly one `ExclusionRecord` for it, not two. Pinned in Task 6 (`test_chain_failing_two_rules_gets_one_record`).
- **A chain whose SEQRES is MSE but whose atoms were mapped to MET:** completeness still reports it. Pinned in Task 3 (`test_completeness_accepts_mapped_modified_residue`).
- **Middle fraction when the whole chain is missing (denominator 0):** 0, no `ZeroDivisionError`. Pinned in Task 6 (`test_middle_fraction_zero_denominator_is_zero`).
- **Method names in policies with odd casing or whitespace:** `" x-ray diffraction "` must match `"X-ray diffraction"`. Pinned in Task 5 (`test_resolution_by_method_strictest_wins`, which uses mixed case).

## Fixture facts (checked while planning)

| Fixture | Fact |
|---|---|
| 22jy | Electron microscopy, `_em_3d_reconstruction.resolution` 2.2, no `_refine` resolution; 5 protein chains of SEQRES 51, fully observed |
| 1a08 | X-ray, resolution 2.2, `r_free` None, `pdbx_Rmerge_I_obs` 0.086, no Rsym; non-standard polymer residues ACE, DIP, FTY; `_pdbx_struct_mod_residue` lists FTY→TYR at (B, 2) and (D, 2) only |
| 1b6w | X-ray; MSE at A 1 and A 35 (both in `_pdbx_struct_mod_residue`, parent MET); MSE rows have **no SE atom**; SEQRES `mon_id` at 1 is `MSE` |
| 10mv | One protein chain A: SEQRES 275, missing N 21, C 18, middle 5 (same under both definitions) |
| 1aui | Protein chains A (SEQRES 521; N 13, C 35, middle 95; 378 observed residues) and B (SEQRES 169; N 4; 165 observed) |
| 1p58 | CA-only trace. Chain A unobserved: (495, 0, 0, 5); incomplete_backbone: (495, 495, 0, 0) |
| 1ayi | X-ray, 704 atoms, one protein chain |

---

### Task 1: Cryo-EM resolution, Rsym and Rmerge in `QualityRecord`

**Files:**
- Modify: `pandora/schemas/metadata.py` (`QualityRecord`)
- Modify: `pandora/metadata/mmcif.py` (`extract_quality`)
- Test: `tests/test_metadata.py`

**Interfaces:**
- Produces: `QualityRecord.r_sym: float | None`, `QualityRecord.r_merge: float | None`; `QualityRecord.resolution` now set for cryo-EM entries.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_metadata.py`)

```python
MMCIF_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"


def _quality(entry_id: str):
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / f"{entry_id}.cif"))
    return collect_metadata(structure).quality


def test_cryo_em_resolution_falls_back_to_em_3d_reconstruction():
    quality = _quality("22jy")
    assert quality.resolution == 2.2
    assert "_em_3d_reconstruction" in quality.provenance.source_category


def test_xray_resolution_comes_from_refine():
    quality = _quality("1a08")
    assert quality.resolution == 2.2
    assert "_em_3d_reconstruction" not in quality.provenance.source_category


def test_rmerge_extracted_when_rsym_missing():
    quality = _quality("1a08")
    assert quality.r_sym is None
    assert quality.r_merge == 0.086
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_metadata.py -k "cryo_em or refine or rmerge" -v`
Expected: `test_cryo_em...` fails (`None != 2.2`); `test_rmerge...` fails (`AttributeError`/no field `r_sym`).

- [ ] **Step 3: Add the fields to `QualityRecord`**

In `pandora/schemas/metadata.py`, add after `r_free: float | None = None`:

```python
    r_sym: float | None = None
    r_merge: float | None = None
```

and to the docstring `Attributes:` after `r_free`:

```
        r_sym: The merging R-factor as Rsym (`_reflns.pdbx_Rsym_value`).
        r_merge: The merging R-factor as Rmerge
            (`_reflns.pdbx_Rmerge_I_obs`); the same statistic as
            `r_sym` under another name.
```

and change `resolution:` in the docstring to:

```
        resolution: The structure's resolution, in angstroms, from
            `_refine`, or `_em_3d_reconstruction` for cryo-EM entries.
```

- [ ] **Step 4: Update `extract_quality`**

Replace the body of `extract_quality` in `pandora/metadata/mmcif.py` from `exptl_rows = ...` to the end with:

```python
    exptl_rows = raw_rows(structure, "_exptl")
    refine = first_row(structure, "_refine")
    reflns = first_row(structure, "_reflns")
    em = first_row(structure, "_em_3d_reconstruction")

    if not exptl_rows and not refine and not reflns and not em:
        return None

    methods = [
        method
        for method in (clean(row.get("method")) for row in exptl_rows)
        if method is not None
    ]

    categories = "_exptl,_refine,_reflns"
    resolution = as_float(refine.get("ls_d_res_high"))
    if resolution is None:
        resolution = as_float(em.get("resolution"))
        if resolution is not None:
            categories += ",_em_3d_reconstruction"

    return QualityRecord(
        experimental_method="; ".join(methods) if methods else None,
        resolution=resolution,
        r_work=as_float(refine.get("ls_R_factor_R_work")),
        r_free=as_float(refine.get("ls_R_factor_R_free")),
        r_sym=as_float(reflns.get("pdbx_Rsym_value")),
        r_merge=as_float(reflns.get("pdbx_Rmerge_I_obs")),
        observed_reflections=as_int(reflns.get("number_obs")),
        percent_possible_observed=as_float(reflns.get("percent_possible_obs")),
        mean_b_factor=as_float(refine.get("B_iso_mean")),
        provenance=provenance(categories),
    )
```

Update the docstring: "Reads the `_exptl`, `_refine`, `_reflns` and `_em_3d_reconstruction` categories ... Resolution comes from `_refine`, falling back to `_em_3d_reconstruction` (cryo-EM)." and the `Returns:` line to list all four categories.

- [ ] **Step 5: Run the tests**

Run: `uv run pytest tests/test_metadata.py tests/test_curation.py -v`
Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add pandora/schemas/metadata.py pandora/metadata/mmcif.py tests/test_metadata.py
git commit -m "Read cryo-EM resolution, Rsym and Rmerge into QualityRecord

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: `modified_residue_rules` in canonicalisation

**Files:**
- Modify: `pandora/schemas/canonicalisation.py`
- Create: `pandora/canonicalisation/modified_residues.py`
- Modify: `pandora/canonicalisation/canonicalise.py`
- Modify: `docs/reference/policies.md`, `docs/policies/canonicalisation.yaml`
- Test: `tests/test_canonicalisation.py`

**Interfaces:**
- Produces: `ModifiedResidueRules(strategy: Literal["preserve", "map_to_parent"] = "preserve", comp_ids: list[str] = ["MSE"])`; `ModifiedResidueMappingItem`, `ModifiedResidueMapping`; `CanonicalMappings.modified_residue_mapping`; `canonicalisationPolicy.modified_residue_rules`.
- Produces: `pandora.canonicalisation.modified_residues.STANDARD_RESIDUES: frozenset[str]` (used by Task 3 and Task 5) and `_map_modified_residues(atoms, mod_rows, rules, diagnostics, entry_id) -> tuple[list[AtomSiteRecord], ModifiedResidueMapping]`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_canonicalisation.py`; add `ModifiedResidueRules` to the existing `from pandora.schemas.canonicalisation import (...)` block)

```python
from pandora.canonicalisation.modified_residues import _map_modified_residues


def test_modified_residues_preserved_by_default():
    structure = _load("1b6w")
    canonical, mappings, _ = canonicalise_structure(
        structure, _lenient_policy()
    )
    assert any(a.label_comp_id == "MSE" for a in canonical.atoms)
    assert mappings.modified_residue_mapping.items == []


def test_mse_mapped_to_met():
    structure = _load("1b6w")
    policy = _lenient_policy(
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent")
    )
    canonical, mappings, provenance = canonicalise_structure(structure, policy)

    assert not any(a.label_comp_id == "MSE" for a in canonical.atoms)
    residue_35 = [
        a
        for a in canonical.atoms
        if a.label_asym_id == "A" and a.label_seq_id == 35
    ]
    assert residue_35
    assert all(
        a.label_comp_id == "MET" and a.auth_comp_id == "MET" for a in residue_35
    )
    assert {
        (i.chain_id, i.seq_id, i.original_comp_id, i.parent_comp_id)
        for i in mappings.modified_residue_mapping.items
    } == {("A", 1, "MSE", "MET"), ("A", 35, "MSE", "MET")}
    assert "modified_residues:map_to_parent" in provenance.transforms
    # input untouched
    assert any(a.label_comp_id == "MSE" for a in structure.atoms)


def test_mse_selenium_becomes_sulfur_delta():
    se = _atom(label_comp_id="MSE", atom_id="SE").model_copy(
        update={"type_symbol": "SE", "auth_atom_id": "SE"}
    )
    atoms, mapping = _map_modified_residues(
        [se],
        [],
        ModifiedResidueRules(strategy="map_to_parent"),
        DiagnosticBundle(),
        "test",
    )
    assert atoms[0].label_comp_id == "MET"
    assert atoms[0].label_atom_id == "SD"
    assert atoms[0].auth_atom_id == "SD"
    assert atoms[0].type_symbol == "S"
    assert len(mapping.items) == 1
    assert se.label_atom_id == "SE"


def test_listed_modified_residue_mapped_from_mod_residue_table():
    structure = _load("1a08")
    diagnostics = DiagnosticBundle()
    atoms, mapping = _map_modified_residues(
        list(structure.atoms),
        structure.raw["_pdbx_struct_mod_residue"],
        ModifiedResidueRules(strategy="map_to_parent", comp_ids=["FTY"]),
        diagnostics,
        structure.entry_id,
    )
    assert not any(a.label_comp_id == "FTY" for a in atoms)
    assert {
        (i.chain_id, i.seq_id, i.parent_comp_id) for i in mapping.items
    } == {
        ("B", 2, "TYR"),
        ("D", 2, "TYR"),
    }
    # ACE and DIP were not asked for: untouched, no warning
    assert any(a.label_comp_id == "ACE" for a in atoms)
    assert diagnostics.warnings == []


def test_unmapped_modified_residue_gets_warning_not_guess():
    structure = _load("1a08")
    diagnostics = DiagnosticBundle()
    atoms, _ = _map_modified_residues(
        list(structure.atoms),
        structure.raw["_pdbx_struct_mod_residue"],
        ModifiedResidueRules(strategy="map_to_parent", comp_ids=[]),
        diagnostics,
        structure.entry_id,
    )
    assert not any(a.label_comp_id == "FTY" for a in atoms)
    assert any(a.label_comp_id == "ACE" for a in atoms)
    unmapped = {
        d.context["comp_id"]
        for d in diagnostics.warnings
        if d.code == "MODIFIED_RESIDUE_UNMAPPED"
    }
    assert unmapped == {"ACE", "DIP"}
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_canonicalisation.py -k "modified or mse" -v`
Expected: collection error, `ModuleNotFoundError: pandora.canonicalisation.modified_residues`.

- [ ] **Step 3: Add the schemas** in `pandora/schemas/canonicalisation.py`

Next to the other strategy literals (after `LigandStrategy`):

```python
ModifiedResidueStrategy = Literal["preserve", "map_to_parent"]
```

After `ResidueNumberMapping`:

```python
class ModifiedResidueMappingItem(BaseModel):
    """One modified residue renamed to its parent.

    Attributes:
        chain_id: The residue's original `label_asym_id`.
        seq_id: The residue's `label_seq_id`.
        auth_seq_id: The residue's `auth_seq_id`.
        original_comp_id: The modified residue's comp_id (e.g. "MSE").
        parent_comp_id: The standard parent it was renamed to (e.g.
            "MET").
    """

    chain_id: str
    seq_id: int | None = None
    auth_seq_id: str | None = None
    original_comp_id: str
    parent_comp_id: str


class ModifiedResidueMapping(BaseModel):
    """Every modified residue renamed during one canonicalisation run.

    Attributes:
        items: One entry per renamed residue.
    """

    items: list[ModifiedResidueMappingItem] = Field(default_factory=list)
```

After `LigandRules`:

```python
class ModifiedResidueRules(BaseModel):
    """Canonicalisation policy for modified polymer residues.

    Attributes:
        strategy: Kept as-is, or renamed to the parent residue given by
            `_pdbx_struct_mod_residue` (MSE always maps to MET).
        comp_ids: Which modified residues to map. Empty means every
            modified residue with a known parent.
    """

    strategy: ModifiedResidueStrategy = "preserve"
    comp_ids: list[str] = Field(default_factory=lambda: ["MSE"])
```

Add to `CanonicalMappings` (field + docstring line `modified_residue_mapping: The modified residues renamed during this run.`):

```python
    modified_residue_mapping: ModifiedResidueMapping = Field(
        default_factory=ModifiedResidueMapping
    )
```

Add to `canonicalisationPolicy` after `ligand_rules` (field + docstring line `modified_residue_rules: The modified-residue mapping rules.`):

```python
    modified_residue_rules: ModifiedResidueRules = Field(
        default_factory=ModifiedResidueRules
    )
```

- [ ] **Step 4: Create `pandora/canonicalisation/modified_residues.py`**

```python
from __future__ import annotations

from pandora.schemas.canonicalisation import (
    ModifiedResidueMapping,
    ModifiedResidueMappingItem,
    ModifiedResidueRules,
)
from pandora.schemas.common import Diagnostic, DiagnosticBundle
from pandora.schemas.structure import AtomSiteRecord

STANDARD_RESIDUES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "UNK",
        "A", "C", "G", "U", "N",
        "DA", "DC", "DG", "DT", "DU", "DN",
    }
)  # fmt: skip

# Parents assumed even when _pdbx_struct_mod_residue doesn't list them.
_BUILTIN_PARENTS = {"MSE": "MET"}


def _map_modified_residues(
    atoms: list[AtomSiteRecord],
    mod_rows: list[dict[str, str | None]],
    rules: ModifiedResidueRules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> tuple[list[AtomSiteRecord], ModifiedResidueMapping]:
    """Rename modified polymer residues to their parent, per the rules."""

    mapping = ModifiedResidueMapping()
    if rules.strategy == "preserve":
        return atoms, mapping

    parents = {
        (row.get("label_asym_id"), row.get("label_seq_id"), comp): parent
        for row in mod_rows
        if (comp := row.get("label_comp_id"))
        and (parent := row.get("parent_comp_id"))
    }
    wanted = set(rules.comp_ids)
    seen: set[tuple[str, str, str]] = set()
    result: list[AtomSiteRecord] = []
    for a in atoms:
        comp = a.label_comp_id
        if (
            a.label_seq_id is None
            or comp in STANDARD_RESIDUES
            or (wanted and comp not in wanted)
        ):
            result.append(a)
            continue

        key = (a.label_asym_id, str(a.label_seq_id), comp)
        parent = parents.get(key) or _BUILTIN_PARENTS.get(comp)
        first = key not in seen
        seen.add(key)
        if parent is None:
            if first:
                diagnostics.warnings.append(
                    Diagnostic(
                        code="MODIFIED_RESIDUE_UNMAPPED",
                        severity="warning",
                        message=(
                            f"No parent known for {comp} {a.label_seq_id} "
                            f"in chain {a.label_asym_id}; left unmapped"
                        ),
                        entry_id=entry_id,
                        context={
                            "chain": a.label_asym_id,
                            "seq_id": a.label_seq_id,
                            "comp_id": comp,
                        },
                    )
                )
            result.append(a)
            continue

        if first:
            mapping.items.append(
                ModifiedResidueMappingItem(
                    chain_id=a.label_asym_id,
                    seq_id=a.label_seq_id,
                    auth_seq_id=a.auth_seq_id,
                    original_comp_id=comp,
                    parent_comp_id=parent,
                )
            )
        update: dict[str, object] = {
            "label_comp_id": parent,
            "auth_comp_id": parent,
            "group_PDB": "ATOM",
        }
        if comp == "MSE" and a.label_atom_id == "SE":
            update |= {
                "label_atom_id": "SD",
                "auth_atom_id": "SD",
                "type_symbol": "S",
            }
        result.append(a.model_copy(update=update))

    return result, mapping
```

Check `Diagnostic`'s fields in `pandora/schemas/common.py` (`code`, `severity`, `message`, `entry_id`, `context`) match; they are used the same way in `canonicalisation/missing_data.py`.

- [ ] **Step 5: Wire it into `canonicalise_structure`** (`pandora/canonicalisation/canonicalise.py`)

Add the import:

```python
from pandora.canonicalisation.modified_residues import _map_modified_residues
```

Add `mrr = policy.modified_residue_rules` next to the other rule aliases. Insert this block **before** `# normalize_chain_ids`, right after `secondary_structure: SSRecord = structure.secondary_structure`:

```python
    # map_modified_residues --------------------------------
    # First, so chain ids and seq ids still match
    # _pdbx_struct_mod_residue.
    atoms, modified_residue_mapping = _map_modified_residues(
        atoms,
        structure.raw.get("_pdbx_struct_mod_residue", []),
        mrr,
        diagnostics,
        structure.entry_id,
    )
    if mrr.strategy != "preserve":
        transforms.append(f"modified_residues:{mrr.strategy}")
```

Add `modified_residue_mapping=modified_residue_mapping,` to the `CanonicalMappings(...)` call. Add "modified-residue mapping" first in the docstring's list of steps.

**Check:** `_apply_chain_map(atoms, asym_units, structure, chain_map)` must use the passed-in `atoms`, not `structure.atoms`. Read `pandora/canonicalisation/chain_ids.py:71` to confirm. If it re-reads `structure.atoms`, the mapping is lost; in that case pass `structure.model_copy(update={"atoms": atoms})` instead.

- [ ] **Step 6: Run the tests**

Run: `uv run pytest tests/test_canonicalisation.py -v`
Expected: all PASS.

- [ ] **Step 7: Document the rule**

In `docs/reference/policies.md`, add after the `### ligand_rules` section (before `### validation_rules`):

```markdown
### `modified_residue_rules`

**`strategy`** (default `preserve`)

| Value | Description |
|---|---|
| `preserve` | Keep modified residues (e.g. `MSE`) as reported. |
| `map_to_parent` | Rename each modified polymer residue in `comp_ids` to its parent from `_pdbx_struct_mod_residue`. `MSE` always maps to `MET`, and its `SE` atom becomes `SD` (element `S`). Each rename is recorded in `CanonicalMappings.modified_residue_mapping`. |

`comp_ids` (default `["MSE"]`): which residues to map. An empty list maps every modified polymer residue with a known parent. A residue with no known parent is left unchanged and gets a `MODIFIED_RESIDUE_UNMAPPED` warning; Pandora doesn't guess a parent.

This step runs first, before chain ids are normalised, so the ids still match `_pdbx_struct_mod_residue`.
```

In `docs/policies/canonicalisation.yaml`, add before `validation_rules:` (same indentation as `ligand_rules:`):

```yaml
  modified_residue_rules:
    strategy: string
    # preserve      — Keep modified residues as reported (default).
    # map_to_parent — Rename listed residues to their parent residue.
    comp_ids: list[string]
    # Residues to map (default ["MSE"]); empty = every one with a parent.
```

- [ ] **Step 8: Lint and commit**

Run: `uv run ruff format . && uv run ruff check .`

```bash
git add pandora/schemas/canonicalisation.py pandora/canonicalisation/modified_residues.py pandora/canonicalisation/canonicalise.py tests/test_canonicalisation.py docs/reference/policies.md docs/policies/canonicalisation.yaml
git commit -m "Map modified residues such as MSE to their parent

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: `chain_completeness`

**Files:**
- Modify: `pandora/schemas/dataset.py`
- Create: `pandora/datasets/completeness.py`
- Modify: `pandora/datasets/__init__.py`
- Test: `tests/test_completeness.py`

**Interfaces:**
- Consumes: `STANDARD_RESIDUES` (Task 2); `_PROTEIN_BACKBONE`, `_backbone_by_entity` from `pandora/canonicalisation/missing_data.py`.
- Produces: `MissingResidueDefinition = Literal["unobserved", "incomplete_backbone"]`; `ChainCompleteness(entry_id, chain_id, seqres_length, missing_n_term, missing_c_term, missing_middle)`; `chain_completeness(structure, definition="incomplete_backbone") -> tuple[list[ChainCompleteness], DiagnosticBundle]`.

- [ ] **Step 1: Write the failing tests** (`tests/test_completeness.py`)

```python
from pathlib import Path

from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import chain_completeness
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    IdentifierRules,
    ModifiedResidueRules,
    ResidueNumberingRules,
    ValidationRules,
    canonicalisationPolicy,
)

MMCIF_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"


def _load(entry_id: str):
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / f"{entry_id}.cif"))
    return structure


def _canonical(structure, **rules):
    policy = canonicalisationPolicy(
        policy_id="t",
        policy_name="t",
        policy_version="1",
        validation_rules=ValidationRules(fail_on_unresolved_issues=False),
        **rules,
    )
    canonical, _, _ = canonicalise_structure(structure, policy)
    return canonical


def _counts(records):
    return {
        r.chain_id: (
            r.seqres_length,
            r.missing_n_term,
            r.missing_c_term,
            r.missing_middle,
        )
        for r in records
    }


def test_completeness_counts_tails_and_middle():
    records, diagnostics = chain_completeness(_load("10mv"))
    assert _counts(records) == {"A": (275, 21, 18, 5)}
    assert records[0].entry_id == "10MV"
    assert diagnostics.warnings == []


def test_completeness_reports_every_polymer_chain():
    records, _ = chain_completeness(_load("1aui"))
    assert _counts(records) == {
        "A": (521, 13, 35, 95),
        "B": (169, 4, 0, 0),
    }


def test_incomplete_backbone_counts_ca_only_residues_as_missing():
    structure = _load("1p58")
    unobserved, _ = chain_completeness(structure, "unobserved")
    backbone, _ = chain_completeness(structure, "incomplete_backbone")
    assert _counts(unobserved)["A"] == (495, 0, 0, 5)
    assert _counts(backbone)["A"] == (495, 495, 0, 0)


def test_chain_without_seqres_gets_diagnostic_not_record():
    structure = _load("10mv")
    raw = {k: v for k, v in structure.raw.items() if k != "_entity_poly_seq"}
    records, diagnostics = chain_completeness(
        structure.model_copy(update={"raw": raw})
    )
    assert records == []
    assert [d.code for d in diagnostics.warnings] == ["NO_SEQRES"]


def test_renumbered_chain_gets_seqres_mismatch():
    canonical = _canonical(
        _load("10mv"),
        identifier_rules=IdentifierRules(
            residue_numbering=ResidueNumberingRules(strategy="renumber")
        ),
    )
    records, diagnostics = chain_completeness(canonical)
    assert records == []
    assert [d.code for d in diagnostics.warnings] == ["SEQRES_MISMATCH"]


def test_completeness_accepts_mapped_modified_residue():
    canonical = _canonical(
        _load("1b6w"),
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
    )
    records, diagnostics = chain_completeness(canonical)
    assert _counts(records) == {"A": (69, 0, 1, 0)}
    assert diagnostics.warnings == []
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_completeness.py -v`
Expected: `ImportError: cannot import name 'chain_completeness'`.

- [ ] **Step 3: Add the schema** to `pandora/schemas/dataset.py` (after `InterfaceRecord`)

```python
MissingResidueDefinition = Literal["unobserved", "incomplete_backbone"]


class ChainCompleteness(BaseModel):
    """How much of one polymer chain's SEQRES has coordinates.

    Attributes:
        entry_id: The structure's entry id.
        chain_id: The chain's `label_asym_id`.
        seqres_length: The number of positions in the chain entity's
            `_entity_poly_seq`.
        missing_n_term: Missing positions before the first present
            residue.
        missing_c_term: Missing positions after the last present
            residue.
        missing_middle: Missing positions between the first and last
            present residue.
    """

    entry_id: str
    chain_id: str
    seqres_length: int
    missing_n_term: int
    missing_c_term: int
    missing_middle: int
```

- [ ] **Step 4: Create `pandora/datasets/completeness.py`**

```python
from __future__ import annotations

from collections import defaultdict

from pandora.canonicalisation.missing_data import (
    _PROTEIN_BACKBONE,
    _backbone_by_entity,
)
from pandora.canonicalisation.modified_residues import STANDARD_RESIDUES
from pandora.schemas.common import Diagnostic, DiagnosticBundle
from pandora.schemas.dataset import ChainCompleteness, MissingResidueDefinition
from pandora.schemas.structure import Structure


def _seqres_by_entity(structure: Structure) -> dict[str, dict[int, set[str]]]:
    """mon_ids per SEQRES position, per entity, from `_entity_poly_seq`."""

    seqres: dict[str, dict[int, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for row in structure.raw.get("_entity_poly_seq", []):
        num, entity_id = row.get("num"), row.get("entity_id")
        if entity_id and num and num.isdigit():
            seqres[entity_id][int(num)].add(row.get("mon_id") or "")
    return seqres


def _matches_seqres(
    residues: dict[int, tuple[str, set[str]]],
    positions: dict[int, set[str]],
) -> bool:
    """Whether every residue sits on a SEQRES position with its comp_id."""

    for seq_id, (comp_id, _) in residues.items():
        mon_ids = positions.get(seq_id)
        if mon_ids is None:
            return False
        if comp_id not in mon_ids and mon_ids <= STANDARD_RESIDUES:
            return False
    return True


def chain_completeness(
    structure: Structure,
    definition: MissingResidueDefinition = "incomplete_backbone",
) -> tuple[list[ChainCompleteness], DiagnosticBundle]:
    """Count the missing residues of each polymer chain against SEQRES.

    SEQRES is `_entity_poly_seq`; a residue lines up with it by
    `label_seq_id`. A chain that has no SEQRES, or whose residues don't
    line up (for example after renumbering), gets no record and a
    `NO_SEQRES` or `SEQRES_MISMATCH` warning instead.

    Args:
        structure: The structure to measure. Its `label_seq_id`s must
            still index SEQRES.
        definition: `"unobserved"` counts a residue as missing when it
            has no atoms; `"incomplete_backbone"` also when any backbone
            atom (N, CA, C, O for protein) is missing.

    Returns:
        `(records, diagnostics)`: one `ChainCompleteness` per measured
        polymer chain, and a warning per chain that could not be
        measured.
    """

    diagnostics = DiagnosticBundle()
    seqres = _seqres_by_entity(structure)
    backbone_by_entity = _backbone_by_entity(structure.entities)
    polymer_entities = {e.id for e in structure.entities if e.type == "polymer"}

    # chain -> seq_id -> (comp_id, atom names)
    residues: dict[str, dict[int, tuple[str, set[str]]]] = defaultdict(dict)
    for a in structure.atoms:
        if a.label_seq_id is None:
            continue
        chain = residues[a.label_asym_id]
        if a.label_seq_id not in chain:
            chain[a.label_seq_id] = (a.label_comp_id, set())
        chain[a.label_seq_id][1].add(a.label_atom_id)

    records: list[ChainCompleteness] = []
    for asym in structure.asym_units:
        if asym.entity_id not in polymer_entities:
            continue
        positions = seqres.get(asym.entity_id, {})
        chain_residues = residues.get(asym.id, {})
        code = None
        if not positions:
            code = "NO_SEQRES"
        elif not _matches_seqres(chain_residues, positions):
            code = "SEQRES_MISMATCH"
        if code is not None:
            diagnostics.warnings.append(
                Diagnostic(
                    code=code,
                    severity="warning",
                    message=f"Chain {asym.id}: completeness not measured "
                    f"({code})",
                    entry_id=structure.entry_id,
                    context={"chain": asym.id},
                )
            )
            continue

        length = len(positions)
        backbone = backbone_by_entity.get(asym.entity_id, _PROTEIN_BACKBONE)
        present = sorted(
            seq_id
            for seq_id, (_, names) in chain_residues.items()
            if definition == "unobserved" or backbone <= names
        )
        if present:
            first, last = present[0], present[-1]
            n_term, c_term = first - 1, length - last
            middle = (last - first + 1) - len(present)
        else:
            n_term, c_term, middle = length, 0, 0
        records.append(
            ChainCompleteness(
                entry_id=structure.entry_id,
                chain_id=asym.id,
                seqres_length=length,
                missing_n_term=n_term,
                missing_c_term=c_term,
                missing_middle=middle,
            )
        )
    return records, diagnostics
```

Note: `_matches_seqres` treats a position whose SEQRES `mon_id`s include any non-standard residue as matching any comp_id (`mon_ids <= STANDARD_RESIDUES` is False), so MSE mapped to MET still lines up. Positions are assumed to be `1..length`; a position outside that set fails `_matches_seqres`.

- [ ] **Step 5: Export it** in `pandora/datasets/__init__.py`

```python
from pandora.datasets.completeness import chain_completeness
```

and add `"chain_completeness",` to `__all__` (alphabetical, before `"curate_structure"`).

- [ ] **Step 6: Run the tests**

Run: `uv run pytest tests/test_completeness.py -v`
Expected: all PASS.

- [ ] **Step 7: Lint and commit**

Run: `uv run ruff format . && uv run ruff check .`

```bash
git add pandora/schemas/dataset.py pandora/datasets/completeness.py pandora/datasets/__init__.py tests/test_completeness.py
git commit -m "Add chain_completeness: missing residues against SEQRES

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: `curate_structure` returns a list of exclusions

No behaviour change: the list holds zero or one record. This task migrates every caller so Tasks 5–6 only add rules.

**Files:**
- Modify: `pandora/schemas/dataset.py` (`ExclusionRecord`)
- Modify: `pandora/datasets/curation.py` (`curate_structure`)
- Modify: `pandora/cli/app.py` (`_cmd_curate`), `pandora/provenance/reproduce.py`, `examples/dataset_pipeline.py`, `examples/ppi_dataset_pipeline.py`
- Modify docs: `docs/usage/datasets.md`, `docs/getting-started/quickstart.md`, `docs/recipes/ppi-01.md`, `docs/contributing/architecture.md`
- Test: `tests/test_curation.py`

**Interfaces:**
- Produces: `ExclusionRecord.chain_id: str | None = None`; `curate_structure(structure, metadata, policy) -> tuple[Structure | None, list[ExclusionRecord], CurationProvenance]`.

- [ ] **Step 1: Update the tests to the new return type**

In `tests/test_curation.py`, every `curated, exclusion, ... = curate_structure(...)` becomes `curated, exclusions, ... = curate_structure(...)`. Replace `assert exclusion is None` with `assert exclusions == []` and `exclusion.reason_code == "X"` with `[e.reason_code for e in exclusions] == ["X"]`. Add:

```python
def test_entry_exclusion_has_no_chain_id():
    structure = _load("1ayi")
    _, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_resolution=0.1)),
    )
    assert [(e.reason_code, e.chain_id) for e in exclusions] == [
        ("NULL_RESOLUTION", None)
    ]
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_curation.py -v`
Expected: FAIL (`exclusions == []` compares `None`; `chain_id` attribute missing).

- [ ] **Step 3: Add `chain_id` to `ExclusionRecord`**

```python
    entry_id: str
    chain_id: str | None = None
    reason_code: ExclusionReason
    message: str
```

Docstring: `chain_id: The excluded chain's label_asym_id, or None when the whole entry was excluded.` Update the class docstring to "Record of why one structure, or one chain of it, was excluded ...".

- [ ] **Step 4: Change `curate_structure`'s return**

```python
    exclusion = _check_quality(
        structure, metadata, policy.quality_rules
    ) or _check_organism(structure, metadata, policy.organism_rules)
    if exclusion is not None:
        return None, [exclusion], provenance
    return (
        _apply_content_rules(structure, policy.content_rules),
        [],
        provenance,
    )
```

Signature return annotation: `tuple[Structure | None, list[ExclusionRecord], CurationProvenance]`. Docstring `Returns:`:

```
        `(curated_structure, exclusions, provenance)`. When the entry
        is kept, `curated_structure` has content rules applied and
        failing chains removed, and `exclusions` holds one record per
        removed chain. When the entry is excluded, `curated_structure`
        is None and the last record in `exclusions` has
        `chain_id=None`. Provenance is always populated.
```

- [ ] **Step 5: Migrate callers**

`pandora/cli/app.py` `_cmd_curate` loop and print:

```python
    for entry_id, structure in structures.items():
        metadata = collect_metadata(structure)
        curated, records, prov = curate_structure(structure, metadata, policy)
        provenance[entry_id] = prov
        exclusions.extend(records)
        if curated is not None:
            retained[entry_id] = curated

    ...
    print(
        f"curated: {len(retained)} retained, "
        f"{len(structures) - len(retained)} excluded -> {output_dir}"
    )
```

`pandora/provenance/reproduce.py` (around line 195):

```python
            curated, records, _ = curate_structure(
                structures[entry_id], metadata, manifest.curation_policy
            )
            excluded.extend(records)
            if curated is None:
                del structures[entry_id]
            else:
                structures[entry_id] = curated
```

`examples/dataset_pipeline.py` and `examples/ppi_dataset_pipeline.py`, same shape:

```python
    curated, records, _ = curate_structure(
        structures[entry_id], None, curation_policy
    )
    exclusions.extend(records)
    if curated is None:
        del structures[entry_id]
    else:
        structures[entry_id] = curated
```

Docs: apply the same shape to the snippets in `docs/getting-started/quickstart.md` (~line 184) and `docs/recipes/ppi-01.md` (~lines 73 and 208). In `docs/usage/datasets.md`:
- lines 25–28: "returns `(curated_structure, exclusions, provenance)`. A kept entry has `exclusions` empty, or one record per chain the chain rules removed; an excluded entry has `curated_structure=None` and ends with an entry-level record (`chain_id=None`)."
- line 39: `curated, exclusions, provenance = curate_structure(` and `print(curated is not None, exclusions)`
- line 56–57: `curated, exclusions, _ = ...`; `print(exclusions[0].reason_code, "-", exclusions[0].message)`
- line 76: unchanged (`curated, _, _`).

Add after the "Curate one structure" examples in `docs/usage/datasets.md`:

```markdown
!!! note "Migrating from the single-exclusion return"
    Before chain-level rules, `curate_structure` returned
    `(structure, exclusion | None, provenance)`. It now returns a list.
    Replace `if curated is None: excluded.append(exclusion)` with
    `excluded.extend(exclusions)`, which also keeps the records of
    chains removed from entries that were kept.
```

`docs/contributing/architecture.md` line 79: return column becomes `` `(Structure | None, list[ExclusionRecord], CurationProvenance)` — `None` structure means excluded; records also cover removed chains ``.

- [ ] **Step 6: Run the full suite**

Run: `uv run pytest`
Expected: all PASS. If `tests/test_cli.py` checks the "excluded" count text, update it to the entry count.

- [ ] **Step 7: Run both examples to check they still work**

Run: `uv run python examples/dataset_pipeline.py > /dev/null && uv run python examples/ppi_dataset_pipeline.py > /dev/null; echo $?`
Expected: `0`. (If either needs network or binaries and was already failing before this change, confirm with `git stash` and note it rather than fixing it here.)

- [ ] **Step 8: Lint and commit**

```bash
uv run ruff format . && uv run ruff check .
git add -A pandora tests examples docs
git commit -m "Return a list of exclusions from curate_structure

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Entry rules: per-method resolution, R-factors, non-standard residues, max_atoms

**Files:**
- Modify: `pandora/schemas/dataset.py` (`QualityRules`, `ExclusionReason`)
- Modify: `pandora/datasets/curation.py`
- Test: `tests/test_curation.py`

**Interfaces:**
- Consumes: `QualityRecord.r_sym`, `r_merge` (Task 1); `STANDARD_RESIDUES` (Task 2); list return (Task 4).
- Produces: `QualityRules.max_resolution_by_method`, `rfactor_methods`, `max_r_free`, `max_r_free_gap`, `max_r_sym`, `null_rfactor_behavior`, `exclude_nonstandard_residues`, `allowed_nonstandard_residues`, `max_atoms`; codes `RFREE_THRESHOLD`, `RFREE_GAP_THRESHOLD`, `RSYM_THRESHOLD`, `NULL_RFACTOR`, `NONSTANDARD_RESIDUE`, `TOO_MANY_ATOMS`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_curation.py`; add `from pandora.canonicalisation import canonicalise_structure`, `from pandora.schemas.canonicalisation import ModifiedResidueRules, ValidationRules, canonicalisationPolicy`, `from pandora.schemas.metadata import MetadataRecord, QualityRecord`)

```python
def _meta(entry_id: str, **quality) -> MetadataRecord:
    return MetadataRecord(entry_id=entry_id, quality=QualityRecord(**quality))


def _codes(exclusions):
    return [e.reason_code for e in exclusions]


def test_cryo_em_entry_passes_resolution_rule_with_real_metadata():
    # Regression: cryo-EM resolution used to be None -> NULL_RESOLUTION.
    structure = _load("22jy")
    metadata = collect_metadata(structure)
    curated, exclusions, _ = curate_structure(
        structure,
        metadata,
        _policy(quality_rules=QualityRules(max_resolution=3.5)),
    )
    assert curated is not None and exclusions == []

    strict = QualityRules(
        max_resolution=3.5,
        max_resolution_by_method={"ELECTRON MICROSCOPY": 2.0},
    )
    curated, exclusions, _ = curate_structure(
        structure, metadata, _policy(quality_rules=strict)
    )
    assert curated is None
    assert _codes(exclusions) == ["RESOLUTION_THRESHOLD"]


def test_resolution_by_method_falls_back_to_max_resolution():
    structure = _load("1ayi")
    rules = QualityRules(
        max_resolution=3.0,
        max_resolution_by_method={"Electron Microscopy": 2.0},
    )
    meta = _meta(
        "1ayi", experimental_method="X-ray diffraction", resolution=2.4
    )
    curated, exclusions, _ = curate_structure(
        structure, meta, _policy(quality_rules=rules)
    )
    assert curated is not None and exclusions == []


def test_resolution_by_method_strictest_wins():
    structure = _load("1ayi")
    rules = QualityRules(
        max_resolution_by_method={
            " x-ray diffraction ": 2.5,
            "NEUTRON DIFFRACTION": 2.0,
        }
    )
    meta = _meta(
        "1ayi",
        experimental_method="X-ray diffraction; Neutron diffraction",
        resolution=2.4,
    )
    _, exclusions, _ = curate_structure(
        structure, meta, _policy(quality_rules=rules)
    )
    assert _codes(exclusions) == ["RESOLUTION_THRESHOLD"]


def test_rfactor_rules():
    structure = _load("1ayi")
    meta = _meta(
        "1ayi",
        experimental_method="X-ray diffraction",
        r_free=0.30,
        r_work=0.20,
        r_merge=0.12,
    )
    cases = [
        (QualityRules(max_r_free=0.25), ["RFREE_THRESHOLD"]),
        (QualityRules(max_r_free_gap=0.07), ["RFREE_GAP_THRESHOLD"]),
        (QualityRules(max_r_sym=0.10), ["RSYM_THRESHOLD"]),  # via r_merge
        (QualityRules(max_r_free=0.35, max_r_free_gap=0.11, max_r_sym=0.2), []),
    ]
    for rules, expected in cases:
        _, exclusions, _ = curate_structure(
            structure, meta, _policy(quality_rules=rules)
        )
        assert _codes(exclusions) == expected, rules


def test_rfactor_rules_skip_cryo_em():
    structure = _load("1ayi")
    meta = _meta(
        "1ayi", experimental_method="Electron Microscopy", resolution=2.0
    )
    rules = QualityRules(max_r_free=0.25, null_rfactor_behavior="exclude")
    curated, exclusions, _ = curate_structure(
        structure, meta, _policy(quality_rules=rules)
    )
    assert curated is not None and exclusions == []


def test_null_rfactor_behavior():
    structure = _load("1ayi")
    meta = _meta("1ayi", experimental_method="X-ray diffraction", r_free=None)
    curated, _, _ = curate_structure(
        structure, meta, _policy(quality_rules=QualityRules(max_r_free=0.25))
    )
    assert curated is not None  # default: include
    _, exclusions, _ = curate_structure(
        structure,
        meta,
        _policy(
            quality_rules=QualityRules(
                max_r_free=0.25, null_rfactor_behavior="exclude"
            )
        ),
    )
    assert _codes(exclusions) == ["NULL_RFACTOR"]


def test_rfactor_rules_skip_entries_without_metadata():
    structure = _load("1ayi")
    rules = QualityRules(max_r_free=0.25, null_rfactor_behavior="exclude")
    curated, exclusions, _ = curate_structure(
        structure, None, _policy(quality_rules=rules)
    )
    assert curated is not None and exclusions == []


def test_nonstandard_residues_exclude_entry_unless_allowed():
    structure = _load("1a08")
    _, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(exclude_nonstandard_residues=True)),
    )
    assert _codes(exclusions) == ["NONSTANDARD_RESIDUE"]
    assert "FTY" in exclusions[0].message

    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(
            quality_rules=QualityRules(
                exclude_nonstandard_residues=True,
                allowed_nonstandard_residues=["ACE", "DIP", "FTY"],
            )
        ),
    )
    assert curated is not None and exclusions == []


def test_mapped_mse_is_standard_for_curation():
    canon_policy = canonicalisationPolicy(
        policy_id="t",
        policy_name="t",
        policy_version="1",
        validation_rules=ValidationRules(fail_on_unresolved_issues=False),
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
    )
    canonical, _, _ = canonicalise_structure(_load("1b6w"), canon_policy)
    curated, exclusions, _ = curate_structure(
        canonical,
        None,
        _policy(quality_rules=QualityRules(exclude_nonstandard_residues=True)),
    )
    assert curated is not None and exclusions == []


def test_max_atoms():
    structure = _load("1ayi")  # 704 atoms
    _, exclusions, _ = curate_structure(
        structure, None, _policy(quality_rules=QualityRules(max_atoms=700))
    )
    assert _codes(exclusions) == ["TOO_MANY_ATOMS"]
    curated, _, _ = curate_structure(
        structure, None, _policy(quality_rules=QualityRules(max_atoms=704))
    )
    assert curated is not None
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_curation.py -v`
Expected: new tests FAIL with pydantic `extra fields`/unexpected keyword errors or wrong codes. `test_cryo_em_entry_passes...` passes its first half already (Task 1) and fails the second.

- [ ] **Step 3: Extend the schema** (`pandora/schemas/dataset.py`)

Add to `ExclusionReason`: `"RFREE_THRESHOLD"`, `"RFREE_GAP_THRESHOLD"`, `"RSYM_THRESHOLD"`, `"NULL_RFACTOR"`, `"NONSTANDARD_RESIDUE"`, `"TOO_MANY_ATOMS"`.

Add to `QualityRules` after `exclude_experimental_methods`:

```python
    max_resolution_by_method: dict[str, float] = Field(default_factory=dict)
    rfactor_methods: list[str] = Field(
        default_factory=lambda: ["X-RAY DIFFRACTION"]
    )
    max_r_free: float | None = None
    max_r_free_gap: float | None = None
    max_r_sym: float | None = None
    null_rfactor_behavior: Literal["exclude", "include"] = "include"
    exclude_nonstandard_residues: bool = False
    allowed_nonstandard_residues: list[str] = Field(default_factory=list)
    max_atoms: int | None = None
```

Docstring lines:

```
        max_resolution_by_method: Resolution limits per experimental
            method, overriding `max_resolution` for entries with that
            method. With several matching methods, the strictest wins.
        rfactor_methods: The methods the R-factor rules apply to.
        max_r_free: The maximum Rfree.
        max_r_free_gap: The maximum |Rfree - Rwork|.
        max_r_sym: The maximum Rsym, or Rmerge when Rsym is missing.
        null_rfactor_behavior: Whether an entry missing a value an
            active R-factor rule needs is excluded or included.
        exclude_nonstandard_residues: Whether an entry with any
            non-standard polymer residue is excluded.
        allowed_nonstandard_residues: Non-standard residues that don't
            trigger `exclude_nonstandard_residues`.
        max_atoms: The maximum number of atoms an entry may have.
```

Update the class summary line to "Curation policy for resolution, R-factors, experimental method, residue content, size, and chain counts."

- [ ] **Step 4: Implement the rules** (`pandora/datasets/curation.py`)

Add imports: `from pandora.canonicalisation.modified_residues import STANDARD_RESIDUES` and `from pandora.schemas.metadata import MetadataRecord, QualityRecord`.

Add helpers above `_check_quality`:

```python
def _upper_set(values: list[str]) -> set[str]:
    """Stripped, upper-cased copies of values."""

    return {v.strip().upper() for v in values}


def _entry_methods(quality: QualityRecord | None) -> set[str]:
    """The entry's experimental methods, split on ';' and normalised."""

    method = quality.experimental_method if quality else None
    return _upper_set([m for m in (method or "").split(";") if m.strip()])


def _resolution_limit(methods: set[str], rules: QualityRules) -> float | None:
    """The strictest per-method limit that applies, else max_resolution."""

    by_method = {
        m.strip().upper(): limit
        for m, limit in rules.max_resolution_by_method.items()
    }
    limits = [by_method[m] for m in methods if m in by_method]
    return min(limits) if limits else rules.max_resolution


def _check_rfactors(
    entry_id: str,
    quality: QualityRecord | None,
    methods: set[str],
    rules: QualityRules,
) -> ExclusionRecord | None:
    """ExclusionRecord if an R-factor rule fails, for rfactor_methods only."""

    if not methods & _upper_set(rules.rfactor_methods):
        return None
    r_free = quality.r_free if quality else None
    r_work = quality.r_work if quality else None
    r_sym = None
    if quality:
        r_sym = quality.r_sym if quality.r_sym is not None else quality.r_merge
    gap = (
        abs(r_free - r_work)
        if r_free is not None and r_work is not None
        else None
    )
    checks = [
        (rules.max_r_free, r_free, "RFREE_THRESHOLD", "r_free"),
        (rules.max_r_free_gap, gap, "RFREE_GAP_THRESHOLD", "|r_free - r_work|"),
        (rules.max_r_sym, r_sym, "RSYM_THRESHOLD", "r_sym"),
    ]
    for limit, value, code, label in checks:
        if limit is None:
            continue
        if value is None:
            if rules.null_rfactor_behavior == "exclude":
                return ExclusionRecord(
                    entry_id=entry_id,
                    reason_code="NULL_RFACTOR",
                    message=f"{label} is null and "
                    "null_rfactor_behavior='exclude'",
                )
            continue
        if value > limit:
            return ExclusionRecord(
                entry_id=entry_id,
                reason_code=code,
                message=f"{label} {value:.3f} exceeds {limit}",
            )
    return None


def _check_composition(
    structure: Structure, rules: QualityRules
) -> ExclusionRecord | None:
    """ExclusionRecord if the entry has non-standard residues or too many
    atoms, else None."""

    if rules.exclude_nonstandard_residues:
        allowed = STANDARD_RESIDUES | set(rules.allowed_nonstandard_residues)
        found = sorted(
            {
                a.label_comp_id
                for a in structure.atoms
                if a.label_seq_id is not None and a.label_comp_id not in allowed
            }
        )
        if found:
            return ExclusionRecord(
                entry_id=structure.entry_id,
                reason_code="NONSTANDARD_RESIDUE",
                message=f"non-standard polymer residues {found}",
            )
    if rules.max_atoms is not None and len(structure.atoms) > rules.max_atoms:
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="TOO_MANY_ATOMS",
            message=f"{len(structure.atoms)} atoms exceed "
            f"max_atoms={rules.max_atoms}",
        )
    return None
```

In `_check_quality`, compute `methods = _entry_methods(quality)` and `limit = _resolution_limit(methods, rules)` right after `resolution = ...`, and replace every `rules.max_resolution` in the resolution block with `limit` (message: `f"resolution {resolution} exceeds limit {limit}"`). Leave the existing include/exclude method check unchanged. After that method check (before `min_chain_length`) add:

```python
    rfactor = _check_rfactors(structure.entry_id, quality, methods, rules)
    if rfactor is not None:
        return rfactor
```

In `curate_structure`:

```python
    exclusion = (
        _check_quality(structure, metadata, policy.quality_rules)
        or _check_organism(structure, metadata, policy.organism_rules)
        or _check_composition(structure, policy.quality_rules)
    )
```

- [ ] **Step 5: Run the tests**

Run: `uv run pytest tests/test_curation.py -v`
Expected: all PASS. The existing `test_curate_structure_pass_and_fail` asserts only the reason code, so the message change is fine; if any test asserts the old message text, update it.

- [ ] **Step 6: Lint and commit**

```bash
uv run ruff format . && uv run ruff check .
git add pandora/schemas/dataset.py pandora/datasets/curation.py tests/test_curation.py
git commit -m "Add per-method resolution, R-factor, residue and size rules

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Chain rules that drop chains

**Files:**
- Modify: `pandora/schemas/dataset.py` (`QualityRules`, `ExclusionReason`)
- Modify: `pandora/datasets/curation.py`
- Test: `tests/test_curation.py`

**Interfaces:**
- Consumes: `chain_completeness`, `MissingResidueDefinition` (Task 3); list return (Task 4).
- Produces: `QualityRules.missing_residue_definition`, `max_missing_tail_fraction`, `max_missing_middle_fraction`, `max_chain_length`; codes `MISSING_TAILS`, `MISSING_MIDDLE`, `CHAIN_TOO_LONG`, `NO_SEQRES`, `NO_CHAINS_LEFT`.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_curation.py`)

```python
def _chain_codes(exclusions):
    return [(e.chain_id, e.reason_code) for e in exclusions]


def test_missing_middle_drops_only_that_chain():
    structure = _load("1aui")  # A: middle 95/473, B: none
    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_missing_middle_fraction=0.1)),
    )
    assert _chain_codes(exclusions) == [("A", "MISSING_MIDDLE")]
    assert {a.label_asym_id for a in curated.atoms if a.label_seq_id} == {"B"}
    assert "A" not in {u.id for u in curated.asym_units}
    assert "A" in {u.id for u in structure.asym_units}  # input untouched


def test_tail_fraction_boundary():
    structure = _load("10mv")  # tails 39/275 = 0.1418
    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_missing_tail_fraction=0.15)),
    )
    assert curated is not None and exclusions == []
    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_missing_tail_fraction=0.14)),
    )
    assert curated is None
    assert _chain_codes(exclusions) == [
        ("A", "MISSING_TAILS"),
        (None, "NO_CHAINS_LEFT"),
    ]


def test_all_chains_dropped_excludes_entry():
    structure = _load("1aui")
    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_missing_tail_fraction=0.01)),
    )
    assert curated is None
    assert _chain_codes(exclusions) == [
        ("A", "MISSING_TAILS"),
        ("B", "MISSING_TAILS"),
        (None, "NO_CHAINS_LEFT"),
    ]


def test_chain_counts_checked_after_drop():
    structure = _load("1aui")
    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(
            quality_rules=QualityRules(
                max_missing_middle_fraction=0.1, min_polymer_chains=2
            )
        ),
    )
    assert curated is None
    assert _chain_codes(exclusions) == [
        ("A", "MISSING_MIDDLE"),
        (None, "TOO_FEW_CHAINS"),
    ]


def test_max_chain_length_drops_long_chain():
    structure = _load("1aui")  # A: 378 observed residues, B: 165
    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_chain_length=200)),
    )
    assert _chain_codes(exclusions) == [("A", "CHAIN_TOO_LONG")]
    assert curated is not None


def test_chain_failing_two_rules_gets_one_record():
    structure = _load("1aui")
    _, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(
            quality_rules=QualityRules(
                max_chain_length=200, max_missing_middle_fraction=0.1
            )
        ),
    )
    assert [e.chain_id for e in exclusions] == ["A"]


def test_chain_without_seqres_dropped_when_completeness_rule_active():
    structure = _load("1aui")
    raw = {k: v for k, v in structure.raw.items() if k != "_entity_poly_seq"}
    structure = structure.model_copy(update={"raw": raw})
    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_chain_length=1000)),
    )
    assert curated is not None and exclusions == []  # no completeness rule
    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_missing_tail_fraction=0.5)),
    )
    assert curated is None
    assert _chain_codes(exclusions) == [
        ("A", "NO_SEQRES"),
        ("B", "NO_SEQRES"),
        (None, "NO_CHAINS_LEFT"),
    ]


def test_middle_fraction_zero_denominator_is_zero():
    # 1p58 is a CA trace: under incomplete_backbone every residue is
    # missing, so all of SEQRES is "tail" and the middle denominator is 0.
    structure = _load("1p58")
    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_missing_middle_fraction=0.0)),
    )
    assert curated is not None and exclusions == []
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_curation.py -v`
Expected: new tests FAIL (unexpected keyword / no chain records).

- [ ] **Step 3: Extend the schema** (`pandora/schemas/dataset.py`)

Add to `ExclusionReason`: `"MISSING_TAILS"`, `"MISSING_MIDDLE"`, `"CHAIN_TOO_LONG"`, `"NO_SEQRES"`, `"NO_CHAINS_LEFT"`.

Add to `QualityRules` (after `max_atoms`). `MissingResidueDefinition` is defined earlier in the same module (Task 3):

```python
missing_residue_definition: MissingResidueDefinition = "incomplete_backbone"
max_missing_tail_fraction: float | None = None
max_missing_middle_fraction: float | None = None
max_chain_length: int | None = None
```

Docstring:

```
        missing_residue_definition: What counts as a missing residue
            for the missing-fraction rules (see `chain_completeness`).
        max_missing_tail_fraction: Chains whose missing N- plus
            C-terminal residues exceed this fraction of SEQRES are
            removed.
        max_missing_middle_fraction: Chains whose missing middle
            residues exceed this fraction of SEQRES minus the missing
            tails are removed.
        max_chain_length: Chains with more observed residues than this
            are removed.
```

Also update `min_chain_length` and `min_polymer_chains` docstrings to say they are checked after chains are removed.

- [ ] **Step 4: Implement** (`pandora/datasets/curation.py`)

Import `from pandora.datasets.completeness import chain_completeness`.

Move the `min_chain_length` and `min_polymer_chains` blocks out of `_check_quality` into a new helper, and add the chain helpers:

```python
def _chain_exclusions(
    structure: Structure, rules: QualityRules
) -> list[ExclusionRecord]:
    """One ExclusionRecord per polymer chain failing a chain rule."""

    failed: dict[str, ExclusionRecord] = {}

    def fail(chain_id: str, code: ExclusionReason, message: str) -> None:
        failed.setdefault(
            chain_id,
            ExclusionRecord(
                entry_id=structure.entry_id,
                chain_id=chain_id,
                reason_code=code,
                message=message,
            ),
        )

    if rules.max_chain_length is not None:
        for chain in extract_chain_records(structure):
            if chain.residue_count > rules.max_chain_length:
                fail(
                    chain.chain_id,
                    "CHAIN_TOO_LONG",
                    f"{chain.residue_count} residues exceed "
                    f"max_chain_length={rules.max_chain_length}",
                )

    max_tail = rules.max_missing_tail_fraction
    max_middle = rules.max_missing_middle_fraction
    if max_tail is not None or max_middle is not None:
        records, _ = chain_completeness(
            structure, rules.missing_residue_definition
        )
        by_chain = {r.chain_id: r for r in records}
        for chain_id in sorted(polymer_asym_ids(structure)):
            c = by_chain.get(chain_id)
            if c is None:
                fail(
                    chain_id,
                    "NO_SEQRES",
                    "completeness could not be measured against SEQRES",
                )
                continue
            tails = c.missing_n_term + c.missing_c_term
            tail_fraction = tails / c.seqres_length
            core = c.seqres_length - tails
            middle_fraction = c.missing_middle / core if core else 0.0
            if max_tail is not None and tail_fraction > max_tail:
                fail(
                    chain_id,
                    "MISSING_TAILS",
                    f"missing tails {tail_fraction:.3f} exceed "
                    f"max_missing_tail_fraction={max_tail}",
                )
            elif max_middle is not None and middle_fraction > max_middle:
                fail(
                    chain_id,
                    "MISSING_MIDDLE",
                    f"missing middle {middle_fraction:.3f} exceeds "
                    f"max_missing_middle_fraction={max_middle}",
                )

    return [failed[c] for c in sorted(failed)]


def _drop_chains(structure: Structure, chain_ids: set[str]) -> Structure:
    """Copy of structure without the atoms and asym units of chain_ids."""

    return structure.model_copy(
        update={
            "atoms": [
                a for a in structure.atoms if a.label_asym_id not in chain_ids
            ],
            "asym_units": [
                u for u in structure.asym_units if u.id not in chain_ids
            ],
        }
    )


def _check_chain_counts(
    structure: Structure, rules: QualityRules, chains_dropped: bool
) -> ExclusionRecord | None:
    """ExclusionRecord if too few or too short chains are left, else None."""

    if chains_dropped and not polymer_asym_ids(structure):
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="NO_CHAINS_LEFT",
            message="every polymer chain was removed by the chain rules",
        )
    # <move the existing min_chain_length block here unchanged>
    # <move the existing min_polymer_chains block here unchanged>
    return None
```

The two `# <move ...>` lines are instructions: cut the `if rules.min_chain_length is not None:` and `if rules.min_polymer_chains is not None:` blocks from `_check_quality` and paste them in their place, unchanged. Import `ExclusionReason` from `pandora.schemas.dataset`.

Sorting: chain records come out sorted by chain id, which is what the tests expect (`A` before `B`).

`curate_structure` body after provenance:

```python
    rules = policy.quality_rules
    exclusion = (
        _check_quality(structure, metadata, rules)
        or _check_organism(structure, metadata, policy.organism_rules)
        or _check_composition(structure, rules)
    )
    if exclusion is not None:
        return None, [exclusion], provenance

    records = _chain_exclusions(structure, rules)
    if records:
        structure = _drop_chains(structure, {r.chain_id for r in records})
    exclusion = _check_chain_counts(structure, rules, bool(records))
    if exclusion is not None:
        return None, [*records, exclusion], provenance

    return (
        _apply_content_rules(structure, policy.content_rules),
        records,
        provenance,
    )
```

Update `_check_quality`'s docstring to "resolution/method/R-factor rules" (chain-length checks moved out).

- [ ] **Step 5: Run the tests**

Run: `uv run pytest tests/test_curation.py tests/test_completeness.py -v`
Expected: all PASS.

- [ ] **Step 6: Full suite, lint, commit**

```bash
uv run pytest && uv run ruff format . && uv run ruff check .
git add pandora/schemas/dataset.py pandora/datasets/curation.py tests/test_curation.py
git commit -m "Drop chains with missing residues or too many residues

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Docs, policy files and diagrams

**Files:**
- Modify: `docs/reference/policies.md`, `docs/policies/dataset.yaml`, `docs/usage/datasets.md`, `docs/usage/canonicalisation.md`, `docs/contributing/architecture.md`
- Regenerate: `docs/assets/diagrams/*`

- [ ] **Step 1: Curation section in `docs/reference/policies.md`**

There is no curation section yet. Add one after the canonicalisation sections, before `## Similarity hit filter`:

```markdown
## Curation

`curate_structure(structure, metadata, policy)` applies a
`DatasetCurationPolicy`. Entry rules exclude the whole entry; chain
rules remove failing chains and keep the rest. An entry left with no
polymer chain is excluded as `NO_CHAINS_LEFT`.

Order: entry rules (resolution, method, R-factors, organism,
non-standard residues, `max_atoms`), then chain rules, then
`min_chain_length` and `min_polymer_chains` on what is left, then
content rules.

### `quality_rules`

| Field | Default | Level | Description |
|---|---|---|---|
| `max_resolution` | `null` | entry | Maximum resolution (Å). Cryo-EM resolution comes from `_em_3d_reconstruction`. |
| `max_resolution_by_method` | `{}` | entry | Per-method limits, e.g. `{"X-RAY DIFFRACTION": 2.5, "ELECTRON MICROSCOPY": 2.0}`. Overrides `max_resolution` for entries with that method; with several matching methods the strictest wins. Names match case- and whitespace-insensitively. |
| `null_resolution_behavior` | `exclude` | entry | What happens to an entry with no resolution when a limit applies. |
| `include_experimental_methods` / `exclude_experimental_methods` | `[]` | entry | Method allow/deny lists. |
| `rfactor_methods` | `["X-RAY DIFFRACTION"]` | entry | Methods the R-factor rules apply to. Other methods (cryo-EM, NMR) are never checked, since they have no R-factors. An entry with no method recorded is not checked. |
| `max_r_free` | `null` | entry | Maximum Rfree. |
| `max_r_free_gap` | `null` | entry | Maximum \|Rfree − Rwork\|. |
| `max_r_sym` | `null` | entry | Maximum Rsym (`_reflns.pdbx_Rsym_value`), falling back to Rmerge (`_reflns.pdbx_Rmerge_I_obs`), the same statistic under another name. |
| `null_rfactor_behavior` | `include` | entry | What happens when an active R-factor rule's value is missing. Defaults to `include` because many X-ray entries lack Rsym. |
| `exclude_nonstandard_residues` | `false` | entry | Exclude an entry with any polymer residue that isn't one of the 20 amino acids, `UNK`, or a standard nucleotide. Runs on the canonical structure, so MSE mapped by `modified_residue_rules` counts as MET. |
| `allowed_nonstandard_residues` | `[]` | entry | Residues that don't trigger the rule above. |
| `max_atoms` | `null` | entry | Maximum atoms in the entry. Checked after parsing, so it protects later steps (contacts, export) but not the parse itself. |
| `missing_residue_definition` | `incomplete_backbone` | chain | `unobserved`: a residue is missing when it has no atoms. `incomplete_backbone`: also when any backbone atom (N, CA, C, O) is missing. |
| `max_missing_tail_fraction` | `null` | chain | Remove a chain whose missing N- plus C-terminal residues exceed this fraction of its SEQRES length. |
| `max_missing_middle_fraction` | `null` | chain | Remove a chain whose missing residues between its first and last present residue exceed this fraction of SEQRES length minus the missing tails. |
| `max_chain_length` | `null` | chain | Remove a chain with more observed residues than this. |
| `min_chain_length` | `null` | entry | Exclude the entry if no remaining chain reaches this many residues. |
| `min_polymer_chains` | `null` | entry | Exclude the entry if fewer polymer chains remain. |

When a missing-fraction rule is set, a chain whose completeness can't be
measured (no `_entity_poly_seq`, or residues that don't line up with it,
for example after renumbering) is removed as `NO_SEQRES`. Use
`chain_completeness()` to inspect the numbers directly.
```

- [ ] **Step 2: `docs/policies/dataset.yaml`**

Under `quality_rules:` after `min_chain_length: int | null`, add:

```yaml
    min_polymer_chains: int | null
    max_resolution_by_method: dict[string, float]
    rfactor_methods: list[string]
    max_r_free: float | null
    max_r_free_gap: float | null
    max_r_sym: float | null
    null_rfactor_behavior: string
    # include (default) | exclude
    exclude_nonstandard_residues: bool
    allowed_nonstandard_residues: list[string]
    max_atoms: int | null
    missing_residue_definition: string
    # incomplete_backbone (default) | unobserved
    max_missing_tail_fraction: float | null
    max_missing_middle_fraction: float | null
    max_chain_length: int | null
```

- [ ] **Step 3: Usage pages**

`docs/usage/datasets.md`: add a `## Chain completeness` section after "Curate one structure":

````markdown
## Chain completeness

`chain_completeness()` counts each polymer chain's missing residues
against SEQRES (`_entity_poly_seq`): missing at the N-terminus, at the
C-terminus, and in the middle. The missing-fraction curation rules use
it; you can also call it directly, for example to pick the chain with
the fewest missing residues.

```python
from pandora.datasets import chain_completeness

records, diagnostics = chain_completeness(canonical)
for r in records:
    print(r.chain_id, r.seqres_length, r.missing_n_term,
          r.missing_c_term, r.missing_middle)
```

Run it on a structure whose `label_seq_id`s still index SEQRES. After
`residue_numbering: renumber` they don't, and the chain is reported as
`SEQRES_MISMATCH` instead of measured.
````

`docs/usage/canonicalisation.md`: add a short `## Modified residues` section pointing to `modified_residue_rules` in the policy reference, with this example:

```python
from pandora.schemas.canonicalisation import ModifiedResidueRules

policy = canonicalisationPolicy(
    policy_id="p",
    policy_name="p",
    policy_version="1.0.0",
    modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
)  # MSE -> MET; mappings in CanonicalMappings.modified_residue_mapping
```

`docs/contributing/architecture.md`: add a row to the datasets table after `curate_structure()`:

```markdown
| `chain_completeness()` | `Structure`, `MissingResidueDefinition` | `(list[ChainCompleteness], DiagnosticBundle)` |
```

and in the prose at line ~140, say curation covers resolution per method, R-factors, residue content, size and missing residues, and removes failing chains.

- [ ] **Step 4: Regenerate the diagrams**

Run: `uv run --extra docs python docs/scripts/generate_erd.py`
Expected: updated SVGs under `docs/assets/diagrams/` (curation policy, canonicalisation policy/mappings).

- [ ] **Step 5: Build the docs**

Use the `docs-build` skill. Run: `uv run zensical build --clean`
Expected: build succeeds with no new broken-link or anchor warnings.

- [ ] **Step 6: Final checks and commit**

```bash
uv run pytest && uv run ruff format --check . && uv run ruff check .
git add docs
git commit -m "Document quality, completeness and modified-residue rules

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
