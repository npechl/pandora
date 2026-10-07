# Contacts: fast, configurable chain-chain interfaces

Date: 2026-10-07. Status: approved design, awaiting spec review.

## Roadmap context

This is sub-project 1 of the plan to let Pandora rebuild the four
benchmark specs in `examples/ppi/benchmarks/`. The plan rests on these
decisions:

- **The four benchmarks are equal peers.** PPI v2.1, Pinder, ATOM3D PIP
  and ProteinFlow have the same weight. Gaps are ranked by how many
  benchmarks they block.
- **Pandora builds datasets.** Auditing, comparing and diffing datasets
  are out of scope for now.
- **Scale target: the whole PDB.** Pandora stays a library of pure
  functions with no batch runner. Users parallelise on their own
  infrastructure, and the examples show a `ProcessPoolExecutor` pattern
  that skips entries already done.
- **Spike 0 result: keep the data model.** Parsing takes 21.5 s and
  6.5 GB on 3j3q (2.44M atoms), about 3× gemmi's memory. The data model
  stays. The scaling work goes into contacts (this spec), similarity
  results as tables (sub-project 2), and a maximum-size exclusion field
  (sub-project 3).
- **External tools are not wrapped.** Pandora owns the data contract:
  it imports tool outputs and maps them to its own residue and chain
  IDs. MMseqs2 and Foldseek keep their thin wrappers and will also
  accept precomputed results.
- **There is no PDB snapshot.** No snapshot can be retrieved. The
  per-entry `retrieved_at` and `revision_date` fields already record
  what was used. Seeds become policy fields wherever randomness is
  added (sub-projects 6 and 7).

Sub-project order:

0. Spike (done)
1. **Contacts (this spec)**
2. Similarity results as tables
3. Quality and completeness
4. External tool outputs
5. Pair layer
6. Splitting
7. Negatives
8. Example scripts (ongoing)

## Problem

`annotate_chain_interfaces` rebuilds a spatial grid for every chain pair
and has no early exit for chains that are far apart. Its cost therefore
grows with the square of the chain count:

| Entry | Chains | Time |
|---|---|---|
| 2ms2 (expanded) | 180 | 47 s |
| 4ug0 | 81 | 29 s |
| 3j3q | 1,356 | about an hour (estimated) |

There are four further problems:

- It counts hydrogens. Structures that model H get more contacts than
  those that don't.
- The atom set can't be chosen. The benchmarks use different ones:
  Pinder uses backbone atoms at 10 Å, ATOM3D uses Cα at 8 Å.
- Chains can't be restricted by polymer type.
- It reports the residues on each side, but not which residue touches
  which. ATOM3D's labels and PPI v2.1's Jaccard need residue pairs.

## Design

### API

```python
annotate_chain_interfaces(
    structure: Structure,
    distance_cutoff: float = 4.0,
    atom_set: ContactAtomSet = "heavy",
    polymer_types: list[str] | None = None,
) -> AnnotationLayer

extract_interface_records(
    structure: Structure,
    distance_cutoff: float = 4.0,
    atom_set: ContactAtomSet = "heavy",
    polymer_types: list[str] | None = None,
) -> list[InterfaceRecord]
```

`ContactAtomSet = Literal["all", "heavy", "backbone", "ca"]`, defined
in `pandora/schemas/annotation.py` (the same pattern as
`ExclusionReason`).

| `atom_set` | Atoms used | Benchmark |
|---|---|---|
| `all` | Every atom, H included | Previous behaviour |
| `heavy` (default) | All atoms except element H or D (`type_symbol`) | ATOM3D heavy-atom option, PRODIGY |
| `backbone` | Protein: N, CA, C, O. Nucleic acid: O5′, C5′, C4′, C3′, O3′ | Pinder 10 Å pair gate |
| `ca` | Atoms named CA | ATOM3D 8 Å Cα default |

How the atom sets are selected:

- **`backbone`** reuses `_backbone_by_entity`, `_PROTEIN_BACKBONE` and
  `_NUCLEIC_BACKBONE` from `canonicalisation/missing_data.py`, so the
  definition isn't duplicated.
- **`ca`** has two edge cases:
  - Calcium ions are also named CA, but they're excluded because only
    polymer chains are considered.
  - Nucleic-acid chains have no CA atoms, so they produce no contacts.

**`polymer_types`:** chains whose entity's `poly.type` is not in the
list are skipped before any contact work. `None` keeps every polymer
chain, as before.

### Output

Each interface in `layer.data["interfaces"]` keeps its existing keys
(`chain_id_1`, `chain_id_2`, `interface_residues_chain_1`,
`interface_residues_chain_2`, `contact_count`) and gains:

- `residue_pairs`: a sorted list of `[residue_id_1, residue_id_2]`
  pairs (e.g. `["A:12", "B:40"]`), one per residue pair with at least
  one selected atom pair within the cutoff.

Layer-level changes:

- **`parameters`** records `distance_cutoff`, `atom_set` and
  `polymer_types`. `reproduce_dataset` calls
  `fn(structure, **layer.parameters)`, so every setting must be a
  JSON-serialisable keyword argument.
- **`method`** becomes `pandora.basic.distance_cutoff_contacts.v2`.

`InterfaceRecord` gains:

- `atom_set: ContactAtomSet`, defaulting to `"all"` so records exported before this change still load (they were all-atom contacts)
- `residue_pairs: list[tuple[str, str]]`, defaulting to an empty list

### Algorithm

1. Select the polymer chains, filtered by `polymer_types`, and their
   atoms, filtered by `atom_set`.
2. For each chain, compute its axis-aligned bounding box and one spatial
   grid with cell size equal to the cutoff, exactly as the existing grid
   does.
3. For each chain pair, skip it if the two boxes are more than the
   cutoff apart on any axis.
4. Otherwise, check only the atoms of chain 1 that lie inside chain 2's
   box grown by the cutoff, against chain 2's grid, in the 27
   neighbouring cells. Record the residue pair for every atom pair
   within the cutoff.
5. Derive the per-side residue lists and `contact_count` from the
   residue pairs.

The box check and the box-restricted atom filter are exact, not
approximate: an atom outside the grown box cannot be within the cutoff
of any atom of the other chain.

**Prototype results** (scratchpad only, not repo code):

| Entry | Before | After |
|---|---|---|
| 1aon | 2.1 s | 0.1 s |
| 2ms2 (expanded) | 47.3 s | 0.9 s |
| 4ug0 | 29.3 s | 1.6 s |
| 3j3q | about 1 h (estimated) | 14.8 s |

With `atom_set="all"` the prototype's output was identical to the
current code on every entry tested.

### Behaviour change

The default changes from all atoms to `heavy`. Pandora is at 0.5.6, so
this is acceptable with a documented note. A manifest written before
this change records only `distance_cutoff`, so `reproduce_dataset`
would rebuild it with `heavy`. `reproduce_dataset` is already
documented as best-effort; the usage docs will state this case.

### Pinder's pair gate

There is no dedicated gate function. A gate is two calls followed by
an intersection of the chain pairs they return:

```python
gate = annotate_chain_interfaces(s, distance_cutoff=10.0, atom_set="backbone")
iface = annotate_chain_interfaces(s, distance_cutoff=4.0)
```

The docs show this pattern.

## Out of scope

- `annotate_ligand_contacts` keeps its all-atom behaviour (it still
  counts H). It can be aligned later.
- No new CLI flags. The CLI uses the defaults.
- No new dependencies. The implementation is pure Python.

## Testing

All tests go in `tests/test_annotations.py`, use local fixtures only,
and are plain pytest functions.

1. **Reference equivalence.** A check-every-atom-pair reference helper
   in the test file must match the function exactly for all four atom
   sets. Inputs: `1a02` (protein and nucleic acid, multi-chain), and a
   synthetic structure with chains far apart to exercise the box skip.
2. **`heavy` vs `all`** on `1a7f` (two chains, 7,300 H atoms):
   - `heavy` has no pair that exists only through an H atom;
   - `all` finds a superset of `heavy`'s pairs.
3. **`backbone` and `ca`.** Synthetic cases where only a side-chain
   atom, or only a non-CA backbone atom, is within the cutoff.
4. **`polymer_types`.** On `1a02` with `["polypeptide(L)"]`, no
   interface involves a nucleic-acid chain.
5. **Consistency.** The sides collected from `residue_pairs` equal the
   per-side lists.
6. **Rebuild.** `annotate_chain_interfaces(s, **layer.parameters)`
   equals the original layer.
7. **`extract_interface_records`.** It passes the arguments through,
   and the records carry `atom_set` and `residue_pairs`.
8. **Existing test.** Update the assertion
   `parameters == {"distance_cutoff": 4.0}` in the existing test.

Speed is not asserted in pytest. The spike script is re-run on
`datasets/large/` and the before/after numbers go in the PR description.

## Docs

- **`docs/usage/annotation.md`:**
  - the new arguments;
  - a table of atom sets with the benchmark that uses each;
  - a warning about the default change and what it means for older
    manifests;
  - the pair-gate recipe.
- **`docs/usage/datasets.md`:** the new `InterfaceRecord` fields.
- **`docs/contributing/architecture.md`:** the updated signatures.
- **ERD diagrams:** regenerate with
  `uv run --extra docs python docs/scripts/generate_erd.py`.
- **Docs site:** build it with the `docs-build` skill.

## Done when

- The CLAUDE.md review checklist passes: `uv run pytest`,
  `ruff format --check` and `ruff check` pass, no input is mutated,
  public functions keep their `__all__` entries (`ContactAtomSet` is
  imported from `pandora.schemas.annotation` directly, like other
  schema types), docstrings are complete, and the docs build.
- On the large entries, the output under `atom_set="all"` matches the
  pre-change output.
