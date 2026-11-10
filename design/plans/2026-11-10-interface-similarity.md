# Interface Similarity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `pandora similarity interfaces` writes an interface-level
`similarity` layer scoring every pair of same-protein-pair interfaces by
per-side Jaccard of their interface residues (issue #134).

**Architecture:** One pure function, `interface_similarity()`, in a new
module `pandora/similarity/interfaces.py`: groups interfaces by their
chain-cluster pair, aligns chain sequences with gemmi, gates on pairwise
identity, scores both chain orders, returns rows plus counts. A thin CLI
handler in `pandora/cli/similarity.py` finds the structures folder the
contact layer came from (the existing `_layer_origin`, moved to
`cli/_io.py`), reads its `chains` table and the `clusters.json`, and
writes a similarity layer folder that `pandora cluster --items interface`
reads unchanged.

**Tech Stack:** Python 3, pydantic v2, gemmi 0.7.5, pyarrow (export
extra, for `layer.parquet`), argparse, pytest.

**Spec:** `design/specs/2026-11-10-interface-similarity-design.md`

## Global Constraints

- `from __future__ import annotations`; type hints; line length 80;
  `uv run ruff format`. Public functions: Google-style docstring,
  exported in `pandora/similarity/__init__.py` `__all__`. Never mutate
  inputs. No new dependency, no new schema model, no policy file.
- Interfaces = exactly `interface_items` (`pandora/similarity/layer.py`):
  rows with both `seq_id`s set and `chain_id_1 != chain_id_2`. Other rows
  are ignored and counted in `rows_ignored` (the spec's
  `non_polymer_rows`, renamed because same-chain rows land there too).
- Item ids from `item_id("interface", entry_id, a, b)` (chains sorted,
  entry id as given in the row). A chain's cluster is looked up by
  `(entry_id.upper(), chain_id)`, parsing cluster components with
  `parse_item_id("chain", ...)`.
- Side sets keyed by `seq_id` only. Sequences keyed by
  `(entry_id.upper(), chain_id)`; position `i` (1-based) = `seq_id i`.
- Sequence guard per interface: for every contact row of the interface
  and each side with `comp_id` set, the expected letter is
  `gemmi.find_tabulated_residue(comp_id).one_letter_code.upper()`
  (`"X"` if not found or blank); mismatch with `sequence[seq_id - 1]`,
  or `seq_id` outside `1..len(sequence)`, skips the interface
  (`interfaces_sequence_mismatch`). `comp_id` None is not checked.
- Alignment: `gemmi.align_string_sequences(list(q), list(t), [])`,
  default scoring. Cigar walk: `M` = aligned column (map `i -> j`, count
  identical), `I` = residue in `q` only, `D` = residue in `t` only.
  Identity = identical / `M` columns (0 if none). Identical sequences:
  no alignment, identity 1.0, map by position. Cache per ordered
  `(q, t)`; `alignments` counts cache misses.
- Identity gate: an order counts only if both its chain pairings have
  identity `>= min_identity` (default 0.9).
- Jaccard for side `Sx` of I against `Sy` of J:
  `|M(Sx) ∩ Sy| / (|Sx| + |Sy| − |M(Sx) ∩ Sy|)`; empty → 0.
- `score = min(jaccard_1, jaccard_2)`; keep the order with the higher
  score, `direct` on a tie.
- Row columns, in this order: `item_id_1, item_id_2, score, jaccard_1,
  jaccard_2, chain_order, identity_1, identity_2`, with
  `item_id_1 < item_id_2`; `jaccard_1`/`identity_1` belong to the
  smaller chain id of `item_id_1`. Sorted by `(item_id_1, item_id_2)`.
- CLI: `pandora similarity interfaces --contacts DIR (dest contacts_dir,
  required) --chain-clusters FILE (required) --output-dir DIR (required)
  [--min-identity FLOAT, default 0.9] [--overwrite]`. Output is always
  `layer.parquet` (no `--table-format`: `cluster` reads only parquet;
  this overrides the spec's sketch).
- `layer.json`: step `similarity_interfaces`; `LayerDescription(
  name="iface_jaccard", template="similarity", items="interface",
  tool="pandora", tool_version=__version__, parameters={"contacts_layer",
  "structures_dir", "chain_clusters", "min_identity"})`;
  `previous_summary = <contacts>/layer.json`; counts (sorted keys):
  `rows_in, rows_ignored, interfaces, interfaces_without_cluster,
  interfaces_without_sequence, interfaces_sequence_mismatch,
  pairs_compared, pairs_below_identity, pairs_swapped, pairs_zero_score,
  alignments`. `pairs_compared` = rows written.
- Printed line: `scored: {pairs_compared} interface pairs
  ({pairs_zero_score} zero, {pairs_below_identity} below identity) from
  {interfaces} interfaces -> {out}`.
- The all-vs-all loop inside each group carries
  `# ponytail: all-vs-all per cluster pair is quadratic; stream rows in batches if a group gets huge.`

## Review Focus

1. Interface given only as B–A rows (chain_id_1 > chain_id_2) → one
   item `E_A_B`, sides assigned to the right chains. Test in Task 1.
2. Homodimer pair where both orders fit clusters and identity → the
   better order wins; tie → `direct`; same-entry pairs (`E_A_B` vs
   `E_C_D`) compared. Tests in Task 1.
3. Selenomethionine (`MSE`) residue against an `M` in the sequence →
   passes the guard (gemmi gives lowercase `m`). Test in Task 1.
4. Contact layer whose origin can't be traced, or whose structures
   folder has no `chains` table → exits before writing, same messages
   as `label-interfaces`. Test in Task 2.
5. No comparable pairs at all → command succeeds, writes `layer.json`
   and no `layer.parquet`; `cluster --items interface` on it gives one
   singleton cluster per interface. Test in Task 2.

---

### Task 1: `interface_similarity` library function

**Files:**
- Create: `pandora/similarity/interfaces.py`
- Modify: `pandora/similarity/__init__.py` (import + `__all__`)
- Test: `tests/test_interface_similarity.py`

**Interfaces:**
- Consumes: `item_id`, `parse_item_id` (`pandora/similarity/layer.py`);
  `ResidueContactRow`, `LayerRow` (`pandora/schemas/annotation.py`);
  `SimilarityCluster` (`pandora/schemas/similarity.py`).
- Produces:
  ```python
  def interface_similarity(
      contacts: Iterable[ResidueContactRow],
      chain_clusters: Iterable[SimilarityCluster],
      sequences: Mapping[tuple[str, str], str],
      min_identity: float = 0.9,
  ) -> tuple[list[LayerRow], dict[str, int]]
  ```
  Rows per Global Constraints; counts dict with every key listed there
  except `rows_in` (the CLI adds it).

Test helpers (in the test file): `_row(entry, a, sa, b, sb, ca=None,
cb=None)` builds a `ResidueContactRow` (`comp_id` None = guard not
checked; only the guard test sets it); `_cl(*chain_items)` builds a
`SimilarityCluster`. Sequences default to a non-repetitive
`_SEQ = ("ACDEFGHIKLMNPQRSTVWY" * 3)[:50]` so alignments are
unambiguous.

- [ ] **Step 1: Write the failing tests**

```python
def test_identical_heterodimers_score_one():
    # 1AAA_A_B and 2BBB_A_B: same residues per side, A~A, B~B clusters
    rows, counts = interface_similarity(rows_1 + rows_2, clusters, seqs)
    assert [(r.item_id_1, r.item_id_2, r.score, r.chain_order) for r in rows] \
        == [("1AAA_A_B", "2BBB_A_B", 1.0, "direct")]
    assert rows[0].identity_1 == rows[0].identity_2 == 1.0
    assert counts["pairs_compared"] == 1 and counts["alignments"] == 0

def test_score_is_min_of_sides():
    # side A identical ({1,2,3}), side B disjoint ({1,2} vs {10,11})
    assert rows[0].jaccard_1 == 1.0 and rows[0].jaccard_2 == 0.0
    assert rows[0].score == 0.0 and counts["pairs_zero_score"] == 1

def test_homodimer_keeps_swapped_order():
    # all four chains in one cluster; 1AAA A:{1,2} B:{20,21};
    # 2BBB A:{20,21} B:{1,2}  -> direct scores 0, swapped scores 1
    assert rows[0].chain_order == "swapped" and rows[0].score == 1.0
    assert counts["pairs_swapped"] == 1

def test_homodimer_tie_is_direct():
    # both orders score 1.0 (each side {1,2})
    assert rows[0].chain_order == "direct"

def test_same_entry_interfaces_compared():
    # 1AAA_A_B and 1AAA_C_D, A,C one cluster; B,D another
    assert [(r.item_id_1, r.item_id_2) for r in rows] == [("1AAA_A_B", "1AAA_C_D")]

def test_reversed_rows_name_one_interface():
    # rows given as (B, A); sides still land on A and B
    assert rows[0].item_id_1 == "1AAA_A_B" and rows[0].jaccard_1 == 1.0
    assert counts["interfaces"] == 2

def test_insertion_maps_through_alignment():
    # 2BBB chain A = 1AAA chain A with one residue inserted after pos 3
    # (len 50 vs 51, identity >= 0.9); 1AAA A:{5,6}, 2BBB A:{6,7}
    assert rows[0].jaccard_1 == 1.0 and counts["alignments"] >= 1
    assert rows[0].identity_1 == 1.0  # identical M columns / M columns

def test_unplaced_residue_counts_against_overlap():
    # residue of Sx lands on a gap: |Sx|=2, |Sy|=1, overlap 1 -> 0.5
    assert rows[0].jaccard_1 == 0.5

def test_below_identity_writes_no_row():
    # same cluster, chain A sequences differ in 10 of 50 positions (80%)
    assert rows == [] and counts["pairs_below_identity"] == 1

def test_min_identity_zero_lets_it_through():
    rows, _ = interface_similarity(..., min_identity=0.0)
    assert len(rows) == 1

def test_different_clusters_not_compared_and_no_cluster_counted():
    # 1AAA_A_B (A~X, B~Y) vs 2BBB_A_B (A~X, B~Z); 3CCC_A_B chains in no cluster
    assert rows == [] and counts["interfaces_without_cluster"] == 1

def test_sequence_guard():
    # comp_id GLY at a position whose sequence letter is A -> skipped
    assert counts["interfaces_sequence_mismatch"] == 1
    # comp_id MSE against "M" passes; seq_id 0 or > len fails
    # missing sequence -> interfaces_without_sequence

def test_ignored_rows():
    # seq_id None row and an A-A row
    assert counts["rows_ignored"] == 2 and counts["interfaces"] == 0

def test_inputs_not_mutated():
    # deep-copy contacts/clusters/sequences before; equal after
```

Each test builds its own rows; the comments give the exact setup and
the asserts are the required values.

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_interface_similarity.py -v`
Expected: FAIL, `ImportError: cannot import name 'interface_similarity'`

- [ ] **Step 3: Implement `interface_similarity` in `pandora/similarity/interfaces.py`**

Order of work: collect side sets and per-interface rows (for the guard)
→ drop interfaces without both clusters / without both sequences /
failing the guard, counting each once in that order → group by
`frozenset({cl(a), cl(b)})` (a homodimer's group key is a one-element
set) → `itertools.combinations` over each group's sorted item ids → for
each order that fits on clusters, compute identities (gate) then
Jaccards → keep best order or count `pairs_below_identity`. Internal
helpers `_align(q, t) -> tuple[dict[int, int], float]` (cached) and
`_jaccard(sx, mapping, sy) -> float`. Full Google docstring on the
public function; `ponytail:` comment from Global Constraints on the
loop. Export in `__init__.py`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_interface_similarity.py -v`
Expected: all PASS

- [ ] **Step 5: Lint, format, commit**

```bash
uv run ruff format . && uv run ruff check .
git add pandora/similarity/interfaces.py pandora/similarity/__init__.py tests/test_interface_similarity.py
git commit -m "feat(similarity): interface_similarity scores interface pairs by per-side Jaccard"
```

### Task 2: `pandora similarity interfaces` command

**Files:**
- Modify: `pandora/cli/_io.py` (add `layer_origin`, moved verbatim from
  `annotate._layer_origin`)
- Modify: `pandora/cli/annotate.py` (use `layer_origin` from `_io`,
  delete `_layer_origin`)
- Modify: `pandora/cli/similarity.py` (`_cmd_interfaces` + parser under
  the `similarity` subparsers, next to `search`/`import`)
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `interface_similarity` (Task 1); `check_contact_layer`,
  `find_table`, `load_json_list`, `prepare_output_dir` (`cli/_io.py`);
  `read_records`, `write_records` (`pandora.export`);
  `DatasetChainRecord`; `_write_summary`, `_require_pyarrow`,
  `_LAYER_OUTPUTS` (already in `cli/similarity.py`).
- Produces: `layer_origin(folder: Path) -> Path | None` in `cli/_io.py`;
  the command and its output folder (Global Constraints).

Test helpers: reuse `_label_inputs(tmp_path)` (1a2c step folder +
`chain_interfaces` layer) but write a `sequence` per chain into
`chains.jsonl`, built from the contact rows so the guard passes: length
= max `seq_id` seen for that chain, one-letter of `comp_id` at each
seen position, `"A"` elsewhere. One `clusters.json` putting every chain
of 1A2C in one cluster.

- [ ] **Step 1: Write the failing tests**

```python
def test_similarity_interfaces_writes_layer(tmp_path):
    # run with --min-identity 0 so 1a2c's different chains get compared
    written = _parquet_rows(out / "layer.parquet")
    expected, counts = interface_similarity(rows, clusters, seqs, 0.0)
    assert written and written == [r.model_dump() for r in expected]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "similarity_interfaces"
    assert summary["policy"]["items"] == "interface"
    assert summary["policy"]["parameters"] == {
        "contacts_layer": "chain_interfaces",
        "structures_dir": str(step.resolve()),
        "chain_clusters": str(clusters_file),
        "min_identity": 0.0,
    }
    assert summary["counts"] == dict(sorted(({"rows_in": len(rows)} | counts).items()))
    assert summary["previous_summary"] == str(contacts / "layer.json")

def test_similarity_interfaces_feeds_cluster(tmp_path):
    # cluster --items interface --contacts <contacts> --min-score 0
    # -> one cluster holding every interface of 1A2C

def test_similarity_interfaces_no_pairs(tmp_path):
    # clusters.json with one cluster per chain -> no row
    assert not (out / "layer.parquet").exists()
    # cluster on it -> one singleton cluster per interface

def test_similarity_interfaces_untraceable_origin_exits(tmp_path):
    # layer.json "input" -> missing dir; SystemExit, nothing written

def test_similarity_interfaces_no_chains_table_exits(tmp_path):
    # delete chains.jsonl; SystemExit matching "pandora curate"

def test_similarity_interfaces_rejects_output_is_contacts(tmp_path):
    # --output-dir == --contacts -> SystemExit
```

The existing `label-interfaces` and `filter-contacts` tests stay green
after the `layer_origin` move.

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_cli.py -k similarity_interfaces -v`
Expected: FAIL (argparse: invalid choice `interfaces`)

- [ ] **Step 3: Move `_layer_origin` to `cli/_io.py` as `layer_origin`; implement `_cmd_interfaces(args: argparse.Namespace) -> None`**

Mirror `annotate._cmd_label`: `check_contact_layer`, output ≠ contacts
check, `layer_origin` + `find_table(origin, "chains")` with the same
error wording (command name `similarity interfaces`). Sequences =
`{(c.entry_id.upper(), c.chain_id): c.sequence}` for chains with a
sequence. Clusters via `load_json_list(SimilarityCluster, path)`.
`prepare_output_dir(out, _LAYER_OUTPUTS, args.overwrite)`; write
`layer.parquet` with `write_records` only when rows exist;
`_write_summary(...)` with the description and counts from Global
Constraints; print the line from Global Constraints. Call
`_require_pyarrow("similarity interfaces")` first.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_cli.py -v`
Expected: all PASS

- [ ] **Step 5: Lint, format, commit**

```bash
uv run ruff format . && uv run ruff check .
git add pandora/cli/_io.py pandora/cli/annotate.py pandora/cli/similarity.py tests/test_cli.py
git commit -m "feat(cli): similarity interfaces writes an interface similarity layer"
```

### Task 3: Docs

**Files:**
- Modify: `docs/usage/similarity.md` (new section "Comparing interfaces")
- Modify: `docs/usage/cli.md` (command entry)
- Modify: `docs/contributing/architecture.md` (function → models row)

- [ ] **Step 1: Write the section**

Explain in plain words: same two proteins (clusters + ≥ 90% pairwise
identity), both chain orders, per-side Jaccard, score = smaller side,
every compared pair written so the cut-off can be chosen from the
distribution, `--min-score` is inclusive. State that structures must
keep `label_seq_id` (`residue_numbering: preserve`) or interfaces are
skipped as `interfaces_sequence_mismatch`. Library tab: full inline
Python calling `interface_similarity` on the 1a2c contact rows; CLI
tab: `annotate compute` → `cluster` (seq clusters) →
`similarity interfaces` → `cluster --items interface --min-score 0.1`
→ `jq -c '.counts'`. The sequence clusters come from a `clusters.json`
written inline with a heredoc (MMseqs2 is not needed to run the
example); one sentence says it normally comes from `similarity search`
+ `cluster` with a hit filter of identity and both-side coverage 0.9. Run both by hand and paste the real output into
` ```text ` blocks.

- [ ] **Step 2: Build the docs**

Run: `uv run zensical build --clean` (use the docs-build skill)
Expected: builds with no new warnings; rendered page contains the
example code, not only the outputs.

- [ ] **Step 3: Full check and commit**

```bash
uv run pytest && uv run ruff format --check .
git add docs/
git commit -m "docs: compare interfaces by per-side Jaccard"
```
