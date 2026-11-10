# Containment dedup and cluster representatives Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `pandora dedup --contained` (a `work/uncontained` step
folder) and `pandora representatives` (one representative chain per
cluster) for PPI v2.1 Step 1 filter 1 and Step 2.4 (#138).

**Architecture:** There are two plain library functions on plain data:
`find_contained` in `pandora/datasets/containment.py` and
`pick_representatives` in `pandora/similarity/representatives.py`. Both
share one ranking helper and a `RankKey` model. The CLI handlers read the
step folder's tables, `clusters.json` and, optionally, similarity layers,
then call the library and write the outputs. The `representatives`
handler lives in `pandora/cli/similarity.py` so it can reuse the
hit-filter and layer-edge helpers that `cluster` already uses.

**Tech Stack:** Python 3, pydantic v2, argparse, pyarrow (Parquet
tables), pytest, Ruff, Zensical.

**Spec:** `docs/superpowers/specs/2026-11-10-contained-dedup-representatives-design.md`

## Global Constraints

- Every module starts with `from __future__ import annotations`. Line
  length is 80, formatted with `ruff format`.
- Public functions get type hints, a Google-style docstring (`Args:`,
  `Returns:`, `Raises:`), and an entry in the package `__all__`.
  Internal helpers start with `_`.
- Never mutate inputs. Return new objects, using `model_copy` for
  pydantic.
- New models go in `pandora/schemas/` and have no behaviour.
- Catch named exceptions only (Ruff BLE001 and S110).
- No new dependencies.
- Tests use local data only: hand-built tables, `clusters.json` files and
  similarity layers, written as JSONL or Parquet under `tmp_path`.
- Item ids come from `pandora.similarity.item_id("chain", entry_id,
  chain_id)` (`"{ENTRY_UPPER}_{chain}"`). Never parse ids.
- The final tie-break is always the lowest item id. It isn't
  configurable. A null rank value ranks last for both `min` and `max`.
- Missing residues are an **absolute** count: `missing_n_term +
  missing_c_term + missing_middle`. If any part is null, the total is
  null.
- Containment **counts chains per cluster**, symmetry copies included.
- `n_hits` counts distinct partners in the same cluster, **after** the
  cluster is limited to the step folder.
- `dedup --contained` ignores `DeduplicationRules.enabled`. With
  `--disable` it writes the `contained` rows but keeps every entry.
- Don't change the output of plain `pandora dedup`.

## Review Focus

- An entry with no rows in `chains` would be contained in every other
  entry, since an empty set is a subset of anything. It must be **kept**
  and never count as contained. The test goes in Task 2.
- Three or more entries with equal counts: exactly one is kept, and every
  removed entry's `contained_in` is that kept winner. The test goes in
  Task 2.
- A chain in a cluster that is strictly contained, where its container is
  in turn contained in a third entry: `contained_in` must name a **kept**
  entry, never a removed one. The test goes in Task 2.
- A cluster whose only members are symmetry copies falls back to every
  member and still returns a representative, without crashing. The test
  goes in Task 4.
- Running `representatives` with the default policy but no
  `--similarity` must exit with a message that names `n_hits`, not raise
  a traceback. The test goes in Task 5.

---

### Task 1: Rank keys and the shared ranking helper

**Files:**
- Modify: `pandora/schemas/common.py` (add `RankKey`)
- Modify: `pandora/_util.py` (add `rank_sort_key`)
- Modify: `pandora/schemas/dataset.py` (`DeduplicationRules.contained_rank`,
  `ContainmentRecord`)
- Modify: `pandora/schemas/similarity.py` (`RepresentativePolicy`,
  `RepresentativeCluster`, `RepresentativesProvenance`)
- Test: `tests/test_ranking.py` (new)

**Interfaces:**
- Produces:
  - `RankKey(BaseModel)` with fields `key: str` and `order:
    Literal["min", "max"] = "min"`.
  - `rank_sort_key(item_id: str, values: Mapping[str, float | int |
    None], rank: Sequence[RankKey]) -> tuple`. A sort key: smaller means
    better. It raises `KeyError` if a rank key is missing from `values`.
  - `DeduplicationRules.contained_rank: list[RankKey] =
    [RankKey(key="resolution", order="min")]`.
  - `ContainmentRecord(entry_id: str, contained_in: str, equal: bool)`.
  - `RepresentativePolicy(items: SimilarityItems = "chain", rank:
    list[RankKey] = [RankKey(key="missing_residues", order="min"),
    RankKey(key="n_hits", order="max")])`.
  - `RepresentativeCluster(representative: str, members: list[str],
    n_members: int)`.
  - `RepresentativesProvenance(selected_at: str, items: SimilarityItems,
    policy: RepresentativePolicy, clusters: str, structures_dir: str,
    layers: list[ClusteredLayer] = [], n_clusters: int, n_members: int,
    n_members_dropped: int, n_clusters_dropped: int)`.

- [ ] **Step 1: Write the failing tests**

```python
from pandora._util import rank_sort_key
from pandora.schemas.common import RankKey

MIN = [RankKey(key="m")]
MAX = [RankKey(key="m", order="max")]


def _best(values, rank):
    return min(values, key=lambda i: rank_sort_key(i, values[i], rank))


def test_min_order_prefers_smaller():
    assert _best({"B": {"m": 1}, "A": {"m": 2}}, MIN) == "B"


def test_max_order_prefers_larger():
    assert _best({"B": {"m": 1}, "A": {"m": 2}}, MAX) == "A"


def test_null_ranks_last_for_both_orders():
    values = {"A": {"m": None}, "B": {"m": 5}}
    assert _best(values, MIN) == "B"
    assert _best(values, MAX) == "B"


def test_ties_fall_back_to_lowest_id():
    assert _best({"B": {"m": 1}, "A": {"m": 1}}, MIN) == "A"


def test_later_keys_break_earlier_ties():
    rank = [RankKey(key="m"), RankKey(key="n", order="max")]
    values = {"A": {"m": 1, "n": 1}, "B": {"m": 1, "n": 3}}
    assert _best(values, rank) == "B"
```

- [ ] **Step 2: Run them to check they fail**

Run: `uv run pytest tests/test_ranking.py -v`
Expected: FAIL with an ImportError for `rank_sort_key` / `RankKey`.

- [ ] **Step 3: Add the models and `rank_sort_key`**

For each key, `rank_sort_key` gives `(value is None, value if order ==
"min" else -value)`. A null becomes `(True, 0)`. After all the keys comes
`item_id`. Give every new model an `Attributes:` docstring in the style
of its neighbours. Export `RankKey` wherever `schemas/common.py` names are
re-exported, if they are.

- [ ] **Step 4: Run the tests and the existing suite**

Run: `uv run pytest tests/test_ranking.py tests/test_curation.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add pandora/_util.py pandora/schemas tests/test_ranking.py
git commit -m "feat(schemas): rank keys and representative/containment models"
```

---

### Task 2: `find_contained`

**Files:**
- Create: `pandora/datasets/containment.py`
- Modify: `pandora/datasets/__init__.py` (export `find_contained`)
- Test: `tests/test_containment.py` (new)

**Interfaces:**
- Consumes: `RankKey`, `rank_sort_key`, `ContainmentRecord` (Task 1).
- Produces: `find_contained(cluster_counts: Mapping[str, Mapping[str,
  int]], rank_values: Mapping[str, Mapping[str, float | int | None]],
  rank: Sequence[RankKey]) -> list[ContainmentRecord]`. Its input is
  entry id → {cluster id: number of the entry's chains in that cluster}.
  It returns one record per removed entry, sorted by `entry_id`.

- [ ] **Step 1: Write the failing tests**

The helper builds counts from strings: `"AA"` means two chains in
cluster `a`, `"AB"` means one in `a` and one in `b`.

```python
from collections import Counter

from pandora.datasets import find_contained
from pandora.schemas.common import RankKey

RES = [RankKey(key="resolution")]


def _counts(**entries):
    return {e: Counter(s.lower()) for e, s in entries.items()}


def _run(counts, res=None):
    res = res or {}
    values = {e: {"resolution": res.get(e)} for e in counts}
    return {
        r.entry_id: (r.contained_in, r.equal)
        for r in find_contained(counts, values, RES)
    }


def test_monomer_contained_in_homodimer():
    assert _run(_counts(E1="A", E2="AA")) == {"E1": ("E2", False)}


def test_heterodimer_contained_in_larger_complex():
    assert _run(_counts(E1="AB", E2="AAB")) == {"E1": ("E2", False)}


def test_unrelated_entries_both_kept():
    assert _run(_counts(E1="AB", E2="AC")) == {}


def test_equal_counts_keep_best_resolution():
    got = _run(_counts(E1="AB", E2="AB"), {"E1": 2.5, "E2": 1.5})
    assert got == {"E1": ("E2", True)}


def test_equal_counts_and_resolution_keep_lowest_id():
    got = _run(_counts(E2="AB", E1="AB"), {"E1": 2.0, "E2": 2.0})
    assert got == {"E2": ("E1", True)}


def test_missing_resolution_loses_a_tie():
    got = _run(_counts(E1="AB", E2="AB"), {"E2": 3.0})
    assert got == {"E1": ("E2", True)}


def test_three_equal_entries_keep_one():
    got = _run(_counts(E1="A", E2="A", E3="A"), {"E1": 3, "E2": 1, "E3": 2})
    assert got == {"E1": ("E2", True), "E3": ("E2", True)}


def test_contained_in_names_a_kept_entry():
    got = _run(_counts(E1="A", E2="AB", E3="ABC"))
    assert got == {"E1": ("E3", False), "E2": ("E3", False)}


def test_entry_without_chains_is_kept():
    assert _run(_counts(E1="", E2="AA")) == {}


def test_inputs_not_mutated():
    counts = _counts(E1="A", E2="AA")
    before = {k: dict(v) for k, v in counts.items()}
    _run(counts)
    assert {k: dict(v) for k, v in counts.items()} == before
```

- [ ] **Step 2: Run them to check they fail**

Run: `uv run pytest tests/test_containment.py -v`
Expected: FAIL with an ImportError for `find_contained`.

- [ ] **Step 3: Implement `find_contained` in `pandora/datasets/containment.py`**

The tests don't fully pin down the algorithm, so here it is:

```text
contains(Y, X) := all(Y.get(c, 0) >= n for c, n in X.items())
index: cluster -> set of entries having it
for each X with a non-empty count:
    candidates = index[X's cluster with the fewest entries] - {X}
    containers[X] = [Y for Y in candidates if contains(Y, X)]
    strict[X] = any(not contains(X, Y) for Y in containers[X])
removed = {X : strict[X]}
for each X not strict with equal-count containers:
    group = {X} + equal containers
    if X is not min(group, by rank_sort_key): removed[X] = equal
contained_in[X] = lowest id among containers[X] that are not removed
```

Strict containment wins over an equal tie (`equal=False`). A kept
container always exists, because the maximal entries keep one per equal
group. Add an `assert` on that, not a fallback. The module docstring
explains the counting rule in one line.

- [ ] **Step 4: Run the tests**

Run: `uv run pytest tests/test_containment.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add pandora/datasets tests/test_containment.py
git commit -m "feat(datasets): find_contained drops biounits contained in another"
```

---

### Task 3: `pandora dedup --contained`

**Files:**
- Modify: `pandora/cli/datasets.py` (new flags and
  `_cmd_dedup_contained`; `_cmd_dedup` calls it when `--contained` is
  set)
- Modify: `pandora/cli/_io.py:29` (add `"contained"` to `_STEP_TABLES`;
  the handler calls `prepare_step_dir(out, "dedup.json", overwrite)`)
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `find_contained` (Task 2), `DeduplicationRules`,
  `ContainmentRecord`, `SimilarityCluster`, `item_id`, `read_columns`,
  `read_records`, `write_table`, `prepare_step_dir`, `load_yaml_model`,
  `load_json_list`, `build_step_summary`, `find_table`.
- Produces: a step folder. `structures/` holds the kept entries' files
  at the same relative `structure_path`, copied with `shutil.copy2`.
  `entries.<fmt>` and `chains.<fmt>` are the input rows (as
  `CurationEntryRecord` / `DatasetChainRecord`) of the kept entries.
  `contained.<fmt>` holds `ContainmentRecord` rows. `failures.jsonl` is
  empty, and `dedup.json` is the `StepSummary`.

New flags on `dedup`:
- `--contained` (store_true).
- `--clusters PATH`.
- `--policy PATH` (a `DeduplicationRules` YAML).
- `--table-format` (`TABLE_FORMATS`, default `parquet`).
- `--overwrite`.

`--contained` without `--clusters` exits with `"dedup: --contained needs
--clusters"`.

Summary counts: `entries_in`, `entries_kept`, `entries_contained`,
`chains_kept`, `chains_unclustered`. `previous_summary` is
`<input>/curation.json` if that file exists. Print one warning line to
stderr when `chains_unclustered > 0`. The stdout line is `dedup: {kept}
retained, {n} contained -> {out}`. A chain that isn't in any cluster gets
its own cluster id, namely its item id. Rank values come from the
`entries` table. The allowed `contained_rank` keys are `resolution` and
`chain_count`; any other key exits naming it.

- [ ] **Step 1: Write the failing tests**

Build a tables-only step folder, as `_cluster_inputs` does, but with
`structure_path` set to `structures/<mid>/<id>.cif` and a dummy file
written there. Use entries `1AAA` (chain A, resolution 2.0), `2BBB`
(chains A and B, resolution 1.0) and `3CCC` (chain A, resolution 1.5).
The `clusters.json` is `[["1AAA_A", "2BBB_A", "2BBB_B"], ["3CCC_A"]]`.

```python
def test_dedup_contained_writes_step_folder(tmp_path):
    step, clusters = _contained_inputs(tmp_path)
    out = tmp_path / "uncontained"
    main(
        [
            "dedup",
            "--contained",
            "--input-dir",
            str(step),
            "--clusters",
            str(clusters),
            "--output-dir",
            str(out),
            "--table-format",
            "jsonl",
        ]
    )
    entries = [
        json.loads(l)["entry_id"]
        for l in (out / "entries.jsonl").read_text().splitlines()
    ]
    assert entries == ["2BBB", "3CCC"]
    chains = (out / "chains.jsonl").read_text().splitlines()
    assert len(chains) == 3
    assert json.loads((out / "contained.jsonl").read_text()) == {
        "entry_id": "1AAA",
        "contained_in": "2BBB",
        "equal": False,
    }
    assert (out / "structures/bb/2bbb.cif").exists()
    assert not (out / "structures/aa/1aaa.cif").exists()
    counts = json.loads((out / "dedup.json").read_text())["counts"]
    assert counts["entries_contained"] == 1
    assert counts["chains_unclustered"] == 0


def test_dedup_contained_counts_unclustered_chains(tmp_path, capsys):
    # clusters.json without 3CCC_A: 3CCC_A becomes its own cluster
    ...
    assert counts["chains_unclustered"] == 1
    assert "3CCC" in entries
    assert "warning" in capsys.readouterr().err


def test_dedup_contained_disable_keeps_everything(tmp_path):
    # --disable: contained.jsonl still has 1AAA, entries keeps all three
    ...


def test_dedup_contained_needs_clusters(tmp_path):
    with pytest.raises(SystemExit, match="--clusters"):
        main(
            [
                "dedup",
                "--contained",
                "--input-dir",
                str(tmp_path),
                "--output-dir",
                str(tmp_path / "o"),
            ]
        )


def test_dedup_contained_rejects_unknown_rank_key(tmp_path):
    # --policy YAML with contained_rank: [{key: n_atoms}] -> SystemExit "n_atoms"
    ...


def test_dedup_contained_refuses_existing_output(tmp_path):
    # second run without --overwrite -> SystemExit; with it -> passes
    ...
```

Complete the `...` tests in the same style, using the comments as the
assertions.

- [ ] **Step 2: Run them to check they fail**

Run: `uv run pytest tests/test_cli.py -k dedup_contained -v`
Expected: FAIL, because `--contained` is an unrecognised argument.

- [ ] **Step 3: Implement the flags and `_cmd_dedup_contained`**

Read the tables with `read_records(find_table(...), Model)`, so that
both `.parquet` and `.jsonl` inputs work. Missing `entries` or `chains`
tables exit with the message `"dedup: --contained needs a curate step
folder (entries and chains tables)"`.

- [ ] **Step 4: Run the CLI tests and the full suite**

Run: `uv run pytest tests/test_cli.py -v`
Expected: PASS, including the existing `test_canonicalise_curate_dedup_export_chain`.

- [ ] **Step 5: Commit**

```bash
git add pandora/cli tests/test_cli.py
git commit -m "feat(cli): dedup --contained writes an uncontained step folder"
```

---

### Task 4: `pick_representatives` and its helpers

**Files:**
- Create: `pandora/similarity/representatives.py`
- Modify: `pandora/similarity/__init__.py` (export the three functions)
- Test: `tests/test_representatives.py` (new)

**Interfaces:**
- Consumes: `RankKey`, `rank_sort_key`, `RepresentativePolicy`,
  `RepresentativeCluster` (Task 1), and `SimilarityCluster`.
- Produces:
  - `restrict_clusters(clusters: Iterable[SimilarityCluster], item_ids:
    Collection[str]) -> list[SimilarityCluster]`. Keeps only members in
    `item_ids` and drops clusters left empty.
  - `cluster_hit_counts(clusters: Iterable[SimilarityCluster], edges:
    Iterable[Sequence[Any]]) -> dict[str, int]`. For every clustered
    item, the number of distinct other items in its cluster that share
    an edge with it (0 if none). Edges to unclustered ids and self-edges
    are ignored.
  - `pick_representatives(clusters: Iterable[SimilarityCluster], values:
    Mapping[str, Mapping[str, float | int | None]], policy:
    RepresentativePolicy, eligible: Collection[str] | None = None) ->
    list[RepresentativeCluster]`. Picks the best eligible member of each
    cluster, falling back to every member if none is eligible. The result
    is sorted by `representative`. It raises `ValueError` naming the key
    when a policy key is missing from a ranked item's values.

- [ ] **Step 1: Write the failing tests**

```python
from pandora.schemas.common import RankKey
from pandora.schemas.similarity import RepresentativePolicy, SimilarityCluster
from pandora.similarity import (
    cluster_hit_counts,
    pick_representatives,
    restrict_clusters,
)


def _c(*ids):
    return SimilarityCluster(components=list(ids), n_components=len(ids))


POLICY = RepresentativePolicy()  # missing_residues min, n_hits max


def _v(missing, hits):
    return {"missing_residues": missing, "n_hits": hits}


def test_fewest_missing_residues_wins():
    reps = pick_representatives(
        [_c("A_A", "B_A")], {"A_A": _v(5, 9), "B_A": _v(1, 0)}, POLICY
    )
    assert reps[0].representative == "B_A"
    assert reps[0].members == ["A_A", "B_A"] and reps[0].n_members == 2


def test_n_hits_then_lowest_id_break_ties():
    values = {"A_A": _v(1, 1), "B_A": _v(1, 2), "C_A": _v(1, 2)}
    reps = pick_representatives([_c("A_A", "B_A", "C_A")], values, POLICY)
    assert reps[0].representative == "B_A"


def test_copy_is_never_picked():
    values = {"A_A": _v(3, 0), "A_A-2": _v(0, 0)}
    reps = pick_representatives(
        [_c("A_A", "A_A-2")], values, POLICY, eligible={"A_A"}
    )
    assert reps[0].representative == "A_A"
    assert reps[0].members == ["A_A", "A_A-2"]


def test_copies_only_cluster_falls_back_to_all_members():
    reps = pick_representatives(
        [_c("A_A-2")], {"A_A-2": _v(0, 0)}, POLICY, eligible=set()
    )
    assert reps[0].representative == "A_A-2"


def test_unknown_key_is_an_error():
    policy = RepresentativePolicy(rank=[RankKey(key="nope")])
    with pytest.raises(ValueError, match="nope"):
        pick_representatives([_c("A_A")], {"A_A": _v(0, 0)}, policy)


def test_restrict_drops_members_and_empty_clusters():
    got = restrict_clusters([_c("A_A", "B_A"), _c("C_A")], {"A_A"})
    assert [c.components for c in got] == [["A_A"]]
    assert got[0].n_components == 1


def test_hit_counts_stay_inside_clusters():
    clusters = [_c("A_A", "B_A", "C_A"), _c("D_A")]
    edges = [
        ("A_A", "B_A"),
        ("B_A", "A_A"),
        ("A_A", "C_A"),
        ("A_A", "D_A"),
        ("A_A", "A_A"),
        ("A_A", "Z_A"),
    ]
    assert cluster_hit_counts(clusters, edges) == {
        "A_A": 2,
        "B_A": 1,
        "C_A": 1,
        "D_A": 0,
    }


def test_inputs_not_mutated():
    clusters = [_c("A_A", "B_A")]
    restrict_clusters(clusters, {"A_A"})
    assert clusters[0].components == ["A_A", "B_A"]
```

- [ ] **Step 2: Run them to check they fail**

Run: `uv run pytest tests/test_representatives.py -v`
Expected: FAIL with an ImportError.

- [ ] **Step 3: Implement the three functions in `pandora/similarity/representatives.py`**

`cluster_hit_counts` maps each item to its cluster index, and collects a
`set` of partners per item for edges whose ends share a cluster index.
`pick_representatives` wraps the `KeyError` from `rank_sort_key` as a
`ValueError(f"rank key {key!r} ...")`.

- [ ] **Step 4: Run the tests**

Run: `uv run pytest tests/test_representatives.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add pandora/similarity tests/test_representatives.py
git commit -m "feat(similarity): pick one representative per cluster"
```

---

### Task 5: `pandora representatives`

**Files:**
- Modify: `pandora/cli/similarity.py` (`_cmd_representatives` and its
  parser in `register()`; `_hit_filters` and `_per_layer` take a
  `command: str` for their messages instead of the hard-coded
  `"cluster:"`)
- Modify: `pandora/cli/app.py` only if subcommands are listed there
  explicitly
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: Task 4's functions, `RepresentativePolicy`,
  `RepresentativesProvenance`, `ClusteredLayer`, `_hit_filters`,
  `_similarity_level`, `_layer_edges`, `_step_table`, `read_columns`,
  `load_json_list`, `load_yaml_model`, `write_records`, `write_json`.
- Produces:
  - `--output` holds a JSON list of `RepresentativeCluster` objects.
  - `representatives_provenance.json` sits next to it.
  - stdout: `{n} representatives from {m} chains -> {output}`.

The flags are:
- `--clusters` (required) and `--structures-dir` (required).
- `--items` (`choices=["chain"]`, default `chain`).
- `--similarity` (append), plus `--hit-filter`, `--min-score` and
  `--min-coverage` (append), with the same help as `cluster`.
- `--policy` (a `RepresentativePolicy` YAML).
- `--output` (required).

Behaviour:
- If `--items` is given, it overrides `policy.items`.
- Values per chain come from the `chains` table:
  - `missing_residues`: the sum of the three missing columns.
  - `n_hits`: from `cluster_hit_counts` over the restricted clusters and
    the union of the layer edges. Don't add copy links.
  - The direct columns `residue_count`, `atom_count`, `seqres_length`,
    `missing_n_term`, `missing_c_term`, `missing_middle`, `resolution`.
- Any other key exits naming it. Only compute `n_hits` when the policy
  uses it. If it's used without `--similarity`, exit with
  `"representatives: rank key n_hits needs --similarity (and a hit
  filter)"`.
- `eligible` = the chain item ids with a null `copy_of`.
- The provenance counts are `n_clusters` and `n_members` (after the
  restriction), plus `n_members_dropped` and `n_clusters_dropped` (what
  the restriction removed).

- [ ] **Step 1: Write the failing tests**

Reuse `_cluster_inputs` and `_similarity_layer`. Add `missing_*` columns
to the chain rows (a new optional argument of `_cluster_inputs`, or a
sibling helper).

```python
def test_representatives_default_policy(tmp_path):
    # chains 104M_A (missing 2), 112M_A (missing 2), 112M_B copy (missing 0),
    # 118L_A (missing 0); one cluster of all four; layer hits 104M_A-112M_A,
    # 104M_A-118L_A, 112M_A-118L_A -> n_hits 2/2/-/2 ; 118L_A wins on missing
    ...
    reps = json.loads(out.read_text())
    assert reps == [
        {
            "representative": "118L_A",
            "members": ["104M_A", "112M_A", "112M_B", "118L_A"],
            "n_members": 4,
        }
    ]
    prov = json.loads(
        (out.parent / "representatives_provenance.json").read_text()
    )
    assert prov["n_clusters"] == 1 and prov["items"] == "chain"


def test_representatives_n_hits_breaks_tie(tmp_path):
    # all missing 0, 112M_A has the most in-cluster hits -> 112M_A
    ...


def test_representatives_drops_chains_outside_step_folder(tmp_path):
    # clusters.json names 999X_A (not in chains table) -> not in members,
    # n_members_dropped == 1; a cluster of only 999X_A disappears,
    # n_clusters_dropped == 1
    ...


def test_representatives_n_hits_needs_similarity(tmp_path):
    with pytest.raises(SystemExit, match="n_hits"):
        main(
            [
                "representatives",
                "--clusters",
                str(clusters),
                "--structures-dir",
                str(step),
                "--output",
                str(out),
            ]
        )


def test_representatives_policy_without_n_hits_needs_no_layer(tmp_path):
    # --policy YAML rank: [{key: missing_residues}] and no --similarity -> runs
    ...


def test_representatives_rejects_interface_items(tmp_path):
    with pytest.raises(SystemExit):
        main([..., "--items", "interface"])


def test_representatives_rejects_unknown_rank_key(tmp_path):
    # policy rank: [{key: bogus}] -> SystemExit naming "bogus"
    ...
```

Fill in the `...` parts following the comments. Also add
`"representatives"` to the expected names in
`test_parser_registers_all_subcommands`.

- [ ] **Step 2: Run them to check they fail**

Run: `uv run pytest tests/test_cli.py -k representatives -v`
Expected: FAIL with an invalid choice `representatives`.

- [ ] **Step 3: Implement `_cmd_representatives` and its parser**

- [ ] **Step 4: Run the CLI tests and the full suite**

Run: `uv run pytest -q`
Expected: PASS. The existing `cluster` tests still pass after the
`_hit_filters` message change. Update any test that matches
`"cluster: set --min-score"` only if its text changed (keep `cluster:`
for `cluster`).

- [ ] **Step 5: Commit**

```bash
git add pandora/cli tests/test_cli.py
git commit -m "feat(cli): representatives picks one chain per cluster"
```

---

### Task 6: Docs and policies

**Files:**
- Create: `docs/policies/deduplication.yaml` (field types of
  `DeduplicationRules`) and `docs/policies/representatives.yaml` (field
  types of `RepresentativePolicy`), in the style of
  `docs/policies/dataset.yaml`
- Modify: `docs/reference/policies.md`. Add a `## Deduplication` section
  for `enabled` (what plain `dedup` does with it, and that `--contained`
  ignores it) and `contained_rank` (allowed keys `resolution`,
  `chain_count`). Add a `## Cluster representatives` section for `items`
  and `rank` (the allowed keys and what each means; the null-last rule;
  the lowest-id final tie-break; missing residues as an absolute count).
- Modify: `docs/usage/datasets.md`. After `## Deduplicate a batch`, add
  `### Contained biounits`, with library and CLI tabs, full code inline,
  and real output pasted from a run.
- Modify: `docs/usage/similarity.md`. After the clustering sections, add
  `### Representatives`, with library and CLI tabs and real output.
- Modify: `docs/usage/cli.md` (`## dedup` gains `--contained`, and a new
  `## representatives`)
- Modify: `pandora/cli/README.md` (add `representatives` to the command
  list, and one paragraph for both commands)
- Modify: `docs/contributing/architecture.md` (add `find_contained` and
  the three similarity functions to the function → model map)
- Regenerate: `uv run --extra docs python docs/scripts/generate_erd.py`

**Interfaces:**
- Consumes: everything above. The docs show exactly the flags and outputs
  from Tasks 3 and 5.

- [ ] **Step 1: Write the docs and policy YAMLs**

Run every example by hand, using the dev fixtures or a small tables-only
folder, and paste its real output. There's no docs-example test in the
repo.

- [ ] **Step 2: Regenerate the ERD**

Run: `uv run --extra docs python docs/scripts/generate_erd.py`
Expected: it exits 0, and the diagram files under `docs/assets/diagrams/`
are updated.

- [ ] **Step 3: Build the docs** (use the `docs-build` skill)

Run: `uv run zensical build --clean`
Expected: no warnings about broken links or anchors. The rendered HTML
contains the new code blocks (grep for `find_contained` and `pandora
representatives`).

- [ ] **Step 4: Lint and format**

Run: `uv run ruff format --check . && uv run ruff check .`
Expected: clean.

- [ ] **Step 5: Commit**

```bash
git add docs pandora/cli/README.md
git commit -m "docs: dedup --contained and representatives"
```

---

### Task 7: PPI v2.1 example and benchmark doc (untracked, no commit)

`examples/` is git-ignored, so these edits are local only.

**Files:**
- Modify: `examples/ppi/v2.1/src/02_sequence_clusters.sh`. Replace the
  "NOT YET IN THE CLI (#138)" block with `pandora dedup --contained
  --input-dir $D/work/curation --clusters $D/clusters/seq90/clusters.json
  --output-dir $D/work/uncontained --overwrite`, then `pandora
  representatives --clusters $D/clusters/seq90/clusters.json
  --structures-dir $D/work/uncontained --similarity $D/layers/seq
  --hit-filter $D/policies/seq90.yaml --policy
  $D/policies/representatives.yaml --output
  $D/clusters/seq90/representatives.json`.
- Create: `examples/ppi/v2.1/dataset/policies/representatives.yaml`
  (`rank: [{key: missing_residues, order: min}, {key: n_hits, order:
  max}]`)
- Modify: `examples/ppi/v2.1/src/03_structure_similarity.sh`,
  `04_dimers.sh`, `05_dimer_clusters.sh`: change `work/curation` to
  `work/uncontained`
- Modify: `examples/ppi/benchmarks/ppi-v2.1-benchmark.md`. In Step 1,
  record that containment counts chains per cluster and runs after
  curation. In Step 2, record absolute missing-residue counts and the
  lowest-id final tie-break. In the "With Pandora" notes, give the two
  commands.

- [ ] **Step 1: Make the edits**
- [ ] **Step 2: Run `bash examples/ppi/v2.1/src/02_sequence_clusters.sh`** if
  `mmseqs` is on `PATH`. Expected: both new commands print their summary
  lines. If `mmseqs` is missing, run only the two new commands against
  the existing `dataset/` outputs.
- [ ] **Step 3: Final check**

Run: `uv run pytest -q && uv run ruff format --check .`
Expected: PASS and clean.
