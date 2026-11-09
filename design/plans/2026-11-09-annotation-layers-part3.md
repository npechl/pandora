# Annotation layers, part 3 (similarity is a layer) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Similarity becomes a `similarity` layer of item pairs (entries,
chains or interfaces). `pandora similarity search` (MMseqs2/Foldseek)
and `pandora similarity import` (any edge list) write a Parquet layer
folder, and `pandora cluster --similarity` clusters any similarity
layer, streaming, at whole-PDBe-archive scale.

**Architecture:** Schema additions in `pandora/schemas/annotation.py`
(`similarity` template, `items`, `SimilarityRow`). A new module
`pandora/similarity/layer.py` holds item IDs, item sets and the
batch-wise matcher. `pandora/similarity/hits.py` keeps the tool TSV
column lists, and gains a TSV/Parquet batch reader and the batch-wise
edge filter (`similarity_edges`), which replaces `iter_edges`.
Everything that streams works on `pyarrow.RecordBatch`es. pyarrow is
imported inside the functions that need it. The CLI reads step-folder
tables column-wise and never parses structures, except for
`min_interface_coverage`.

**Tech Stack:** Python 3, pydantic v2, pyarrow (the `export` extra,
already a dependency there), argparse, pytest. No new dependencies.

**Spec:** `design/specs/2026-11-06-annotation-layers-design.md`, section
"Part 3: similarity is a layer".

## Global Constraints

- PEP 8, `ruff format`, line length 80. `from __future__ import
  annotations` at the top of new modules.
- Type hints everywhere. Public functions get a Google-style docstring
  (`Args:`, `Returns:`, `Raises:`). Internal helpers start with `_` and
  get a one-line docstring.
- Public names go in the package `__init__.py` `__all__`.
- Never mutate inputs.
- Catch named exception types only (Ruff BLE001, S110).
- No new dependency. pyarrow imported lazily inside functions.
- Tests are offline. Step folders for similarity tests can be just
  hand-written `entries.jsonl` / `chains.jsonl` (+ `curation.json`):
  similarity code reads only tables. Fake MMseqs2/Foldseek binaries
  are shell scripts, as in `tests/test_sequence_similarity.py` and
  `tests/test_structure_similarity.py`.
- No transition period for removed CLI flags, models or outputs
  (Pandora is 0.x).
- Similarity layers are **Parquet only**. Without pyarrow, the commands
  stop with: `similarity layers need pyarrow (the export extra)`.
- `examples/` is untracked: update the scripts there, but they get no
  tests and are not committed.
- Before calling the work done: `uv run pytest`, `uv run ruff check .`,
  `uv run ruff format --check .`, `uv run zensical build --clean`,
  `uv run --extra docs python docs/scripts/generate_erd.py`.

## Review Focus

1. A user TSV whose extra (uninterpreted) column holds numbers for the
   first rows and `n/a` later: the import stops with a message naming
   the file and the column, not a pyarrow traceback (Task 4 test).
2. A similarity layer where every row was unmatched (no
   `layer.parquet`): `cluster` gives all items as single-item clusters
   and does not crash (Task 5 test).
3. Chain IDs are case-sensitive while entry IDs are not: `104m_A`
   matches `104M_A`, but `104M_a` is `chain_not_found` (Task 3 test).
4. A Parquet table whose `score` column is integer (a count) is
   accepted and stored as float (Task 3 test).
5. `min_coverage` set and a row with null `coverage_1`: the row fails
   the threshold and doesn't crash (Task 5 test).

## Design decisions made in this plan

- The item list a layer was matched against is recorded as the
  summary's `previous_summary`: the step folder's
  `curation.json`/`canonicalisation.json`, or the contacts layer's
  `layer.json`. `StepSummary` does not change.
- `cluster_edges` keeps its dict-based union-find (string IDs). At
  about 1M items the dict is a few hundred MB, and the per-edge Python
  loop costs the same with integers. Mark it
  `# ponytail: dict union-find; integer arrays if memory bites at archive scale`.
  This departs from the spec's "integer indices". Say so at handoff.
- Old `HitFilter` JSON with `coverage_of: query` no longer loads (no
  transition). `tm_normalisation` in old JSON is ignored (pydantic's
  default `extra="ignore"`).
- Symmetry-copy links (copy → source) are fed to `cluster_edges` as
  extra edges, so `ClusteringProvenance.n_edges` counts them. Docstring
  says so.

---

### Task 1: Remove pairwise identity and the similarity replay in `reproduce`

**Files:**
- Delete: `pandora/annotations/pairwise.py`
- Modify: `pandora/annotations/__init__.py`, `pandora/metadata/__init__.py`
  (drop `annotate_pairwise_sequence_identity` from imports and
  `__all__`)
- Modify: `pandora/provenance/reproduce.py`
- Test: `tests/test_provenance.py`, and any test importing the removed
  function (`grep -rn pairwise_sequence_identity tests`)

**Interfaces:**
- Produces: `reproduce_dataset(manifest, output_dir, *, fetch_options)`
  with the same signature. It no longer replays clustering, partition
  or annotations: the returned manifest has `clustering=None`,
  `partition=None`, `splits=None`, and bundles with no annotations.
  `_reproduce_similarity`, `_reproduce_clustering`,
  `_reproduce_annotations` and `PAIRWISE_ANNOTATION_DISPATCH` are
  deleted.

- [ ] **Step 1: Replace the clustering reproduce tests**

Delete `test_reproduce_reruns_computed_search`,
`test_reproduce_reuses_existing_precomputed_hits`,
`test_reproduce_precomputed_hits_missing_raises`,
`test_reproduce_old_manifest_clustering_raises`,
`test_reproduce_interface_filter_raises` and their helpers. Add:

```python
def test_reproduce_does_not_replay_clustering(tmp_path, monkeypatch):
    # same fetch monkeypatching as test_reproduce_dataset
    manifest = _manifest_with_clustering(
        ClusteringProvenance(clustered_at="2026-01-01T00:00:00Z")
    )
    _, new_manifest = reproduce_dataset(manifest, tmp_path / "out")
    assert new_manifest.clustering is None
    assert new_manifest.partition is None
```

- [ ] **Step 2: Run, expect FAIL** (clustering replay raises without a
  search)

Run: `uv run pytest tests/test_provenance.py -v`

- [ ] **Step 3: Delete the function, its exports and the replay code.**
  Update the `reproduce_dataset` docstring: clustering, partition and
  annotations are not replayed, and the `Raises:` entries about
  searches go.

- [ ] **Step 4: Run all tests**

Run: `uv run pytest`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git commit -am "refactor!: remove pairwise sequence identity and reproduce's similarity replay"
```

---

### Task 2: `similarity` template, `items`, `SimilarityRow`

**Files:**
- Modify: `pandora/schemas/annotation.py`
- Modify: `pandora/annotations/layers.py` (`check_layer_columns`)
- Test: `tests/test_annotations.py`

**Interfaces:**
- Produces (all in `pandora/schemas/annotation.py`):
  - `LayerTemplate` gains `"similarity"`.
  - `SimilarityItems = Literal["entry", "chain", "interface"]`.
  - `LayerDescription.items: SimilarityItems | None = None`. A
    `model_validator(mode="after")` raises `ValueError` when
    `template == "similarity"` and `items is None` ("similarity layers
    need `items`: entry, chain or interface"), and when `items` is set
    for any other template.
  - `UnmatchedReason` gains `"interface_not_found"`, `"bad_item_id"`.
  - `class SimilarityRow(BaseModel)`, `extra="allow"`: `item_id_1: str`,
    `item_id_2: str`, `score: float`, `identity: float | None = None`,
    `coverage_1`/`coverage_2: float | None = None`,
    `start_1`/`end_1`/`start_2`/`end_2: int | None = None`. Docstring
    says that bulk reading never builds one per row.
- `check_layer_columns(columns, "similarity")` raises `ValueError`:
  `similarity layers are imported with \`pandora similarity import\``.

- [ ] **Step 1: Write the failing tests**

```python
def test_similarity_description_needs_items():
    with pytest.raises(ValidationError, match="items"):
        LayerDescription(name="g", template="similarity", tool="R")
    d = LayerDescription(
        name="g", template="similarity", tool="R", items="chain"
    )
    assert d.items == "chain"


def test_items_rejected_for_other_templates():
    with pytest.raises(ValidationError, match="items"):
        LayerDescription(
            name="c", template="residue_contacts", tool="x", items="chain"
        )


def test_similarity_row_coerces_int_score():
    row = SimilarityRow(
        item_id_1="1ABC_A", item_id_2="2XYZ_B", score=3, extra_col="k"
    )
    assert row.score == 3.0 and isinstance(row.score, float)
    assert row.model_dump()["extra_col"] == "k"


def test_annotate_import_refuses_similarity_template():
    with pytest.raises(ValueError, match="pandora similarity import"):
        check_layer_columns({"item_id_1", "item_id_2", "score"}, "similarity")
```

- [ ] **Step 2: Run, expect FAIL**

Run: `uv run pytest tests/test_annotations.py -k "similarity or items" -v`

- [ ] **Step 3: Implement the schema additions and the
  `check_layer_columns` guard** (first line of the function).

- [ ] **Step 4: Run all tests, expect PASS** (`uv run pytest`)

- [ ] **Step 5: Commit**

```bash
git commit -am "feat(schemas): similarity template, items level and SimilarityRow"
```

---

### Task 3: Item IDs, item sets and the batch matcher

**Files:**
- Create: `pandora/similarity/layer.py`
- Modify: `pandora/similarity/__init__.py`
- Test: `tests/test_similarity_layer.py` (new)

**Interfaces:**
- Consumes: `SimilarityItems`, `UnmatchedReason` (Task 2).
- Produces, in `pandora/similarity/layer.py`:
  - `item_id(level: SimilarityItems, entry_id: str, *chain_ids: str) -> str`.
    `entry` → `E`, `chain` → `E_C`, `interface` → `E_C1_C2` with the
    two chains sorted. Wrong number of chains raises `ValueError`.
  - `parse_item_id(level: SimilarityItems, item: str) -> tuple[str, ...]`:
    `(entry,)`, `(entry, chain)`, `(entry, chain_1, chain_2)`. Splits
    from the right (`rsplit("_", n)`), so `pdb_00001abc_A` →
    `("pdb_00001abc", "A")`. Raises `ValueError` on too few parts or
    an empty part.
  - `@dataclass(frozen=True) class ItemSet`: `level: SimilarityItems`;
    `ids: Mapping[str, str]` (lookup key → stored ID);
    `entries: frozenset[str]` (upper-cased entry IDs);
    `chains: frozenset[tuple[str, str]]` ((upper entry, chain); empty
    for `entry` and `interface`). The lookup key is the item ID with
    its entry part upper-cased and interface chains sorted. The stored
    ID uses the entry ID as the source table spells it.
  - `entry_items(entry_ids: Iterable[str]) -> ItemSet`
  - `chain_items(chains: Iterable[tuple[str, str]]) -> ItemSet`
    (entry_id, chain_id) pairs
  - `interface_items(rows: Iterable[Mapping[str, Any]]) -> ItemSet`:
    rows with keys `entry_id`, `chain_id_1`, `chain_id_2`, `seq_id_1`,
    `seq_id_2` (contact-layer columns as dicts, so the CLI can stream
    them without models). Keeps rows with both `seq_id`s present and
    different chains (the same filter as `extract_interface_records`).
    Unordered.
  - `SIMILARITY_REQUIRED = ("item_id_1", "item_id_2", "score")`,
    `SIMILARITY_FLOATS = ("score", "identity", "coverage_1", "coverage_2")`,
    `SIMILARITY_INTS = ("start_1", "end_1", "start_2", "end_2")`.
  - `check_similarity_columns(columns: Collection[str]) -> None`:
    raises `ValueError` naming the missing required columns.
  - `match_similarity_batches(batches: Iterable[pa.RecordBatch], items: ItemSet) -> Iterator[tuple[pa.RecordBatch, pa.RecordBatch]]`:
    yields one `(matched, unmatched)` pair per input batch. In
    `matched`, the item columns hold stored IDs, recognised float
    columns are `float64` and int columns `int64`, and every other
    column is kept as given. In `unmatched`, every input column is cast
    to string, plus a `reason` column. The reason is the first that
    applies, checking side 1 then side 2: `bad_item_id`,
    `entry_not_found`, `chain_not_found` (chain level only),
    `interface_not_found` (interface level, entry known), then
    `bad_value` (a recognised column present but not a number, or a
    null `score`).

  Approach for `match_similarity_batches`: a fast path using
  `pyarrow.compute.is_in` of both item columns against the stored IDs
  (a `pa.array` built once per call), plus a successful
  `pc.cast(..., float64/int64)` of the recognised columns. Only rows
  that fail it go through a per-row Python path (parse, look up the
  key, find the reason). Built-in search output takes the fast path
  almost entirely. Self-hits and repeated pairs are kept.

- [ ] **Step 1: Write the failing tests** (`tests/test_similarity_layer.py`)

```python
def test_item_id_round_trip_and_interface_sorted():
    assert item_id("interface", "1ABC", "B", "A") == "1ABC_A_B"
    assert parse_item_id("chain", "pdb_00001abc_A") == ("pdb_00001abc", "A")
    with pytest.raises(ValueError):
        parse_item_id("interface", "1ABC_A")


def _batch(**cols):
    return pa.RecordBatch.from_pydict(cols)


def test_chain_match_is_case_insensitive_on_entry_only():
    items = chain_items([("104M", "A"), ("112M", "A")])
    batch = _batch(
        item_id_1=["104m_A", "104M_a", "999X_A", "112M"],
        item_id_2=["112M_A"] * 4,
        score=[1, 1, 1, 1],
    )
    matched, unmatched = next(match_similarity_batches([batch], items))
    assert matched.column("item_id_1").to_pylist() == ["104M_A"]
    assert matched.schema.field("score").type == pa.float64()
    assert unmatched.column("reason").to_pylist() == [
        "chain_not_found",
        "entry_not_found",
        "bad_item_id",
    ]


def test_interface_items_come_from_contacts_only():
    rows = [
        dict(
            entry_id="1A2C",
            chain_id_1="B",
            chain_id_2="A",
            seq_id_1=5,
            seq_id_2=7,
        ),
        dict(
            entry_id="1A2C",
            chain_id_1="B",
            chain_id_2="E",
            seq_id_1=5,
            seq_id_2=None,
        ),
    ]
    items = interface_items(rows)
    batch = _batch(
        item_id_1=["1A2C_B_A", "1A2C_B_C"],
        item_id_2=["1A2C_A_B", "1A2C_A_B"],
        score=[0.9, 0.9],
    )
    matched, unmatched = next(match_similarity_batches([batch], items))
    assert matched.column("item_id_1").to_pylist() == ["1A2C_A_B"]
    assert unmatched.column("reason").to_pylist() == ["interface_not_found"]


def test_entry_items_and_bad_value():
    items = entry_items(["104M", "112M"])
    batch = _batch(
        item_id_1=["104M", "104M"],
        item_id_2=["112M", "112M"],
        score=["0.5", "x"],
        keep=["a", "b"],
    )
    matched, unmatched = next(match_similarity_batches([batch], items))
    assert matched.column("score").to_pylist() == [0.5]
    assert matched.column("keep").to_pylist() == ["a"]
    assert unmatched.column("reason").to_pylist() == ["bad_value"]


def test_self_hits_and_repeats_are_kept():
    items = chain_items([("104M", "A")])
    batch = _batch(
        item_id_1=["104M_A"] * 2, item_id_2=["104M_A"] * 2, score=[1.0, 1.0]
    )
    matched, _ = next(match_similarity_batches([batch], items))
    assert matched.num_rows == 2


def test_check_similarity_columns_names_missing():
    with pytest.raises(ValueError, match="score"):
        check_similarity_columns({"item_id_1", "item_id_2"})
```

- [ ] **Step 2: Run, expect FAIL** (`ModuleNotFoundError`)

Run: `uv run pytest tests/test_similarity_layer.py -v`

- [ ] **Step 3: Implement `pandora/similarity/layer.py` and export its
  public names**

- [ ] **Step 4: Run, expect PASS**, then `uv run pytest`

- [ ] **Step 5: Commit**

```bash
git add pandora/similarity tests/test_similarity_layer.py
git commit -m "feat(similarity): item IDs, item sets and batch matching for similarity layers"
```

---

### Task 4: `pandora similarity search` and `pandora similarity import`

**Files:**
- Modify: `pandora/similarity/hits.py` (add the readers. `iter_edges` and
  `load_similarity_search` stay until Task 5)
- Modify: `pandora/similarity/sequence.py`, `pandora/similarity/structure.py`
  (return type)
- Modify: `pandora/cli/_io.py` (`read_columns`, `ParquetSink`)
- Modify: `pandora/cli/similarity.py` (`similarity` becomes a group.
  `cluster`/`partition` untouched here)
- Test: `tests/test_hits.py`, `tests/test_sequence_similarity.py`,
  `tests/test_structure_similarity.py`, `tests/test_cli.py`

**Interfaces:**
- Consumes: Task 3's `chain_items`, `entry_items`, `interface_items`,
  `check_similarity_columns`, `match_similarity_batches`;
  `LayerDescription` (Task 2).
- Produces:
  - `pandora/similarity/hits.py`:
    - `read_hit_tsv(path: str | Path, engine: str, *, batch_size: int = 1 << 20) -> Iterator[pa.RecordBatch]`:
      a headerless tool TSV with `HIT_COLUMNS[engine]`, renamed to
      the template. `query`→`item_id_1`, `target`→`item_id_2`,
      `fident`→`identity`, `qcov`→`coverage_1`, `tcov`→`coverage_2`,
      `qstart`/`qend`/`tstart`/`tend`→`start_1`/`end_1`/`start_2`/`end_2`.
      Adds `score` = `fident` (MMseqs2) or `alntmscore` (Foldseek).
      `alnlen`, `alntmscore`, `qtmscore`, `ttmscore` are kept under
      their own names. Raises `ValueError` for an unknown engine.
    - `read_similarity_table(path: str | Path) -> Iterator[pa.RecordBatch]`:
      `.parquet` via `pq.ParquetFile(path).iter_batches()`. `.tsv` via
      `pyarrow.csv.open_csv` (tab delimiter, header row), with the item
      columns and the recognised value columns read as strings, so that
      the matcher decides `bad_value`. Any other suffix raises
      `ValueError`. A pyarrow conversion error while streaming is
      re-raised as `ValueError` naming the file and pyarrow's message
      (which names the column).
  - `compute_sequence_similarity(...)` and `compute_structure_similarity(...)`
    keep their arguments and still write the TSV to `hits_path`. They
    now return `LayerDescription(name="mmseqs2"|"foldseek",
    template="similarity", items="chain", tool="MMseqs2"|"Foldseek",
    tool_version=<version>, parameters=<as today's parameters>)`.
  - `pandora/cli/_io.py`:
    - `read_columns(path: Path, columns: Sequence[str]) -> Iterator[dict[str, Any]]`:
      rows as dicts with only `columns`. Parquet is read batch-wise
      with pyarrow, JSON Lines line by line. A column absent from the
      file gives `None`.
    - `class ParquetSink`: `ParquetSink(path)`, `.write(batch)` opens a
      `pq.ParquetWriter` on the first non-empty batch, `.close()`,
      `.rows: int`. Usable as a context manager. Writes no file if
      nothing was written.
  - CLI: `pandora similarity search` and `pandora similarity import`,
    flags as in the spec. The search also takes `--structures-dir`
    (required) and `--input-dir` (Foldseek only, required there). Both
    write `layer.parquet`, `unmatched.parquet` and `layer.json`, and
    `import` also writes `source/`. `prepare_output_dir` covers these
    names plus `hits.tsv`. The search writes the tool TSV to
    `<output-dir>/hits.tsv`, converts it, then deletes it. Steps:
    `"similarity_search"` / `"similarity_import"`. `previous_summary`
    is the item source's summary (see "Design decisions"). Counts:
    `rows_in`, `rows_matched`, `rows_unmatched`, `unmatched:<reason>`.
  - Item list helper (in `pandora/cli/similarity.py`):
    `_item_set(command: str, level: SimilarityItems, structures_dir: Path | None, contacts: Path | None) -> tuple[ItemSet, Path | None]`
    returns the item set and the previous-summary path. `entry` reads
    `entries.<fmt>` (`entry_id`, `structure_path`; rows with a
    structure). `chain` reads `chains.<fmt>` (`entry_id`, `chain_id`).
    `interface` reads the contacts layer's `layer.<fmt>`, whose
    `layer.json` template must be `residue_contacts`. A missing table,
    or the wrong flag for the level, is a `SystemExit` naming the step
    that writes it (`chains` → `pandora curate`). Uses
    `find_table` and `read_columns`, never `read_structures_dir`.
  - MMseqs2 search input: `{item_id("chain", e, c): seq}` from
    `chains.<fmt>` (`entry_id`, `chain_id`, `sequence`, `copy_of`), for
    rows with a sequence and **no `copy_of`**.

- [ ] **Step 1: Write the failing tests**

`tests/test_hits.py`:
```python
def test_read_hit_tsv_maps_foldseek_columns(tmp_path):
    # one Foldseek row: q t fident alnlen qcov tcov alntm qtm ttm qs qe ts te
    batch = next(read_hit_tsv(_foldseek_file(tmp_path), "Foldseek"))
    row = batch.to_pylist()[0]
    assert row["item_id_1"] == "a_A" and row["score"] == row["alntmscore"]
    assert {"coverage_1", "start_2", "qtmscore"} <= set(row)


def test_read_similarity_table_bad_extra_column_names_file(tmp_path):
    path = tmp_path / "edges.tsv"
    lines = ["item_id_1\titem_id_2\tscore\tnote_num"]
    lines += [f"104M_A\t112M_A\t0.5\t{i}" for i in range(5000)] + [
        "104M_A\t112M_A\t0.5\tn/a"
    ]
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="edges.tsv"):
        list(read_similarity_table(path))
```
(If pyarrow infers types from a block bigger than 5000 rows, make the
file large enough that `n/a` falls in a later block. Pass a small
`block_size` in `ReadOptions` only in the test, through a keyword
argument, if needed.)

`tests/test_sequence_similarity.py`: update the existing tests so they
assert the returned `LayerDescription` (`items == "chain"`,
`tool == "MMseqs2"`, `parameters["max_seqs"] == 300`). Same for
Foldseek in `tests/test_structure_similarity.py`.

`tests/test_cli.py`. Replace the old `similarity` tests
(`test_similarity_writes_hits_and_search_record`,
`test_similarity_wraps_precomputed_hits`,
`test_similarity_needs_input_dir_without_precomputed_hits`,
`test_similarity_records_absolute_hits_path`,
`test_similarity_rejects_search_options_with_precomputed_hits`) with
tests that use a hand-written step folder:

```python
def _similarity_step(tmp_path):
    """Step folder with tables only: entries 104M, 112M; chains 104M_A, 112M_A, 112M_B (copy of A)."""
    ...  # entries.jsonl, chains.jsonl (with sequence, copy_of), curation.json


def test_similarity_search_mmseqs_writes_layer_without_copies(tmp_path):
    # fake mmseqs: records the query FASTA, writes rows 104M_A->112M_A and 104M_A->999X_A
    ...
    assert ">112M_B" not in query_fasta
    rows = pq.read_table(out / "layer.parquet").to_pylist()
    assert [(r["item_id_1"], r["item_id_2"]) for r in rows] == [
        ("104M_A", "112M_A")
    ]
    assert pq.read_table(out / "unmatched.parquet").column(
        "reason"
    ).to_pylist() == ["entry_not_found"]
    assert not (out / "hits.tsv").exists()
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "similarity_search"
    assert summary["policy"]["items"] == "chain"
    assert summary["previous_summary"].endswith("curation.json")


def test_similarity_import_tsv_entry_level(tmp_path):
    ...
    # description items: entry; rows 104m->112M, 104M->999X; layer has 104M/112M, unmatched entry_not_found


def test_similarity_import_interface_level_needs_contacts(tmp_path):
    ...
    # --structures-dir with items: interface -> SystemExit mentioning --contacts


def test_similarity_import_chain_level_on_canonicalise_folder(tmp_path):
    ...
    # folder without chains table -> SystemExit mentioning "pandora curate"


def test_similarity_import_missing_items_stops_before_writing(tmp_path):
    ...
    # description without items -> SystemExit, output dir has no layer.json


def test_similarity_search_refuses_existing_layer(tmp_path):
    ...
    # second run without --overwrite -> SystemExit
```

- [ ] **Step 2: Run, expect FAIL**

Run: `uv run pytest tests/test_hits.py tests/test_sequence_similarity.py tests/test_structure_similarity.py tests/test_cli.py -k "similarity or read_hit or read_similarity or search" -v`

- [ ] **Step 3: Implement the readers, the return-type change,
  `read_columns`/`ParquetSink`, and the `similarity` group.** Keep
  `cluster` and `partition` registering as they are. The old cluster
  tests build `SimilaritySearch` by hand and keep passing until
  Task 5.

- [ ] **Step 4: Run all tests, expect PASS** (`uv run pytest`)

- [ ] **Step 5: Commit**

```bash
git commit -am "feat(cli)!: similarity search/import write similarity layers"
```

---

### Task 5: Cluster any similarity layer

**Files:**
- Modify: `pandora/schemas/similarity.py`
- Modify: `pandora/similarity/hits.py`, `pandora/similarity/clustering.py`,
  `pandora/similarity/__init__.py`
- Modify: `pandora/datasets/records.py`, `pandora/datasets/__init__.py`
  (remove `entry_sequences` if `grep` shows no other user)
- Modify: `pandora/cli/similarity.py` (`cluster`), `pandora/cli/_io.py`
  if needed
- Modify: `pandora/provenance/manifest.py`, `pandora/schemas/provenance.py`,
  `pandora/export/json.py`, `pandora/schemas/metadata.py` (only where
  they import removed names)
- Test: `tests/test_hits.py`, `tests/test_clustering.py`, `tests/test_cli.py`

**Interfaces:**
- Consumes: Task 3's `parse_item_id`, `ItemSet`; Task 4's
  `read_columns`, `_item_set`, `read_similarity_table` (for
  `layer.parquet`).
- Produces:
  - `HitFilter`: `tm_normalisation` removed;
    `score_column: str = "score"`;
    `coverage_of: Literal["both", "side_1", "side_2", "either"] = "both"`.
    The other fields are unchanged. Update the docstring.
  - `Edge`: `identity: float | None`, `coverage: float | None`.
  - `ClusteringProvenance`: `search` removed;
    `similarity_layer: str | None = None` (path of the layer's
    `layer.json`).
  - Removed: `SimilaritySearch`, `SimilarityEngine`,
    `SimilarityRelationship`, `SimilarityRelationshipProvenance`,
    `SimilarityMethod`, `SimilarityType` (`grep` first; keep any still
    used); `iter_edges`, `load_similarity_search`,
    `cluster_similar_items`; `entry_sequences`.
  - `pandora/similarity/hits.py`:
    `similarity_edges(batches: Iterable[pa.RecordBatch], hit_filter: HitFilter, *, interface_residues: Mapping[str, set[int]] | None = None) -> Iterator[Edge]`.
    Before yielding anything, it checks the first batch's schema and
    raises `ValueError` for:
    - a set threshold whose column is absent (`score_column`;
      `identity`; `coverage_1`/`coverage_2` as `coverage_of` needs);
    - `min_interface_coverage` without `interface_residues`;
    - `interface_residues` without the `start_*`/`end_*` columns.

    It filters each batch with `pyarrow.compute`. Thresholds are
    inclusive. Null or NaN fails (`fill_null(False)`; NaN `>=` is
    false). Coverage is min, `side_1`, `side_2` or max of
    `coverage_1`/`coverage_2`. Interface coverage is computed per row,
    only for rows that pass the other thresholds and whose two items
    are both in `interface_residues`. It is the min over the two sides
    of the fraction of interface positions inside `[start, end]`,
    reusing today's `_interface_coverage`. Self-hits are skipped, and
    `Edge` IDs are sorted.
  - `pandora/similarity/clustering.py`:
    - `entry_edges(edges: Iterable[Sequence[Any]]) -> Iterator[tuple[str, str]]`:
      maps each chain-item edge to its two entries (`parse_item_id`)
      and skips pairs within one entry.
    - `cluster_edges` unchanged apart from the `ponytail:` comment.
  - CLI `pandora cluster`: `--similarity` (required), `--structures-dir`,
    `--contacts`, `--items {entry,chain,interface}` (default: the
    layer's level; only `entry` on a chain layer may differ, anything
    else is `SystemExit`), `--hit-filter`, `--min-score`,
    `--min-coverage`, `--pairs`, `--output`. `--search` and
    `--input-dir` are removed. The layer folder's `layer.json`
    template must be `similarity`. Item list:
    - `entry` and `chain` from `--structures-dir`;
    - `interface` from `--contacts`;
    - `chain` with `--contacts` additionally builds interface residues
      with `interface_residues_from_annotation(read_structures_dir(structures_dir), contact rows)`.
      This is the known ceiling.

    Edges = `similarity_edges(read_similarity_table(layer.parquet), ...)`,
    or nothing if the file is absent; wrapped in `entry_edges` for the
    roll-up. For chain items, chained with copy links
    `(item_id("chain", e, c), item_id("chain", e, copy_of))` from
    `chains.<fmt>`. `ClusteringProvenance.similarity_layer` is the
    `layer.json` path.

- [ ] **Step 1: Write the failing tests**

`tests/test_hits.py`: rewrite the `iter_edges` tests over batches. Keep
their cases: inclusive thresholds, any passing row makes an edge, NaN
fails, interface coverage, self-hits skipped. Add:

```python
def test_similarity_edges_score_column():
    ...
    # HitFilter(min_score=0.6, score_column="qtmscore") keeps a row with score 0.1, qtmscore 0.7


@pytest.mark.parametrize(
    "coverage_of, expected",
    [("both", 0.4), ("side_1", 0.9), ("side_2", 0.4), ("either", 0.9)],
)
def test_similarity_edges_coverage_of(coverage_of, expected): ...


def test_similarity_edges_null_coverage_fails_threshold():
    ...
    # coverage_1=None with min_coverage=0.1 -> no edge, no error


def test_similarity_edges_threshold_on_absent_column_raises():
    ...
    # min_identity set, no identity column -> ValueError before iteration yields


def test_similarity_edges_batching_does_not_change_edges():
    ...
    # same rows as one batch vs three batches -> identical edge lists
```

`tests/test_clustering.py`:
```python
def test_entry_edges_links_entries_through_any_chain():
    edges = [("1AAA_B", "2BBB_C"), ("1AAA_A", "1AAA_B")]
    assert list(entry_edges(edges)) == [("1AAA", "2BBB")]
```
Delete `test_cluster_similar_items_filters_hits_and_records_search`.
Keep `test_old_clustering_provenance_still_loads`; an old `search` key
is ignored.

`tests/test_cli.py`. Replace `_cluster_inputs` and its users with a
helper that writes a similarity layer (`layer.parquet` +
`layer.json`) next to `_similarity_step`:

```python
def test_cluster_chain_layer_with_copies_and_singletons(tmp_path):
    # layer: 104M_A-112M_A score 0.9; items 104M_A, 112M_A, 112M_B (copy of A), 118L_A
    # -> clusters [[104M_A, 112M_A, 112M_B], [118L_A]]

def test_cluster_rolls_chain_layer_up_to_entries(tmp_path):
    # chains 1AAA_A, 1AAA_B, 2BBB_A; edge only 1AAA_B-2BBB_A (B is not the longest chain)
    # --items entry -> one cluster [1AAA, 2BBB]

def test_cluster_layer_with_no_matched_rows_gives_singletons(tmp_path):
    # layer folder with layer.json but no layer.parquet -> one cluster per item

def test_cluster_drops_edges_outside_smaller_item_list(tmp_path):
    # layer has 104M_A-999X_A; step folder lacks 999X -> cluster_provenance n_edges_unknown_ids == 1

def test_cluster_rejects_non_similarity_layer(tmp_path): ...
def test_cluster_rejects_items_other_than_entry_rollup(tmp_path): ...
    # chain layer, --items interface -> SystemExit
```
Update `test_cluster_with_pairs`, `test_cluster_reads_hit_filter_yaml_and_records_it`,
`test_cluster_flag_overrides_hit_filter_yaml` and
`test_cluster_requires_a_threshold` to the new flags. Their assertions
stay. The record assertion checks `similarity_layer` instead of
`search`. Update the existing `cluster --contacts` test from part 2 so
it uses a Foldseek-shaped similarity layer (with `start_*`/`end_*`).

- [ ] **Step 2: Run, expect FAIL**

Run: `uv run pytest tests/test_hits.py tests/test_clustering.py tests/test_cli.py -k "edges or cluster" -v`

- [ ] **Step 3: Implement the schema changes, `similarity_edges`,
  `entry_edges`, the removals and the `cluster` command.** `grep -rn`
  each removed name across `pandora/ tests/ docs/scripts/` and fix
  every hit.

- [ ] **Step 4: Run all tests and lint**

Run: `uv run pytest && uv run ruff check .`
Expected: PASS, no errors.

- [ ] **Step 5: Commit**

```bash
git commit -am "feat(similarity)!: cluster any similarity layer; entry roll-up; drop SimilaritySearch"
```

---

### Task 6: Docs, recipe, ERD, examples

**Files:**
- Modify: `docs/usage/similarity.md`, `docs/usage/cli.md`,
  `docs/usage/annotation.md`, `docs/usage/export.md`,
  `docs/usage/provenance.md`, `docs/getting-started/glossary.md`,
  `docs/contributing/architecture.md`, `docs/contributing/contribute.md`,
  `docs/reference/policies.md` (`HitFilter`: `score_column`,
  `coverage_of` sides, no `tm_normalisation`), `docs/recipes/ppi-01.md`,
  `pandora/cli/README.md`, `README.md`
- Regenerate: `docs/assets/diagrams/*` via `docs/scripts/generate_erd.py`
  (drop the removed models from its list first)
- Update (not committed): `examples/ppi_dataset_pipeline.py`,
  `examples/dataset_pipeline.py`, `examples/overview.py`,
  `examples/ppi/v2.1/b_sequence_clusters.py`

**Interfaces:**
- Consumes: the commands and functions of Tasks 4–5.

- [ ] **Step 1: Rewrite the similarity usage page** around the three
  commands (`similarity search`, `similarity import`, `cluster
  --similarity`): the template table, the item ID table, where the
  item list comes from, the entry roll-up, symmetry copies, and
  Parquet only. Follow the docs-inline convention (full code inline,
  library and CLI tabs).

- [ ] **Step 2: Move the PPI recipe** (`docs/recipes/ppi-01.md`) to a
  chain-level MMseqs2 search plus the `--items entry` roll-up (library:
  `entry_edges`). Run the recipe's code and replace the documented
  cluster count with the new measured number.

- [ ] **Step 3: Update the remaining pages, READMEs and examples**
  listed above. Run `grep -rn "precomputed-hits\|SimilaritySearch\|iter_edges\|cluster_similar_items\|entry_sequences\|pairwise_sequence_identity\|tm_normalisation\|--search" docs README.md pandora examples`
  and expect no hits outside historical design docs.

- [ ] **Step 4: Regenerate the ERD and build the docs**

Run: `uv run --extra docs python docs/scripts/generate_erd.py && uv run zensical build --clean`
Expected: no errors, no broken-link warnings (use the docs-build skill).

- [ ] **Step 5: Full check and commit**

Run: `uv run pytest && uv run ruff check . && uv run ruff format --check .`

```bash
git add docs README.md pandora/cli/README.md
git commit -m "docs: similarity layers, similarity search/import, cluster --similarity"
```
