# Similarity results as streamed hit files

Date: 2026-10-07. Status: approved design, awaiting spec review.

## Context

This is sub-project 2 of the benchmark roadmap. The roadmap decisions
are in `design/specs/2026-10-07-contacts-design.md` under "Roadmap
context". The ones that matter here:

- Pandora builds datasets at whole-PDB scale.
- It stays a library of pure functions.
- It adds no new dependency without agreement.
- MMseqs2 and Foldseek keep their thin wrappers and also accept
  precomputed results.

## Problem

`compute_sequence_similarity` and `compute_structure_similarity` read
the whole result file into memory, keep a dict entry for the best hit of
each pair, and return one `SimilarityRelationship` per pair, each with
its own nested `SimilarityMethod`.

Measured memory is about 2.4 KB per pair (211 B for the dict entry,
2,208 B for the object):

| Pairs | Memory |
|---|---|
| 10⁷ | ~24 GB |
| 10⁸ | ~240 GB |

An all-vs-all search over the whole PDB lands in that range.

There are three further problems:

- Clustering ignores coverage, so the specs' combined rules ("identity
  ≥ 0.9 **and** coverage ≥ 0.9", "TM ≥ 0.5 **and** coverage ≥ 0.8") have
  to be pre-filtered by hand.
- Search caps (`--max-seqs`) are neither exposed nor recorded. This is
  PPI v2.1 limitation D3: capped searches hide leakage.
- TM-score normalisation is fixed to the alignment length (limitation
  D4).

## Decisions (agreed in brainstorming)

1. **The full hit file stays on disk and is streamed.** Memory is
   bounded by the edges kept, not by the hits. Re-thresholding means
   re-reading the file, not re-running the search. The stdlib is enough,
   so no pyarrow or numpy in the core dependencies.
2. **The stored file is the tool's raw TSV, with columns Pandora
   fixes.** Everything else is derived at read time: pair orientation,
   interface coverage and TM normalisation.
3. **Thresholds are a policy model, `HitFilter`,** loadable from YAML
   and recorded in clustering provenance.
4. **The list-returning API is replaced, not kept alongside.** This is a
   breaking change at 0.5.x, with a migration note in the docs.
   `SimilarityRelationship` stays as a schema for small, hand-built
   networks.
5. **A pair is an edge if any single hit row passes every threshold.**
   Edges are produced row by row, with each pair's IDs sorted. A pair
   may be yielded more than once (once per passing alignment or
   direction). Clustering is unaffected by duplicates. Callers that need
   unique pairs deduplicate with a `set`.

## Design

### Fixed output columns

| Engine | `--format-output` columns |
|---|---|
| MMseqs2 | `query,target,fident,alnlen,qcov,tcov` (unchanged) |
| Foldseek | `query,target,fident,alnlen,qcov,tcov,alntmscore,qtmscore,ttmscore,qstart,qend,tstart,tend` (adds `qtmscore`, `ttmscore`) |

Coverage and identity are fractions in [0, 1], which is the default for
both tools.

### Schemas (`pandora/schemas/similarity.py`)

```python
class SimilaritySearch(BaseModel):
    engine: Literal["MMseqs2", "Foldseek"]
    version: str | None = None
    hits_path: str
    columns: list[str]
    parameters: dict[str, Any] = {}
    origin: Literal["computed", "precomputed"] = "computed"
    searched_at: str | None = None


class HitFilter(BaseModel):
    min_score: float | None = None  # TM-score (Foldseek) / identity (MMseqs2)
    min_identity: float | None = None
    min_coverage: float | None = None
    coverage_of: Literal["both", "query", "target", "either"] = "both"
    min_interface_coverage: float | None = None
    tm_normalisation: Literal["alignment", "query", "target", "max"] = (
        "alignment"
    )
```

How the filter reads a hit:

- **Score:**
  - MMseqs2: `fident`.
  - Foldseek: `alntmscore`, `qtmscore`, `ttmscore` or
    `max(qtmscore, ttmscore)`, chosen by `tm_normalisation`.
- **Coverage**, chosen by `coverage_of`:
  - `both`: `min(qcov, tcov)`;
  - `query`: `qcov`;
  - `target`: `tcov`;
  - `either`: `max(qcov, tcov)`.
- **Thresholds:** a `None` threshold is not applied, and every other
  threshold is inclusive (`>=`).

`ClusteringProvenance` becomes:

```python
class ClusteringProvenance(BaseModel):
    clustered_at: str
    hit_filter: HitFilter | None = None
    search: SimilaritySearch | None = None
    n_edges: int
    n_edges_unknown_ids: int = 0
    n_clusters: int
```

It replaces `threshold`, `n_relationships` and `similarity_method`.
`hit_filter` and `search` are `None` when clustering hand-built edges
through `cluster_edges`.

`Edge` is a `NamedTuple` in `pandora/schemas/similarity.py`, so it
isn't a per-edge pydantic object:

```python
class Edge(NamedTuple):
    source_id: str
    target_id: str
    score: float
    identity: float
    coverage: float
    interface_coverage: float | None
```

`source_id < target_id` always holds, and `coverage` is the value
selected by `coverage_of`.

### Functions (`pandora/similarity/`)

| Function | Behaviour |
|---|---|
| `compute_sequence_similarity(sequences, hits_path, *, mmseqs_bin="mmseqs", sensitivity=5.7, max_seqs=300, tmp_dir=None, mmseqs_options=None) -> SimilaritySearch` | Runs `easy-search` with the fixed columns, `-s` and `--max-seqs`, writes the TSV to `hits_path` and returns the record. `parameters` holds `sensitivity`, `max_seqs` and `mmseqs_options`. `sequences` accepts the same input types as today. |
| `compute_structure_similarity(structures, hits_path, *, foldseek_bin="foldseek", sensitivity=9.5, alignment_type=2, max_seqs=1000, exhaustive_search=False, tmp_dir=None, foldseek_options=None) -> SimilaritySearch` | The same for Foldseek, adding `--alignment-type`, `--max-seqs` and `--exhaustive-search 1` when set. `interface_residues` is removed from this function. |
| `load_similarity_search(hits_path, engine, *, version=None, parameters=None) -> SimilaritySearch` | Wraps an existing TSV (`origin="precomputed"`). It checks that the first non-empty row has the engine's column count, so an obviously wrong file fails early. |
| `iter_edges(search, hit_filter, *, interface_residues=None) -> Iterator[Edge]` | Opens `hits_path` and reads it line by line. It skips self-hits, validates each row's column count, computes interface coverage when `interface_residues` is given (the existing `_interface_coverage`, min over both sides), and yields an `Edge` for every row passing all thresholds. |
| `cluster_edges(item_ids, edges) -> tuple[list[SimilarityCluster], ClusteringProvenance]` | The existing union-find, taking any iterable of `(source, target, ...)` tuples (an `Edge` or a plain pair). Edges whose IDs aren't in `item_ids` are ignored and counted in `n_edges_unknown_ids`. |
| `cluster_similar_items(item_ids, search, hit_filter, *, interface_residues=None) -> tuple[list[SimilarityCluster], ClusteringProvenance]` | Runs `cluster_edges(item_ids, iter_edges(...))`, then records `hit_filter` and `search`. |

Default `max_seqs` values are each tool's own default. They are always
recorded, and the docs warn that large families exceed them.

**Unchanged:** `pair_cluster_keys`, `partition_dataset`,
`interface_residues_from_annotation`, `residue_positions`,
`chain_item_id`, and the `SimilarityRelationship` schema.

`HitFilter` is the implemented form of the edge rules sketched in
`docs/policies/leakage.yaml` (`similarity_rules.*.threshold` and
`coverage_threshold`). That design file is updated to reference
`HitFilter` instead of the two loose fields, so there is one definition.

### Errors

- **Malformed rows:** a row with the wrong column count raises
  `ValueError("<path>:<line>: expected N columns (<names>), got M")`
  while streaming.
- **TM normalisation on MMseqs2:** `tm_normalisation != "alignment"` on
  an MMseqs2 search raises `ValueError` when `iter_edges` is called.
- **Interface filter without interfaces:** `min_interface_coverage` set
  without `interface_residues` raises `ValueError`.
- **One side without interface residues:** the pair has
  `interface_coverage=None`, and it fails `min_interface_coverage` when
  that threshold is set.
- **Missing binaries:** a missing `hits_path` or binary raises the same
  `RuntimeError` / `FileNotFoundError` as today.
- **Empty file:** an empty hit file yields no edges. That is not an
  error.

### CLI (`pandora/cli/app.py`)

**`similarity`:**

- `--output` is the hit TSV, and `<output>.search.json` holds the
  `SimilaritySearch` record.
- New flags: `--max-seqs`, `--exhaustive-search` (Foldseek), and
  `--precomputed-hits PATH`, which calls `load_similarity_search` with
  `--engine` instead of running a search.
- `--interface-residues` moves to `cluster`.

**`cluster`:**

- `--search` (a `*.search.json` file) replaces `--relationships`.
- The filter comes from `--hit-filter` (YAML file) or the shortcut flags
  `--min-score` (replaces `--threshold`) and `--min-coverage`. Shortcut
  flags override the same fields in the YAML.
- `--interface-residues` and `--pairs` are accepted.
- Outputs are unchanged: `clusters.json`, `cluster_provenance.json`,
  `cluster_pairs.json`.

`pandora/cli/README.md` is updated to match.

### `reproduce_dataset`

When `manifest.clustering.search` is set:

- **`origin="computed"`:** re-run the search with the recorded engine
  and parameters, writing to `output_dir / "hits.tsv"`.
- **`origin="precomputed"`:** re-use `hits_path` if it exists.
  Otherwise raise `ValueError` saying the search must be redone and
  loaded with `load_similarity_search`.

Then cluster with the recorded `hit_filter`. A manifest from before this
change (no `search`) raises the existing "did not record which
similarity method" `ValueError`.

## Testing

All tests are offline. MMseqs2 and Foldseek are replaced by shell
scripts that record their arguments and write canned TSVs, following
`tests/test_sequence_similarity.py`.

1. **Searches:** each search keeps the hit file at `hits_path`, returns
   the right `SimilaritySearch` fields, and passes `--max-seqs`,
   `--alignment-type` and `--exhaustive-search` to the binary.
2. **`iter_edges`:**
   - each threshold alone;
   - all four `coverage_of` values and all four `tm_normalisation`
     values;
   - inclusive bounds;
   - a pair whose best-scoring row fails coverage while a second row
     passes is yielded;
   - self-hits are skipped;
   - `source_id < target_id`;
   - interface coverage is computed at read time from
     `interface_residues`.
3. **Errors:** the malformed row (message has the line number), TM
   normalisation on MMseqs2, and the interface filter without residues.
4. **`load_similarity_search`:** accepts a well-formed file and rejects
   a wrong column count.
5. **Clustering:**
   - `cluster_edges` on plain pairs, ignoring and counting unknown IDs;
   - `cluster_similar_items` records `hit_filter` and `search`;
   - duplicate edges don't change clusters;
   - the existing clustering and pair-key tests are ported.
6. **CLI:** `similarity` writes the TSV and search record, including the
   `--precomputed-hits` path; `cluster` works with `--search` and
   `--hit-filter` / `--min-score`.
7. **`reproduce_dataset`:** a computed search is re-run; a precomputed
   one is re-used when the file exists and raises `ValueError` when it
   doesn't.

**Scale verification** (not in pytest): generate a synthetic Foldseek
TSV of about 20M rows over about 1M IDs, then run
`cluster_similar_items` with a filter that keeps about 1% of rows.
Report peak RSS and time. Peak RSS must stay well under 1 GB (union-find
over about 1M IDs plus line buffers). Compare against the pre-change
code on a 1M-row file, scaled by row count.

## Docs

- **`docs/usage/similarity.md`:**
  - the new flow (search → `SimilaritySearch` → `HitFilter` → cluster);
  - a migration note from the list API;
  - precomputed searches;
  - a warning about the `max_seqs` defaults;
  - the TM normalisation and coverage choices.
- **`docs/reference/policies.md`:** a `HitFilter` section.
- **`docs/policies/leakage.yaml`:** reference `HitFilter`.
- **`pandora/cli/README.md`:** the new flags.
- **Recipes and examples:** `docs/recipes/ppi-01.md`,
  `examples/ppi_dataset_pipeline.py` and `examples/dataset_pipeline.py`
  (the last one via `cluster_edges`).
- **Diagrams and site:** regenerate the ERD diagrams and run a strict
  docs build.

## Out of scope

- Community detection and bounded-depth deleaking (sub-project 6).
- Pair-level (dimer) similarity (sub-project 5).
- A unique-pairs or best-hit helper. Use a `set` until one is needed.

## Done when

- The CLAUDE.md review checklist passes.
- The scale verification meets its memory bound.
- The examples run against the new API (MMseqs2/Foldseek-dependent
  examples are checked with the fake binaries).
