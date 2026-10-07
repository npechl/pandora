# Similarity

`pandora.similarity` finds similar items and turns them into leakage-safe dataset splits. `compute_sequence_similarity()` and `compute_structure_similarity()` run MMseqs2 and Foldseek (external binaries, installed separately — see `CLAUDE.md`), write the tool's hits to a file on disk, and return a small `SimilaritySearch` record pointing at it. A `HitFilter` policy decides which hits count as similar, and `cluster_similar_items()` streams the hit file through it, so memory doesn't grow with the number of hits — an all-vs-all search over the whole PDB can be clustered on a workstation. `cluster_edges()` clusters any pairs you build yourself, and `partition_dataset()` turns clusters into splits. See [Functions](../reference/functions.md#pandora.similarity) for full signatures.

```python
from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.schemas.canonicalisation import canonicalisationPolicy

policy = canonicalisationPolicy(
    policy_id="p", policy_name="p", policy_version="1.0.0"
)
entry_ids = ["104m", "112m", "118l", "138l", "1ayi"]
structures = {}
for entry_id in entry_ids:
    structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
    canonical, _, _ = canonicalise_structure(structure, policy)
    structures[canonical.entry_id] = canonical
```

## Sequence similarity (MMseqs2)

`compute_sequence_similarity()` runs an all-vs-all `mmseqs easy-search` and keeps its output at the path you give. It accepts a `{id: sequence}` mapping — `pandora.datasets.entry_sequences()` builds exactly that from a batch of structures (one representative sequence per entry, its longest polymer chain).

=== "`library`"

    ```python
    from pandora.datasets import entry_sequences
    from pandora.similarity import compute_sequence_similarity

    sequences = entry_sequences(structures)
    search = compute_sequence_similarity(
        sequences, "datasets/output/hits/mmseqs.tsv", max_seqs=300
    )
    print(search.engine, search.hits_path, search.parameters["max_seqs"])
    # MMseqs2 datasets/output/hits/mmseqs.tsv 300
    ```

    The hit file is MMseqs2's own TSV, one row per hit, with the columns `query, target, fident, alnlen, qcov, tcov`. Only related pairs appear, in both directions, plus self-hits:

    ```text
    104M  104M  1.000  153  1.000  1.000
    104M  112M  0.993  153  1.000  0.994
    112M  104M  0.993  153  0.994  1.000
    118L  138L  0.987  164  1.000  1.000
    ...
    ```

=== "`cli`"

    ```bash
    pandora similarity --input-dir deduped/ --engine mmseqs2 --max-seqs 300 --output mmseqs.tsv
    # hits -> mmseqs.tsv; search record -> mmseqs.tsv.search.json
    ```

    `entry_sequences()` is called internally over every `*.cif` in `--input-dir`. `mmseqs.tsv.search.json` is the `SimilaritySearch` record `pandora cluster` reads.

!!! warning "`max_seqs` caps how many hits each query keeps"
    Both tools stop after a fixed number of hits per query (MMseqs2 300, Foldseek 1000 by default) and skip candidates in a fast prefilter. In large families (antibodies, kinases) similar pairs beyond the cap are silently missing, and a missing pair counts as "not similar" — that's leakage across splits. For split-building, raise `max_seqs`, or pass `exhaustive_search=True` to Foldseek. Whatever you use is recorded in `search.parameters`.

## Structural similarity (Foldseek)

`compute_structure_similarity()` runs an all-vs-all `foldseek easy-search` over structure files on disk — pass a `{id: path}` mapping, or a directory of files to use their filenames as ids directly.

=== "`library`"

    ```python
    from pathlib import Path
    from pandora.export import structure_to_mmcif
    from pandora.similarity import compute_structure_similarity

    output_dir = Path("./datasets/output/struct/")
    paths = {
        entry_id: structure_to_mmcif(structure, output_dir / f"{entry_id}.cif")
        for entry_id, structure in structures.items()
    }
    fsearch = compute_structure_similarity(
        paths, "datasets/output/hits/foldseek.tsv"
    )
    print(fsearch.engine, fsearch.version)
    # Foldseek 10.941cd33
    ```

    Structural similarity finds distant relationships sequence similarity misses (21 hit rows here against MMseqs2's 9), since a fold can be conserved long after sequence identity drops. The Foldseek hit file has these columns:

    | Column | Meaning |
    |---|---|
    | `query`, `target` | The two item ids |
    | `fident` | Fraction of identical aligned residues |
    | `alnlen` | Alignment length |
    | `qcov`, `tcov` | Fraction of the query / target covered by the alignment |
    | `alntmscore`, `qtmscore`, `ttmscore` | TM-score normalised by the alignment, the query, or the target |
    | `qstart`, `qend`, `tstart`, `tend` | Alignment range on each side (used for interface coverage) |

=== "`cli`"

    ```bash
    pandora similarity --input-dir deduped/ --engine foldseek --output foldseek.tsv
    # hits -> foldseek.tsv; search record -> foldseek.tsv.search.json
    ```

    `--input-dir` is passed straight to `compute_structure_similarity()` as a directory — the `*.cif` filenames (uppercased stem) become the ids, which is why [keeping ids consistent between stages](#keep-ids-consistent-between-stages) is automatic here. `--exhaustive-search` and `--max-seqs` map to the function's arguments.

### Keep ids consistent between stages

!!! warning
    Both functions take their ids from whatever you give them. If you write structures with `structure_to_mmcif()` yourself (rather than the CLI's `pandora similarity`/`pandora cluster`, which already keep this consistent) and pass a *directory* instead of an explicit `{id: path}` mapping, make sure the filenames you feed `compute_structure_similarity()` use the same casing as the ids you pass to `cluster_similar_items()` below. Edges naming an id outside `item_ids` are ignored; `ClusteringProvenance.n_edges_unknown_ids` counts them, so check it is 0.

## Filtering hits: `HitFilter`

A `HitFilter` decides which hit rows become similarity edges. A pair of items is an edge if **any single hit row** passes every threshold you set — so a pair whose best-scoring alignment covers too little of the chain still counts if another alignment of the pair passes. `None` thresholds aren't applied; the others are inclusive (`>=`). `iter_edges()` streams the passing rows:

```python
from pandora.schemas.similarity import HitFilter
from pandora.similarity import iter_edges

hit_filter = HitFilter(min_score=0.5, min_coverage=0.8)
for edge in iter_edges(fsearch, hit_filter):
    print(edge.source_id, edge.target_id, round(edge.score, 3), edge.coverage)
# 118L 138L 1.006 1.0
# 104M 112M 0.999 0.994
# 118L 138L 1.006 1.0
# ...
```

Each edge has `source_id < target_id`. A pair can appear more than once (once per passing alignment or direction); clustering doesn't mind, and `{(e.source_id, e.target_id) for e in edges}` gives unique pairs.

| Field | Default | Meaning |
|---|---|---|
| `min_score` | `None` | TM-score for Foldseek (see `tm_normalisation`), identity for MMseqs2 |
| `min_identity` | `None` | Fraction of identical aligned residues |
| `min_coverage` | `None` | Alignment coverage (see `coverage_of`) |
| `coverage_of` | `"both"` | `both` (the smaller of query and target coverage), `query`, `target`, or `either` (the larger) |
| `min_interface_coverage` | `None` | Fraction of each side's interface residues inside the alignment (Foldseek only, see below) |
| `tm_normalisation` | `"alignment"` | Foldseek score: `alignment` (`alntmscore`), `query` (`qtmscore`), `target` (`ttmscore`) or `max` (the larger of the two) |

`HitFilter` is a policy, so it can live in YAML (`docs/reference/policies.md`) and is recorded in `ClusteringProvenance`. The choices of `coverage_of` and `tm_normalisation` change results: state them explicitly when you describe a dataset.

### Interface-restricted coverage (PPI pairs)

For PPI work, whole-chain coverage can hide that two complexes only resemble each other away from the interface (or vice versa). Pass `interface_residues={id: {positions...}}` when reading a Foldseek search, and filter on `min_interface_coverage`:

```python
edges = iter_edges(
    fsearch,
    HitFilter(min_interface_coverage=0.5),
    interface_residues={"104M_A": {4, 5, 12}, "112M_A": {4, 6, 13}},
)
```

Interface coverage is computed when the hit file is read, from each row's alignment range, so you can try different interface definitions without re-running the search. A pair gets `interface_coverage` only when both its items appear in `interface_residues`; otherwise it's `None`, and the pair fails `min_interface_coverage`.

!!! warning
    `interface_residues` positions must match Foldseek's own 1-indexed residue numbering for that item's structure file — the order residues appear *in the file*, not `label_seq_id`. These only coincide for a single-chain, gap-free file; a missing loop shifts every residue after it, and a multi-chain file adds chain-order ambiguity on top. Don't hand-build this mapping — derive it, below.

#### Deriving `interface_residues` correctly

`annotate_chain_interfaces()` reports residues as `label_asym_id:label_seq_id` strings (mmCIF's own numbering, gaps and all) — not Foldseek positions. Bridge the two with `export_chain_mmcif()` (writes one chain per file, removing the multi-chain ambiguity) and `interface_residues_from_annotation()` (converts `label_seq_id` to that file's actual residue order):

```python
from pandora.annotations import annotate_chain_interfaces
from pandora.export import export_chain_mmcif
from pandora.similarity import (
    chain_item_id,
    cluster_similar_items,
    compute_structure_similarity,
    interface_residues_from_annotation,
)

interface_layers = {
    entry_id: annotate_chain_interfaces(structure)
    for entry_id, structure in structures.items()
}
interface_residues = interface_residues_from_annotation(
    structures, interface_layers
)

paths = {}
for entry_id, structure in structures.items():
    for chain_id in {a.label_asym_id for a in structure.atoms}:
        item_id = chain_item_id(entry_id, chain_id)
        paths[item_id] = export_chain_mmcif(
            structure, chain_id, output_dir / f"{item_id}.cif"
        )

chain_search = compute_structure_similarity(
    paths, "datasets/output/hits/chains.tsv"
)
clusters, provenance = cluster_similar_items(
    sorted(paths),
    chain_search,
    HitFilter(min_interface_coverage=0.5),
    interface_residues=interface_residues,
)
```

The CLI's `pandora cluster --interface-residues` takes the same mapping as a JSON file — write `interface_residues` from the recipe above as `{item_id: [position, ...]}` (`json.dumps` needs lists, not sets), run `pandora similarity` with `--input-dir` pointed at the `chains/` directory `paths` was written into, then:

```python
import json

output_dir.joinpath("interfaces.json").write_text(
    json.dumps({k: sorted(v) for k, v in interface_residues.items()})
)
```

```bash
pandora similarity --input-dir chains/ --engine foldseek --output chains.tsv
pandora cluster --input-dir chains/ --search chains.tsv.search.json --hit-filter filter.yaml --interface-residues interfaces.json --output clusters.json
# filter.yaml: min_interface_coverage: 0.5
```

## Precomputed searches

A whole-PDB search is often run on a cluster, outside Pandora. Run the tool with Pandora's columns — `--format-output` set to `",".join(HIT_COLUMNS[engine])` from `pandora.similarity.hits`:

| Engine | `--format-output` |
|---|---|
| MMseqs2 | `query,target,fident,alnlen,qcov,tcov` |
| Foldseek | `query,target,fident,alnlen,qcov,tcov,alntmscore,qtmscore,ttmscore,qstart,qend,tstart,tend` |

then wrap the file with `load_similarity_search()`, recording how it was run:

```python
from pandora.similarity import load_similarity_search

pre = load_similarity_search(
    "datasets/output/hits/foldseek.tsv",
    "Foldseek",
    version="10.941cd33",
    parameters={"max_seqs": 5000, "exhaustive_search": True},
)
print(pre.origin)
# precomputed
```

The CLI form is `pandora similarity --engine foldseek --precomputed-hits foldseek.tsv --output foldseek.tsv`, which writes `foldseek.tsv.search.json`. A wrong column count raises `ValueError` naming the file and line. `reproduce_dataset()` re-uses a precomputed hit file if it still exists at the recorded path, and otherwise asks you to redo the search.

## Clustering

`cluster_similar_items()` groups ids into connected-component clusters: two items land in the same cluster iff connected through a chain of edges that pass the `HitFilter`. Items with no edges become their own singleton cluster.

=== "`library`"

    ```python
    from pandora.similarity import cluster_similar_items

    clusters, provenance = cluster_similar_items(
        list(structures), search, HitFilter(min_score=0.9)
    )
    for cluster in clusters:
        print(cluster.components)
    # ['104M', '112M']
    # ['118L', '138L']
    # ['1AYI']
    ```

    `provenance` records the `hit_filter` and `search`, plus `n_edges` (4 here: each pair once per direction).

    For pairs you build yourself, `cluster_edges()` takes any iterable of `(id, id)` pairs:

    ```python
    from pandora.similarity import cluster_edges

    own_clusters, _ = cluster_edges(["a", "b", "c"], [("a", "b")])
    # [['a', 'b'], ['c']]
    ```

=== "`cli`"

    ```bash
    pandora cluster --input-dir deduped/ --search mmseqs.tsv.search.json --min-score 0.9 --output clusters.json
    # 3 clusters from 4 edges -> clusters.json
    ```

    `--hit-filter filter.yaml` loads a full `HitFilter`; `--min-score` and `--min-coverage` override its fields. At least one of the three is required: with no threshold every reported hit would become an edge. Item ids come from `--input-dir`'s `*.cif` filenames (uppercased) — this is the ordering [Keep ids consistent between stages](#keep-ids-consistent-between-stages) warns about.

### Paired cluster keys (PPI pairs)

`pair_cluster_keys()` looks up each side of a PPI pair in item-level clusters and returns a `(cluster_id_1, cluster_id_2)` key per pair (Pinder-style `{cluster_id_R, cluster_id_L}`) — useful for deduplicating PPI pairs by which fold-pair they represent, not just which structure they came from:

```python
from pandora.similarity import pair_cluster_keys

keys = pair_cluster_keys([("104M", "112M")], clusters)
# {('104M', '112M'): ('104M', '104M')}
```

Each cluster's own lexicographically-smallest member id is used as its key, so it's stable regardless of clustering order.

The CLI's `cluster` subcommand takes the same pairs as a JSON file and writes the keys alongside its usual output:

```bash
pandora cluster --input-dir deduped/ --search mmseqs.tsv.search.json --min-score 0.9 --pairs pairs.json --output clusters.json
# pairs.json: [["104M", "112M"]]
# 3 clusters from 4 edges -> clusters.json
# 1 paired cluster keys -> cluster_pairs.json
```

## Migrating from the list API (≤ 0.5.6)

Up to 0.5.6 the searches returned a list of `SimilarityRelationship` objects held in memory, which runs out of memory at whole-PDB scale. The searches now keep their hit file on disk:

| Before | After |
|---|---|
| `rels = compute_sequence_similarity(seqs)` | `search = compute_sequence_similarity(seqs, "hits.tsv")` |
| `compute_structure_similarity(paths, interface_residues=ir)` | `compute_structure_similarity(paths, "hits.tsv")`, then `interface_residues=ir` on `iter_edges` / `cluster_similar_items` |
| `cluster_similar_items(ids, rels, threshold=0.9)` | `cluster_similar_items(ids, search, HitFilter(min_score=0.9))` |
| hand-built `SimilarityRelationship` list + threshold | `cluster_edges(ids, ((r.source_id, r.target_id) for r in rels if r.score >= 0.9))` |
| `ClusteringProvenance.threshold` / `.similarity_method` | `.hit_filter` / `.search` |
| CLI `similarity --output relationships.json` | `--output hits.tsv` (+ `hits.tsv.search.json`) |
| CLI `cluster --relationships r.json --threshold 0.9` | `cluster --search hits.tsv.search.json --min-score 0.9` |
| CLI `cluster` with no `--threshold` (defaulted to 0.9) | no default: `cluster` refuses to run without `--min-score`, `--min-coverage` or `--hit-filter` |

Manifests written before this change still load, but `reproduce_dataset()` can't rebuild their clustering: it raises a `ValueError` saying the search and hit filter weren't recorded.

## Leakage-safe partitioning

`partition_dataset()` assigns whole clusters to train/val/test —
similar structures never end up split across partitions, since a
whole cluster moves together to whichever split is furthest from its
target share.

=== "`library`"

    ```python
    from pandora.similarity import partition_dataset

    splits, provenance = partition_dataset(
        clusters, pct_train=0.6, pct_val=0.2, pct_test=0.2
    )
    print(splits)
    # {'train': ['104M', '112M', '118L', '138L'], 'val': ['1AYI'], 'test': []}
    ```

    Pass `keep_similar_items=False` to instead divide each cluster's
    members proportionally across splits — only do this if leakage
    between splits genuinely doesn't matter for your use case.

=== "`cli`"

    ```bash
    pandora partition --clusters clusters.json --pct-train 0.6 --pct-val 0.2 --pct-test 0.2 --output splits.json
    # split sizes: {'train': 4, 'val': 1, 'test': 0} -> splits.json
    ```

    `--no-keep-similar-items` is the CLI form of `keep_similar_items=False`.
