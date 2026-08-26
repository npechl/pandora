# Similarity

`pandora.similarity` builds pairwise similarity networks and turns them
into leakage-safe dataset splits. `compute_sequence_similarity()` and
`compute_structure_similarity()` shell out to external binaries
(MMseqs2 and Foldseek respectively — install separately, see
`CLAUDE.md`); `cluster_similar_items()` and `partition_dataset()` are
pure Python and only depend on the resulting `SimilarityRelationship`
objects, not how they were computed. See
[Functions](../reference/functions.md#pandora.similarity) for full
signatures.

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

`compute_sequence_similarity()` runs an all-vs-all `mmseqs easy-search`.
It accepts a `{id: sequence}` mapping — `pandora.datasets.entry_sequences()`
builds exactly that from a batch of structures (one representative
sequence per entry, its longest polymer chain).

```python
from pandora.datasets import entry_sequences
from pandora.similarity import compute_sequence_similarity

sequences = entry_sequences(structures)
relationships = compute_sequence_similarity(sequences, sensitivity=5.7)

for r in relationships:
    print(r.source_id, r.target_id, round(r.score, 3), r.method.engine)
# 104M 112M 0.993 MMseqs2
# 118L 138L 0.987 MMseqs2
```

Only related pairs come back — unrelated pairs simply don't appear
(`104M`/`118L`, `1AYI`/anything, etc. are absent above).

## Structural similarity (Foldseek)

`compute_structure_similarity()` runs an all-vs-all `foldseek
easy-search` over structure files on disk — pass a `{id: path}` mapping,
or a directory of files to use their filenames as ids directly.

```python
from pathlib import Path
from pandora.export import structure_to_mmcif
from pandora.similarity import compute_structure_similarity

output_dir = Path("./datasets/output/struct/")
paths = {
    entry_id: structure_to_mmcif(structure, output_dir / f"{entry_id}.cif")
    for entry_id, structure in structures.items()
}
relationships = compute_structure_similarity(paths, sensitivity=9.5)

for r in relationships:
    print(r.source_id, r.target_id, round(r.score, 3))
# 104M 112M 0.999
# 104M 118L 0.265
# 104M 138L 0.257
# ... (9 total — structural similarity finds distant relationships
#      sequence similarity misses, since a fold can be conserved
#      long after sequence identity drops)
```

`score` is the alignment's TM-score; `identity` is its fraction of
identical aligned residues.

### Keep ids consistent between stages

!!! warning
    Both functions derive `source_id`/`target_id` from whatever id you
    give them. If you write structures with `structure_to_mmcif()`
    yourself (rather than the CLI's `pandora similarity`/`pandora
    cluster`, which already keep this consistent) and pass a *directory*
    instead of an explicit `{id: path}` mapping, make sure the filenames
    you feed `compute_structure_similarity()` use the same casing as the
    ids you pass to `cluster_similar_items()` below — a mismatch means
    `cluster_similar_items()` silently drops every relationship instead
    of erroring.

## Clustering

`cluster_similar_items()` groups ids into connected-component clusters:
two items land in the same cluster iff connected through a chain of
relationships scoring at or above `threshold`. Items with no qualifying
edges become their own singleton cluster.

```python
from pandora.similarity import cluster_similar_items

clusters, provenance = cluster_similar_items(
    list(structures), relationships, threshold=0.9
)
for cluster in clusters:
    print(cluster.components)
# ['104M', '112M']
# ['118L', '138L']
# ['1AYI']
```

## Leakage-safe partitioning

`partition_dataset()` assigns whole clusters to train/val/test —
similar structures never end up split across partitions, since a
whole cluster moves together to whichever split is furthest from its
target share.

```python
from pandora.similarity import partition_dataset

splits, provenance = partition_dataset(
    clusters, pct_train=0.6, pct_val=0.2, pct_test=0.2
)
print(splits)
# {'train': ['104M', '112M', '118L', '138L'], 'val': ['1AYI'], 'test': []}
```

Pass `keep_similar_items=False` to instead divide each cluster's
members proportionally across splits — only do this if leakage between
splits genuinely doesn't matter for your use case.
