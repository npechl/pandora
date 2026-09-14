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

=== "`library`"

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

=== "`cli`"

    ```bash
    pandora similarity --input-dir deduped/ --engine mmseqs2 --sensitivity 5.7 --output relationships.json
    # computed 2 relationships -> relationships.json
    ```

    `entry_sequences()` is called internally over every `*.cif` in
    `--input-dir`.

## Structural similarity (Foldseek)

`compute_structure_similarity()` runs an all-vs-all `foldseek
easy-search` over structure files on disk — pass a `{id: path}` mapping,
or a directory of files to use their filenames as ids directly.

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

=== "`cli`"

    ```bash
    pandora similarity --input-dir deduped/ --engine foldseek --sensitivity 9.5 --output relationships.json
    # computed 9 relationships -> relationships.json
    ```

    `--input-dir` is passed straight to `compute_structure_similarity()`
    as a directory — the `*.cif` filenames (uppercased stem) become the
    ids, which is why [keeping ids consistent between
    stages](#keep-ids-consistent-between-stages) is automatic here.

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

### Interface-restricted coverage (PPI pairs)

For PPI work, whole-chain coverage can hide that two complexes only
resemble each other away from the interface (or vice versa). Pass
`interface_residues={id: {positions...}}` to restrict coverage to each
item's interface residues instead:

```python
relationships = compute_structure_similarity(
    paths,
    interface_residues={
        "104M_A": {4, 5, 12},
        "112M_A": {4, 6, 13},
    },
)
for r in relationships:
    print(r.source_id, r.target_id, r.coverage, r.interface_coverage)
```

`interface_coverage` is only set on a relationship when both its items
appear in `interface_residues` — otherwise it's `None` and `coverage`
(whole-chain) is unaffected either way.

!!! warning
    `interface_residues` positions must match Foldseek's own 1-indexed
    residue numbering for that item's structure file — the order residues
    appear *in the file*, not `label_seq_id`. These only coincide for a
    single-chain, gap-free file; a missing loop shifts every residue after
    it, and a multi-chain file adds chain-order ambiguity on top. Don't
    hand-build this mapping — derive it, below.

#### Deriving `interface_residues` correctly

`annotate_chain_interfaces()` reports residues as `label_asym_id:label_seq_id`
strings (mmCIF's own numbering, gaps and all) — not Foldseek positions.
Bridge the two with `export_chain_mmcif()` (writes one chain per file,
removing the multi-chain ambiguity) and `interface_residues_from_annotation()`
(converts `label_seq_id` to that file's actual residue order):

```python
from pandora.annotations import annotate_chain_interfaces
from pandora.export import export_chain_mmcif
from pandora.similarity import (
    chain_item_id,
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

relationships = compute_structure_similarity(
    paths, interface_residues=interface_residues
)
```

The CLI's `--interface-residues` takes the same mapping as a JSON file —
write `interface_residues` from the recipe above as `{item_id: [position,
...]}` (`json.dumps` needs lists, not sets), point `--input-dir` at the
`chains/` directory `paths` was written into, then:

```python
import json

output_dir.joinpath("interfaces.json").write_text(
    json.dumps({k: sorted(v) for k, v in interface_residues.items()})
)
```

```bash
pandora similarity --input-dir chains/ --engine foldseek --interface-residues interfaces.json --output relationships.json
```

## Clustering

`cluster_similar_items()` groups ids into connected-component clusters:
two items land in the same cluster iff connected through a chain of
relationships scoring at or above `threshold`. Items with no qualifying
edges become their own singleton cluster.

=== "`library`"

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

=== "`cli`"

    ```bash
    pandora cluster --input-dir deduped/ --relationships relationships.json --threshold 0.9 --output clusters.json
    # 3 clusters at threshold=0.9 -> clusters.json
    ```

    Item ids come from `--input-dir`'s `*.cif` filenames (uppercased),
    not from `relationships.json` — this is the ordering
    [Keep ids consistent between stages](#keep-ids-consistent-between-stages)
    warns about.

### Paired cluster keys (PPI pairs)

`pair_cluster_keys()` looks up each side of a PPI pair in item-level
clusters and returns a `(cluster_id_1, cluster_id_2)` key per pair
(Pinder-style `{cluster_id_R, cluster_id_L}`) — useful for deduplicating
PPI pairs by which fold-pair they represent, not just which structure
they came from:

```python
from pandora.similarity import pair_cluster_keys

keys = pair_cluster_keys([("104M", "112M")], clusters)
# {('104M', '112M'): ('104M', '112M')}
```

Each cluster's own lexicographically-smallest member id is used as its
key, so it's stable regardless of clustering order.

The CLI's `cluster` subcommand takes the same pairs as a JSON file and
writes the keys alongside its usual output:

```bash
pandora cluster --input-dir deduped/ --relationships relationships.json --threshold 0.9 --pairs pairs.json --output clusters.json
# pairs.json: [["104M", "112M"]]
# 1 clusters at threshold=0.9 -> clusters.json
# 1 paired cluster keys -> cluster_pairs.json
```

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
