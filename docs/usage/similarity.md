# Similarity

`pandora.similarity` finds similar items and turns them into leakage-safe dataset splits. Similarity is a [layer](annotation.md): a `similarity` table of hits between two items. Pandora's own MMseqs2 and Foldseek searches write one, and `pandora similarity import` turns anyone else's edge list (another tool, or your own script) into the same layer. `pandora cluster` clusters any similarity layer, and `partition_dataset()` turns clusters into splits. See [Functions](../reference/functions.md#pandora.similarity) for full signatures.

Everything here is built to run on the whole PDB archive: layers are Parquet, read batch by batch, thresholds are applied to whole columns at once, and the list of valid items comes from a step folder's tables, never from parsing structures again. Similarity layers need pyarrow, the `export` extra.

The examples below start from a `curated/` step folder holding 104M, 112M, 118L, 138L and 1AYI (one polymer chain each):

```bash
pandora canonicalise --input-dir raw/ --policy datasets/canonicalisation.yaml --output-dir canonical/
pandora curate --input-dir canonical/ --policy curate.yaml --output-dir curated/
```

## The `similarity` layer

One row per hit between two items:

| Column | Type | Required |
|---|---|---|
| `item_id_1`, `item_id_2` | string | yes |
| `score` | number | yes |
| `identity` | number | recognised if present |
| `coverage_1`, `coverage_2` | number | recognised if present (coverage of item 1 / item 2) |
| `start_1`, `end_1`, `start_2`, `end_2` | integer | recognised if present (1-based alignment range on each item, for [interface coverage](#interface-restricted-coverage-ppi-pairs)) |

Any other column is kept as given. Rows are not deduplicated or reordered: A→B and B→A may both be there, and self-hits are kept. Clustering is unaffected by either.

### Items {: #items }

Each layer holds one level of item, set by `items` in its description:

| `items` | ID form | Example | Valid items come from |
|---|---|---|---|
| `entry` | `{entry}` | `1ABC` | the step folder's `entries` table (entries with a structure), `--structures-dir` |
| `chain` | `{entry}_{chain}` | `1ABC_A` | the step folder's `chains` table, written by `pandora curate`, `--structures-dir` |
| `interface` | `{entry}_{chain1}_{chain2}` | `1ABC_A_B` | a [`residue_contacts` layer](annotation.md#contacts-are-one-layer): its pairs of polymer chains in contact, `--contacts` |

Chains are Pandora's chain IDs (`label_asym_id`). The entry part matches case-insensitively and is stored as Pandora spells it; chain IDs are case-sensitive. Interface IDs are stored with their two chains sorted, so `1ABC_B_A` becomes `1ABC_A_B`. An interface is an item only if a contact layer found it, so similarity never recomputes contacts. In Python, `item_id()` and `parse_item_id()` build and split these IDs.

A row naming something that isn't an item goes to `unmatched.parquet` with a `reason`: `bad_item_id` (doesn't fit the level's form), `entry_not_found`, `chain_not_found`, `interface_not_found`, or `bad_value` (a recognised column that isn't a number, or no `score`).

## Sequence similarity (MMseqs2)

`pandora similarity search --engine mmseqs2` runs an all-vs-all `mmseqs easy-search` over the step folder's chains, with sequences taken from its `chains` table, and writes the hits as a chain-level layer folder: `layer.parquet`, `unmatched.parquet` (if any) and `layer.json`. Symmetry copies (`copy_of` set in the `chains` table) are left out of the search: identical copies would take up the per-query hit slots and push out real homologues. They are still items, and `pandora cluster` puts each copy in its source chain's cluster.

=== "`cli`"

    ```bash
    pandora similarity search --engine mmseqs2 --structures-dir curated/ --max-seqs 300 --output-dir similarity/seq
    ```

    ```text
    searched: 9 hits, 0 unmatched -> similarity/seq
    ```

    `layer.parquet` holds the hits with the template's column names. `score` is MMseqs2's identity (`fident`):

    ```text
    item_id_1  item_id_2  identity  alnlen  coverage_1  coverage_2  score
    104M_A     104M_A     1.0       153     1.0         1.0         1.0
    104M_A     112M_A     0.993     153     1.0         0.994       0.993
    112M_A     104M_A     0.993     153     0.994       1.0         0.993
    118L_A     138L_A     0.987     164     1.0         1.0         0.987
    ...
    ```

    `layer.json` records the search as the layer's description: `{"name": "mmseqs2", "template": "similarity", "items": "chain", "tool": "MMseqs2", "tool_version": "18.8cc5c", "parameters": {"mmseqs_bin": "mmseqs", "sensitivity": 5.7, "max_seqs": 300, "mmseqs_options": []}}`, plus the counts and the step folder's summary as `previous_summary`.

=== "`library`"

    ```python
    from pandora.export import read_records
    from pandora.schemas.dataset import DatasetChainRecord
    from pandora.similarity import (
        chain_items,
        compute_sequence_similarity,
        item_id,
        match_similarity_batches,
        read_hit_tsv,
    )

    chains = read_records("curated/chains.parquet", DatasetChainRecord)
    items = chain_items((c.entry_id, c.chain_id) for c in chains)
    sequences = {
        item_id("chain", c.entry_id, c.chain_id): c.sequence
        for c in chains
        if c.sequence
        and not c.copy_of  # symmetry copies would crowd out homologues
    }

    description = compute_sequence_similarity(
        sequences, "datasets/output/hits/mmseqs.tsv", max_seqs=300
    )
    print(description.tool, description.tool_version, description.items)

    for matched, unmatched in match_similarity_batches(
        read_hit_tsv("datasets/output/hits/mmseqs.tsv", "MMseqs2"), items
    ):
        print(matched.num_rows, "matched,", unmatched.num_rows, "unmatched")
        print(matched.slice(0, 2).to_pylist())
    ```

    ```text
    MMseqs2 18.8cc5c chain
    9 matched, 0 unmatched
    [{'item_id_1': '104M_A', 'item_id_2': '104M_A', 'identity': 1.0, 'alnlen': 153, 'coverage_1': 1.0, 'coverage_2': 1.0, 'score': 1.0}, {'item_id_1': '104M_A', 'item_id_2': '112M_A', 'identity': 0.993, 'alnlen': 153, 'coverage_1': 1.0, 'coverage_2': 0.994, 'score': 0.993}]
    ```

    `compute_sequence_similarity()` writes MMseqs2's own TSV and returns the layer's `LayerDescription`. `read_hit_tsv()` streams that TSV as record batches in the template's columns, and `match_similarity_batches()` links each batch to the items. Write the matched batches with `pyarrow.parquet.ParquetWriter` to get a layer.

!!! warning "`--max-seqs` caps how many hits each query keeps"
    Both tools stop after a fixed number of hits per query (MMseqs2 300, Foldseek 1000 by default) and skip candidates in a fast prefilter. In large families (antibodies, kinases) similar pairs beyond the cap are silently missing, and a missing pair counts as "not similar": that's leakage across splits. For split-building, raise `--max-seqs`, or pass `--exhaustive-search` to Foldseek. Whatever you use is recorded in the layer's `parameters`.

## Structural similarity (Foldseek)

`pandora similarity search --engine foldseek` runs an all-vs-all `foldseek easy-search` over a folder of per-chain structure files, named `{entry}_{chain}.cif` so that the file names are chain item IDs. No command writes these files yet; `export_chain_mmcif()` does, for the chains in the step folder's `chains` table:

```python
from pathlib import Path

from pandora.cli._io import read_structures_dir
from pandora.export import export_chain_mmcif, read_records
from pandora.schemas.dataset import DatasetChainRecord
from pandora.similarity import item_id

structures = read_structures_dir(Path("curated"))
for chain in read_records("curated/chains.parquet", DatasetChainRecord):
    name = item_id("chain", chain.entry_id, chain.chain_id)
    export_chain_mmcif(
        structures[chain.entry_id],
        chain.chain_id,
        Path("chains") / f"{name}.cif",
    )
print(sorted(p.name for p in Path("chains").iterdir()))
```

```text
['104M_A.cif', '112M_A.cif', '118L_A.cif', '138L_A.cif', '1AYI_A.cif']
```

```bash
pandora similarity search --engine foldseek --input-dir chains/ --structures-dir curated/ --output-dir similarity/struct
```

```text
searched: 21 hits, 0 unmatched -> similarity/struct
```

Structural similarity finds distant relationships sequence similarity misses (21 hits here against MMseqs2's 9), since a fold can be conserved long after sequence identity drops. `score` is Foldseek's TM-score normalised by the alignment (`alntmscore`); the layer also keeps `alntmscore`, `qtmscore` and `ttmscore` (normalised by item 1 and item 2), `alnlen`, and the alignment ranges `start_1`, `end_1`, `start_2`, `end_2`:

```text
item_id_1  item_id_2  identity  coverage_1  coverage_2  qtmscore  ttmscore  start_1  end_1  start_2  end_2  score
118L_A     138L_A     0.987     1.0         1.0         0.9995    0.9995    1        162    1        162    1.006
118L_A     1AYI_A     0.161     0.148       0.36        0.1223    0.2083    112      135    13       43     0.4387
...
```

In Python, `compute_structure_similarity(paths, hits_path)` runs the search and `read_hit_tsv(hits_path, "Foldseek")` reads it, as for MMseqs2 above.

## Importing an edge list

`pandora similarity import` brings in similarity computed any other way: an R or Python script, a graph exported from igraph, or MMseqs2/Foldseek run on a cluster outside Pandora. Convert it to the template's columns (Parquet, or a TSV with a header row) and describe it once in YAML. A TSV from MMseqs2 or Foldseek needs a header with the template's names: `item_id_1`, `item_id_2`, `score`, and `identity`, `coverage_1`, ... for the columns you want recognised.

=== "`cli`"

    ```yaml
    # r_graph.yaml
    name: r_graph
    template: similarity
    items: entry
    tool: R script
    tool_version: "4.4.1"
    parameters:
      graph: igraph
    ```

    ```text
    item_id_1	item_id_2	score	method
    104m	112M	0.97	graph_r
    118L	138L	0.95	graph_r
    138L	9XYZ	0.40	graph_r
    ```

    ```bash
    pandora similarity import --table r_graph_edges.tsv --description r_graph.yaml --structures-dir curated/ --output-dir similarity/r_graph
    ```

    ```text
    imported: 2 rows, 1 unmatched -> similarity/r_graph
    ```

    `layer.parquet` holds `104M`–`112M` (the entry ID as Pandora spells it) and `118L`–`138L`, with `method` kept as given. `unmatched.parquet` holds the `9XYZ` row with reason `entry_not_found`. `--source-files` copies the original output into `source/`, as in [`annotate import`](annotation.md). A description without `items`, or a table without `item_id_1`, `item_id_2` or `score`, stops the command before anything is written.

=== "`library`"

    ```python
    from pandora.similarity import (
        entry_items,
        match_similarity_batches,
        read_similarity_table,
    )

    items = entry_items(["104M", "112M", "118L", "138L", "1AYI"])
    for matched, unmatched in match_similarity_batches(
        read_similarity_table("r_graph_edges.tsv"), items
    ):
        print(
            matched.column("item_id_1").to_pylist(),
            unmatched.column("reason").to_pylist(),
        )
    ```

    ```text
    ['104M', '118L'] ['entry_not_found']
    ```

### Interface items

For interfaces, the valid items are the chain pairs of a contact layer, given with `--contacts`. Here `annotations/interfaces` is the `chain_interfaces` layer of 1A2C, whose interfaces are A–B, B–C and B–D:

```bash
pandora annotate compute --layer chain_interfaces --structures-dir curated2/ --output-dir annotations/interfaces
pandora similarity import --table iface_edges.tsv --description iface.yaml --contacts annotations/interfaces --output-dir similarity/iface
```

```text
computed: 120 rows for 1 entries -> annotations/interfaces
imported: 1 rows, 1 unmatched -> similarity/iface
```

With `iface.yaml` setting `items: interface`, the row `1A2C_B_A`–`1A2C_B_C` is stored as `1A2C_A_B`–`1A2C_B_C`, and `1A2C_A_C`–`1A2C_B_D` is unmatched with reason `interface_not_found`: chains A and C are both in the entry but never touch.

## Filtering hits: `HitFilter`

A `HitFilter` decides which rows of a similarity layer become edges. A pair of items is an edge if **any single row** passes every threshold you set, so a pair whose best-scoring alignment covers too little of the chain still counts if another alignment of the pair passes. `None` thresholds aren't applied; the others are inclusive (`>=`). A null or NaN value fails. A threshold on a column the layer doesn't have is an error, raised before any row is read.

| Field | Default | Meaning |
|---|---|---|
| `min_score` | `None` | Minimum value of `score_column` |
| `score_column` | `"score"` | The column `min_score` tests: `score` (identity for Pandora's MMseqs2 layers, alignment TM-score for its Foldseek layers), or any number column, e.g. `qtmscore` |
| `min_identity` | `None` | Minimum `identity` |
| `min_coverage` | `None` | Minimum alignment coverage (see `coverage_of`) |
| `coverage_of` | `"both"` | `both` (the smaller of `coverage_1` and `coverage_2`), `side_1`, `side_2`, or `either` (the larger) |
| `min_interface_coverage` | `None` | Fraction of each side's interface residues inside the alignment (see below) |

`HitFilter` is a policy, so it can live in YAML ([Policies](../reference/policies.md)) and is recorded in `ClusteringProvenance`. The choices of `coverage_of` and `score_column` change results: state them explicitly when you describe a dataset.

`similarity_edges()` streams the passing rows of a layer as `Edge`s, each with `source_id < target_id`. A pair can appear more than once (once per passing row); clustering doesn't mind. Self-hits are skipped.

### Interface-restricted coverage (PPI pairs)

For PPI work, whole-chain coverage can hide that two complexes only resemble each other away from the interface (or vice versa). On a chain layer with alignment ranges (a Foldseek layer), `pandora cluster --contacts` takes a `residue_contacts` layer folder and, with `--structures-dir`, works out each chain's interface residues in the per-chain files' own numbering (`interface_residues_from_annotation()`), so `min_interface_coverage` can be tested:

```bash
pandora cluster --similarity similarity/struct2 --structures-dir curated2/ --contacts annotations/interfaces --hit-filter filter.yaml --output clusters/interface_cov.json
# filter.yaml: min_interface_coverage: 0.5
```

```text
4 clusters from 0 edges -> clusters/interface_cov.json
```

(1A2C's four chains share no fold, so the only hits are self-hits.) A pair gets an interface coverage only when both its chains have interface residues; otherwise it fails `min_interface_coverage`. Coverage is computed when the layer is read, so you can try different interface definitions without re-running the search. This is the one step that parses the structures again, so it doesn't reach whole-archive scale.

!!! warning
    Interface positions must match Foldseek's own 1-indexed residue numbering for each per-chain file: the order residues appear *in the file*, not `label_seq_id`. These only coincide for a gap-free chain. Don't hand-build the mapping; `interface_residues_from_annotation()` derives it from the structures `export_chain_mmcif()` wrote.

## Clustering

`pandora cluster` groups items into connected-component clusters: two items land in the same cluster if they are connected through a chain of edges that pass the `HitFilter`. Every valid item appears exactly once, so items with no hits are single-item clusters. Edges naming an item outside the list are dropped and counted (`n_edges_unknown_ids`), so you can search the whole archive once and cluster any curated subset of it.

=== "`cli`"

    ```bash
    pandora cluster --similarity similarity/seq --structures-dir curated/ --min-score 0.9 --output clusters/chains.json
    ```

    ```text
    3 clusters from 4 edges -> clusters/chains.json
    ```

    ```text
    ['104M_A', '112M_A']
    ['118L_A', '138L_A']
    ['1AYI_A']
    ```

    `--hit-filter filter.yaml` loads a full `HitFilter`; `--min-score` and `--min-coverage` override its fields. At least one of the three is required: with no threshold every hit would become an edge. `cluster_provenance.json`, next to `--output`, records each layer under `layers` (its `layer.json`, item level, hit filter and edge count) and the total edge and cluster counts.

=== "`library`"

    ```python
    from pandora.export import read_records
    from pandora.schemas.dataset import DatasetChainRecord
    from pandora.schemas.similarity import HitFilter
    from pandora.similarity import (
        chain_items,
        cluster_edges,
        read_similarity_table,
        similarity_edges,
    )

    chains = read_records("curated/chains.parquet", DatasetChainRecord)
    items = chain_items((c.entry_id, c.chain_id) for c in chains)
    edges = similarity_edges(
        read_similarity_table("similarity/seq/layer.parquet"),
        HitFilter(min_score=0.9),
    )
    clusters, provenance = cluster_edges(sorted(items.ids.values()), edges)
    print([c.components for c in clusters], provenance.n_edges)
    ```

    ```text
    [['104M_A', '112M_A'], ['118L_A', '138L_A'], ['1AYI_A']] 4
    ```

    `cluster_edges()` takes any iterable of `(id, id)` pairs, so it also clusters pairs you build yourself. The CLI additionally links each symmetry copy to its source chain.

### Clustering entries from a chain layer

Splits are usually made by entry, so that chains of one entry never land in different splits. `--items entry` on a chain layer rolls the chains up to their entries: two entries share a cluster when **any** of their chains are similar. This is stricter against leakage than comparing one representative chain per entry, and needs only the one chain search.

=== "`cli`"

    ```bash
    pandora cluster --similarity similarity/seq --structures-dir curated/ --items entry --min-score 0.9 --output clusters/entries.json
    ```

    ```text
    3 clusters from 4 edges -> clusters/entries.json
    ```

    ```text
    ['104M', '112M']
    ['118L', '138L']
    ['1AYI']
    ```

    `--items` may differ from the layer's own level only on a chain layer, which rolls up to `entry` or `interface` (below); anything else is an error.

=== "`library`"

    ```python
    from pandora.similarity import entry_edges, entry_items

    entries = entry_items(sorted({c.entry_id for c in chains}))
    edges = similarity_edges(
        read_similarity_table("similarity/seq/layer.parquet"),
        HitFilter(min_score=0.9),
    )
    clusters, _ = cluster_edges(sorted(entries.ids.values()), entry_edges(edges))
    print([c.components for c in clusters])
    ```

    ```text
    [['104M', '112M'], ['118L', '138L'], ['1AYI']]
    ```

An interface layer clusters interfaces, with the items from `--contacts`:

```bash
pandora cluster --similarity similarity/iface --contacts annotations/interfaces --min-score 0.5 --output clusters/iface.json
```

```text
2 clusters from 1 edges -> clusters/iface.json
```

```text
['1A2C_A_B', '1A2C_B_C']
['1A2C_B_D']
```

### Clustering interfaces from a chain layer, and several layers

`--items interface` on a chain layer rolls the chains up to the interfaces of `--contacts`: two interfaces share a cluster when **any** chain of one is similar to any chain of the other. Interfaces that share a chain are always joined, since a chain matches itself, so every interface of one assembly whose chains connect lands in one cluster. A hit involving a chain that is in no interface is dropped: a chain that is not part of any dimer must not join two dimers it happens to resemble.

`--similarity` can be repeated, and an edge from any layer links two items. This is how a dimer similarity built from several checks (say chain fold *or* interface similarity) is clustered in one go. `--min-score`, `--min-coverage` and `--hit-filter` are each given once for every layer, or once per `--similarity` in the same order. Each layer must be at the `--items` level or be a chain layer that rolls up to it.

The examples use `1a3n` (deoxyhemoglobin) from [Comparing interfaces](#comparing-interfaces), with its interface layer `similarity/iface_jaccard` and a Foldseek layer `similarity/struct` over its four chains.

=== "`cli`"

    ```bash
    pandora cluster --items interface --contacts annotations/interfaces \
      --similarity similarity/iface_jaccard --min-score 0.1 \
      --similarity similarity/struct --min-score 0.5 \
      --structures-dir curated/ --output clusters/dimers.json
    jq -c '.[].components' clusters/dimers.json
    jq -c '.layers[] | [.similarity_layer, .items, .n_edges]' clusters/cluster_provenance.json
    ```

    ```text
    1 clusters from 18 edges -> clusters/dimers.json
    ["1A3N_A_B","1A3N_A_C","1A3N_A_D","1A3N_B_C","1A3N_B_D"]
    ["similarity/iface_jaccard/layer.json","interface",2]
    ["similarity/struct/layer.json","chain",12]
    ```

    The interface layer alone gives three clusters; α and β are both globins (TM-score above 0.9), so the chain layer joins them. The 18 edges are the interface layer's 2, the 10 chain hits between chains of different interfaces, and 6 links between interfaces sharing a chain. `--structures-dir` is optional here; when given, symmetry copies are linked to their source chain as for chain clustering.

=== "`library`"

    ```python
    import pandas as pd

    from pandora.schemas.similarity import HitFilter
    from pandora.similarity import (
        cluster_edges,
        interface_edges,
        interface_items,
        read_similarity_table,
        similarity_edges,
    )

    contacts = pd.read_parquet("annotations/interfaces/layer.parquet")
    interfaces = sorted(interface_items(contacts.to_dict("records")).ids.values())
    edges = similarity_edges(
        read_similarity_table("similarity/struct/layer.parquet"),
        HitFilter(min_score=0.5),
    )
    clusters, _ = cluster_edges(interfaces, interface_edges(edges, interfaces))
    print([c.components for c in clusters])
    ```

    ```text
    [['1A3N_A_B', '1A3N_A_C', '1A3N_A_D', '1A3N_B_C', '1A3N_B_D']]
    ```

    To join several layers, chain their edge streams (`itertools.chain`) before `cluster_edges()`, rolling up only the chain-level ones.

### Paired cluster keys (PPI pairs)

`pair_cluster_keys()` looks up each side of a PPI pair in item-level clusters and returns a `(cluster_id_1, cluster_id_2)` key per pair (Pinder-style `{cluster_id_R, cluster_id_L}`), useful for deduplicating PPI pairs by which fold-pair they represent, not just which structure they came from:

```python
from pandora.similarity import pair_cluster_keys

keys = pair_cluster_keys([("104M_A", "112M_A")], clusters)
# {('104M_A', '112M_A'): ('104M_A', '104M_A')}
```

Each cluster's own lexicographically smallest member is its key, so keys are stable regardless of clustering order. On the command line, `pandora cluster --pairs pairs.json` (a JSON array of `[item_id_1, item_id_2]`) writes the keys to `cluster_pairs.json` next to `--output`.

## Comparing interfaces

Two dimers are redundant when they are the same two proteins touching at the same spot. `pandora similarity interfaces` scores that for every pair of interfaces in a contact layer, and writes the scores as a similarity layer with `items: interface`, ready for `pandora cluster --items interface`.

- **Same two proteins.** Two interfaces are compared only when their chains fall in the same chain clusters, for example 90% sequence clusters from `similarity search` + `cluster` (a hit filter with `min_identity: 0.9`, `min_coverage: 0.9` and `coverage_of: both`). Clusters only decide which pairs are looked at. Each chain must also be at least `--min-identity` (default 0.9) identical to its counterpart on a pairwise alignment, since connected components can join chains less similar than that.
- **Both chain orders.** For a homomer, (A, B) lines up with (A′, B′) and with (B′, A′). Both are scored and the better one is kept; `chain_order` says which (`direct` or `swapped`).
- **Per side.** `jaccard_1` compares the interface residues of the first chain of `item_id_1` with those of its counterpart; `jaccard_2` the other chain. Residues are matched through the alignment, and a residue it can't place counts as not shared.
- **`score` is the smaller side.** Both proteins must touch with the same residues. One protein reusing a patch to bind its partner at two different places scores 0.
- **Every compared pair is written**, zeros included, so you can look at the distribution before choosing a cut-off. `cluster --min-score` keeps scores at or above it.

The residue numbers must be positions in the chain's sequence: canonicalise with `residue_numbering: preserve` (renumbering moves them). Each standard residue in the contact rows is checked against the sequence (modified residues such as `MSE` are not, since the sequence may give their parent), and an interface that disagrees is skipped and counted as `interfaces_sequence_mismatch`, never matched by guess.

The examples use `1a3n` (deoxyhemoglobin): chains A and B are α, C and D are β. α1β1 and α2β2 (`A_C`, `B_D`) are the same interface; so are the two α1β2 contacts (`A_D`, `B_C`).

=== "`cli`"

    The chain clusters are written by hand here; normally they come from `pandora cluster` on a sequence layer.

    ```bash
    pandora annotate compute --layer chain_interfaces --structures-dir curated/ \
      --output-dir annotations/interfaces
    cat > clusters/seq90.json <<'EOF'
    [
      {"components": ["1A3N_A", "1A3N_B"], "n_components": 2},
      {"components": ["1A3N_C", "1A3N_D"], "n_components": 2}
    ]
    EOF
    pandora similarity interfaces --contacts annotations/interfaces \
      --chain-clusters clusters/seq90.json --output-dir similarity/iface_jaccard
    pandora cluster --similarity similarity/iface_jaccard --items interface \
      --contacts annotations/interfaces --min-score 0.1 --output clusters/dimers.json
    jq -c '.counts' similarity/iface_jaccard/layer.json
    jq -c '.[].components' clusters/dimers.json
    ```

    ```text
    computed: 135 rows for 1 entries -> annotations/interfaces
    scored: 6 interface pairs (4 zero, 0 below identity) from 5 interfaces -> similarity/iface_jaccard
    3 clusters from 2 edges -> clusters/dimers.json
    {"alignments":0,"interfaces":5,"interfaces_sequence_mismatch":0,"interfaces_without_cluster":0,"interfaces_without_sequence":0,"pairs_below_identity":0,"pairs_compared":6,"pairs_swapped":0,"pairs_zero_score":4,"rows_ignored":0,"rows_in":135}
    ["1A3N_A_B"]
    ["1A3N_A_C","1A3N_B_D"]
    ["1A3N_A_D","1A3N_B_C"]
    ```

    The structures folder and its `chains` table (sequences) are found from the contact layer's `layer.json`, as for `annotate label-interfaces`. `alignments` is 0 because identical sequences match by position without aligning. `layer.json` records the contact layer, the structures folder, the clusters file and `--min-identity` under `policy.parameters`.

=== "`library`"

    ```python
    import json

    import pandas as pd

    from pandora.export import read_records
    from pandora.schemas.annotation import ResidueContactRow
    from pandora.schemas.similarity import SimilarityCluster
    from pandora.similarity import interface_similarity

    contacts = read_records(
        "annotations/interfaces/layer.parquet", ResidueContactRow
    )
    with open("clusters/seq90.json") as f:
        clusters = [SimilarityCluster.model_validate(c) for c in json.load(f)]
    chains = pd.read_parquet("curated/chains.parquet")
    sequences = {
        (row.entry_id.upper(), row.chain_id): row.sequence
        for row in chains.itertuples()
    }
    rows, counts = interface_similarity(contacts, clusters, sequences)
    columns = ["item_id_1", "item_id_2", "jaccard_1", "jaccard_2", "score"]
    print(pd.DataFrame([row.model_dump() for row in rows])[columns].round(2))
    ```

    ```text
      item_id_1 item_id_2  jaccard_1  jaccard_2  score
    0  1A3N_A_C  1A3N_A_D       0.00       0.07   0.00
    1  1A3N_A_C  1A3N_B_C       0.00       0.07   0.00
    2  1A3N_A_C  1A3N_B_D       0.84       1.00   0.84
    3  1A3N_A_D  1A3N_B_C       1.00       0.93   0.93
    4  1A3N_A_D  1A3N_B_D       0.00       0.07   0.00
    5  1A3N_B_C  1A3N_B_D       0.00       0.07   0.00
    ```

    `1A3N_A_C` against `1A3N_A_D` shares a little of β's surface (`jaccard_2` 0.07) but none of α's, so its score is 0: α1β1 and α1β2 are different interfaces.

## Migrating from hit-file searches

Up to this release the searches kept a raw hit TSV with a `SimilaritySearch` record, and `pandora cluster` read only that. Similarity is now a layer:

| Before | After |
|---|---|
| `pandora similarity --engine mmseqs2 --input-dir d/ --output hits.tsv` | `pandora similarity search --engine mmseqs2 --structures-dir curated/ --output-dir similarity/seq` (searches chains) |
| `pandora similarity --precomputed-hits hits.tsv ...` | `pandora similarity import --table hits.tsv ...`, with a header row in the template's names |
| `pandora cluster --input-dir d/ --search hits.tsv.search.json` | `pandora cluster --similarity similarity/seq --structures-dir curated/` |
| entry-level search (each entry's longest chain, `entry_sequences()`) | chain search plus `--items entry` / `entry_edges()` |
| `compute_*_similarity()` returning `SimilaritySearch` | returning the layer's `LayerDescription`; read the TSV with `read_hit_tsv()` |
| `iter_edges()`, `cluster_similar_items()` | `similarity_edges()` + `cluster_edges()` |
| `HitFilter.tm_normalisation` | `HitFilter.score_column` (`qtmscore`, `ttmscore`; no `max`) |
| `HitFilter.coverage_of: query / target` | `side_1` / `side_2` |
| `ClusteringProvenance.search` | `ClusteringProvenance.layers[].similarity_layer` |
| `ClusteringProvenance.hit_filter`, `.similarity_layer` (one layer) | `ClusteringProvenance.layers` (one entry per `--similarity`) |
| `annotate_pairwise_sequence_identity()` | removed: use an MMseqs2 search |

`reproduce_dataset()` no longer replays similarity or clustering.

## Leakage-safe partitioning

`partition_dataset()` assigns whole clusters to train/val/test: similar items never end up split across partitions, since a whole cluster moves together to whichever split is furthest from its target share.

=== "`library`"

    ```python
    from pandora.similarity import partition_dataset

    splits, provenance = partition_dataset(
        clusters, pct_train=0.6, pct_val=0.2, pct_test=0.2
    )
    print(splits)
    ```

    ```text
    {'train': ['104M', '112M', '118L', '138L'], 'val': ['1AYI'], 'test': []}
    ```

    Pass `keep_similar_items=False` to instead divide each cluster's members proportionally across splits; only do this if leakage between splits genuinely doesn't matter for your use case.

=== "`cli`"

    ```bash
    pandora partition --clusters clusters/entries.json --pct-train 0.6 --pct-val 0.2 --pct-test 0.2 --output splits.json
    ```

    ```text
    split sizes: {'train': 4, 'val': 1, 'test': 0} -> splits.json
    ```

    `--no-keep-similar-items` is the CLI form of `keep_similar_items=False`.

### Every kind in every split

Splitting by size alone can leave a kind of item (heteromer interfaces, say) out of a split entirely, and a test set without heteromers tells you nothing about them. Give the split a label per item, from a `chain_labels` or `interface_labels` [layer](annotation.md), and clusters are placed so that every label reaches every split: the rarest labels go first, each using its smallest clusters, then the remaining clusters are placed by size as before. Clusters are still never split, so this is best effort: a label in fewer clusters than there are splits can't reach them all. The partition goes ahead anyway and reports which labels are missing from which splits. Items with no label count as the label `unlabelled`.

Here a labels table, written in pandas, sorts the five chains into two made-up kinds:

```python
import pandas as pd

labels = pd.DataFrame(
    {
        "entry_id": ["104M", "112M", "118L", "138L", "1AYI"],
        "chain_id": ["A", "A", "A", "A", "A"],
        "label": ["kind_a", "kind_b", "kind_a", "kind_b", "kind_a"],
    }
)
labels.to_parquet("kinds.parquet")
```

```yaml
# kinds.yaml
name: kinds
template: chain_labels
tool: pandas
tool_version: "2.2"
parameters: {}
```

Symmetry copies sit in their source chain's cluster as items of their own, so label them too (e.g. copy the source chain's label, using the `copy_of` column of the step folder's `chains` table); an unlabelled copy counts as `unlabelled`.

=== "`cli`"

    ```bash
    pandora annotate import --structures-dir curated/ --table kinds.parquet \
      --description kinds.yaml --output-dir annotations/kinds
    pandora partition --clusters clusters/chains.json \
      --labels annotations/kinds --output splits.json
    jq -c '.label_counts, .missing_labels' partition_provenance.json
    ```

    ```text
    imported: 5 rows, 0 unmatched -> annotations/kinds
    split sizes: {'train': 2, 'val': 2, 'test': 1} -> splits.json
    unlabelled items: 0 of 5
    partition: warning: test has no item labelled kind_b
    {"train":{"kind_a":1,"kind_b":1},"val":{"kind_a":1,"kind_b":1},"test":{"kind_a":1}}
    {"test":["kind_b"]}
    ```

    Without `--labels`, the same clusters split as `{'train': 4, 'val': 1, 'test': 0}`: test gets nothing. With labels, `kind_b` sits in two clusters, so it reaches only two splits; the warning and `missing_labels` say which split it misses. `--labels` can't be combined with `--no-keep-similar-items`. A clusters file of another level than the labels (chain clusters, interface labels) runs with every item unlabelled, and the `unlabelled items` line shows it.

=== "`library`"

    ```python
    import json

    from pandora.export import read_records
    from pandora.schemas.annotation import LayerRow
    from pandora.schemas.similarity import SimilarityCluster
    from pandora.similarity import item_labels, partition_dataset

    clusters = [
        SimilarityCluster(**c) for c in json.load(open("clusters/chains.json"))
    ]
    rows = read_records("annotations/kinds/layer.parquet", LayerRow)
    labels = item_labels(rows, "chain_labels")
    splits, provenance = partition_dataset(clusters, labels=labels)
    print(splits)
    print(provenance.missing_labels)
    ```

    ```text
    {'train': ['104M_A', '112M_A'], 'val': ['118L_A', '138L_A'], 'test': ['1AYI_A']}
    {'test': ['kind_b']}
    ```

    `item_labels()` turns a `chain_labels` or `interface_labels` layer into item IDs (`1ABC_A`, `1ABC_A_B`; a label on B–A labels `1ABC_A_B`) and raises `ValueError` if one item has two different labels. `partition_dataset(labels=...)` raises `ValueError` with `keep_similar_items=False`.
