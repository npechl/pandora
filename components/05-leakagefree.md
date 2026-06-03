# Component 05 — Similarity Analysis & Leakage-Safe Dataset Splitting Layer

## Purpose

The Similarity Analysis & Leakage-Safe Dataset Splitting Layer takes a
curated `Dataset` and turns it into a reproducible, similarity-aware,
leakage-safe dataset partition.

This component is responsible for:

* computing structural similarity relationships,
* building similarity networks,
* clustering related structures,
* creating leakage-safe train/validation/test partitions,
* recording similarity provenance,
* and exposing relationship artifacts for auditability.

This component does **not** compute new biological annotations or curate
the dataset itself. It operates on the output of Component 4.

---

# 1. Architectural Role

```text
PandoraDataset (Dataset | ChainDataset | InterfaceDataset | ResidueDataset)
  → SimilarityRelationship
  → SimilarityNetwork
  → SimilarityClusters
  → LeakageSafeDataset
```

---

# 2. Core Design Principles

## Similarity-aware splitting

Pandora prevents leakage by grouping related structures into clusters before
partitioning. Structures within the same cluster are always assigned to the
same split.

---

## Policy-driven splitting

Users explicitly define:

* similarity engines and thresholds,
* clustering strategy,
* split ratios,
* and leakage constraints.

---

## External similarity engines

Pandora wraps external similarity engines rather than reimplementing them.

Supported in V1:

* MMseqs2 — sequence similarity,
* Foldseek — structure similarity.

---

## Relationship artifacts are first-class outputs

Component 5 exposes:

* pairwise similarity records (`SimilarityRelationship`),
* similarity networks (`SimilarityNetwork`),
* clusters (`SimilarityClusters`),
* leakage-safe partitions (`LeakageSafeDataset`).

All artifacts are inspectable and reproducible.

---

# 3. Input Schemas

## 3.0 PandoraDataset — Union Type

Component 05 accepts any dataset produced by Component 04, regardless of
granularity level.

```yaml
PandoraDataset:
  # Dataset | ChainDataset | InterfaceDataset | ResidueDataset
  # Discriminated by the granularity field:
  #   Dataset:           granularity = "structure"
  #   ChainDataset:      granularity = "chain"
  #   InterfaceDataset:  granularity = "interface"
  #   ResidueDataset:    granularity = "residue"
  # All four types share: dataset_id, dataset_version,
  # applied_policy, provenance, and granularity.
```

---

## 3.1 Leakage analysis input

```yaml
LeakageAnalysisInput:
  dataset: PandoraDataset

  leakage_policy: LeakagePolicy
  # See Section 8 for the full policy schema and strategy definitions.
```

---

## 3.2 Batch leakage analysis input

```yaml
LeakageAnalysisBatchInput:
  datasets:
    - PandoraDataset
  # Each element is an independent splitting job producing one
  # LeakageSafeDataset.

  leakage_policy: LeakagePolicy

  mode: string
  # sequential | parallel

  parallel_options:
    max_workers: int | null
    # Number of concurrent jobs in parallel mode.
    # null uses the system default.
    # Ignored in sequential mode.

    fail_fast: bool
    # If true, abort remaining jobs on the first failure.
    # If false (default), isolate failures and continue.
```

---

# 4. Output Schemas

## 4.1 SimilarityRelationship

A single pairwise similarity record between two structures.

```yaml
SimilarityRelationship:
  source_id: string
  # Item identifier of the first item. Format depends on granularity:
  #   structure:  entry_id (e.g. "1abc")
  #   chain:      "{entry_id}_{chain_id}" (e.g. "1abc_A")
  #   interface:  "{entry_id}_{chain_id_1}_{chain_id_2}" (e.g. "1abc_A_B")
  # source_id and target_id are always ordered lexicographically
  # (source_id < target_id) to avoid duplicate pairs.

  target_id: string
  # Item identifier of the second item. Same format as source_id.

  similarity_type: string
  # sequence_similarity | structure_similarity | custom

  score: float
  # Primary similarity score in range [0.0, 1.0].
  # For sequence similarity: sequence identity (MMseqs2 default output).
  # For structure similarity: TM-score or equivalent (Foldseek default output).

  coverage: float | null
  # Alignment coverage in range [0.0, 1.0]. null if not reported by engine.

  identity: float | null
  # Sequence identity in range [0.0, 1.0].
  # Populated for sequence_similarity; null for structure_similarity.

  method:
    engine: string
    # Supported values: MMseqs2 | Foldseek | custom
    # Custom engines must be registered and must produce normalised outputs.

    version: string | null
    # Engine version string. null if not available from engine output.

    parameters: object | null
    # Key engine parameters used (e.g. sensitivity, e-value cutoff).
    # null if not recorded.

  provenance:
    computed_at: string | null
    # ISO 8601 timestamp. null if not recorded.
    source_dataset_id: string | null
    # dataset_id of the Dataset this pair was drawn from.
```

---

## 4.2 SimilarityNetwork

A graph where nodes are structures and edges are similarity relationships
at or above the policy threshold.

```yaml
SimilarityNetwork:
  network_id: string
  dataset_id: string

  relationships: list[SimilarityRelationship]
  # All pairwise relationships computed, including those below the threshold.
  # The threshold is applied when constructing edges (see below).

  nodes: list[string]
  # All item identifiers in the source dataset, including isolates (no edges).
  # Format matches the SimilarityRelationship identifier format for the
  # dataset's granularity.

  edges: list[SimilarityEdge]
  # Relationships with score >= policy threshold that form graph edges.

  graph_statistics:
    node_count: int
    edge_count: int
    connected_components: int | null
    # Number of connected components in the graph.
    # null if clustering has not yet been applied.
    largest_component_size: int | null

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  provenance:
    built_at: string | null
    # ISO 8601 timestamp.

SimilarityEdge:
  source_id: string
  target_id: string
  score: float
  similarity_type: string
```

---

## 4.3 SimilarityClusters

```yaml
SimilarityClusters:
  clustering_id: string
  dataset_id: string

  clusters: list[Cluster]

  clustering_summary:
    total_clusters: int
    singleton_clusters: int
    # Clusters with exactly one member (no similar neighbours).
    largest_cluster_size: int
    mean_cluster_size: float

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  provenance:
    clustered_at: string | null
    # ISO 8601 timestamp.

Cluster:
  cluster_id: string

  members: list[string]
  # Item identifiers of all items in this cluster.
  # Format matches the SimilarityRelationship identifier format for the
  # dataset's granularity.

  cluster_size: int

  representative_id: string | null
  # The entry_id selected as the cluster representative.
  # Selection criterion: the member with the highest mean similarity score
  # to all other cluster members (most central node).
  # null for singleton clusters (no pairwise scores available).
```

---

## 4.4 LeakageSafeDataset

```yaml
LeakageSafeDataset:
  dataset_id: string
  dataset_name: string
  dataset_version: string

  granularity: string
  # Granularity of items in the partition lists.
  # Values: structure | chain | interface | residue
  # Matches source_dataset.granularity.

  source_dataset: PandoraDataset
  # The curated dataset from Component 4. Retrieve items by identifier:
  #   Dataset:           source_dataset.structures filtered by entry_id
  #   ChainDataset:      source_dataset.chains filtered by chain_id
  #   InterfaceDataset:  source_dataset.interfaces filtered by interface_id
  #   ResidueDataset:    source_dataset.residues filtered by residue_id

  similarity_network: SimilarityNetwork
  similarity_clusters: SimilarityClusters

  partitions:
    train: list[string]
    validation: list[string]
    test: list[string]
    # All lists contain item identifier strings.
    # Format depends on granularity:
    #   structure:  entry_id (e.g. "1abc")
    #   chain:      "{entry_id}_{chain_id}" (e.g. "1abc_A")
    #   interface:  "{entry_id}_{chain_id_1}_{chain_id_2}" (e.g. "1abc_A_B")
    # To retrieve items: look up each identifier in source_dataset.

  partition_summary:
    train_count: int
    validation_count: int
    test_count: int
    train_fraction_achieved: float
    validation_fraction_achieved: float
    test_fraction_achieved: float
    # Achieved fractions may differ slightly from requested fractions due
    # to discrete cluster sizes. See Section 5.4 for the assignment algorithm.

  leakage_summary:
    max_cross_split_similarity: float | null
    # Highest similarity score observed between any pair of structures
    # assigned to different splits.
    # Measured on the same [0.0, 1.0] scale as SimilarityRelationship.score.
    # null if no cross-split edges exist in the similarity network.

    leakage_detected: bool
    # True if max_cross_split_similarity > leakage_rules.max_allowed_cross_split_similarity,
    # or if a cluster was split across splits when enforce_cluster_isolation: true.

    leakage_diagnostics: list[Diagnostic]
    # Detailed per-pair or per-cluster leakage records.

  applied_policy:
    policy_id: string
    policy_name: string
    policy_version: string

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]

  provenance:
    split_at: string | null
    # ISO 8601 timestamp.
    similarity_engines: list[string]
    # Engine names used (e.g. ["MMseqs2", "Foldseek"]).
```

---

## 4.5 Batch output

```yaml
LeakageAnalysisBatchResult:
  mode: string
  # sequential | parallel

  summary:
    total: int
    success: int
    warning: int
    failed: int

  results:
    - dataset_id: string
      status: string
      # success | warning | failed

      leakage_safe_dataset: LeakageSafeDataset | null
      # null when status == "failed".

      diagnostics:
        warnings: list[Diagnostic]
        errors: list[Diagnostic]
```

---

# 5. Public Functions

## 5.1 `compute_similarity_relationships()`

### Responsibility

Compute pairwise similarity relationships for all structures in a dataset
by delegating to external similarity engines.

### Internal Workflow

```text
if sequence_similarity.enabled:
    run_sequence_similarity_engine(dataset, policy.similarity_rules.sequence_similarity)
    → list[SimilarityRelationship] with similarity_type="sequence_similarity"

if structure_similarity.enabled:
    run_structure_similarity_engine(dataset, policy.similarity_rules.structure_similarity)
    → list[SimilarityRelationship] with similarity_type="structure_similarity"

if both are enabled:
    merge both relationship lists (no deduplication — each type produces
    independent records for the same pairs)

Normalise all engine outputs to SimilarityRelationship schema:
    - map engine-specific score fields to score, coverage, identity
    - record engine name and version in method
    - enforce source_id < target_id lexicographic ordering

return merged list of SimilarityRelationship records
```

### Input Schema

```yaml
compute_similarity_relationships:
  dataset: PandoraDataset
  leakage_policy: LeakagePolicy
```

### Output Schema

```yaml
compute_similarity_relationships_result:
  relationships: list[SimilarityRelationship]

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

---

## 5.2 `build_similarity_network()`

### Responsibility

Convert pairwise similarity relationships into a `SimilarityNetwork` by
applying the similarity threshold to determine graph edges.

### Edge creation rule

A `SimilarityRelationship` becomes a `SimilarityEdge` if and only if
`score >= policy.similarity_rules.{type}.threshold` for its
`similarity_type`. Relationships below the threshold are retained in
`network.relationships` for auditability but do not form edges.

### Input Schema

```yaml
build_similarity_network:
  relationships: list[SimilarityRelationship]
  dataset_id: string
  leakage_policy: LeakagePolicy
```

### Output Schema

```yaml
build_similarity_network_result:
  network: SimilarityNetwork
```

---

## 5.3 `cluster_similar_items()`

### Responsibility

Group structures into `SimilarityClusters` using the network's connected
structure and the clustering strategy from the policy.

### Input Schema

```yaml
cluster_similar_items:
  network: SimilarityNetwork

  clustering_rules:
    enabled: bool
    strategy: string
    # connected_components | single_linkage | custom
    # See Section 8 for strategy definitions.
    min_cluster_size: int | null
    custom_clustering_fn: callable | null
    # Required when strategy == "custom". Must accept SimilarityNetwork
    # and return list[Cluster].
```

### Output Schema

```yaml
cluster_similar_items_result:
  clusters: SimilarityClusters
```

### Notes

When `clustering_rules.enabled: false`, each structure is placed in its
own singleton cluster before partitioning. This disables similarity-based
grouping but still produces a valid `SimilarityClusters` object.

---

## 5.4 `partition_dataset()`

### Responsibility

Assign clusters to train/validation/test splits while enforcing leakage
constraints.

### Partition assignment algorithm

```text
Precondition: train_fraction + validation_fraction + test_fraction == 1.0
              (within floating-point tolerance 0.001). Raise configuration
              error if violated.

1. Sort clusters by size descending (largest clusters first).
   This ensures large clusters are placed early when splits have more
   remaining capacity.

2. Initialise split buckets: train=[], validation=[], test=[]
   with targets: n_train = round(total * train_fraction), etc.

3. For each cluster (largest to smallest):
   a. Compute current fill fractions for each split.
   b. Assign the cluster to the split that is most under-represented
      relative to its target fraction (greedy bin-packing).
   c. If keep_similar_items_together: true (default), all members of the
      cluster go to the same split. This is the primary leakage-prevention
      mechanism.

4. Record achieved fractions in partition_summary.

Note: Achieved fractions may deviate from targets when clusters are large
relative to dataset size. This is expected and reported as
UNBALANCED_PARTITION_WARNING if deviation > 5%.
```

### Input Schema

```yaml
partition_dataset:
  dataset: PandoraDataset
  clusters: SimilarityClusters
  leakage_policy: LeakagePolicy
```

### Output Schema

```yaml
partition_dataset_result:
  leakage_safe_dataset: LeakageSafeDataset
```

---

## 5.5 `analyze_leakage()`

### Responsibility

Assess whether the partition violates similarity-based leakage constraints
and populate `leakage_summary`.

### V1 Leakage Analysis Rules

```yaml
error_rules:
  FORBIDDEN_CROSS_SPLIT_SIMILARITY:
    condition: "A SimilarityEdge connects structures in different splits AND
                leakage_rules.forbid_cross_split_similarity: true."
    result_status: failed

  CLUSTER_ISOLATION_VIOLATED:
    condition: "Members of the same cluster are assigned to different splits AND
                leakage_rules.enforce_cluster_isolation: true."
    result_status: failed

warning_rules:
  CROSS_SPLIT_SIMILARITY_ABOVE_THRESHOLD:
    condition: "max_cross_split_similarity > leakage_rules.max_allowed_cross_split_similarity."
    result_status: warning

  UNBALANCED_PARTITION:
    condition: "Achieved split fraction deviates by more than 5% from the
                requested fraction for any split."
    result_status: warning

  LARGE_CLUSTER:
    condition: "A single cluster contains more than 10% of total dataset structures."
    result_status: warning

  SINGLETON_DOMINATED:
    condition: "More than 80% of clusters are singletons (no similar neighbours
                found). May indicate threshold is too strict or dataset is highly
                diverse."
    result_status: warning
```

### Status determination

```yaml
status_rules:
  failed:  Any error_rule fires.
  warning: No error_rules fire but one or more warning_rules fire.
  clean:   No rules fire.
```

### Remediation guidance

```yaml
remediation:
  FORBIDDEN_CROSS_SPLIT_SIMILARITY:
    - Lower the similarity threshold to reduce the number of edges.
    - Use a stricter clustering strategy to merge more structures.
    - Accept a higher leakage tolerance via max_allowed_cross_split_similarity.

  CLUSTER_ISOLATION_VIOLATED:
    - Check that keep_similar_items_together: true is set.
    - Investigate oversized clusters using similarity_clusters output.

  UNBALANCED_PARTITION:
    - Expected when clusters are large relative to dataset size.
    - Accept the imbalance or pre-filter the dataset to reduce cluster sizes.
```

### Input Schema

```yaml
analyze_leakage:
  leakage_safe_dataset: LeakageSafeDataset
  leakage_policy: LeakagePolicy
```

### Output Schema

```yaml
analyze_leakage_result:
  leakage_detected: bool

  diagnostics:
    warnings: list[Diagnostic]
    errors: list[Diagnostic]
```

---

## 5.6 `build_leakage_safe_dataset()`

### Responsibility

Run the complete leakage-safe splitting workflow for a single dataset.
This is the main orchestrator for Component 05.

### Internal Workflow

```text
1. Validate partition fractions:
   train_fraction + validation_fraction + test_fraction must == 1.0
   (within 0.001). Raise configuration error if violated.

2. Compute similarity relationships:
   compute_similarity_relationships(dataset, policy)
   → list[SimilarityRelationship]

3. Build similarity network:
   build_similarity_network(relationships, dataset_id, policy)
   → SimilarityNetwork

4. Cluster similar items:
   cluster_similar_items(network, policy.clustering_rules)
   → SimilarityClusters

5. Partition dataset:
   partition_dataset(dataset, clusters, policy)
   → LeakageSafeDataset (partitions populated, leakage_summary empty)

6. Analyze leakage:
   analyze_leakage(leakage_safe_dataset, policy)
   → Populate leakage_summary and diagnostics

7. If leakage_detected AND forbid_cross_split_similarity: true:
   return LeakageSafeDataset with status "failed" and leakage diagnostics.
   The partitions are still populated for inspection.

8. Return LeakageSafeDataset with status "success" or "warning".
```

### Input Schema

```yaml
build_leakage_safe_dataset:
  input: LeakageAnalysisInput
```

### Output Schema

```yaml
build_leakage_safe_dataset_result:
  leakage_safe_dataset: LeakageSafeDataset
```

---

## 5.7 `build_leakage_safe_dataset_many()`

### Responsibility

Run leakage-safe splitting for multiple independent datasets.

### Input Schema

```yaml
build_leakage_safe_dataset_many:
  input: LeakageAnalysisBatchInput
```

### Output Schema

```yaml
build_leakage_safe_dataset_many_result:
  result: LeakageAnalysisBatchResult
```

---

# 6. Internal Helper Functions

## 6.1 `run_sequence_similarity_engine()`

### Responsibility

Invoke MMseqs2 (or a registered sequence similarity engine) and normalise
its output to `SimilarityRelationship` records.

### Input

```yaml
run_sequence_similarity_engine:
  dataset: PandoraDataset
  # Sequence extraction per granularity:
  #   Dataset:      representative chain sequences from AnnotatedStructureWithPlugins
  #   ChainDataset: ChainRecord.sequence fed directly to MMseqs2
  # InterfaceDataset/ResidueDataset: use structure_similarity instead.
  sequence_similarity_rules:
    engine: string
    threshold: float | null
    coverage_threshold: float | null
```

### Output

```yaml
run_sequence_similarity_engine_result:
  relationships: list[SimilarityRelationship]
  # All pairs with score >= 0.0. Threshold filtering occurs in
  # build_similarity_network(), not here.
  engine_version: string | null
  diagnostics: list[Diagnostic]
```

---

## 6.2 `run_structure_similarity_engine()`

### Responsibility

Invoke Foldseek (or a registered structure similarity engine) and normalise
its output to `SimilarityRelationship` records.

### Input

```yaml
run_structure_similarity_engine:
  dataset: PandoraDataset
  # Coordinate extraction per granularity:
  #   Dataset:           all-atom coordinates from AnnotatedStructureWithPlugins
  #   ChainDataset:      per-chain backbone coordinates from ChainRecord.residues
  #   InterfaceDataset:  interface partner coordinates from InterfaceRecord
  # ResidueDataset: structure_similarity not supported in V1.
  structure_similarity_rules:
    engine: string
    threshold: float | null
    coverage_threshold: float | null
```

### Output

```yaml
run_structure_similarity_engine_result:
  relationships: list[SimilarityRelationship]
  engine_version: string | null
  diagnostics: list[Diagnostic]
```

---

## 6.3 `build_clusters_from_network()`

### Responsibility

Apply the clustering strategy to the `SimilarityNetwork` edges and produce
a list of `Cluster` records.

### Input

```yaml
build_clusters_from_network:
  network: SimilarityNetwork
  strategy: string
  min_cluster_size: int | null
  custom_clustering_fn: callable | null
```

### Output

```yaml
build_clusters_from_network_result:
  clusters: list[Cluster]
  diagnostics: list[Diagnostic]
```

---

## 6.4 `assign_partitions_by_cluster()`

### Responsibility

Execute the greedy bin-packing assignment algorithm to assign clusters to
train/validation/test splits.

### Input

```yaml
assign_partitions_by_cluster:
  clusters: list[Cluster]
  partition_rules:
    train_fraction: float
    validation_fraction: float
    test_fraction: float
    keep_similar_items_together: bool
    stratify_by_cluster: bool
```

### Output

```yaml
assign_partitions_by_cluster_result:
  train: list[string]
  validation: list[string]
  test: list[string]
  partition_summary:
    train_count: int
    validation_count: int
    test_count: int
    train_fraction_achieved: float
    validation_fraction_achieved: float
    test_fraction_achieved: float
```

---

## 6.5 `summarize_leakage_risks()`

### Responsibility

Scan the similarity network for cross-split edges and compute
`leakage_summary` fields.

### Input

```yaml
summarize_leakage_risks:
  network: SimilarityNetwork
  partitions:
    train: list[string]
    validation: list[string]
    test: list[string]
  leakage_rules:
    forbid_cross_split_similarity: bool
    max_allowed_cross_split_similarity: float | null
    enforce_cluster_isolation: bool
```

### Output

```yaml
summarize_leakage_risks_result:
  max_cross_split_similarity: float | null
  leakage_detected: bool
  leakage_diagnostics: list[Diagnostic]
```

---

## 6.6 `record_split_provenance()`

### Responsibility

Construct the provenance record for the `LeakageSafeDataset`.

### Input

```yaml
record_split_provenance:
  dataset: PandoraDataset
  policy: LeakagePolicy
  engine_versions: list[string]
```

### Output

```yaml
record_split_provenance_result:
  provenance:
    split_at: string
    # ISO 8601 timestamp.
    similarity_engines: list[string]
```

---

# 7. V1 Similarity Support

## Supported in V1

```yaml
supported_similarity_types:
  sequence_similarity:
    engine: MMseqs2
    score_metric: sequence_identity
    score_range: [0.0, 1.0]

  structure_similarity:
    engine: Foldseek
    score_metric: TM-score
    score_range: [0.0, 1.0]
```

## Recommended similarity type per granularity

| Granularity | Recommended similarity | Engine | Notes |
|-------------|----------------------|--------|-------|
| structure | sequence_similarity | MMseqs2 | Fast; representative chain sequence per entry |
| structure | structure_similarity | Foldseek | Use when fold-level leakage matters |
| chain | sequence_similarity | MMseqs2 | Primary use case; `ChainRecord.sequence` fed directly |
| chain | structure_similarity | Foldseek | Optional; use for fold-level cluster separation |
| interface | structure_similarity | Foldseek | Interface geometry comparison |
| residue | — | — | Not recommended in V1; inherit splits from chain level |

## Future support

```yaml
future_similarity_types:
  - fold_similarity
  - assembly_similarity
  - interface_similarity
  - custom_similarity_plugins
```

---

# 8. Policy Schema

## 8.1 LeakagePolicy

```yaml
LeakagePolicy:
  policy_id: string
  policy_name: string
  policy_version: string
  description: string

  similarity_rules:
    sequence_similarity:
      enabled: bool

      engine: string
      # Supported: MMseqs2
      # Custom engines must be registered and produce normalised outputs.

      threshold: float | null
      # Minimum score to form a graph edge. Range [0.0, 1.0].
      # null = include all pairs as edges regardless of score.

      coverage_threshold: float | null
      # Minimum alignment coverage. Range [0.0, 1.0].
      # null = no coverage filter.

    structure_similarity:
      enabled: bool

      engine: string
      # Supported: Foldseek

      threshold: float | null
      # Minimum TM-score to form a graph edge. Range [0.0, 1.0].

      coverage_threshold: float | null

  clustering_rules:
    enabled: bool
    # If false, each structure is its own singleton cluster.
    # Partitioning still proceeds but no similarity grouping is applied.

    strategy: string
    # connected_components — All structures transitively connected through
    #   edges above the similarity threshold form one cluster. Uses a standard
    #   graph connected-components algorithm (e.g. union-find). No distance
    #   metric is needed beyond the edge threshold.
    #
    # single_linkage       — Hierarchical agglomerative clustering using
    #   single linkage (merge two groups if any pair across them has score >=
    #   threshold). The dendrogram is cut at the similarity threshold.
    #
    # custom               — Use a user-supplied clustering function provided
    #   as custom_clustering_fn. Must accept SimilarityNetwork and return
    #   list[Cluster].

    min_cluster_size: int | null
    # Clusters smaller than this value are merged into a residual "small
    # clusters" group before partitioning.
    # null = no minimum (singleton clusters are allowed).

    custom_clustering_fn: callable | null
    # Required when strategy == "custom". Ignored otherwise.

  partition_rules:
    train_fraction: float
    validation_fraction: float
    test_fraction: float
    # Must sum to 1.0 (within tolerance 0.001). Configuration error if not.

    keep_similar_items_together: bool
    # If true (default and recommended), all members of a cluster are
    # assigned to the same split. This is the primary leakage-prevention
    # mechanism. Setting to false allows intra-cluster splitting and is
    # not recommended for leakage-safe datasets.

    stratify_by_cluster: bool
    # If true, the greedy bin-packing assignment attempts to distribute
    # clusters of varying sizes across all splits, so each split receives
    # a mix of small and large clusters rather than all large clusters
    # going to train.
    # If false, clusters are assigned purely to minimise fraction deviation.
    # This is independent of keep_similar_items_together.

  leakage_rules:
    forbid_cross_split_similarity: bool
    # If true, any cross-split similarity edge causes status "failed".

    max_allowed_cross_split_similarity: float | null
    # Maximum tolerated cross-split similarity score. Range [0.0, 1.0].
    # Measured on the same scale as SimilarityRelationship.score.
    # null = no tolerance limit (only forbid_cross_split_similarity applies).

    enforce_cluster_isolation: bool
    # If true, any cluster whose members are split across different
    # partitions causes status "failed".

  provenance_rules:
    record_similarity_method: bool
    record_engine_versions: bool
    record_thresholds: bool
    record_partition_history: bool
```

---

# 9. Non-Responsibilities

Component 05 is not responsible for:
  - dataset_construction
  - metadata_integration
  - canonicalization
  - graph_generation
  - embeddings
  - model_training
  - benchmark_reporting
  - task_specific_annotation_creation
  - biological_annotation_computation

---

# 10. Component Definition

The Similarity Analysis & Leakage-Safe Dataset Splitting Layer computes structural similarity relationships, clusters related structures, and generates reproducible leakage-safe partitions using explicit similarity and leakage policies.
