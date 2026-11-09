from __future__ import annotations

from typing import Literal, NamedTuple

from pydantic import BaseModel, Field


class SimilarityCluster(BaseModel):
    """One connected-component cluster of similar items.

    Attributes:
        components: Every item id in this cluster.
        n_components: The number of items in this cluster.
    """

    components: list[str]
    n_components: int


class HitFilter(BaseModel):
    """Policy deciding which rows of a similarity layer become edges.

    A pair of items is an edge if any single row passes every
    threshold. `None` thresholds are not applied; the others are
    inclusive (`>=`). A null or NaN value fails. A threshold on a
    column the layer doesn't have is an error.

    Attributes:
        min_score: Minimum value of `score_column`.
        score_column: The column `min_score` tests: "score" (identity
            for Pandora's MMseqs2 layers, alignment TM-score for its
            Foldseek layers), or any other number column, e.g.
            "qtmscore" / "ttmscore" in a Foldseek layer.
        min_identity: Minimum `identity` (fraction of identical aligned
            residues).
        min_coverage: Minimum alignment coverage (see `coverage_of`).
        coverage_of: Which coverage is tested: "both" (the smaller of
            `coverage_1` and `coverage_2`), "side_1", "side_2", or
            "either" (the larger).
        min_interface_coverage: Minimum fraction of each side's
            interface residues inside the alignment (needs the
            alignment ranges and `interface_residues`).
    """

    min_score: float | None = None
    score_column: str = "score"
    min_identity: float | None = None
    min_coverage: float | None = None
    coverage_of: Literal["both", "side_1", "side_2", "either"] = "both"
    min_interface_coverage: float | None = None


class Edge(NamedTuple):
    """One layer row that passed a `HitFilter` (source_id < target_id).

    A plain tuple, not a pydantic model, so streaming millions of edges
    stays cheap.

    Attributes:
        source_id: The lexicographically smaller item id.
        target_id: The lexicographically larger item id.
        score: The row's `score_column` value.
        identity: The row's fraction of identical aligned residues, or
            None if the layer has no `identity`.
        coverage: The row's coverage, as selected by `coverage_of`, or
            None if the layer lacks those columns.
        interface_coverage: Interface-restricted coverage, or None if
            not computed for this pair.
    """

    source_id: str
    target_id: str
    score: float
    identity: float | None
    coverage: float | None
    interface_coverage: float | None


class ClusteringProvenance(BaseModel):
    """Record of one clustering run: which similarity layer and filter
    produced the edges, and how many edges and clusters resulted.

    Attributes:
        clustered_at: When this clustering run completed, as an ISO
            8601 timestamp.
        hit_filter: The filter that turned hits into edges, or None
            for hand-built edges (`cluster_edges`).
        similarity_layer: The `layer.json` of the similarity layer that
            was clustered, or None for hand-built edges.
        n_edges: How many edges were read (duplicates included; with
            chain items, the links from symmetry copies to their source
            chain count too).
        n_edges_unknown_ids: How many of those named an id outside
            `item_ids` and were ignored.
        n_clusters: How many clusters resulted.
    """

    clustered_at: str
    hit_filter: HitFilter | None = None
    similarity_layer: str | None = None
    n_edges: int = 0
    n_edges_unknown_ids: int = 0
    n_clusters: int


class PartitionProvenance(BaseModel):
    """Record of one train/val/test partition run: target fractions,
    mode, and resulting split sizes.

    Attributes:
        partitioned_at: When this partition run completed, as an ISO
            8601 timestamp.
        pct_train: The target fraction of items assigned to train.
        pct_val: The target fraction of items assigned to val.
        pct_test: The target fraction of items assigned to test.
        keep_similar_items: Whether each cluster was kept whole in
            one split (leakage-safe) or divided proportionally.
        split_sizes: The resulting number of items per split.
        label_counts: With labels: split -> label -> number of items,
            listing only labels present in that split. Empty otherwise.
        missing_labels: With labels: split -> sorted labels with no item
            in that split. Splits missing none, and splits with a 0
            fraction, are left out. Empty otherwise.
        labels_layer: The labels layer's `layer.json` (it records the
            layer's name), when the CLI was given one.
    """

    partitioned_at: str
    pct_train: float
    pct_val: float
    pct_test: float
    keep_similar_items: bool
    split_sizes: dict[str, int] = Field(default_factory=dict)
    label_counts: dict[str, dict[str, int]] = Field(default_factory=dict)
    missing_labels: dict[str, list[str]] = Field(default_factory=dict)
    labels_layer: str | None = None
