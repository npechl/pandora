from __future__ import annotations

from typing import Any, Literal, NamedTuple

from pydantic import BaseModel, Field

SimilarityType = Literal[
    "sequence_similarity", "structure_similarity", "custom"
]


class SimilarityMethod(BaseModel):
    """Which similarity engine (and version/parameters) produced a
    similarity network.

    Attributes:
        engine: The similarity engine used (e.g. "MMseqs2",
            "Foldseek").
        version: The engine's version string, if determined.
        parameters: The parameters the engine was run with.
    """

    engine: str
    version: str | None = None
    parameters: dict[str, Any] | None = None


class SimilarityRelationshipProvenance(BaseModel):
    """When (and from which source dataset) one similarity relationship
    was computed.

    Attributes:
        computed_at: When this relationship was computed, as an ISO
            8601 timestamp.
        source_dataset_id: The dataset this relationship was computed
            for, if known.
    """

    computed_at: str | None = None
    source_dataset_id: str | None = None


class SimilarityRelationship(BaseModel):
    """One pairwise similarity score between two items
    (source_id < target_id).

    Attributes:
        source_id: The lexicographically smaller of the two item ids.
        target_id: The lexicographically larger of the two item ids.
        similarity_type: What kind of similarity this measures.
        score: The relationship's primary similarity score.
        coverage: The alignment coverage between the two items, if
            reported.
        interface_coverage: The alignment coverage restricted to each
            item's interface residues, if computed (see
            `compute_structure_similarity(interface_residues=...)`).
        identity: The sequence/structural identity between the two
            items, if reported.
        method: Which engine/parameters produced this relationship.
        provenance: When and for which dataset this relationship was
            computed.
    """

    source_id: str
    target_id: str
    # source_id < target_id lexicographically, one record per unordered pair.
    similarity_type: SimilarityType
    score: float
    coverage: float | None = None
    interface_coverage: float | None = None
    identity: float | None = None
    method: SimilarityMethod
    provenance: SimilarityRelationshipProvenance = Field(
        default_factory=SimilarityRelationshipProvenance
    )


class SimilarityCluster(BaseModel):
    """One connected-component cluster of similar items.

    Attributes:
        components: Every item id in this cluster.
        n_components: The number of items in this cluster.
    """

    components: list[str]
    n_components: int


SimilarityEngine = Literal["MMseqs2", "Foldseek"]


class SimilaritySearch(BaseModel):
    """One all-vs-all search: where its hit file is and how it was run.

    The hit file is the tool's raw TSV with Pandora's fixed
    `--format-output` columns; it stays on disk and is streamed by
    `iter_edges()`, never loaded whole.

    Attributes:
        engine: Which tool produced the hits.
        version: The tool's version string, if determined.
        hits_path: Path to the hit TSV.
        columns: The TSV's columns, in order.
        parameters: The search settings (binary, sensitivity,
            max_seqs, extra options, ...), as passed to the search
            function, so the search can be re-run.
        origin: "computed" if Pandora ran the search, "precomputed" if
            an existing file was wrapped with `load_similarity_search()`.
        searched_at: When the search ran, as an ISO 8601 timestamp.
    """

    engine: SimilarityEngine
    version: str | None = None
    hits_path: str
    columns: list[str]
    parameters: dict[str, Any] = Field(default_factory=dict)
    origin: Literal["computed", "precomputed"] = "computed"
    searched_at: str | None = None


class HitFilter(BaseModel):
    """Policy deciding which hit rows become similarity edges.

    A pair of items is an edge if any single hit row passes every
    threshold. `None` thresholds are not applied; the others are
    inclusive (`>=`).

    Attributes:
        min_score: Minimum score: TM-score for Foldseek (see
            `tm_normalisation`), identity for MMseqs2.
        min_identity: Minimum fraction of identical aligned residues.
        min_coverage: Minimum alignment coverage (see `coverage_of`).
        coverage_of: Which coverage is tested: "both" (the smaller of
            query and target), "query", "target", or "either" (the
            larger).
        min_interface_coverage: Minimum fraction of each side's
            interface residues inside the alignment (Foldseek only;
            needs `interface_residues`).
        tm_normalisation: Foldseek TM-score used as the score:
            "alignment" (alntmscore), "query" (qtmscore), "target"
            (ttmscore), or "max" (the larger of the two).
    """

    min_score: float | None = None
    min_identity: float | None = None
    min_coverage: float | None = None
    coverage_of: Literal["both", "query", "target", "either"] = "both"
    min_interface_coverage: float | None = None
    tm_normalisation: Literal["alignment", "query", "target", "max"] = (
        "alignment"
    )


class Edge(NamedTuple):
    """One hit row that passed a `HitFilter` (source_id < target_id).

    A plain tuple, not a pydantic model, so streaming millions of edges
    stays cheap.

    Attributes:
        source_id: The lexicographically smaller item id.
        target_id: The lexicographically larger item id.
        score: The row's score, as selected by the filter.
        identity: The row's fraction of identical aligned residues.
        coverage: The row's coverage, as selected by `coverage_of`.
        interface_coverage: Interface-restricted coverage, or None if
            not computed for this pair.
    """

    source_id: str
    target_id: str
    score: float
    identity: float
    coverage: float
    interface_coverage: float | None


class ClusteringProvenance(BaseModel):
    """Record of one clustering run: threshold, relationship/cluster
    counts, and the similarity method used.

    Attributes:
        clustered_at: When this clustering run completed, as an ISO
            8601 timestamp.
        threshold: The similarity score threshold used to connect
            items.
        n_relationships: How many relationships were considered.
        n_clusters: How many clusters resulted.
        similarity_method: Which engine/parameters produced the
            clustered relationships.
    """

    clustered_at: str
    threshold: float
    n_relationships: int
    n_clusters: int
    similarity_method: SimilarityMethod | None = None


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
    """

    partitioned_at: str
    pct_train: float
    pct_val: float
    pct_test: float
    keep_similar_items: bool
    split_sizes: dict[str, int] = Field(default_factory=dict)
