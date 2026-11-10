from __future__ import annotations

from typing import Literal, NamedTuple

from pydantic import BaseModel, Field

from pandora.schemas.annotation import SimilarityItems
from pandora.schemas.common import RankKey


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


class ClusteredLayer(BaseModel):
    """One similarity layer read by a clustering run.

    Attributes:
        similarity_layer: The layer's `layer.json`, or None when the
            edges came from elsewhere (e.g. a hit TSV).
        items: The layer's item level; chain layers can be rolled up
            to entries or interfaces.
        hit_filter: The filter that turned the layer's rows into edges.
        n_edges: How many edges the layer gave, before any roll-up.
    """

    similarity_layer: str | None = None
    items: SimilarityItems
    hit_filter: HitFilter
    n_edges: int = 0


class ClusteringProvenance(BaseModel):
    """Record of one clustering run: which similarity layers and filters
    produced the edges, and how many edges and clusters resulted.

    Attributes:
        clustered_at: When this clustering run completed, as an ISO
            8601 timestamp.
        layers: The similarity layers whose edges were joined (an edge
            from any layer links two items), with their filters. Empty
            for hand-built edges (`cluster_edges`).
        n_edges: How many edges were read (duplicates included; with
            chain items, the links from symmetry copies to their source
            chain count too, and with a roll-up to interfaces, the
            links between interfaces sharing a chain).
        n_edges_unknown_ids: How many of those named an id outside
            `item_ids` and were ignored.
        n_clusters: How many clusters resulted.
        chain_representatives: The representatives JSON that mapped
            chains to their representative before the roll-up to
            interfaces, if any.
    """

    clustered_at: str
    layers: list[ClusteredLayer] = Field(default_factory=list)
    n_edges: int = 0
    n_edges_unknown_ids: int = 0
    n_clusters: int
    chain_representatives: str | None = None


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


class RepresentativePolicy(BaseModel):
    """Policy for picking one representative item per cluster.

    Attributes:
        items: The cluster items' level: "chain" or "interface" (entry
            is not supported).
        rank: Ranking keys, best first; the lowest item id breaks the
            last tie. Chain keys: "missing_residues" (absolute count),
            "n_hits" (partners in the same cluster), or a number column
            of the `chains` table. Interface keys: "interface_residues"
            (residues in contact, both sides), "resolution",
            "missing_residues" (both chains) and "n_hits". The default
            depends on `items`: chains rank by fewest missing residues,
            then most hits; interfaces by most interface residues, then
            best resolution, then fewest missing residues.
    """

    items: SimilarityItems = "chain"
    rank: list[RankKey] = Field(
        default_factory=lambda data: _DEFAULT_RANK[data["items"]]()
    )


_DEFAULT_RANK = {
    "chain": lambda: [
        RankKey(key="missing_residues", order="min"),
        RankKey(key="n_hits", order="max"),
    ],
    "interface": lambda: [
        RankKey(key="interface_residues", order="max"),
        RankKey(key="resolution", order="min"),
        RankKey(key="missing_residues", order="min"),
    ],
    "entry": list,
}


class RepresentativeCluster(BaseModel):
    """One cluster with its chosen representative.

    Attributes:
        representative: The item that stands for the cluster.
        members: Every item of the cluster, the representative included.
        n_members: The number of members.
    """

    representative: str
    members: list[str]
    n_members: int


class RepresentativesProvenance(BaseModel):
    """Record of one `representatives` run.

    Attributes:
        selected_at: When the run completed, as an ISO 8601 timestamp.
        items: The items' level.
        policy: The ranking policy used.
        clusters: The clusters file read.
        structures_dir: Chains: the step folder the clusters were limited
            to. Interfaces: the structures folder `contacts` traces to,
            whose `chains` table gave the chain values.
        contacts: Interfaces: the contact layer folder the clusters were
            limited to.
        layers: The similarity layers read for `n_hits`, if any.
        n_clusters: Clusters written.
        n_members: Members written, over all clusters.
        n_members_dropped: Members not in the step folder (chains) or
            the contact layer (interfaces).
        n_clusters_dropped: Clusters left empty by that, so not written.
    """

    selected_at: str
    items: SimilarityItems
    policy: RepresentativePolicy
    clusters: str
    structures_dir: str
    contacts: str | None = None
    layers: list[ClusteredLayer] = Field(default_factory=list)
    n_clusters: int
    n_members: int
    n_members_dropped: int = 0
    n_clusters_dropped: int = 0


class PositiveDimer(BaseModel):
    """One dimer kept as a positive by `rebalance`.

    Attributes:
        dimer: The interface item id (`{entry}_{chain}_{chain}`).
        entry_id: The dimer's entry.
        chain_id_1: The first chain (sorted order of the item id).
        chain_id_2: The second chain.
        group: The dimer group it came from (the group's first
            component).
        kind: "homomer" (both chains in one 90% cluster) or "heteromer".
        degree: Its number of neighbours when it was kept.
    """

    dimer: str
    entry_id: str
    chain_id_1: str
    chain_id_2: str
    group: str
    kind: Literal["homomer", "heteromer"]
    degree: int


class DroppedDimer(BaseModel):
    """One dimer `rebalance` did not keep, and why.

    Attributes:
        dimer: The interface item id.
        group: The dimer group it came from.
        reason: "homomer_in_mixed_group" (a homomer in a group with
            heteromers), "similar_to_kept" (a neighbour of a kept
            heteromer) or "not_most_connected" (a homomer-only group's
            other dimers).
        kept: The kept dimer that caused the drop, for
            "similar_to_kept".
    """

    dimer: str
    group: str
    reason: Literal[
        "homomer_in_mixed_group", "similar_to_kept", "not_most_connected"
    ]
    kept: str | None = None


class RebalanceSettings(BaseModel):
    """What one `rebalance` run used, recorded as its summary's policy.

    Attributes:
        rank: The ranking that breaks ties between equally connected
            dimers.
        layers: The similarity layers and hit filters of the dimer
            graph (the same as the Step 5 `cluster` run).
        clusters: The dimer groups file.
        chain_clusters: The 90% chain clusters file (homomer test).
        chain_representatives: The representatives JSON mapping chains
            to their representative, if any.
    """

    rank: RepresentativePolicy
    layers: list[ClusteredLayer] = Field(default_factory=list)
    clusters: str
    chain_clusters: str
    chain_representatives: str | None = None


class NegativePair(BaseModel):
    """One negative chain pair drawn for a positive by `negatives`.

    Attributes:
        chain_1: The positive's chain the negative was drawn for (A), a
            chain item id.
        chain_2: The drawn partner (C), a representative chain item id.
        positive: The positive dimer it was drawn for.
    """

    chain_1: str
    chain_2: str
    positive: str
