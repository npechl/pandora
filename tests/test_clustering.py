from __future__ import annotations

from pandora.schemas.similarity import (
    ClusteringProvenance,
    HitFilter,
    SimilaritySearch,
)
from pandora.similarity.clustering import (
    cluster_edges,
    cluster_similar_items,
    pair_cluster_keys,
)
from pandora.similarity.hits import HIT_COLUMNS

ITEM_IDS = ["a", "b", "c", "d", "e"]
EDGES = [("a", "b"), ("b", "c")]  # a-b-c joined transitively, no a-c edge


def test_transitive_merge_and_isolate() -> None:
    clusters, provenance = cluster_edges(ITEM_IDS, EDGES)

    assert [c.components for c in clusters] == [["a", "b", "c"], ["d"], ["e"]]
    assert provenance.n_edges == 2
    assert provenance.n_edges_unknown_ids == 0
    assert provenance.n_clusters == 3
    assert provenance.hit_filter is None
    assert provenance.search is None


def test_every_item_appears_exactly_once() -> None:
    clusters, _ = cluster_edges(ITEM_IDS, EDGES)

    seen = [item for cluster in clusters for item in cluster.components]
    assert sorted(seen) == sorted(ITEM_IDS)


def test_unknown_ids_are_ignored_and_counted() -> None:
    clusters, provenance = cluster_edges(["a", "b"], [("a", "b"), ("a", "z")])

    assert [c.components for c in clusters] == [["a", "b"]]
    assert provenance.n_edges == 2
    assert provenance.n_edges_unknown_ids == 1


def test_duplicate_edges_do_not_change_clusters() -> None:
    once, _ = cluster_edges(ITEM_IDS, EDGES)
    repeated, _ = cluster_edges(ITEM_IDS, EDGES + EDGES + [("b", "a")])

    assert repeated == once


def test_cluster_similar_items_filters_hits_and_records_search(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text(
        "a\tb\t0.95\t100\t0.9\t0.9\n"
        "b\tc\t0.95\t100\t0.9\t0.9\n"
        "d\te\t0.10\t100\t0.9\t0.9\n"
    )
    search = SimilaritySearch(
        engine="MMseqs2",
        hits_path=str(path),
        columns=list(HIT_COLUMNS["MMseqs2"]),
    )
    hit_filter = HitFilter(min_score=0.5)

    clusters, provenance = cluster_similar_items(ITEM_IDS, search, hit_filter)

    assert [c.components for c in clusters] == [["a", "b", "c"], ["d"], ["e"]]
    assert provenance.hit_filter == hit_filter
    assert provenance.search == search
    assert provenance.n_edges == 2


def test_old_clustering_provenance_still_loads() -> None:
    old = {
        "clustered_at": "2026-01-01T00:00:00+00:00",
        "threshold": 0.9,
        "n_relationships": 3,
        "n_clusters": 2,
    }

    provenance = ClusteringProvenance.model_validate(old)

    assert provenance.search is None
    assert provenance.hit_filter is None
    assert provenance.n_clusters == 2


def test_pair_cluster_keys_uses_each_clusters_min_component_id() -> None:
    clusters, _ = cluster_edges(ITEM_IDS, EDGES)

    keys = pair_cluster_keys([("a", "d"), ("c", "e")], clusters)

    # a, c share the {a,b,c} cluster (keyed "a"); d, e are singletons.
    assert keys == {("a", "d"): ("a", "d"), ("c", "e"): ("a", "e")}
