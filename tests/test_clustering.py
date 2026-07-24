from __future__ import annotations

from pandora.schemas.similarity import SimilarityMethod, SimilarityRelationship
from pandora.similarity.clustering import cluster_similar_items


def _rel(
    source_id: str, target_id: str, score: float
) -> SimilarityRelationship:
    return SimilarityRelationship(
        source_id=source_id,
        target_id=target_id,
        similarity_type="sequence_similarity",
        score=score,
        method=SimilarityMethod(engine="MMseqs2"),
    )


ITEM_IDS = ["a", "b", "c", "d", "e"]
RELATIONSHIPS = [
    _rel("a", "b", 0.9),
    _rel("b", "c", 0.9),  # transitively joins a-b-c despite no direct a-c edge
    _rel("d", "e", 0.1),  # below threshold, should not merge d and e
]


def test_transitive_merge_and_isolate() -> None:
    clusters = cluster_similar_items(ITEM_IDS, RELATIONSHIPS, threshold=0.5)
    by_size = sorted(clusters, key=lambda c: c.n_components)

    singletons = [c for c in by_size if c.n_components == 1]
    assert {c.components[0] for c in singletons} == {"d", "e"}

    merged = [c for c in by_size if c.n_components == 3][0]
    assert merged.components == ["a", "b", "c"]


def test_every_item_appears_exactly_once() -> None:
    clusters = cluster_similar_items(ITEM_IDS, RELATIONSHIPS, threshold=0.5)
    seen = [item for cluster in clusters for item in cluster.components]
    assert sorted(seen) == sorted(ITEM_IDS)


def test_higher_threshold_splits_the_cluster() -> None:
    clusters = cluster_similar_items(ITEM_IDS, RELATIONSHIPS, threshold=0.95)
    assert all(c.n_components == 1 for c in clusters)


if __name__ == "__main__":
    test_transitive_merge_and_isolate()
    test_every_item_appears_exactly_once()
    test_higher_threshold_splits_the_cluster()
    print("ok")
