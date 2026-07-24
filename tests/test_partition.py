from __future__ import annotations

from pandora.schemas.similarity import SimilarityCluster
from pandora.similarity.partition import partition_dataset

CLUSTERS = [
    SimilarityCluster(components=[f"a{i}" for i in range(10)], n_components=10),
    SimilarityCluster(components=[f"b{i}" for i in range(8)], n_components=8),
    SimilarityCluster(components=[f"c{i}" for i in range(6)], n_components=6),
    SimilarityCluster(components=[f"d{i}" for i in range(4)], n_components=4),
    SimilarityCluster(components=[f"e{i}" for i in range(2)], n_components=2),
]


def test_keep_similar_items_hits_target_sizes() -> None:
    result = partition_dataset(CLUSTERS, 0.6, 0.2, 0.2, keep_similar_items=True)
    assert [len(result[k]) for k in ("train", "val", "test")] == [18, 6, 6]


def test_keep_similar_items_never_splits_a_cluster() -> None:
    result = partition_dataset(CLUSTERS, 0.6, 0.2, 0.2, keep_similar_items=True)
    for cluster in CLUSTERS:
        hits = [
            k for k, v in result.items() if set(cluster.components) & set(v)
        ]
        assert len(hits) == 1


def test_no_keep_similar_items_splits_proportionally_and_flat() -> None:
    result = partition_dataset(
        CLUSTERS, 0.6, 0.2, 0.2, keep_similar_items=False
    )
    assert all(isinstance(item, str) for item in result["train"])
    assert sum(len(v) for v in result.values()) == sum(
        c.n_components for c in CLUSTERS
    )


def test_bad_fractions_raise() -> None:
    try:
        partition_dataset(CLUSTERS, 0.6, 0.2, 0.3)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for fractions != 1.0")


if __name__ == "__main__":
    test_keep_similar_items_hits_target_sizes()
    test_keep_similar_items_never_splits_a_cluster()
    test_no_keep_similar_items_splits_proportionally_and_flat()
    test_bad_fractions_raise()
    print("ok")
