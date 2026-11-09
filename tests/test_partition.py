from __future__ import annotations

import pytest

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
    result, provenance = partition_dataset(
        CLUSTERS, 0.6, 0.2, 0.2, keep_similar_items=True
    )
    assert [len(result[k]) for k in ("train", "val", "test")] == [18, 6, 6]
    assert provenance.split_sizes == {"train": 18, "val": 6, "test": 6}
    assert provenance.keep_similar_items is True


def test_keep_similar_items_never_splits_a_cluster() -> None:
    result, _ = partition_dataset(
        CLUSTERS, 0.6, 0.2, 0.2, keep_similar_items=True
    )
    for cluster in CLUSTERS:
        hits = [
            k for k, v in result.items() if set(cluster.components) & set(v)
        ]
        assert len(hits) == 1


def test_no_keep_similar_items_splits_proportionally_and_flat() -> None:
    result, _ = partition_dataset(
        CLUSTERS, 0.6, 0.2, 0.2, keep_similar_items=False
    )
    assert all(isinstance(item, str) for item in result["train"])
    assert sum(len(v) for v in result.values()) == sum(
        c.n_components for c in CLUSTERS
    )


def test_no_keep_similar_items_hits_targets_with_small_clusters() -> None:
    for size in (1, 2):
        clusters = [
            SimilarityCluster(
                components=[f"i{i}_{j}" for j in range(size)],
                n_components=size,
            )
            for i in range(100 // size)
        ]
        result, _ = partition_dataset(
            clusters, 0.6, 0.2, 0.2, keep_similar_items=False
        )
        assert [len(result[k]) for k in ("train", "val", "test")] == [
            60,
            20,
            20,
        ]


def test_bad_fractions_raise() -> None:
    try:
        partition_dataset(CLUSTERS, 0.6, 0.2, 0.3)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for fractions != 1.0")


KINDS = [
    [f"a{i}" for i in range(10)],
    [f"b{i}" for i in range(8)],
    [f"c{i}" for i in range(6)],
    ["h0", "h1"],
    ["h2"],
    ["h3"],
    ["s0"],
    ["s1"],
    ["s2"],
    ["s3"],
]
KIND_LABELS = {i: "homomer" for c in KINDS for i in c} | {
    h: "heteromer" for h in ("h0", "h1", "h2", "h3")
}


def _clusters(groups: list[list[str]]) -> list[SimilarityCluster]:
    return [
        SimilarityCluster(components=g, n_components=len(g)) for g in groups
    ]


def _kinds_in(items: list[str], labels: dict[str, str]) -> set[str]:
    return {labels.get(i, "unlabelled") for i in items}


def test_labels_reach_every_split() -> None:
    clusters = _clusters(KINDS)
    plain, _ = partition_dataset(clusters)
    assert any(
        "heteromer" not in _kinds_in(v, KIND_LABELS) for v in plain.values()
    )

    result, provenance = partition_dataset(clusters, labels=KIND_LABELS)
    for items in result.values():
        assert _kinds_in(items, KIND_LABELS) == {"homomer", "heteromer"}
    for cluster in clusters:
        hits = [
            k for k, v in result.items() if set(cluster.components) & set(v)
        ]
        assert len(hits) == 1
    assert sum(len(v) for v in result.values()) == 32
    assert provenance.missing_labels == {}
    assert sum(c["heteromer"] for c in provenance.label_counts.values()) == 4


def test_rare_label_missing_is_reported() -> None:
    groups = [*KINDS[:3], ["h0", "h1"], *KINDS[6:]]
    result, provenance = partition_dataset(
        _clusters(groups), labels=KIND_LABELS
    )
    assert len(provenance.missing_labels) == 2
    assert all(v == ["heteromer"] for v in provenance.missing_labels.values())


def test_unlabelled_items_are_a_kind() -> None:
    labels = {"a0": "homomer", "h0": "unlabelled"}
    _, provenance = partition_dataset(_clusters(KINDS), labels=labels)
    counts = provenance.label_counts
    assert sum(c.get("unlabelled", 0) for c in counts.values()) == 31
    assert sum(c.get("homomer", 0) for c in counts.values()) == 1


def test_no_labels_is_unchanged() -> None:
    result, provenance = partition_dataset(CLUSTERS)
    assert [len(result[k]) for k in ("train", "val", "test")] == [18, 6, 6]
    assert provenance.label_counts == {}
    assert provenance.missing_labels == {}
    assert provenance.labels_layer is None


def test_zero_fraction_split_skipped() -> None:
    result, provenance = partition_dataset(
        _clusters(KINDS), 0.8, 0.0, 0.2, labels=KIND_LABELS
    )
    assert provenance.split_sizes["val"] == 0
    assert "val" not in provenance.missing_labels
    for split in ("train", "test"):
        assert "heteromer" in _kinds_in(result[split], KIND_LABELS)


def test_labels_need_keep_similar_items() -> None:
    with pytest.raises(ValueError, match="keep_similar_items"):
        partition_dataset(
            _clusters(KINDS), labels=KIND_LABELS, keep_similar_items=False
        )


def test_labels_input_not_mutated() -> None:
    clusters = _clusters(KINDS)
    before = [c.model_dump() for c in clusters]
    labels = dict(KIND_LABELS)
    partition_dataset(clusters, labels=labels)
    assert [c.model_dump() for c in clusters] == before
    assert labels == KIND_LABELS


if __name__ == "__main__":
    test_keep_similar_items_hits_target_sizes()
    test_keep_similar_items_never_splits_a_cluster()
    test_no_keep_similar_items_splits_proportionally_and_flat()
    test_bad_fractions_raise()
    print("ok")
