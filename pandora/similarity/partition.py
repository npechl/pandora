from __future__ import annotations

from pandora.schemas.similarity import SimilarityCluster

_FRACTION_TOLERANCE = 0.001


def partition_dataset(
    clusters: list[SimilarityCluster],
    pct_train: float = 0.6,
    pct_val: float = 0.2,
    pct_test: float = 0.2,
    keep_similar_items: bool = True,
) -> dict[str, list[str]]:
    """Assign cluster members to train/val/test splits.

    keep_similar_items=True (default): each cluster goes to a single split
    (whichever is most under its target share) — the leakage-safe mode.
    keep_similar_items=False: each cluster's members are divided
    proportionally across all three splits.
    """

    if abs(pct_train + pct_val + pct_test - 1.0) > _FRACTION_TOLERANCE:
        raise ValueError(
            f"pct_train + pct_val + pct_test must equal 1.0, got "
            f"{pct_train + pct_val + pct_test}"
        )

    ordered = sorted(clusters, key=lambda c: c.n_components, reverse=True)
    total = sum(c.n_components for c in ordered)

    targets = {
        "train": round(total * pct_train),
        "val": round(total * pct_val),
        "test": round(total * pct_test),
    }
    splits: dict[str, list[str]] = {"train": [], "val": [], "test": []}

    if keep_similar_items:
        for cluster in ordered:
            split = max(
                splits, key=lambda name: targets[name] - len(splits[name])
            )
            splits[split].extend(cluster.components)
        return splits

    for cluster in ordered:
        n_train = round(cluster.n_components * pct_train)
        n_val = round(cluster.n_components * pct_val)
        splits["train"].extend(cluster.components[:n_train])
        splits["val"].extend(cluster.components[n_train : n_train + n_val])
        splits["test"].extend(cluster.components[n_train + n_val :])

    return splits
