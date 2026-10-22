from __future__ import annotations

from pandora._util import now_iso
from pandora.schemas.similarity import PartitionProvenance, SimilarityCluster

_FRACTION_TOLERANCE = 0.001


def partition_dataset(
    clusters: list[SimilarityCluster],
    pct_train: float = 0.6,
    pct_val: float = 0.2,
    pct_test: float = 0.2,
    keep_similar_items: bool = True,
) -> tuple[dict[str, list[str]], PartitionProvenance]:
    """Assign cluster members to train/val/test splits.

    keep_similar_items=True (default): each cluster goes to a single split
    (whichever is most under its target share) — the leakage-safe mode.
    keep_similar_items=False: each cluster's members are divided
    proportionally across all three splits.

    Args:
        clusters: The similarity clusters to assign, as returned by
            `cluster_similar_items`.
        pct_train: Target fraction of items assigned to the train split.
        pct_val: Target fraction of items assigned to the val split.
        pct_test: Target fraction of items assigned to the test split.
            `pct_train + pct_val + pct_test` must equal 1.0.
        keep_similar_items: If True (default), each cluster is assigned
            whole to a single split (leakage-safe). If False, each
            cluster's members are divided proportionally across splits.

    Returns:
        `(splits, provenance)` — the split assignment, and a record of
        the fractions/mode applied and the resulting split sizes (for
        reporting how a dataset's split was built).

    Raises:
        ValueError: `pct_train + pct_val + pct_test` does not equal 1.0.
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
    else:
        # Each item goes to the split furthest below its share of the items
        # seen so far; rounding per cluster sent every singleton to train.
        fractions = {"train": pct_train, "val": pct_val, "test": pct_test}
        seen = 0
        for cluster in ordered:
            for item_id in cluster.components:
                seen += 1
                split = max(
                    splits,
                    key=lambda name: fractions[name] * seen - len(splits[name]),
                )
                splits[split].append(item_id)

    provenance = PartitionProvenance(
        partitioned_at=now_iso(),
        pct_train=pct_train,
        pct_val=pct_val,
        pct_test=pct_test,
        keep_similar_items=keep_similar_items,
        split_sizes={name: len(members) for name, members in splits.items()},
    )
    return splits, provenance
