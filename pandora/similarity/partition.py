from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Mapping

from pandora._util import now_iso
from pandora.schemas.similarity import PartitionProvenance, SimilarityCluster

_FRACTION_TOLERANCE = 0.001
_UNLABELLED = "unlabelled"


def partition_dataset(
    clusters: list[SimilarityCluster],
    pct_train: float = 0.6,
    pct_val: float = 0.2,
    pct_test: float = 0.2,
    keep_similar_items: bool = True,
    labels: Mapping[str, str] | None = None,
) -> tuple[dict[str, list[str]], PartitionProvenance]:
    """Assign cluster members to train/val/test splits.

    keep_similar_items=True (default): each cluster goes to a single split
    (whichever is most under its target share) — the leakage-safe mode.
    keep_similar_items=False: each cluster's members are divided
    proportionally across all three splits.

    With `labels`, clusters are first placed so that every label appears
    in every split with a non-zero fraction, where keeping clusters whole
    allows: the rarest labels go first, each using its smallest clusters.
    The rest are then placed by size as without labels. Labels that
    could not reach a split are reported in the provenance, not raised.

    Args:
        clusters: The similarity clusters to assign, as returned by
            `cluster_edges`.
        pct_train: Target fraction of items assigned to the train split.
        pct_val: Target fraction of items assigned to the val split.
        pct_test: Target fraction of items assigned to the test split.
            `pct_train + pct_val + pct_test` must equal 1.0.
        keep_similar_items: If True (default), each cluster is assigned
            whole to a single split (leakage-safe). If False, each
            cluster's members are divided proportionally across splits.
        labels: Item ID -> label (e.g. from `item_labels`). Items not in
            it count as the label `"unlabelled"`.

    Returns:
        `(splits, provenance)` — the split assignment, and a record of
        the fractions/mode applied, the resulting split sizes and, with
        `labels`, the labels per split and those missing from a split.

    Raises:
        ValueError: `pct_train + pct_val + pct_test` does not equal 1.0,
            or `labels` is given with `keep_similar_items=False`.
    """

    if abs(pct_train + pct_val + pct_test - 1.0) > _FRACTION_TOLERANCE:
        raise ValueError(
            f"pct_train + pct_val + pct_test must equal 1.0, got "
            f"{pct_train + pct_val + pct_test}"
        )
    if labels is not None and not keep_similar_items:
        raise ValueError("labels need keep_similar_items=True")

    ordered = sorted(clusters, key=lambda c: c.n_components, reverse=True)
    total = sum(c.n_components for c in ordered)

    fractions = {"train": pct_train, "val": pct_val, "test": pct_test}
    targets = {name: round(total * f) for name, f in fractions.items()}
    splits: dict[str, list[str]] = {"train": [], "val": [], "test": []}

    def deficit(name: str) -> int:
        return targets[name] - len(splits[name])

    active = [name for name, f in fractions.items() if f > 0]
    if keep_similar_items:
        placed: set[int] = set()
        if labels is not None:
            placed = _cover_labels(clusters, labels, active, splits, deficit)
        for cluster in sorted(
            (c for i, c in enumerate(clusters) if i not in placed),
            key=lambda c: c.n_components,
            reverse=True,
        ):
            splits[max(splits, key=deficit)].extend(cluster.components)
    else:
        # Each item goes to the split furthest below its share of the items
        # seen so far; rounding per cluster sent every singleton to train.
        seen = 0
        for cluster in ordered:
            for item_id in cluster.components:
                seen += 1
                split = max(
                    splits,
                    key=lambda name: fractions[name] * seen - len(splits[name]),
                )
                splits[split].append(item_id)

    label_counts: dict[str, dict[str, int]] = {}
    missing_labels: dict[str, list[str]] = {}
    if labels is not None:
        label_counts, missing_labels = _label_report(
            clusters, labels, splits, active
        )

    provenance = PartitionProvenance(
        partitioned_at=now_iso(),
        pct_train=pct_train,
        pct_val=pct_val,
        pct_test=pct_test,
        keep_similar_items=keep_similar_items,
        split_sizes={name: len(members) for name, members in splits.items()},
        label_counts=label_counts,
        missing_labels=missing_labels,
    )
    return splits, provenance


def _label(labels: Mapping[str, str], item_id: str) -> str:
    """An item's label, `unlabelled` if it has none."""

    return labels.get(item_id, _UNLABELLED)


def _cover_labels(
    clusters: list[SimilarityCluster],
    labels: Mapping[str, str],
    active: list[str],
    splits: dict[str, list[str]],
    deficit: Callable[[str], int],
) -> set[int]:
    """Place clusters so each label reaches each active split; return
    the indexes of the placed clusters."""

    kinds = [{_label(labels, i) for i in c.components} for c in clusters]
    by_label: dict[str, list[int]] = {}
    for i in sorted(
        range(len(clusters)),
        key=lambda i: (clusters[i].n_components, clusters[i].components[0]),
    ):
        for label in kinds[i]:
            by_label.setdefault(label, []).append(i)
    covered: dict[str, set[str]] = {name: set() for name in splits}
    placed: set[int] = set()
    for label in sorted(by_label, key=lambda x: (len(by_label[x]), x)):
        candidates = (i for i in by_label[label] if i not in placed)
        for i in candidates:
            lacking = [n for n in active if label not in covered[n]]
            if not lacking:
                break
            split = max(lacking, key=deficit)
            splits[split].extend(clusters[i].components)
            covered[split] |= kinds[i]
            placed.add(i)
    return placed


def _label_report(
    clusters: list[SimilarityCluster],
    labels: Mapping[str, str],
    splits: dict[str, list[str]],
    active: list[str],
) -> tuple[dict[str, dict[str, int]], dict[str, list[str]]]:
    """Labels per split, and the labels missing from each active split."""

    counts = {
        name: dict(Counter(_label(labels, i) for i in members))
        for name, members in splits.items()
    }
    every = {_label(labels, i) for c in clusters for i in c.components}
    missing = {
        name: sorted(every - counts[name].keys())
        for name in active
        if every - counts[name].keys()
    }
    return counts, missing
