from __future__ import annotations

from pandora.schemas.similarity import SimilarityCluster, SimilarityRelationship


def cluster_similar_items(
    item_ids: list[str],
    relationships: list[SimilarityRelationship],
    threshold: float,
) -> list[SimilarityCluster]:
    """Group item ids into clusters via connected components.

    Two items land in the same cluster iff connected through a chain of
    relationships each scoring >= threshold. Items with no such edges
    (isolates) form their own singleton cluster, so every id in
    `item_ids` ends up in exactly one cluster.
    """

    parent = {item_id: item_id for item_id in item_ids}

    def find(item_id: str) -> str:
        while parent[item_id] != item_id:
            parent[item_id] = parent[parent[item_id]]
            item_id = parent[item_id]
        return item_id

    for rel in relationships:
        if (
            rel.score >= threshold
            and rel.source_id in parent
            and rel.target_id in parent
        ):
            root_source, root_target = find(rel.source_id), find(rel.target_id)
            if root_source != root_target:
                parent[root_target] = root_source

    groups: dict[str, list[str]] = {}
    for item_id in item_ids:
        groups.setdefault(find(item_id), []).append(item_id)

    return [
        SimilarityCluster(components=sorted(members), n_components=len(members))
        for members in sorted(groups.values(), key=lambda members: members[0])
    ]
