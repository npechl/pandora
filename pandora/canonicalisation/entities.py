from __future__ import annotations

from collections import defaultdict

from pandora.schemas.canonicalisation import (
    EntityMapping,
    EntityMappingItem,
)
from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
    EntityRecord,
)


def _normalize_entities(
    entities: list[EntityRecord],
    asym_units: list[AsymRecord],
    atoms: list[AtomSiteRecord],
    rules,
    record: bool,
) -> tuple[
    list[EntityRecord],
    list[AsymRecord],
    list[AtomSiteRecord],
    EntityMapping,
    dict[str, str],
]:
    """Apply the entity id/merge strategy to entities, propagating ids to
    asym_units/atoms, per rules. Also returns the original -> canonical id
    map (empty when ids are kept)."""

    if rules.strategy == "preserve":
        mapping = EntityMapping()
        if record:
            for ent in entities:
                mapping.items.append(
                    EntityMappingItem(
                        canonical_entity_id=ent.id,
                        original_entity_ids=[ent.id],
                    )
                )
        return entities, asym_units, atoms, mapping, {}

    if rules.strategy == "standardize":
        new_entities, entity_id_map, mapping = _standardize_entities(
            entities, rules.preserve_original_entity_ids, record
        )
    elif rules.strategy == "merge_equivalent_entities":
        new_entities, entity_id_map, mapping = _merge_equivalent_entities(
            entities, record
        )
    else:
        return entities, asym_units, atoms, EntityMapping(), {}

    new_asyms, new_atoms = _propagate_entity_ids(
        asym_units, atoms, entity_id_map
    )
    return new_entities, new_asyms, new_atoms, mapping, entity_id_map


def _standardize_entities(
    entities: list[EntityRecord], preserve_original_ids: bool, record: bool
) -> tuple[list[EntityRecord], dict[str, str], EntityMapping]:
    """Entities renumbered 1, 2, ..., the id map, and the mapping."""

    mapping = EntityMapping()
    entity_id_map: dict[str, str] = {}
    new_entities = []
    for i, ent in enumerate(entities, 1):
        new_id = str(i)
        entity_id_map[ent.id] = new_id
        new_entities.append(ent.model_copy(update={"id": new_id}))
        if record:
            orig = ent.id if preserve_original_ids else new_id
            mapping.items.append(
                EntityMappingItem(
                    canonical_entity_id=new_id,
                    original_entity_ids=[orig],
                )
            )
    return new_entities, entity_id_map, mapping


def _merge_equivalent_entities(
    entities: list[EntityRecord], record: bool
) -> tuple[list[EntityRecord], dict[str, str], EntityMapping]:
    """Entities with the same canonical sequence merged into the first
    one, renumbered 1, 2, ..., the id map, and the mapping."""

    entity_id_map: dict[str, str] = {}
    seq_to_id: dict[str, str] = {}
    merged_orig: dict[str, list[str]] = defaultdict(list)
    new_entities = []
    for ent in entities:
        seq_key = (
            ent.poly.pdbx_seq_one_letter_code_can
            if ent.poly and ent.poly.pdbx_seq_one_letter_code_can
            else None
        )
        if seq_key and seq_key in seq_to_id:
            canonical_id = seq_to_id[seq_key]
            entity_id_map[ent.id] = canonical_id
            merged_orig[canonical_id].append(ent.id)
        else:
            canonical_id = str(len(new_entities) + 1)
            entity_id_map[ent.id] = canonical_id
            if seq_key:
                seq_to_id[seq_key] = canonical_id
            merged_orig[canonical_id] = [ent.id]
            new_entities.append(ent.model_copy(update={"id": canonical_id}))

    mapping = EntityMapping()
    if record:
        for cid, orig_ids in merged_orig.items():
            mapping.items.append(
                EntityMappingItem(
                    canonical_entity_id=cid,
                    original_entity_ids=orig_ids,
                )
            )
    return new_entities, entity_id_map, mapping


def _propagate_entity_ids(
    asym_units: list[AsymRecord],
    atoms: list[AtomSiteRecord],
    entity_id_map: dict[str, str],
) -> tuple[list[AsymRecord], list[AtomSiteRecord]]:
    """asym_units and atoms with entity ids rewritten through
    entity_id_map."""

    new_asyms = [
        a.model_copy(
            update={"entity_id": entity_id_map.get(a.entity_id, a.entity_id)}
        )
        for a in asym_units
    ]
    new_atoms = [
        (
            a.model_copy(
                update={
                    "label_entity_id": entity_id_map.get(
                        a.label_entity_id, a.label_entity_id
                    )
                }
            )
            if entity_id_map.get(a.label_entity_id, a.label_entity_id)
            != a.label_entity_id
            else a
        )
        for a in atoms
    ]
    return new_asyms, new_atoms


def _remap_raw_entity_ids(
    raw: dict[str, list[dict[str, str | None]]], entity_id_map: dict[str, str]
) -> dict[str, list[dict[str, str | None]]]:
    """Raw with every `entity_id` column in canonical ids. Rows of an entity
    merged into an earlier one are dropped: they would repeat its rows."""

    if not entity_id_map:
        return raw
    first_original = {}
    for original, canonical in entity_id_map.items():
        first_original.setdefault(canonical, original)
    remapped = {}
    for category, rows in raw.items():
        if not rows or "entity_id" not in rows[0]:
            remapped[category] = rows
            continue
        remapped[category] = [
            {**row, "entity_id": entity_id_map[row["entity_id"]]}
            if row["entity_id"] in entity_id_map
            else row
            for row in rows
            if row["entity_id"] not in entity_id_map
            or first_original[entity_id_map[row["entity_id"]]]
            == row["entity_id"]
        ]
    return remapped
