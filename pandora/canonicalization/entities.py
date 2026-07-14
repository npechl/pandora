from __future__ import annotations

from collections import defaultdict

from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
    EntityRecord,
)
from pandora.schemas.canonicalization import (
    EntityMapping,
    EntityMappingItem,
)


def _normalize_entities(
    entities: list[EntityRecord],
    asym_units: list[AsymRecord],
    atoms: list[AtomSiteRecord],
    rules,
    record: bool,
) -> tuple[
    list[EntityRecord], list[AsymRecord], list[AtomSiteRecord], EntityMapping
]:
    mapping = EntityMapping()

    if rules.strategy == "preserve":
        if record:
            for ent in entities:
                mapping.items.append(
                    EntityMappingItem(
                        canonical_entity_id=ent.id,
                        original_entity_ids=[ent.id],
                    )
                )
        return entities, asym_units, atoms, mapping

    entity_id_map: dict[str, str] = {}

    if rules.strategy == "standardize":
        new_entities = []
        for i, ent in enumerate(entities, 1):
            new_id = str(i)
            entity_id_map[ent.id] = new_id
            new_entities.append(ent.model_copy(update={"id": new_id}))
            if record:
                orig = ent.id if rules.preserve_original_entity_ids else new_id
                mapping.items.append(
                    EntityMappingItem(
                        canonical_entity_id=new_id,
                        original_entity_ids=[orig],
                    )
                )

    elif rules.strategy == "merge_equivalent_entities":
        seq_to_id: dict[str, str] = {}
        merged_orig: dict[str, list[str]] = defaultdict(list)
        new_entities = []
        counter = 1
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
                canonical_id = str(counter)
                counter += 1
                entity_id_map[ent.id] = canonical_id
                if seq_key:
                    seq_to_id[seq_key] = canonical_id
                merged_orig[canonical_id] = [ent.id]
                new_entities.append(ent.model_copy(update={"id": canonical_id}))

        if record:
            for cid, orig_ids in merged_orig.items():
                mapping.items.append(
                    EntityMappingItem(
                        canonical_entity_id=cid,
                        original_entity_ids=orig_ids,
                    )
                )
    else:
        return entities, asym_units, atoms, mapping

    new_asyms = [
        a.model_copy(
            update={"entity_id": entity_id_map.get(a.entity_id, a.entity_id)}
        )
        for a in asym_units
    ]
    new_atoms = [
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
        for a in atoms
    ]
    return new_entities, new_asyms, new_atoms, mapping
