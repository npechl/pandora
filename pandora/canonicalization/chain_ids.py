from __future__ import annotations

import string

from typing import Iterator

from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
    Structure,
)
from pandora.schemas.canonicalization import (
    ChainIdMapping,
    ChainIdMappingItem,
)

# dump function for remap
def _sequential_chain_ids() -> Iterator[str]:
    letters = string.ascii_uppercase
    for c in letters:
        yield c
    for c1 in letters:
        for c2 in letters:
            yield c1 + c2


def _normalize_chain_ids(
    asym_units: list[AsymRecord],
    strategy: str,
    record: bool,
) -> tuple[dict[str, str], ChainIdMapping]:
    chain_map: dict[str, str] = {}
    mapping = ChainIdMapping()
    id_gen = _sequential_chain_ids()

    for asym in asym_units:
        if strategy == "preserve":
            canonical = asym.id
        elif strategy == "use_auth_chain_id":
            canonical = asym.auth_id or asym.id
        else:  # remap
            canonical = next(id_gen)

        chain_map[asym.id] = canonical
        if record:
            mapping.items.append(
                ChainIdMappingItem(
                    canonical_chain_id=canonical,
                    original_chain_id=asym.id,
                    original_auth_chain_id=asym.auth_id or asym.id,
                )
            )

    return chain_map, mapping


def _apply_chain_map(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    structure: Structure,
    chain_map: dict[str, str],
) -> tuple[list[AtomSiteRecord], list[AsymRecord], Structure]:

    new_atoms = [
        a.model_copy(
            update={
                "label_asym_id": chain_map.get(a.label_asym_id, a.label_asym_id)
            }
        )
        if chain_map.get(a.label_asym_id, a.label_asym_id) != a.label_asym_id
        else a
        for a in atoms
    ]

    new_asyms = [
        AsymRecord(
            id=chain_map.get(a.id, a.id),
            entity_id=a.entity_id,
            auth_id=a.auth_id,
        )
        for a in asym_units
    ]

    # Update assembly generators
    new_assemblies = []
    for asm in structure.assemblies:
        new_gens = []
        for gen in asm.generators:
            new_asym_list = [
                chain_map.get(aid, aid) for aid in gen.asym_id_list
            ]
            new_gens.append(
                gen.model_copy(update={"asym_id_list": new_asym_list})
            )
        new_assemblies.append(asm.model_copy(update={"generators": new_gens}))

    # Update connections
    new_conns = []
    for conn in structure.connections:
        p1 = conn.ptnr1.model_copy(
            update={
                "label_asym_id": chain_map.get(
                    conn.ptnr1.label_asym_id, conn.ptnr1.label_asym_id
                )
            }
        )
        p2 = conn.ptnr2.model_copy(
            update={
                "label_asym_id": chain_map.get(
                    conn.ptnr2.label_asym_id, conn.ptnr2.label_asym_id
                )
            }
        )
        new_conns.append(conn.model_copy(update={"ptnr1": p1, "ptnr2": p2}))

    # Update secondary structure
    new_conf = [
        r.model_copy(
            update={
                "beg_label_asym_id": chain_map.get(
                    r.beg_label_asym_id, r.beg_label_asym_id
                ),
                "end_label_asym_id": chain_map.get(
                    r.end_label_asym_id, r.end_label_asym_id
                ),
            }
        )
        for r in structure.secondary_structure.conf_records
    ]

    new_strands = [
        r.model_copy(
            update={
                "beg_label_asym_id": chain_map.get(
                    r.beg_label_asym_id, r.beg_label_asym_id
                ),
                "end_label_asym_id": chain_map.get(
                    r.end_label_asym_id, r.end_label_asym_id
                ),
            }
        )
        for r in structure.secondary_structure.sheet_strands
    ]

    new_ss = structure.secondary_structure.model_copy(
        update={"conf_records": new_conf, "sheet_strands": new_strands}
    )

    new_structure = structure.model_copy(
        update={
            "assemblies": new_assemblies,
            "connections": new_conns,
            "secondary_structure": new_ss,
        }
    )

    return new_atoms, new_asyms, new_structure
