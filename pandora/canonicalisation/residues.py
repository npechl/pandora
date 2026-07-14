from __future__ import annotations
from collections import defaultdict
from pandora.schemas.structure import AtomSiteRecord
from pandora.schemas.canonicalisation import (
    ResidueNumberMapping,
    ResidueNumberMappingItem,
)

def _normalize_residue_numbering(
    atoms: list[AtomSiteRecord],
    chain_map: dict[str, str],
    strategy: str,
    preserve_insertion_codes: bool,
    record: bool,
) -> tuple[list[AtomSiteRecord], ResidueNumberMapping]:
    mapping = ResidueNumberMapping()

    if strategy == "preserve":
        if preserve_insertion_codes:
            return atoms, mapping
        result = []
        for a in atoms:
            ins = a.pdbx_PDB_ins_code
            if ins and ins not in (".", "?"):
                a = a.model_copy(update={"pdbx_PDB_ins_code": None})
            result.append(a)
        return result, mapping

    seen: set[tuple] = set()

    if strategy == "use_auth_seq":
        result = []
        for a in atoms:
            raw = a.auth_seq_id or ""
            digits = "".join(c for c in raw if c.isdigit() or c == "-")
            try:
                new_seq = int(digits)
            except ValueError:
                new_seq = a.label_seq_id

            ins = a.pdbx_PDB_ins_code if preserve_insertion_codes else None

            key = (a.label_asym_id, a.label_seq_id, a.auth_seq_id)
            if record and key not in seen:
                seen.add(key)
                mapping.items.append(
                    ResidueNumberMappingItem(
                        canonical_chain_id=chain_map.get(
                            a.label_asym_id, a.label_asym_id
                        ),
                        canonical_seq_id=new_seq or 0,
                        original_chain_id=a.label_asym_id,
                        original_seq_id=a.label_seq_id,
                        original_auth_seq_id=a.auth_seq_id or "",
                        original_insertion_code=a.pdbx_PDB_ins_code,
                    )
                )
            result.append(
                a.model_copy(
                    update={"label_seq_id": new_seq, "pdbx_PDB_ins_code": ins}
                )
            )
        return result, mapping

    # strategy == "renumber"
    chain_res_order: dict[str, list[tuple]] = defaultdict(list)
    chain_res_seen: dict[str, set] = defaultdict(set)
    for a in atoms:
        rk = (a.label_seq_id, a.label_comp_id, a.pdbx_PDB_ins_code)
        if rk not in chain_res_seen[a.label_asym_id]:
            chain_res_seen[a.label_asym_id].add(rk)
            chain_res_order[a.label_asym_id].append(rk)

    renumber: dict[str, dict[tuple, int]] = {
        ch: {rk: i + 1 for i, rk in enumerate(rks)}
        for ch, rks in chain_res_order.items()
    }

    result = []
    for a in atoms:
        rk = (a.label_seq_id, a.label_comp_id, a.pdbx_PDB_ins_code)
        new_seq = renumber[a.label_asym_id].get(rk, a.label_seq_id)

        key = (a.label_asym_id, a.label_seq_id, a.auth_seq_id)
        if record and key not in seen:
            seen.add(key)
            mapping.items.append(
                ResidueNumberMappingItem(
                    canonical_chain_id=chain_map.get(
                        a.label_asym_id, a.label_asym_id
                    ),
                    canonical_seq_id=new_seq or 0,
                    original_chain_id=a.label_asym_id,
                    original_seq_id=a.label_seq_id,
                    original_auth_seq_id=a.auth_seq_id or "",
                    original_insertion_code=a.pdbx_PDB_ins_code,
                )
            )
        result.append(
            a.model_copy(
                update={"label_seq_id": new_seq, "pdbx_PDB_ins_code": None}
            )
        )

    return result, mapping

