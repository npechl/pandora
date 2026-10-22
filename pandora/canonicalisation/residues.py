from __future__ import annotations
from collections import defaultdict
from collections.abc import Callable
from pandora.schemas.structure import AtomSiteRecord
from pandora.schemas.canonicalisation import (
    ResidueNumberMapping,
    ResidueNumberMappingItem,
)


def _residue_key(a: AtomSiteRecord) -> tuple:
    """What tells a's residue apart from the others in its chain."""

    return (
        a.label_seq_id,
        a.auth_seq_id,
        a.label_comp_id,
        a.pdbx_PDB_ins_code,
    )


def _auth_seq_number(a: AtomSiteRecord) -> int | None:
    """a's auth_seq_id as an int, or its label_seq_id if that fails."""

    raw = a.auth_seq_id or ""
    digits = "".join(c for c in raw if c.isdigit() or c == "-")
    try:
        return int(digits)
    except ValueError:
        return a.label_seq_id


def _sequential_numbering(
    atoms: list[AtomSiteRecord],
) -> Callable[[AtomSiteRecord], int | None]:
    """Number for each polymer residue: 1, 2, ... per chain, in order of
    first appearance."""

    chain_res_order: dict[str, list[tuple]] = defaultdict(list)
    chain_res_seen: dict[str, set] = defaultdict(set)
    for a in atoms:
        if a.label_seq_id is None:
            continue
        rk = _residue_key(a)
        if rk not in chain_res_seen[a.label_asym_id]:
            chain_res_seen[a.label_asym_id].add(rk)
            chain_res_order[a.label_asym_id].append(rk)

    renumber: dict[str, dict[tuple, int]] = {
        ch: {rk: i + 1 for i, rk in enumerate(rks)}
        for ch, rks in chain_res_order.items()
    }
    return lambda a: renumber[a.label_asym_id].get(
        _residue_key(a), a.label_seq_id
    )


def _apply_numbering(
    atoms: list[AtomSiteRecord],
    new_seq_id: Callable[[AtomSiteRecord], int | None],
    keep_insertion_codes: bool,
    chain_map: dict[str, str],
    record: bool,
) -> tuple[list[AtomSiteRecord], ResidueNumberMapping]:
    """Give each polymer atom label_seq_id new_seq_id(atom), drop insertion
    codes unless keep_insertion_codes, and record each residue's original
    numbering if record."""

    mapping = ResidueNumberMapping()
    seen: set[tuple] = set()
    # Atoms already carry canonical chain ids (the chain map ran first).
    original_chain = {new: old for old, new in chain_map.items()}

    result = []
    for a in atoms:
        ins = a.pdbx_PDB_ins_code if keep_insertion_codes else None
        if a.label_seq_id is None:
            result.append(a.model_copy(update={"pdbx_PDB_ins_code": ins}))
            continue
        new_seq = new_seq_id(a)

        key = (a.label_asym_id, a.label_seq_id, a.auth_seq_id)
        if record and key not in seen:
            seen.add(key)
            mapping.items.append(
                ResidueNumberMappingItem(
                    canonical_chain_id=a.label_asym_id,
                    canonical_seq_id=new_seq or 0,
                    original_chain_id=original_chain.get(
                        a.label_asym_id, a.label_asym_id
                    ),
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


def _normalize_residue_numbering(
    atoms: list[AtomSiteRecord],
    chain_map: dict[str, str],
    strategy: str,
    preserve_insertion_codes: bool,
    record: bool,
) -> tuple[list[AtomSiteRecord], ResidueNumberMapping]:
    """Renumber residues (or fall back to auth numbering) per the
    residue-numbering strategy."""

    if strategy == "preserve":
        if preserve_insertion_codes:
            return atoms, ResidueNumberMapping()
        return [
            a.model_copy(update={"pdbx_PDB_ins_code": None})
            if a.pdbx_PDB_ins_code not in (None, "", ".", "?")
            else a
            for a in atoms
        ], ResidueNumberMapping()

    # Only polymer residues (label_seq_id set) are numbered; ligands and
    # waters keep label_seq_id None, as mmCIF requires.
    if strategy == "use_auth_seq":
        return _apply_numbering(
            atoms, _auth_seq_number, preserve_insertion_codes, chain_map, record
        )
    # renumber
    return _apply_numbering(
        atoms, _sequential_numbering(atoms), False, chain_map, record
    )
