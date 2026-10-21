from __future__ import annotations

from pandora.schemas.canonicalisation import (
    ModifiedResidueMapping,
    ModifiedResidueMappingItem,
    ModifiedResidueRules,
)
from pandora.schemas.common import Diagnostic, DiagnosticBundle
from pandora.schemas.structure import AsymRecord, AtomSiteRecord, ConnRecord

STANDARD_RESIDUES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "UNK",
        "A", "C", "G", "U", "N",
        "DA", "DC", "DG", "DT", "DU", "DN",
    }
)  # fmt: skip

# Parents assumed even when _pdbx_struct_mod_residue doesn't list them.
_BUILTIN_PARENTS = {"MSE": "MET"}


def _map_modified_residues(
    atoms: list[AtomSiteRecord],
    mod_rows: list[dict[str, str | None]],
    rules: ModifiedResidueRules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> tuple[list[AtomSiteRecord], ModifiedResidueMapping]:
    """Rename modified polymer residues to their parent, per the rules."""

    mapping = ModifiedResidueMapping()
    if rules.strategy == "preserve":
        return atoms, mapping

    parents = {
        (row.get("label_asym_id"), row.get("label_seq_id"), comp): parent
        for row in mod_rows
        if (comp := row.get("label_comp_id"))
        and (parent := row.get("parent_comp_id"))
    }
    wanted = set(rules.comp_ids)
    seen: set[tuple[str, str, str]] = set()
    result: list[AtomSiteRecord] = []
    for a in atoms:
        comp = a.label_comp_id
        if (
            a.label_seq_id is None
            or comp in STANDARD_RESIDUES
            or (wanted and comp not in wanted)
        ):
            result.append(a)
            continue

        key = (a.label_asym_id, str(a.label_seq_id), comp)
        parent = parents.get(key) or _BUILTIN_PARENTS.get(comp)
        first = key not in seen
        seen.add(key)
        if parent is None:
            if first:
                diagnostics.warnings.append(
                    Diagnostic(
                        code="MODIFIED_RESIDUE_UNMAPPED",
                        severity="warning",
                        message=(
                            f"No parent known for {comp} {a.label_seq_id} "
                            f"in chain {a.label_asym_id}; left unmapped"
                        ),
                        entry_id=entry_id,
                        context={
                            "chain": a.label_asym_id,
                            "seq_id": a.label_seq_id,
                            "comp_id": comp,
                        },
                    )
                )
            result.append(a)
            continue

        if first:
            mapping.items.append(
                ModifiedResidueMappingItem(
                    chain_id=a.label_asym_id,
                    seq_id=a.label_seq_id,
                    auth_seq_id=a.auth_seq_id,
                    original_comp_id=comp,
                    parent_comp_id=parent,
                )
            )
        update: dict[str, object] = {
            "label_comp_id": parent,
            "auth_comp_id": parent,
            "group_PDB": "ATOM",
        }
        if comp == "MSE" and a.label_atom_id == "SE":
            update |= {
                "label_atom_id": "SD",
                "auth_atom_id": "SD",
                "type_symbol": "S",
            }
        result.append(a.model_copy(update=update))

    return result, mapping


def _renamed_residues(
    mapping: ModifiedResidueMapping,
) -> dict[tuple[str, int | None], ModifiedResidueMappingItem]:
    """Mapped residues by (chain, seq_id)."""

    return {(i.chain_id, i.seq_id): i for i in mapping.items}


def _map_connections(
    connections: list[ConnRecord], mapping: ModifiedResidueMapping
) -> list[ConnRecord]:
    """Connections renamed like the atoms: comp_id to its parent, and an
    MSE's SE atom to SD."""

    renamed = _renamed_residues(mapping)
    if not renamed:
        return connections

    def side(partner):
        item = renamed.get((partner.label_asym_id, partner.label_seq_id))
        if item is None or partner.label_comp_id != item.original_comp_id:
            return partner
        update = {"label_comp_id": item.parent_comp_id}
        if item.original_comp_id == "MSE" and partner.label_atom_id == "SE":
            update["label_atom_id"] = "SD"
        return partner.model_copy(update=update)

    return [
        c.model_copy(update={"ptnr1": side(c.ptnr1), "ptnr2": side(c.ptnr2)})
        for c in connections
    ]


def _map_entity_poly_seq(
    raw: dict[str, list[dict[str, str | None]]],
    mapping: ModifiedResidueMapping,
    asym_units: list[AsymRecord],
) -> dict[str, list[dict[str, str | None]]]:
    """Raw with `_entity_poly_seq` mon_ids renamed at mapped positions."""

    rows = raw.get("_entity_poly_seq")
    if not rows or not mapping.items:
        return raw
    entity_of = {a.id: a.entity_id for a in asym_units}
    parent_at = {
        (entity_of.get(i.chain_id), str(i.seq_id), i.original_comp_id): (
            i.parent_comp_id
        )
        for i in mapping.items
    }
    mapped = [
        {**row, "mon_id": parent}
        if (
            parent := parent_at.get(
                (row.get("entity_id"), row.get("num"), row.get("mon_id"))
            )
        )
        else row
        for row in rows
    ]
    return {**raw, "_entity_poly_seq": mapped}
