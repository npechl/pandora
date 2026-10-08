from __future__ import annotations

from pandora.schemas.canonicalisation import (
    ModifiedResidueMapping,
    ModifiedResidueMappingItem,
    ModifiedResidueRules,
)
from pandora.schemas.common import Diagnostic, DiagnosticBundle
from pandora.schemas.structure import AtomSiteRecord

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
