from __future__ import annotations

from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
    EntityRecord,
)

from pandora.schemas.common import Diagnostic, DiagnosticBundle


def _filter_ligands(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    entities: list[EntityRecord],
    rules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> tuple[list[AtomSiteRecord], list[AsymRecord]]:
    if rules.strategy not in ("filter", "annotate_only"):
        return atoms, asym_units

    entity_type: dict[str, str] = {e.id: e.type for e in entities}
    entity_desc: dict[str, str] = {
        e.id: (e.pdbx_description or "").upper() for e in entities
    }
    asym_entity: dict[str, str] = {a.id: a.entity_id for a in asym_units}

    _ION_KEYWORDS = frozenset(
        {
            "ION",
            "ZINC",
            "CALCIUM",
            "MAGNESIUM",
            "SODIUM",
            "POTASSIUM",
            "IRON",
            "COPPER",
            "MANGANESE",
            "COBALT",
            "NICKEL",
            "CHLORIDE",
            "SULFATE",
            "PHOSPHATE",
        }
    )

    def _keep(asym_id: str) -> bool:
        eid = asym_entity.get(asym_id)
        if eid is None:
            return True
        etype = entity_type.get(eid, "polymer")
        if etype == "polymer":
            return True
        if rules.strategy == "annotate_only":
            return False
        if etype == "water":
            return rules.keep_waters
        desc = entity_desc.get(eid, "")
        is_ion = any(kw in desc for kw in _ION_KEYWORDS)
        if is_ion:
            return rules.keep_ions
        return rules.keep_nonpolymer_ligands

    keep = {a.id for a in asym_units if _keep(a.id)}

    if rules.strategy == "annotate_only":
        for a in asym_units:
            if a.id in keep:
                continue
            eid = asym_entity.get(a.id)
            diagnostics.warnings.append(
                Diagnostic(
                    code="LIGAND_ANNOTATED_ONLY",
                    severity="warning",
                    message=(
                        f"Ligand in asym {a.id} excluded from canonical "
                        "structure (ligand_rules.strategy=annotate_only)"
                    ),
                    entry_id=entry_id,
                    context={
                        "asym_id": a.id,
                        "entity_id": eid,
                        "entity_type": entity_type.get(eid, "polymer"),
                        "description": entity_desc.get(eid, ""),
                    },
                )
            )

    return (
        [a for a in atoms if a.label_asym_id in keep],
        [a for a in asym_units if a.id in keep],
    )
