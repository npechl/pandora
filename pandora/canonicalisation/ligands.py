from __future__ import annotations

import re

from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
    EntityRecord,
)

from pandora.schemas.canonicalisation import LigandRules
from pandora.schemas.common import Diagnostic, DiagnosticBundle


def _keep(entity_type: str, description: str, rules: LigandRules) -> bool:
    """Whether an asym unit of this entity type and (upper-cased)
    description is kept under the ligand rules."""

    if entity_type == "polymer":
        return True
    if rules.strategy == "annotate_only":
        return False
    if entity_type == "water":
        return rules.keep_waters
    # PDB names ions "<ELEMENT/ANION> ION" (ZINC ION, SULFATE ION).
    # Whole-word match only: "PHOSPHATE" alone would also catch
    # GMP/ATP/NADP, and "ION" as a substring catches glutathione.
    if "ION" in re.findall(r"[A-Z0-9]+", description):
        return rules.keep_ions
    return rules.keep_nonpolymer_ligands


def filter_ligands(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    entities: list[EntityRecord],
    rules: LigandRules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> tuple[list[AtomSiteRecord], list[AsymRecord]]:
    """Filter non-polymer ligands/waters/ions out of atoms/asym_units per
    the ligand rules.

    Categorizes each non-polymer asym unit as water, ion (by keyword
    match on the entity description), or other non-polymer ligand, and
    keeps or drops it according to `rules.keep_waters`/`keep_ions`/
    `keep_nonpolymer_ligands`. With `rules.strategy="annotate_only"`,
    every non-polymer asym unit is dropped from the returned atoms/
    asym_units (but recorded as a diagnostic) rather than selectively
    kept.

    Args:
        atoms: The structure's atom records.
        asym_units: The structure's asym unit records.
        entities: The structure's entity records, used to classify each
            asym unit's entity type and description.
        rules: The `LigandRules` governing which categories are kept.
            No-op (returns `atoms, asym_units` unchanged) unless
            `rules.strategy` is "filter" or "annotate_only".
        diagnostics: Bundle to append a warning to for each asym unit
            excluded under `annotate_only`.
        entry_id: The structure's entry id, attached to any diagnostics
            raised.

    Returns:
        `(atoms, asym_units)` with excluded ligand atoms/asym units
        removed; polymer atoms/asym units are always kept.
    """

    if rules.strategy not in ("filter", "annotate_only"):
        return atoms, asym_units

    entity_type: dict[str, str] = {e.id: e.type for e in entities}
    entity_desc: dict[str, str] = {
        e.id: (e.pdbx_description or "").upper() for e in entities
    }
    asym_entity: dict[str, str] = {a.id: a.entity_id for a in asym_units}

    keep = {
        a.id
        for a in asym_units
        if asym_entity.get(a.id) is None
        or _keep(
            entity_type.get(asym_entity[a.id], "polymer"),
            entity_desc.get(asym_entity[a.id], ""),
            rules,
        )
    }

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
