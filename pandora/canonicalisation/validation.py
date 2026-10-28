from __future__ import annotations

from collections import defaultdict

from pandora.canonicalisation.altlocs import _has_altloc
from pandora.schemas.common import Diagnostic, DiagnosticBundle
from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
)


def _validate(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    rules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> str:
    """Check for chain-id/residue-number collisions and compute a
    "failed"/"warning"/"success" status per validation_rules."""

    chain_ids = [a.id for a in asym_units]
    if len(chain_ids) != len(set(chain_ids)):
        diagnostics.errors.append(
            Diagnostic(
                code="CANONICAL_CHAIN_ID_COLLISION",
                severity="error",
                message="Duplicate canonical chain IDs",
                entry_id=entry_id,
            )
        )

    _check_residue_collisions(atoms, diagnostics, entry_id)

    if rules.strictness == "strict":
        for w in list(diagnostics.warnings):
            diagnostics.errors.append(
                w.model_copy(update={"severity": "error"})
            )

    has_errors = bool(diagnostics.errors)
    has_warnings = bool(diagnostics.warnings)

    if has_errors and rules.fail_on_unresolved_issues:
        return "failed"
    if rules.warnings_as_errors and has_warnings:
        return "failed"
    if has_errors or has_warnings:
        return "warning"
    return "success"


def _check_residue_collisions(
    atoms: list[AtomSiteRecord],
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> None:
    """Add an error to diagnostics for each residue-number collision."""

    # Distinct named altlocs can carry different comp_ids at one slot.
    # Without an explicit altloc, an atom is shared: conservatively check
    # its comp_id against every variant. Repeated atoms are not collisions.
    chain_res_comp: dict[str, dict[tuple, dict[str | None, str]]] = defaultdict(
        dict
    )
    # Two distinct representatives suffice: any incoming shared atom must
    # differ from at least one of them. Never scan every named variant.
    chain_res_representatives: dict[str, dict[tuple, list[str]]] = defaultdict(
        dict
    )
    flagged: set[tuple[str, tuple]] = set()
    for a in atoms:
        if a.group_PDB == "ATOM" and a.label_seq_id is not None:
            rk = (a.label_seq_id, a.pdbx_PDB_ins_code)
            seen = chain_res_comp[a.label_asym_id].setdefault(rk, {})
            representatives = chain_res_representatives[
                a.label_asym_id
            ].setdefault(rk, [])
            alt_id = a.label_alt_id if _has_altloc(a) else None
            candidates = (
                representatives
                if alt_id is None
                else (seen.get(None), seen.get(alt_id))
            )
            prior_comp_id = next(
                (
                    comp_id
                    for comp_id in candidates
                    if comp_id is not None and comp_id != a.label_comp_id
                ),
                None,
            )
            seen.setdefault(alt_id, a.label_comp_id)
            if (
                len(representatives) < 2
                and a.label_comp_id not in representatives
            ):
                representatives.append(a.label_comp_id)
            if (
                prior_comp_id is not None
                and (a.label_asym_id, rk) not in flagged
            ):
                flagged.add((a.label_asym_id, rk))
                diagnostics.errors.append(
                    Diagnostic(
                        code="RESIDUE_NUMBER_COLLISION",
                        severity="error",
                        message=(
                            f"Residue number collision in chain "
                            f"{a.label_asym_id} at seq_id "
                            f"{a.label_seq_id}: {prior_comp_id} vs "
                            f"{a.label_comp_id}"
                        ),
                        entry_id=entry_id,
                    )
                )
