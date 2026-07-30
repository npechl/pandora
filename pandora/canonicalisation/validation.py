from __future__ import annotations

from collections import defaultdict

from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
)
from pandora.schemas.common import Diagnostic, DiagnosticBundle


def _validate(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    rules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> str:
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

    # A collision is two *different* residues (distinct comp_id) sharing
    # a (seq_id, ins_code) slot in the same chain — not merely a residue
    # having more than one atom, which is the normal case.
    chain_res_comp: dict[str, dict[tuple, str]] = defaultdict(dict)
    flagged: set[tuple[str, tuple]] = set()
    for a in atoms:
        if a.group_PDB == "ATOM" and a.label_seq_id is not None:
            rk = (a.label_seq_id, a.pdbx_PDB_ins_code)
            seen = chain_res_comp[a.label_asym_id]
            prior_comp_id = seen.get(rk)
            if prior_comp_id is None:
                seen[rk] = a.label_comp_id
            elif (
                prior_comp_id != a.label_comp_id
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
