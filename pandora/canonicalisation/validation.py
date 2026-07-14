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

    chain_res_keys: dict[str, set] = defaultdict(set)
    for a in atoms:
        if a.group_PDB == "ATOM" and a.label_seq_id is not None:
            rk = (a.label_seq_id, a.pdbx_PDB_ins_code)
            if rk in chain_res_keys[a.label_asym_id]:
                diagnostics.errors.append(
                    Diagnostic(
                        code="RESIDUE_NUMBER_COLLISION",
                        severity="error",
                        message=(
                            f"Residue number collision in chain "
                            f"{a.label_asym_id} at seq_id "
                            f"{a.label_seq_id}"
                        ),
                        entry_id=entry_id,
                    )
                )
            chain_res_keys[a.label_asym_id].add(rk)

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
