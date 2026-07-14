from __future__ import annotations

from collections import defaultdict

from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
)

from pandora.schemas.common import Diagnostic, DiagnosticBundle

_BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})


def _handle_missing_atoms(
    atoms: list[AtomSiteRecord],
    rules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> list[AtomSiteRecord]:
    strategy = rules.strategy
    if strategy in ("preserve", "impute"):
        return atoms

    polymer_by_residue: dict[tuple, list[AtomSiteRecord]] = defaultdict(list)
    hetatm: list[AtomSiteRecord] = []
    for a in atoms:
        if a.group_PDB == "ATOM":
            polymer_by_residue[
                (a.label_asym_id, a.label_seq_id, a.label_comp_id)
            ].append(a)
        else:
            hetatm.append(a)

    result: list[AtomSiteRecord] = []
    for key, res_atoms in polymer_by_residue.items():
        present = {a.label_atom_id for a in res_atoms}
        missing = _BACKBONE_ATOMS - present
        if missing:
            asym_id, seq_id, comp_id = key
            if rules.record_missingness:
                diagnostics.warnings.append(
                    Diagnostic(
                        code="MISSING_ATOMS",
                        severity="warning",
                        message=(
                            f"Residue {comp_id} {seq_id} in chain "
                            f"{asym_id} missing backbone atoms"
                        ),
                        entry_id=entry_id,
                        context={
                            "chain": asym_id,
                            "seq_id": seq_id,
                            "comp_id": comp_id,
                            "missing": sorted(missing),
                        },
                    )
                )
            if strategy == "drop_partial_residue":
                continue
        result.extend(res_atoms)

    result.extend(hetatm)
    return result


def _handle_missing_residues(
    atoms: list[AtomSiteRecord],
    rules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> list[AtomSiteRecord]:
    strategy = rules.strategy
    record = rules.record_gaps

    if strategy == "preserve" and not record:
        return atoms

    chain_seq_ids: dict[str, list[int]] = defaultdict(list)
    chain_seen: dict[str, set] = defaultdict(set)
    for a in atoms:
        if a.group_PDB == "ATOM" and a.label_seq_id is not None:
            if a.label_seq_id not in chain_seen[a.label_asym_id]:
                chain_seen[a.label_asym_id].add(a.label_seq_id)
                chain_seq_ids[a.label_asym_id].append(a.label_seq_id)

    gap_chains: set[str] = set()
    for chain, seq_ids in chain_seq_ids.items():
        sorted_ids = sorted(seq_ids)
        for i in range(len(sorted_ids) - 1):
            if sorted_ids[i + 1] - sorted_ids[i] > 1:
                gap_chains.add(chain)
                if record:
                    diagnostics.warnings.append(
                        Diagnostic(
                            code="SEQUENCE_GAP",
                            severity="warning",
                            message=(
                                f"Sequence gap in chain {chain} between "
                                f"residues {sorted_ids[i]} and "
                                f"{sorted_ids[i + 1]}"
                            ),
                            entry_id=entry_id,
                            context={
                                "chain": chain,
                                "gap_start": sorted_ids[i],
                                "gap_end": sorted_ids[i + 1],
                            },
                        )
                    )

    if strategy == "drop_chain_segment" and gap_chains:
        return [a for a in atoms if a.label_asym_id not in gap_chains]

    return atoms


def _handle_incomplete_chains(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    rules,
) -> tuple[list[AtomSiteRecord], list[AsymRecord]]:
    strategy = rules.strategy
    if strategy == "preserve":
        return atoms, asym_units

    chain_seq_ids: dict[str, list[int]] = defaultdict(list)
    chain_seen: dict[str, set] = defaultdict(set)
    for a in atoms:
        if a.group_PDB == "ATOM" and a.label_seq_id is not None:
            if a.label_seq_id not in chain_seen[a.label_asym_id]:
                chain_seen[a.label_asym_id].add(a.label_seq_id)
                chain_seq_ids[a.label_asym_id].append(a.label_seq_id)

    incomplete: set[str] = set()
    for chain, seq_ids in chain_seq_ids.items():
        s = sorted(seq_ids)
        if s[-1] - s[0] + 1 > len(s):
            incomplete.add(chain)

    if strategy == "exclude":
        keep = {a.id for a in asym_units} - incomplete
        return (
            [a for a in atoms if a.label_asym_id in keep],
            [a for a in asym_units if a.id in keep],
        )

    # truncate_to_complete_regions: keep longest contiguous run per
    # incomplete chain
    keep_ids: dict[str, set[int]] = {}
    for chain in incomplete:
        s = sorted(chain_seq_ids[chain])
        best: list[int] = []
        run: list[int] = [s[0]]
        for i in range(1, len(s)):
            if s[i] == s[i - 1] + 1:
                run.append(s[i])
            else:
                if len(run) > len(best):
                    best = run
                run = [s[i]]
        if len(run) > len(best):
            best = run
        keep_ids[chain] = set(best)

    result = []
    for a in atoms:
        if a.label_asym_id in incomplete:
            if a.label_seq_id in keep_ids.get(a.label_asym_id, set()):
                result.append(a)
        else:
            result.append(a)

    return result, asym_units
