"""Estimate an entry's atom count before parsing it into a `Structure`."""

from __future__ import annotations

from pathlib import Path

import gemmi

from pandora.canonicalisation.assemblies import _select_biological_assembly
from pandora.schemas.canonicalisation import canonicalisationPolicy
from pandora.schemas.common import DiagnosticBundle


def estimate_atom_count(
    path: str | Path, policy: canonicalisationPolicy
) -> int:
    """Estimate how many atoms canonicalising an mmCIF file would produce.

    Reads the file with gemmi only (no pydantic records), so it is cheap
    enough to run before `mmcif_to_structure` on very large entries.
    Counts model 1. When the policy expands a biological assembly, each
    generator's chains are counted once per operator. Waters, ligands
    and every altloc are included, so the result is an upper bound on
    the canonicalised structure's atoms.

    Args:
        path: The mmCIF file (plain or gzipped).
        policy: The canonicalisation policy whose assembly rules apply.

    Returns:
        The estimated atom count.

    Raises:
        RuntimeError: If gemmi cannot read the file.
        ValueError: If gemmi cannot read the file.
    """

    st = gemmi.read_structure(str(path), format=gemmi.CoorFormat.Mmcif)
    if len(st) == 0:
        return 0
    per_subchain: dict[str, int] = {}
    for chain in st[0]:
        for sub in chain.subchains():
            sid = sub.subchain_id()  # a subchain can come in several runs
            per_subchain[sid] = per_subchain.get(sid, 0) + sum(
                len(res) for res in sub
            )
    total = sum(per_subchain.values())

    rules = policy.assembly_rules
    if rules.strategy != "standardize_biological_assembly":
        return total
    # gemmi assemblies carry the author/software flags the selector reads.
    selected = _select_biological_assembly(
        list(st.assemblies),
        rules.preferred_assembly_source,
        DiagnosticBundle(),
        st.name,
    )
    if selected is None:
        return total
    return sum(
        len(gen.operators) * sum(per_subchain.get(s, 0) for s in gen.subchains)
        for gen in selected.generators
    )
