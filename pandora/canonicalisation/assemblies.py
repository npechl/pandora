from pandora.schemas.structure import AssemblyRecord, AsymRecord, AtomSiteRecord
from pandora.schemas.canonicalisation import (
    AssemblyChainCopy,
    AssemblyMapping,
    AssemblyMappingItem,
)
from pandora.schemas.common import Diagnostic, DiagnosticBundle

_IDENTITY_MATRIX = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
_TOLERANCE = 1e-6


def _is_identity_transform(
    matrix: list[list[float]] | None, vector: list[float] | None
) -> bool:
    """Whether a symmetry operator is the identity (no-op) transform."""

    if matrix is None or vector is None:
        return True
    for row, ident_row in zip(matrix, _IDENTITY_MATRIX):
        for value, ident_value in zip(row, ident_row):
            if abs(value - ident_value) > _TOLERANCE:
                return False
    for value in vector:
        if abs(value) > _TOLERANCE:
            return False
    return True


def _apply_transform(
    matrix: list[list[float]], vector: list[float], x: float, y: float, z: float
) -> tuple[float, float, float]:
    """Apply a 3x3 rotation + translation to one coordinate."""

    return (
        matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + vector[0],
        matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + vector[1],
        matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + vector[2],
    )


def _select_biological_assembly(
    assemblies: list[AssemblyRecord],
    preferred_assembly_source: str,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> AssemblyRecord | None:
    """Pick the assembly to standardize to, per `preferred_assembly_source`."""

    if not assemblies:
        return None

    if preferred_assembly_source == "pdbe":
        diagnostics.warnings.append(
            Diagnostic(
                code="ASSEMBLY_SOURCE_UNAVAILABLE",
                severity="warning",
                message=(
                    "preferred_assembly_source='pdbe' has no external "
                    "metadata source implemented; falling back to "
                    "author/software-determined, then first assembly."
                ),
                entry_id=entry_id,
            )
        )

    for asm in assemblies:
        if asm.author_determined:
            return asm
    for asm in assemblies:
        if asm.software_determined:
            return asm
    return assemblies[0]


def _expand_biological_assembly(
    selected: AssemblyRecord,
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    record: bool,
) -> tuple[
    AssemblyRecord,
    list[AtomSiteRecord],
    list[AsymRecord],
    list[AssemblyChainCopy],
]:
    """Materialize `selected`'s symmetry operators into concrete chains."""

    atoms_by_chain: dict[str, list[AtomSiteRecord]] = {}
    for atom in atoms:
        atoms_by_chain.setdefault(atom.label_asym_id, []).append(atom)
    asym_by_id = {asym.id: asym for asym in asym_units}
    operators_by_id = {op.id: op for op in selected.operators}

    new_atoms: list[AtomSiteRecord] = []
    new_asym_units: list[AsymRecord] = []
    existing_chain_ids = {asym.id for asym in asym_units}
    seen_identity_chains: set[str] = set()
    chain_copies: list[AssemblyChainCopy] = []

    for gen in selected.generators:
        oper_ids = [op_id for op_id in gen.oper_expression.split(",") if op_id]
        for source_chain in gen.asym_id_list:
            if source_chain not in atoms_by_chain:
                continue
            for oper_id in oper_ids:
                op = operators_by_id.get(oper_id)
                matrix = op.matrix if op else None
                vector = op.vector if op else None

                if _is_identity_transform(matrix, vector):
                    if source_chain in seen_identity_chains:
                        continue
                    seen_identity_chains.add(source_chain)
                    continue

                new_id = f"{source_chain}_{oper_id}"
                while new_id in existing_chain_ids:
                    new_id = f"{new_id}_"
                existing_chain_ids.add(new_id)

                for atom in atoms_by_chain[source_chain]:
                    x, y, z = _apply_transform(
                        matrix, vector, atom.Cartn_x, atom.Cartn_y, atom.Cartn_z
                    )
                    new_atoms.append(
                        atom.model_copy(
                            update={
                                "label_asym_id": new_id,
                                "auth_asym_id": new_id,
                                "Cartn_x": x,
                                "Cartn_y": y,
                                "Cartn_z": z,
                            }
                        )
                    )

                source_asym = asym_by_id.get(source_chain)
                new_asym_units.append(
                    AsymRecord(
                        id=new_id,
                        entity_id=source_asym.entity_id if source_asym else "",
                        auth_id=None,
                    )
                )

                if record:
                    chain_copies.append(
                        AssemblyChainCopy(
                            canonical_chain_id=new_id,
                            source_chain_id=source_chain,
                            operator_id=oper_id,
                        )
                    )

    all_atoms = atoms + new_atoms
    all_asym_units = asym_units + new_asym_units
    new_oligomeric_count = (
        selected.oligomeric_count + len(chain_copies)
        if selected.oligomeric_count is not None
        else None
    )
    materialized = selected.model_copy(
        update={
            "generators": [],
            "operators": [],
            "oligomeric_count": new_oligomeric_count,
        }
    )
    return materialized, all_atoms, all_asym_units, chain_copies


def _normalize_assemblies(
    assemblies: list[AssemblyRecord],
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    assembly_rules,
    id_strategy: str,
    record: bool,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> tuple[
    list[AssemblyRecord],
    list[AtomSiteRecord],
    list[AsymRecord],
    AssemblyMapping,
]:
    """Apply the assembly selection/standardization/id strategy, per
    assembly_rules."""

    mapping = AssemblyMapping()
    result = list(assemblies)
    chain_copies_by_assembly_id: dict[str, list[AssemblyChainCopy]] = {}

    if assembly_rules.strategy == "select_first_assembly" and result:
        result = [result[0]]

    elif assembly_rules.strategy == "standardize_biological_assembly":
        selected = _select_biological_assembly(
            result,
            assembly_rules.preferred_assembly_source,
            diagnostics,
            entry_id,
        )
        if selected is not None:
            materialized, atoms, asym_units, chain_copies = (
                _expand_biological_assembly(selected, atoms, asym_units, record)
            )
            result = [materialized]
            chain_copies_by_assembly_id[materialized.id] = chain_copies

    if id_strategy == "preserve":
        if record:
            for asm in result:
                mapping.items.append(
                    AssemblyMappingItem(
                        canonical_assembly_id=asm.id,
                        original_assembly_id=asm.id,
                        chain_copies=chain_copies_by_assembly_id.get(
                            asm.id, []
                        ),
                    )
                )
        return result, atoms, asym_units, mapping

    # remap or standardize → sequential integers
    new_result = []
    for i, asm in enumerate(result, 1):
        new_id = str(i)
        if record:
            mapping.items.append(
                AssemblyMappingItem(
                    canonical_assembly_id=new_id,
                    original_assembly_id=asm.id,
                    chain_copies=chain_copies_by_assembly_id.get(asm.id, []),
                )
            )
        new_result.append(asm.model_copy(update={"id": new_id}))

    return new_result, atoms, asym_units, mapping
