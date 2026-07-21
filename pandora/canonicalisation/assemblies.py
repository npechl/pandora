from pandora.schemas.structure import AssemblyRecord
from pandora.schemas.canonicalisation import (
    AssemblyMapping,
    AssemblyMappingItem,
)


def _normalize_assemblies(
    assemblies: list[AssemblyRecord],
    assembly_rules,
    id_strategy: str,
    record: bool,
) -> tuple[list[AssemblyRecord], AssemblyMapping]:
    mapping = AssemblyMapping()
    result = list(assemblies)

    if assembly_rules.strategy == "select_first_assembly" and result:
        result = [result[0]]

    # standardize_biological_assembly: without external metadata we keep
    # order as-is; the preferred assembly (if annotated) would be
    # identified from pdbx flags in raw data.

    if id_strategy == "preserve":
        if record:
            for asm in result:
                mapping.items.append(
                    AssemblyMappingItem(
                        canonical_assembly_id=asm.id,
                        original_assembly_id=asm.id,
                    )
                )
        return result, mapping

    # remap or standardize → sequential integers
    new_result = []
    for i, asm in enumerate(result, 1):
        new_id = str(i)
        if record:
            mapping.items.append(
                AssemblyMappingItem(
                    canonical_assembly_id=new_id,
                    original_assembly_id=asm.id,
                )
            )
        new_result.append(asm.model_copy(update={"id": new_id}))

    return new_result, mapping
