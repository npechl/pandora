from pandora.canonicalisation import canonicalise_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    AssemblyRules,
    canonicalisationPolicy,
)

# 13dg: assembly 1 is chains A-D under operators 1, 2 and 3.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/13dg.cif")
print("asymmetric unit:", [u.id for u in structure.asym_units])

policy = canonicalisationPolicy(
    policy_id="p",
    policy_name="p",
    policy_version="1.0.0",
    assembly_rules=AssemblyRules(strategy="standardize_biological_assembly"),
)
canonical, mappings, _ = canonicalise_structure(structure, policy)
print("after expansion:", [u.id for u in canonical.asym_units])
for item in mappings.assembly_mapping.items:
    for copy in item.chain_copies[:2]:
        print(
            f"  {copy.canonical_chain_id} = chain {copy.source_chain_id} "
            f"under operator {copy.operator_id}"
        )
