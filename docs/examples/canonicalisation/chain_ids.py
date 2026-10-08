from pandora.canonicalisation import canonicalise_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    ChainIdRules,
    IdentifierRules,
    canonicalisationPolicy,
)


def chain_map(path: str, strategy: str) -> str:
    structure, _, _ = mmcif_to_structure(path)
    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        identifier_rules=IdentifierRules(
            chain_id=ChainIdRules(strategy=strategy)
        ),
    )
    _, mappings, _ = canonicalise_structure(structure, policy)
    return ", ".join(
        f"{item.original_chain_id}->{item.canonical_chain_id}"
        for item in mappings.chain_id_mapping.items
    )


# 13dg lists its chains as A, E, B, F, ...; remap relabels them in order.
print("remap:", chain_map("datasets/dev/mmcif/13dg.cif", "remap"))
# 1a08: label chain C is author chain B, and ligands share author ids.
print(
    "use_auth_chain_id:",
    chain_map("datasets/dev/mmcif/1a08.cif", "use_auth_chain_id"),
)
