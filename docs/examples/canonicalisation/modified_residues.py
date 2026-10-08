from pandora.canonicalisation import canonicalise_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    ModifiedResidueRules,
    canonicalisationPolicy,
)


def mapped(path: str, rules: ModifiedResidueRules) -> None:
    structure, _, _ = mmcif_to_structure(path)
    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        modified_residue_rules=rules,
    )
    _, mappings, provenance = canonicalise_structure(structure, policy)
    print(f"{structure.entry_id} {provenance.transforms[0]}")
    for item in mappings.modified_residue_mapping.items:
        print(
            f"  chain {item.chain_id} {item.seq_id}: "
            f"{item.original_comp_id} -> {item.parent_comp_id}"
        )


# Selenomethionine, mapped by default once the strategy is on.
mapped(
    "datasets/dev/mmcif/1b6w.cif",
    ModifiedResidueRules(strategy="map_to_parent"),
)
# Any residue named in comp_ids whose parent the file lists.
mapped(
    "datasets/dev/mmcif/1a08.cif",
    ModifiedResidueRules(strategy="map_to_parent", comp_ids=["FTY"]),
)
