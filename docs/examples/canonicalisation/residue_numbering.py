from pandora.canonicalisation import canonicalise_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    IdentifierRules,
    ResidueNumberingRules,
    canonicalisationPolicy,
)

# 10mv: chain A starts at SEQRES position 22 and has a gap in the middle.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/10mv.cif")

for strategy in ["preserve", "use_auth_seq", "renumber"]:
    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        identifier_rules=IdentifierRules(
            residue_numbering=ResidueNumberingRules(strategy=strategy)
        ),
    )
    canonical, mappings, _ = canonicalise_structure(structure, policy)
    seq_ids = sorted(
        {a.label_seq_id for a in canonical.atoms if a.label_asym_id == "A"}
    )
    print(f"{strategy}: first {seq_ids[:3]}, last {seq_ids[-1]}")

item = mappings.residue_number_mapping.items[0]
print(
    f"renumber mapping: chain {item.original_chain_id} "
    f"{item.original_seq_id} (auth {item.original_auth_seq_id}) "
    f"-> {item.canonical_seq_id}"
)
