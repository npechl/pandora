from pandora.canonicalisation import canonicalise_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    AltlocRules,
    canonicalisationPolicy,
)

# 1b6w: residue 35 has two conformers, A and B, at occupancy 0.5 each.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1b6w.cif")

rules = [
    AltlocRules(strategy="preserve"),
    AltlocRules(
        strategy="select_best_occupancy", tie_breaker="lowest_b_factor"
    ),
    AltlocRules(strategy="select_user_defined", user_defined_altloc="B"),
]
for altloc_rules in rules:
    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        altloc_rules=altloc_rules,
    )
    canonical, mappings, _ = canonicalise_structure(structure, policy)
    atoms = [a for a in canonical.atoms if a.label_seq_id == 35]
    print(f"{altloc_rules.strategy}: {len(atoms)} atoms in residue 35")
    for item in mappings.altloc_selection_mapping.items:
        print(
            f"  {item.residue_id}: kept {item.selected_altloc} of "
            f"{item.available_altlocs} ({item.selection_reason})"
        )
