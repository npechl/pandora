from pandora.canonicalisation import canonicalise_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    IncompleteChainRules,
    MissingAtomsRules,
    MissingDataRules,
    canonicalisationPolicy,
)


def residues(structure, chain: str) -> int:
    return len(
        {a.label_seq_id for a in structure.atoms if a.label_asym_id == chain}
    )


def canonicalise(structure, rules: MissingDataRules):
    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        missing_data_rules=rules,
    )
    canonical, _, _ = canonicalise_structure(structure, policy)
    return canonical


# 1a08 chain B: two residues lack part of their backbone.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1a08.cif")
dropped = canonicalise(
    structure,
    MissingDataRules(
        missing_atoms=MissingAtomsRules(strategy="drop_partial_residue")
    ),
)
print(
    "drop_partial_residue, 1a08 chain B:",
    residues(structure, "B"),
    "->",
    residues(dropped, "B"),
    "residues",
)

# 10mv chain A has a gap; keep only its longest unbroken stretch.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/10mv.cif")
truncated = canonicalise(
    structure,
    MissingDataRules(
        incomplete_chains=IncompleteChainRules(
            strategy="truncate_to_complete_regions"
        )
    ),
)
print(
    "truncate_to_complete_regions, 10mv chain A:",
    residues(structure, "A"),
    "->",
    residues(truncated, "A"),
    "residues",
)
