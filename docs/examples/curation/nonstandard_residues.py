from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import curate_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    ModifiedResidueRules,
    canonicalisationPolicy,
)
from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules


def curate(structure, quality_rules: QualityRules) -> str:
    policy = DatasetCurationPolicy(
        policy_id="c",
        policy_name="c",
        policy_version="1.0.0",
        quality_rules=quality_rules,
    )
    curated, exclusions, _ = curate_structure(structure, None, policy)
    return "kept" if curated else exclusions[-1].message


reject = QualityRules(exclude_nonstandard_residues=True)

# 1a08 contains ACE and DIP caps and FTY, a modified tyrosine.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1a08.cif")
print("1a08:", curate(structure, reject))
allow = reject.model_copy(
    update={"allowed_nonstandard_residues": ["ACE", "DIP", "FTY"]}
)
print("1a08, allow-list:", curate(structure, allow))

# 1b6w contains MSE. Mapping it to MET in canonicalisation makes it
# standard by the time curation runs.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1b6w.cif")
print("1b6w as deposited:", curate(structure, reject))
canonical, _, _ = canonicalise_structure(
    structure,
    canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
    ),
)
print("1b6w after MSE -> MET:", curate(canonical, reject))
