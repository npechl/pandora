from pandora.datasets import curate_structure
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules

# 22jy: a cryo-EM structure at 2.2 A.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/22jy.cif")
metadata = collect_metadata(structure)
print(metadata.quality.experimental_method, metadata.quality.resolution)

rules = {
    "one limit, 3.5": QualityRules(max_resolution=3.5),
    "X-ray 2.5, cryo-EM 2.0": QualityRules(
        max_resolution_by_method={
            "X-RAY DIFFRACTION": 2.5,
            "ELECTRON MICROSCOPY": 2.0,
        }
    ),
}
for label, quality_rules in rules.items():
    policy = DatasetCurationPolicy(
        policy_id="c",
        policy_name="c",
        policy_version="1.0.0",
        quality_rules=quality_rules,
    )
    curated, exclusions, _ = curate_structure(structure, metadata, policy)
    print(f"{label}: kept={curated is not None}", end="")
    print("".join(f" ({e.reason_code}: {e.message})" for e in exclusions))
