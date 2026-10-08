from pandora.datasets import curate_structure
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules

structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
metadata = collect_metadata(structure)

policy = DatasetCurationPolicy(
    policy_id="c1", policy_name="Default", policy_version="1.0.0"
)
curated, exclusions, provenance = curate_structure(structure, metadata, policy)
print("default:", curated is not None, exclusions)

strict = policy.model_copy(
    update={"quality_rules": QualityRules(max_resolution=1.0)}
)
curated, exclusions, provenance = curate_structure(structure, metadata, strict)
print("strict:", curated is not None)
for e in exclusions:
    print(f"  {e.reason_code} - {e.message}")
