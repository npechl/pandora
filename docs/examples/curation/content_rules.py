from pandora.datasets import curate_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import ContentRules, DatasetCurationPolicy

structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")

policy = DatasetCurationPolicy(
    policy_id="c",
    policy_name="Protein only",
    policy_version="1.0.0",
    content_rules=ContentRules(
        keep_ligands=False, keep_waters=False, keep_ions=False
    ),
)
curated, exclusions, _ = curate_structure(structure, None, policy)
print(len(structure.atoms), "->", len(curated.atoms), "atoms;", exclusions)
