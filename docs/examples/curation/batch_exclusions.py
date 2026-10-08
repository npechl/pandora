from pandora.datasets import curate_structure
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules

policy = DatasetCurationPolicy(
    policy_id="c",
    policy_name="Sharp, complete chains",
    policy_version="1.0.0",
    quality_rules=QualityRules(
        max_resolution=2.5,
        max_missing_middle_fraction=0.1,
    ),
)

kept = {}
exclusions = []
for entry_id in ["104m", "1aui", "1a08", "1p58"]:
    structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
    curated, records, _ = curate_structure(
        structure, collect_metadata(structure), policy
    )
    exclusions.extend(records)  # chain records from kept entries too
    if curated is not None:
        kept[entry_id] = curated

entries_out = [e for e in exclusions if e.chain_id is None]
chains_out = [e for e in exclusions if e.chain_id is not None]
print("kept:", sorted(kept))
print("entries excluded:", [(e.entry_id, e.reason_code) for e in entries_out])
print(
    "chains removed:",
    [(e.entry_id, e.chain_id, e.reason_code) for e in chains_out],
)
