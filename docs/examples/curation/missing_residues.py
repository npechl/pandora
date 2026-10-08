from pandora.datasets import curate_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules

# 1aui chain A misses 95 residues in its middle; chain B is nearly whole.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1aui.cif")


def curate(quality_rules: QualityRules) -> None:
    policy = DatasetCurationPolicy(
        policy_id="c",
        policy_name="c",
        policy_version="1.0.0",
        quality_rules=quality_rules,
    )
    curated, exclusions, _ = curate_structure(structure, None, policy)
    for e in exclusions:
        print(f"  {e.chain_id or 'entry'}: {e.reason_code} - {e.message}")
    if curated is not None:
        chains = sorted(
            {a.label_asym_id for a in curated.atoms if a.label_seq_id}
        )
        print(f"  kept polymer chains: {chains}")


print("max 30% missing tails, 10% missing middle:")
curate(
    QualityRules(max_missing_tail_fraction=0.3, max_missing_middle_fraction=0.1)
)
print("the same, but the entry needs two chains:")
curate(
    QualityRules(
        max_missing_tail_fraction=0.3,
        max_missing_middle_fraction=0.1,
        min_polymer_chains=2,
    )
)
