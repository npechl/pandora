from pandora.datasets import curate_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules


def curate(entry_id: str, quality_rules: QualityRules) -> None:
    structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
    policy = DatasetCurationPolicy(
        policy_id="c",
        policy_name="c",
        policy_version="1.0.0",
        quality_rules=quality_rules,
    )
    curated, exclusions, _ = curate_structure(structure, None, policy)
    print(f"{entry_id}: kept={curated is not None}")
    for e in exclusions:
        print(f"  {e.chain_id or 'entry'}: {e.reason_code} - {e.message}")


# Whole entries: 1aui has 4,832 atoms.
curate("1aui", QualityRules(max_atoms=4000))
# Single chains: 1aui chain A has 378 observed residues, chain B 165.
curate("1aui", QualityRules(max_chain_length=300))
