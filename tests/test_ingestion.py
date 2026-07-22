from pathlib import Path

from pandora.ingestion.policy import load_policy

POLICY_PATH = Path(__file__).parent.parent / "examples" / "canonicalisation.yaml"


def test_load_policy():
    policy = load_policy(str(POLICY_PATH))

    assert policy.policy_id == "overview-remap"
    assert policy.identifier_rules.chain_id.strategy == "remap"
    assert policy.altloc_rules.strategy == "select_best_occupancy"
    assert policy.ligand_rules.keep_waters is False
