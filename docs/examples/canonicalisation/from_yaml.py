from pandora.canonicalisation import canonicalise_structure
from pandora.ingestion.policy import load_policy
from pandora.parsing import mmcif_to_structure

policy = load_policy("datasets/canonicalisation.yaml")
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
canonical, mappings, provenance = canonicalise_structure(structure, policy)

print(policy.policy_id, policy.policy_version)
print(provenance.transforms)
