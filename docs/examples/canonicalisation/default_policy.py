from pandora.canonicalisation import canonicalise_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import canonicalisationPolicy

structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")

policy = canonicalisationPolicy(
    policy_id="p1", policy_name="Default", policy_version="1.0.0"
)
canonical, mappings, provenance = canonicalise_structure(structure, policy)

print(provenance.transforms)
