from pandora.canonicalisation import canonicalise_structure
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    ChainIdRules,
    IdentifierRules,
    ModifiedResidueRules,
    ResidueNumberingRules,
    canonicalisationPolicy,
    canonicalisationProvenanceRules,
)

structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1b6w.cif")

policy = canonicalisationPolicy(
    policy_id="p",
    policy_name="p",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        chain_id=ChainIdRules(strategy="use_auth_chain_id"),
        residue_numbering=ResidueNumberingRules(strategy="renumber"),
    ),
    modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
    provenance_rules=canonicalisationProvenanceRules(
        emit_canonicalisation_report=True
    ),
)
canonical, mappings, provenance = canonicalise_structure(structure, policy)

# Which rule groups changed something, in the order they ran.
print("transforms:", provenance.transforms)
# How many diagnostics the run raised (the diagnostics themselves are
# not returned yet).
print("report:", provenance.report)
# One mapping list per kind of change, to trace any id back.
for name, mapping in mappings:
    print(f"{name}: {len(mapping.items)} items")
