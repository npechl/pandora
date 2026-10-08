from pandora.canonicalisation import canonicalise_structure, filter_ligands
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    LigandRules,
    canonicalisationPolicy,
)
from pandora.schemas.common import DiagnosticBundle

# 104m: myoglobin, a sulfate ion, a heme, n-butyl isocyanide and waters.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
names = {e.id: e.pdbx_description for e in structure.entities}

policy = canonicalisationPolicy(
    policy_id="p",
    policy_name="p",
    policy_version="1.0.0",
    ligand_rules=LigandRules(
        strategy="filter", keep_waters=False, keep_ions=False
    ),
)
canonical, _, _ = canonicalise_structure(structure, policy)
print("in the policy:", len(structure.atoms), "->", len(canonical.atoms))
kept = sorted({names[u.entity_id] for u in canonical.asym_units})
print("  kept:", kept)

# The same filter on its own, without the rest of canonicalisation.
atoms, asym_units = filter_ligands(
    list(structure.atoms),
    list(structure.asym_units),
    structure.entities,
    LigandRules(strategy="filter", keep_waters=False, keep_ions=True),
    DiagnosticBundle(),
    structure.entry_id,
)
print("filter_ligands, keep ions:", len(structure.atoms), "->", len(atoms))
