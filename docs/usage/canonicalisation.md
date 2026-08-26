# Canonicalisation

`pandora.canonicalisation` normalizes a parsed `Structure` according to a
policy — chain IDs, residue numbering, assemblies, entities, missing
data, altlocs, ligands, then validation, always in that order. See
[Architecture](../architecture.md#canonicalisation-policy) for how the
policy's rule groups relate to each other, and
[Functions](../reference/functions.md#pandora.canonicalisation) for full
signatures.

## Canonicalise with the default policy

An empty `canonicalisationPolicy` isn't a no-op: several rule groups
default to something other than `"preserve"` (missing-atom/residue
flagging, best-occupancy altloc selection), so even the default policy
records transforms.

```python
from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.schemas.canonicalisation import canonicalisationPolicy

structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")

policy = canonicalisationPolicy(
    policy_id="p1", policy_name="Default", policy_version="1.0.0"
)
canonical, mappings, provenance = canonicalise_structure(structure, policy)

print(provenance.transforms)
# ['missing_atoms:annotate', 'missing_residues:annotate', 'altloc:select_best_occupancy']
```

## Canonicalise with a custom policy

Remap chain IDs, renumber residues, merge equivalent entities, and drop
waters/ions — the policy `examples/overview.py` and
`datasets/canonicalisation.yaml` both use:

```python
from pandora.schemas.canonicalisation import (
    canonicalisationPolicy,
    IdentifierRules,
    ChainIdRules,
    ResidueNumberingRules,
    AltlocRules,
    EntityRules,
    LigandRules,
)

policy = canonicalisationPolicy(
    policy_id="remap-1",
    policy_name="Remap",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        chain_id=ChainIdRules(strategy="remap"),
        residue_numbering=ResidueNumberingRules(strategy="renumber"),
    ),
    altloc_rules=AltlocRules(
        strategy="select_best_occupancy", tie_breaker="lowest_b_factor"
    ),
    entity_rules=EntityRules(strategy="merge_equivalent_entities"),
    ligand_rules=LigandRules(
        strategy="filter", keep_waters=False, keep_ions=False
    ),
)
canonical, mappings, provenance = canonicalise_structure(structure, policy)

print(provenance.transforms)
# ['chain_id:remap', 'residue_numbering:renumber', 'missing_atoms:annotate',
#  'missing_residues:annotate', 'altloc:select_best_occupancy',
#  'entity:merge_equivalent_entities', 'ligands:filter']
print(len(mappings.chain_id_mapping.items))
# 5 — one entry per original chain, mapping it to its new remapped id
```

You can load the same policy from YAML instead of building it in code —
see [`load_policy()`](ingestion.md#load-a-canonicalisation-policy).

## Filter ligands directly

`filter_ligands()` is what `ligand_rules` calls internally, but it's
also exported standalone — useful if you want ligand filtering without
running the rest of canonicalisation (this is exactly how
`pandora.datasets.curate_structure()`'s content rules reuse it).

```python
from pandora.canonicalisation import filter_ligands
from pandora.schemas.canonicalisation import LigandRules
from pandora.schemas.common import DiagnosticBundle

atoms, asym_units = filter_ligands(
    list(structure.atoms),
    list(structure.asym_units),
    structure.entities,
    LigandRules(strategy="filter", keep_waters=False, keep_ions=True),
    DiagnosticBundle(),
    structure.entry_id,
)
print(len(structure.atoms), "->", len(atoms))
# 1450 -> 1271
```
