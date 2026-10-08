from pathlib import Path

import yaml

from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import curate_structure
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    ModifiedResidueRules,
    canonicalisationPolicy,
)
from pandora.schemas.dataset import DatasetCurationPolicy

# Map MSE to MET first, so the non-standard residue rule doesn't drop
# every selenomethionine crystal.
canonical_policy = canonicalisationPolicy(
    policy_id="mse",
    policy_name="Map MSE",
    policy_version="1.0.0",
    modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
)

entries = ["104m", "1a08", "1aui", "1b6w", "10mv", "10tm", "22jy", "1p58"]
structures = {}
for entry_id in entries:
    structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
    canonical, _, _ = canonicalise_structure(structure, canonical_policy)
    structures[entry_id] = (canonical, collect_metadata(canonical))

for path in sorted(Path("docs/examples/recipes").glob("*.yaml")):
    policy = DatasetCurationPolicy.model_validate(
        yaml.safe_load(path.read_text())
    )
    print(policy.policy_name)
    for entry_id, (structure, metadata) in structures.items():
        curated, exclusions, _ = curate_structure(structure, metadata, policy)
        line = f"  {entry_id}: {'kept' if curated else 'excluded'}"
        if exclusions:
            line += " - " + "; ".join(
                f"{e.chain_id or 'entry'} {e.reason_code}" for e in exclusions
            )
        print(line)
