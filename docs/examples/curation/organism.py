from pandora.datasets import curate_structure
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import DatasetCurationPolicy, OrganismRules

policy = DatasetCurationPolicy(
    policy_id="c",
    policy_name="Human only",
    policy_version="1.0.0",
    organism_rules=OrganismRules(include_taxa=["9606"]),
)
for entry_id in ["1aui", "1ayi"]:
    structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
    metadata = collect_metadata(structure)
    curated, exclusions, _ = curate_structure(structure, metadata, policy)
    organisms = sorted({t.organism_scientific for t in metadata.taxonomies})
    outcome = "kept" if curated else exclusions[-1].reason_code
    print(f"{entry_id} {organisms}: {outcome}")
