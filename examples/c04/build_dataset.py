"""Build a curated dataset from a list of structures (C01 → C02 → C03 → C04)."""
import sys
sys.path.insert(0, "../..")

from pandora.c01_ingestion import ingest_mmcif
from pandora.c02_canonicalization import canonicalize_structure
from pandora.c03_metadata import attach_metadata, attach_plugins
from pandora.c04_curation import build_dataset, extract_chains, extract_residues
from pandora.schemas.c01_ingestion import MmCIFIngestionInput
from pandora.schemas.c02_canonicalization import CanonicalizationPolicy
from pandora.schemas.c03_metadata import AnnotationPluginPolicy, MetadataIntegrationPolicy
from pandora.schemas.c04_curation import DatasetCurationPolicy

entry_ids = sys.argv[1:] if len(sys.argv) > 1 else ["1cbs", "4hhb", "1tup"]

canon_policy = CanonicalizationPolicy(
    policy_id="canon-v1", policy_name="Default Canonicalization", policy_version="1.0.0",
)
meta_policy = MetadataIntegrationPolicy(
    policy_id="meta-v1", policy_name="Default Metadata", policy_version="1.0.0",
)
plugin_policy = AnnotationPluginPolicy(
    policy_id="plugin-v1", policy_name="Default Plugins", policy_version="1.0.0",
)
curation_policy = DatasetCurationPolicy(
    policy_id="curation-v1", policy_name="Default Curation", policy_version="1.0.0",
)

# C01 → C02 → C03 for each entry
annotated_structures = []
for eid in entry_ids:
    ing = ingest_mmcif(MmCIFIngestionInput(entry_id=eid, provider="pdbe"))
    can = canonicalize_structure(ing, canon_policy)
    meta = attach_metadata(can, meta_policy)
    ann = attach_plugins(meta, plugin_policy, plugins=[])
    annotated_structures.append(ann)
    print(f"  prepared {eid}: ingestion={ing.status}")

# C04 — build structure-level dataset
dataset = build_dataset(
    annotated_structures,
    curation_policy,
    dataset_id="example-dataset",
    dataset_name="Example Dataset",
    dataset_version="1.0.0",
)
print(f"\n[C04] Dataset:    {dataset.dataset_id}")
print(f"      structures: {dataset.counts.total_selected}")

# Extract chains
chain_dataset = extract_chains(dataset, curation_policy)
print(f"      chains:     {chain_dataset.counts.total_chains_extracted} (stub — implement extraction)")

# Extract residues from chain dataset
residue_dataset = extract_residues(chain_dataset, curation_policy)
print(f"      residues:   {residue_dataset.counts.total_residues_extracted} (stub — implement extraction)")
print()
print("NOTE: filtering, deduplication, and extraction rules are stubs.")
