"""Create leakage-safe train/validation/test splits (C01–C04 → C05)."""
import sys
sys.path.insert(0, "../..")

from pandora.c01_ingestion import ingest_mmcif
from pandora.c02_canonicalization import canonicalize_structure
from pandora.c03_metadata import attach_metadata, attach_plugins
from pandora.c04_curation import build_dataset
from pandora.c05_splitting import create_leakage_safe_dataset
from pandora.schemas.c01_ingestion import MmCIFIngestionInput
from pandora.schemas.c02_canonicalization import CanonicalizationPolicy
from pandora.schemas.c03_metadata import AnnotationPluginPolicy, MetadataIntegrationPolicy
from pandora.schemas.c04_curation import DatasetCurationPolicy
from pandora.schemas.c05_splitting import LeakagePolicy

entry_ids = sys.argv[1:] if len(sys.argv) > 1 else ["1cbs", "4hhb", "1tup", "2hho", "1a4w"]

canon_policy   = CanonicalizationPolicy(policy_id="canon-v1",   policy_name="Default Canonicalization", policy_version="1.0.0")
meta_policy    = MetadataIntegrationPolicy(policy_id="meta-v1", policy_name="Default Metadata",         policy_version="1.0.0")
plugin_policy  = AnnotationPluginPolicy(policy_id="plugin-v1",  policy_name="Default Plugins",          policy_version="1.0.0")
curation_policy = DatasetCurationPolicy(policy_id="cur-v1",     policy_name="Default Curation",         policy_version="1.0.0")
leakage_policy  = LeakagePolicy(policy_id="split-v1",           policy_name="Default Splitting",        policy_version="1.0.0")

# C01 → C02 → C03
annotated = []
for eid in entry_ids:
    ing = ingest_mmcif(MmCIFIngestionInput(entry_id=eid, provider="pdbe"))
    can = canonicalize_structure(ing, canon_policy)
    meta = attach_metadata(can, meta_policy)
    ann = attach_plugins(meta, plugin_policy, plugins=[])
    annotated.append(ann)
    print(f"  prepared {eid}")

# C04 — curated dataset
dataset = build_dataset(annotated, curation_policy,
    dataset_id="split-example", dataset_name="Split Example", dataset_version="1.0.0")

# C05 — leakage-safe split
split = create_leakage_safe_dataset(dataset, leakage_policy)

ps = split.partition_summary
print(f"\n[C05] Split: {split.dataset_id}")
print(f"      train:      {ps.train_count}  ({ps.train_fraction_achieved:.0%})")
print(f"      validation: {ps.validation_count}  ({ps.validation_fraction_achieved:.0%})")
print(f"      test:       {ps.test_count}  ({ps.test_fraction_achieved:.0%})")
print(f"      clusters:   {split.similarity_clusters.clustering_summary.total_clusters}")
print()
print("NOTE: MMseqs2/Foldseek similarity and cluster-aware assignment are stubs.")
