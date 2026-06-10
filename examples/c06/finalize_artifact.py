"""Assemble provenance and finalize a Pandora artifact (C01–C05 → C06)."""
import sys
sys.path.insert(0, "../..")

from pandora.c01_ingestion import ingest_mmcif
from pandora.c02_canonicalization import canonicalize_structure
from pandora.c03_metadata import attach_metadata, attach_plugins
from pandora.c04_curation import build_dataset
from pandora.c05_splitting import create_leakage_safe_dataset
from pandora.c06_provenance import finalize_artifact
from pandora.schemas.c01_ingestion import MmCIFIngestionInput
from pandora.schemas.c02_canonicalization import CanonicalizationPolicy
from pandora.schemas.c03_metadata import AnnotationPluginPolicy, MetadataIntegrationPolicy
from pandora.schemas.c04_curation import DatasetCurationPolicy
from pandora.schemas.c05_splitting import LeakagePolicy
from pandora.schemas.c06_provenance import ProvenancePolicy

entry_ids = sys.argv[1:] if len(sys.argv) > 1 else ["1cbs", "4hhb", "1tup", "2hho", "1a4w"]

canon_policy    = CanonicalizationPolicy(policy_id="canon-v1",  policy_name="Default Canonicalization", policy_version="1.0.0")
meta_policy     = MetadataIntegrationPolicy(policy_id="meta-v1",policy_name="Default Metadata",         policy_version="1.0.0")
plugin_policy   = AnnotationPluginPolicy(policy_id="plugin-v1", policy_name="Default Plugins",          policy_version="1.0.0")
curation_policy = DatasetCurationPolicy(policy_id="cur-v1",     policy_name="Default Curation",         policy_version="1.0.0")
leakage_policy  = LeakagePolicy(policy_id="split-v1",           policy_name="Default Splitting",        policy_version="1.0.0")
prov_policy     = ProvenancePolicy(policy_id="prov-v1",         policy_name="Default Provenance",       policy_version="1.0.0")

# C01 → C02 → C03
annotated = []
for eid in entry_ids:
    ing  = ingest_mmcif(MmCIFIngestionInput(entry_id=eid, provider="pdbe"))
    can  = canonicalize_structure(ing, canon_policy)
    meta = attach_metadata(can, meta_policy)
    ann  = attach_plugins(meta, plugin_policy, plugins=[])
    annotated.append(ann)
    print(f"  prepared {eid}")

# C04
dataset = build_dataset(annotated, curation_policy,
    dataset_id="artifact-example", dataset_name="Artifact Example", dataset_version="1.0.0")

# C05
split = create_leakage_safe_dataset(dataset, leakage_policy)

# C06
artifact = finalize_artifact(
    artifact_id="pandora-artifact-001",
    leakage_safe_dataset=split,
    provenance_policy=prov_policy,
    artifact_name="Example Artifact",
)

print(f"\n[C06] Artifact:  {artifact.artifact_id}")
print(f"      name:      {artifact.artifact_name}")
print(f"      manifest:  {artifact.manifest.manifest_id}")
print(f"      version:   {artifact.manifest.pandora_version}")
print(f"      generated: {artifact.provenance.generated_at}")
print(f"      train:     {artifact.manifest.dataset_summary.train_count}")
print(f"      val:       {artifact.manifest.dataset_summary.validation_count}")
print(f"      test:      {artifact.manifest.dataset_summary.test_count}")
print()
print("NOTE: SHA-256 checksums, lineage export, and full provenance traversal are stubs.")
