"""Attach metadata and annotation to a canonical structure (C01 → C02 → C03)."""
import sys
sys.path.insert(0, "../..")

from pandora.c01_ingestion import ingest_mmcif
from pandora.c02_canonicalization import canonicalize_structure
from pandora.c03_metadata import attach_metadata
from pandora.schemas.c01_ingestion import MmCIFIngestionInput
from pandora.schemas.c02_canonicalization import CanonicalizationPolicy
from pandora.schemas.c03_metadata import MetadataIntegrationPolicy

entry_id = sys.argv[1] if len(sys.argv) > 1 else "1cbs"

# C01
ingestion = ingest_mmcif(MmCIFIngestionInput(entry_id=entry_id, provider="pdbe"))
print(f"[C01] ingestion:       {ingestion.status}")

# C02
canon_policy = CanonicalizationPolicy(
    policy_id="canon-v1", policy_name="Default Canonicalization", policy_version="1.0.0",
)
canonical = canonicalize_structure(ingestion, canon_policy)
print(f"[C02] canonicalization: {canonical.status}")

# C03
meta_policy = MetadataIntegrationPolicy(
    policy_id="meta-v1", policy_name="Default Metadata", policy_version="1.0.0",
)
annotated = attach_metadata(canonical, meta_policy)

print(f"[C03] metadata attached")
print(f"      entry_id:  {annotated.canonical_structure_result.entry_id}")
print(f"      retrieved: {annotated.provenance.retrieved_at}")
print()
print("NOTE: C03 API calls (PDBe, SIFTS, UniProt) are stubs — metadata fields are empty until implemented.")
