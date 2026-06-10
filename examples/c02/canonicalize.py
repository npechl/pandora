"""Canonicalize a parsed structure (C01 → C02)."""
import sys
sys.path.insert(0, "../..")

from pandora.c01_ingestion import ingest_mmcif
from pandora.c02_canonicalization import canonicalize_structure
from pandora.schemas.c01_ingestion import MmCIFIngestionInput
from pandora.schemas.c02_canonicalization import CanonicalizationPolicy

entry_id = sys.argv[1] if len(sys.argv) > 1 else "1cbs"

# C01 — ingest
ingestion_result = ingest_mmcif(MmCIFIngestionInput(
    entry_id=entry_id,
    provider="pdbe",
))
print(f"[C01] {entry_id} ingestion status: {ingestion_result.status}")

# C02 — canonicalize
policy = CanonicalizationPolicy(
    policy_id="canon-v1",
    policy_name="Default Canonicalization",
    policy_version="1.0.0",
)

result = canonicalize_structure(ingestion_result, policy)

print(f"[C02] canonicalization status: {result.status}")
if result.canonical_structure:
    cs = result.canonical_structure
    print(f"      chains:    {len(cs.chains)}")
    print(f"      residues:  {len(cs.residues)}")
    print(f"      atoms:     {len(cs.atoms)}")
print(f"      policy:    {result.applied_policy.policy_id} v{result.applied_policy.policy_version}")
print()
print("NOTE: C02 normalization steps are stubs — mappings are empty until implemented.")
