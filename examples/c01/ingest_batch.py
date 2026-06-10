"""Ingest a batch of PDB entries sequentially."""
import sys
sys.path.insert(0, "../..")

from pandora.c01_ingestion import ingest_list_mmcif
from pandora.schemas.c01_ingestion import MmCIFBatchInput, MmCIFIngestionInput

entry_ids = sys.argv[1:] if len(sys.argv) > 1 else ["1cbs", "4hhb", "1tup"]

batch = MmCIFBatchInput(entries=[
    MmCIFIngestionInput(entry_id=eid, provider="pdbe")
    for eid in entry_ids
])

result = ingest_list_mmcif(batch)

print(f"Total:   {result.summary.total}")
print(f"Success: {result.summary.success}")
print(f"Warning: {result.summary.warning}")
print(f"Failed:  {result.summary.failed}")
print()

for r in result.results:
    if r.parsed_structure:
        ps = r.parsed_structure
        cached = r.provenance.from_cache
        print(f"  {r.entry_id:6s}  {r.status:8s}  chains={len(ps.chains)}  atoms={len(ps.atoms):5d}  cached={cached}")
    else:
        msg = r.diagnostics.errors[0].message if r.diagnostics.errors else "unknown"
        print(f"  {r.entry_id:6s}  FAILED  {msg}")
