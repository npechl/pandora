"""Fetch and parse a single mmCIF entry from PDBe."""
import sys
sys.path.insert(0, "../..")

from pandora.c01_ingestion import ingest_mmcif
from pandora.schemas.c01_ingestion import MmCIFIngestionInput

entry_id = sys.argv[1] if len(sys.argv) > 1 else "1cbs"

result = ingest_mmcif(MmCIFIngestionInput(
    entry_id=entry_id,
    provider="pdbe",
))

print(f"Entry:      {result.entry_id}")
print(f"Status:     {result.status}")
print(f"From cache: {result.provenance.from_cache}")

if result.parsed_structure:
    ps = result.parsed_structure
    print(f"Chains:     {len(ps.chains)}")
    print(f"Residues:   {len(ps.residues)}")
    print(f"Atoms:      {len(ps.atoms)}")
    print(f"Entities:   {len(ps.entities)}")
    print(f"Assemblies: {len(ps.assemblies)}")
    print(f"Ligands:    {len(ps.ligands)}")
    print()
    for ch in ps.chains:
        print(f"  chain {ch.chain_id!r:4s}  type={ch.chain_type:12s}  entity={ch.entity_id!r}")
    print()
    for ent in ps.entities:
        seq_len = len(ent.sequence) if ent.sequence else 0
        print(f"  entity {ent.entity_id!r}  type={ent.entity_type:12s}  seq_len={seq_len}")
else:
    print("Errors:")
    for err in result.diagnostics.errors:
        print(f"  [{err.code}] {err.message}")
