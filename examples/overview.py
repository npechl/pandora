from pathlib import Path

from pandora.ingestion import fetch_mmcif
from pandora.parsing import mmcif_to_structure
from pandora.canonicalization import canonicalize_structure
from pandora.schemas.canonicalization import (
    CanonicalizationPolicy,
    IdentifierRules,
    ChainIdRules,
    ResidueNumberingRules,
    AltlocRules,
    LigandRules,
)

ENTRY_ID = "1cbs"
MMCIF_DIR = Path("./examples/mmcif")
MMCIF_PATH = MMCIF_DIR / f"{ENTRY_ID}.cif"

ingestion_prov = fetch_mmcif(
    entry_id=ENTRY_ID,
    provider="pdbe",
    source_uri=None,
    output_dir=MMCIF_DIR,
)

print(f"Fetched : {ingestion_prov.source_uri}")
print(f"Written : {MMCIF_PATH}")
print(
    f"Provider: {ingestion_prov.provider}  cached={ingestion_prov.from_cache}"
)

structure, diag, status = mmcif_to_structure(str(MMCIF_PATH))

if structure is None:
    print(f"Parse failed ({status}):")
    for e in diag.errors:
        print(f"  ERROR {e.code}: {e.message}")
    raise SystemExit(1)

print(f"\nParsed  : {structure.entry_id}  (status={status})")
print(f"  chains   : {len(structure.asym_units)}")
print(f"  entities : {len(structure.entities)}")
print(f"  atoms    : {len(structure.atoms)}")
print(f"  assemblies: {len(structure.assemblies)}")

if diag.warnings:
    for w in diag.warnings:
        print(f"  WARN {w.code}: {w.message}")


policy = CanonicalizationPolicy(
    policy_id="overview-v1",
    policy_name="Overview Example Policy",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        chain_id=ChainIdRules(strategy="preserve"),
        residue_numbering=ResidueNumberingRules(strategy="preserve"),
    ),
    altloc_rules=AltlocRules(strategy="select_best_occupancy"),
    ligand_rules=LigandRules(strategy="preserve", keep_waters=False),
)

canonical, mappings, canon_prov = canonicalize_structure(structure, policy)

print(f"\nCanonical: {canonical.entry_id}  (policy={canon_prov.policy_id})")
print(f"  atoms    : {len(canonical.atoms)}")
print(f"  transforms: {canon_prov.transforms or ['none']}")

if mappings.altloc_selection_mapping.items:
    print(
        f"  altloc selections: {len(mappings.altloc_selection_mapping.items)}"
    )
    for item in mappings.altloc_selection_mapping.items[:3]:
        print(
            f"    chain={item.canonical_chain_id} res={item.residue_id} "
            f"→ altloc={item.selected_altloc} "
            f"({item.selection_reason})"
        )
