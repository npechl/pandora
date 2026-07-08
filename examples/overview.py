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
    AssemblyRules,
    EntityRules,
    LigandRules,
)

ENTRY_ID = "1cbs"
MMCIF_DIR = Path("./examples/mmcif")

ingestion_prov = fetch_mmcif(
    entry_id=ENTRY_ID,
    provider="pdbe",
    source_uri=None,
    output_dir=MMCIF_DIR,
)

print(f"Fetched : {ingestion_prov.source_uri}")
print(
    f"Provider: {ingestion_prov.provider}  cached={ingestion_prov.from_cache}"
)

structure, diag, status = mmcif_to_structure(str(MMCIF_DIR / f"{ENTRY_ID}.cif"))


# Example 1: preserve everything as reported -----------------------------

preserve_policy = CanonicalizationPolicy(
    policy_id="overview-preserve",
    policy_name="Preserve As Reported",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        chain_id=ChainIdRules(strategy="preserve"),
        residue_numbering=ResidueNumberingRules(strategy="preserve"),
    ),
    altloc_rules=AltlocRules(strategy="select_best_occupancy"),
    ligand_rules=LigandRules(strategy="preserve", keep_waters=False),
)

canonical, mappings, canon_prov = canonicalize_structure(
    structure, preserve_policy
)

# Example 2: remap identifiers and renumber residues ---------------------

remap_policy = CanonicalizationPolicy(
    policy_id="overview-remap",
    policy_name="Remap Chains And Renumber",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        chain_id=ChainIdRules(strategy="remap"),
        residue_numbering=ResidueNumberingRules(strategy="renumber"),
    ),
    altloc_rules=AltlocRules(
        strategy="select_best_occupancy", tie_breaker="lowest_b_factor"
    ),
    entity_rules=EntityRules(strategy="merge_equivalent_entities"),
    ligand_rules=LigandRules(
        strategy="filter", keep_waters=False, keep_ions=False
    ),
)
canonical, mappings, canon_prov = canonicalize_structure(
    structure, remap_policy
)


# Example 3: author-centric identifiers, first altloc, strict validation --

author_policy = CanonicalizationPolicy(
    policy_id="overview-author",
    policy_name="Author Identifiers",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        chain_id=ChainIdRules(strategy="use_auth_chain_id"),
        residue_numbering=ResidueNumberingRules(strategy="use_auth_seq"),
    ),
    altloc_rules=AltlocRules(strategy="select_first"),
    assembly_rules=AssemblyRules(
        strategy="standardize_biological_assembly",
        preferred_assembly_source="author",
    ),
    ligand_rules=LigandRules(strategy="annotate_only"),
)

canonical, mappings, canon_prov = canonicalize_structure(
    structure, author_policy
)
