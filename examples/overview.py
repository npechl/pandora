"""End-to-end walkthrough of every implemented Pandora component, run
against two local entries from datasets/dev/mmcif/:

    1ayi — apo monomer (no ligand)
    104m — ligand-bound monomer (myoglobin + heme)

Ingestion (pandora.ingestion.fetch_mmcif) is skipped: it only fetches
from PDBe/PDB over HTTP, and these files are already on disk.
"""

from pathlib import Path

from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.metadata import collect_metadata
from pandora.annotations import (
    annotate_ligand_contacts,
    annotate_pairwise_sequence_identity,
    annotate_structure_counts,
)
from pandora.schemas.canonicalisation import (
    canonicalisationPolicy,
    IdentifierRules,
    ChainIdRules,
    ResidueNumberingRules,
    AltlocRules,
    EntityRules,
    LigandRules,
)

MMCIF_DIR = Path("./datasets/dev/mmcif")
ENTRY_IDS = ["1ayi", "104m"]

policy = canonicalisationPolicy(
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

canonical_structures = {}

for entry_id in ENTRY_IDS:
    # 1. Parsing — raw mmCIF -> Pandora's typed Structure.
    structure, diagnostics, status = mmcif_to_structure(
        str(MMCIF_DIR / f"{entry_id}.cif")
    )
    print(f"[{entry_id}] parsed: status={status} atoms={len(structure.atoms)}")

    # 2. Canonicalisation — apply the policy (chain IDs, altlocs, ligands, ...).
    canonical, mappings, canon_prov = canonicalise_structure(structure, policy)
    canonical_structures[entry_id] = canonical
    print(f"[{entry_id}] canonicalised: transforms={canon_prov.transforms}")

    # 3. Metadata — source-backed entry/entity/quality/taxonomy records.
    metadata = collect_metadata(canonical)
    print(f"[{entry_id}] metadata: title={metadata.entry.title!r}")

    # 4. Annotations — derived, per-entry structural summaries.
    counts = annotate_structure_counts(canonical)
    print(f"[{entry_id}] structure_counts: {counts.data}")

    contacts = annotate_ligand_contacts(canonical, distance_cutoff=4.0)
    n_ligands = len(contacts.data["ligands"])
    print(f"[{entry_id}] ligand_contacts: {n_ligands} ligand(s)")

# 4b. Pairwise annotation — compares entities across two canonical structures.
identity = annotate_pairwise_sequence_identity(
    canonical_structures["1ayi"], canonical_structures["104m"]
)
best_identity = identity.data["best_identity"]
print(f"\n1ayi vs 104m best sequence identity: {best_identity:.3f}")
