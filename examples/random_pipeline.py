"""This is a random pipeline.

No biologial knowledge/meaning is extracted.
"""

from pandora.schemas.canonicalisation import (
    canonicalisationPolicy,
    IdentifierRules,
    IncompleteChainRules,
    MissingResiduesRules,
    MissingDataRules,
    MissingAtomsRules,
    AssemblyIdRules,
    ResidueNumberingRules,
    AltlocRules,
    AssemblyRules,
    LigandRules,
    EntityRules,
)

from pathlib import Path

from pandora.canonicalisation import canonicalise_structure
from pandora.parsing import mmcif_to_structure

from pandora.datasets.records import (
    extract_chain_records,
    extract_residue_records,
)

from os import listdir

MMCIF_DIR = Path("./datasets/dev/mmcif")
ENTRY_IDS = [file[0:-4] for file in listdir(MMCIF_DIR)]


policy = canonicalisationPolicy(
    policy_id="random-1",
    policy_name="Random example 1",
    policy_version="1.0.0",
    identifier_rules=IdentifierRules(
        residue_numbering=ResidueNumberingRules(strategy="renumber"),
        assembly_id=AssemblyIdRules(strategy="standardize"),
    ),
    missing_data_rules=MissingDataRules(
        missing_atoms=MissingAtomsRules(strategy="annotate"),
        missing_residues=MissingResiduesRules(strategy="drop_chain_segment"),
        incomplete_chains=IncompleteChainRules(
            strategy="truncate_to_complete_regions"
        ),
    ),
    altloc_rules=AltlocRules(strategy="select_best_occupancy"),
    assembly_rules=AssemblyRules(strategy="standardize_biological_assembly"),
    entity_rules=EntityRules(strategy="merge_equivalent_entities"),
    ligand_rules=LigandRules(strategy="annotate_only"),
)


structures = {}

for entry_id in ENTRY_IDS:
    structure, _, status = mmcif_to_structure(
        str(MMCIF_DIR / f"{entry_id}.cif")
    )
    canonical, _, _ = canonicalise_structure(structure, policy)

    structures[entry_id] = canonical

    print(
        f"[{structure.entry_id}] parsed+canonicalised: status={status} "
        f"chains={len(canonical.asym_units)}"
    )


residue_records = {}
chain_records = {}

for entry_id in ENTRY_IDS:
    residue_records[entry_id] = extract_residue_records(structures[entry_id])
    chain_records[entry_id] = extract_chain_records(structures[entry_id])
