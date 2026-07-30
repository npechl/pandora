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

from pandora.datasets.records import extract_chain_records, extract_residue_records

# from pandora.similarity.sequence import compute_sequence_similarity
# from pandora.similarity.clustering import cluster_similar_items
# from pandora.similarity.partition import partition_dataset

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
        incomplete_chains=IncompleteChainRules(strategy="truncate_to_complete_regions")
    ),
    altloc_rules=AltlocRules(strategy="select_best_occupancy"),
    assembly_rules=AssemblyRules(strategy="standardize_biological_assembly"),
    entity_rules=EntityRules(strategy="merge_equivalent_entities"),
    ligand_rules=LigandRules(strategy="annotate_only"),
)


structures = {}

for entry_id in ENTRY_IDS:
    structure, _, status = mmcif_to_structure(str(MMCIF_DIR / f"{entry_id}.cif"))
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



# # 2. All-vs-all sequence identity -> a SimilarityRelationship network -------
# print("\npairwise sequence identity:")
# relationships = []
# computed_at = datetime.now(timezone.utc).isoformat()


# for left_id, right_id in combinations(ENTRY_IDS, 2):
#     layer = annotate_pairwise_sequence_identity(
#         structures[left_id], structures[right_id]
#     )
#     score = layer.data["best_identity"] or 0.0
#     source_id, target_id = sorted((left_id, right_id))
#     relationships.append(
#         SimilarityRelationship(
#             source_id=source_id,
#             target_id=target_id,
#             similarity_type="sequence_similarity",
#             score=score,
#             identity=score,
#             method=SimilarityMethod(
#                 engine="pandora.basic.ungapped_entity_identity.v1"
#             ),
#             provenance=SimilarityRelationshipProvenance(
#                 computed_at=computed_at
#             ),
#         )
#     )
#     print(f"  {source_id} vs {target_id}: identity={score:.3f}")

# # 3. Connected-component clustering ----------------------------------------
# clusters = cluster_similar_items(ENTRY_IDS, relationships, IDENTITY_THRESHOLD)
# print(f"\n{len(clusters)} cluster(s) at threshold={IDENTITY_THRESHOLD}:")
# for cluster in clusters:
#     print(f"  {cluster.components}")

# # 4. Leakage-safe train/val/test split — whole clusters move together, so
# #    near-duplicate structures never end up split across partitions.
# splits = partition_dataset(clusters, pct_train=0.6, pct_val=0.2, pct_test=0.2)
# print("\nleakage-safe split:")
# for split_name, entry_ids in splits.items():
#     print(f"  {split_name}: {entry_ids}")

# # 5. Dataset records for the whole batch ------------------------------------
# all_chain_records = [
#     record
#     for entry_id in ENTRY_IDS
#     for record in extract_chain_records(structures[entry_id])
# ]
# print(
#     f"\n{len(all_chain_records)} chain records across {len(ENTRY_IDS)} entries"
# )
