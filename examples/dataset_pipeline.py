"""Multi-structure dataset workflow: batch canonicalisation, an
all-vs-all sequence-similarity network, connected-component clustering,
and a leakage-safe train/val/test split.

pandora.similarity.sequence.compute_sequence_similarity() is the "real"
way to build the similarity network (MMseqs2 easy-search), but that
shells out to an external binary this repo doesn't bundle. This swaps
in annotate_pairwise_sequence_identity() instead — a simple ungapped,
position-by-position comparison with no alignment step. That's enough
to catch 118l/138l and 1a08/1a1b (same reading frame throughout), but
it's genuinely fooled by 104m/112m: one has an uncleaved N-terminal
Met the other lacks, so every position after residue 1 is off by one
and the ungapped score collapses even though the two are otherwise
near-identical — a real limitation of this annotation, not a bug in
this example. 

Everything downstream (cluster_similar_items(),
partition_dataset()) only cares about the resulting
SimilarityRelationship objects, not how they were computed, so swapping
in a real MMseqs2 network later (which aligns first) is a one-function
change and would cluster 104m/112m too.
"""

from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

from pandora.annotations import annotate_pairwise_sequence_identity
from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import extract_chain_records
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import canonicalisationPolicy
from pandora.schemas.similarity import (
    SimilarityMethod,
    SimilarityRelationship,
    SimilarityRelationshipProvenance,
)
from pandora.similarity import cluster_similar_items, partition_dataset

from os import listdir

MMCIF_DIR = Path("./datasets/dev/mmcif")
ENTRY_IDS = [file[0:-4] for file in listdir(MMCIF_DIR)]
# ENTRY_IDS = ["104m", "112m", "118l", "138l", "134l", "1ayi", "1a08", "1a1b"]
IDENTITY_THRESHOLD = 0.9

# No identifier remapping needed for this workflow — the default
# "preserve" policy is enough to resolve altlocs consistently before
# extracting per-chain sequences.
policy = canonicalisationPolicy(
    policy_id="dataset-default",
    policy_name="Default",
    policy_version="1.0.0",
)

# 1. Parse + canonicalise every entry --------------------------------------

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

# 2. All-vs-all sequence identity -> a SimilarityRelationship network -------
print("\npairwise sequence identity:")
relationships = []
computed_at = datetime.now(timezone.utc).isoformat()


for left_id, right_id in combinations(ENTRY_IDS, 2):
    layer = annotate_pairwise_sequence_identity(
        structures[left_id], structures[right_id]
    )
    score = layer.data["best_identity"] or 0.0
    source_id, target_id = sorted((left_id, right_id))
    relationships.append(
        SimilarityRelationship(
            source_id=source_id,
            target_id=target_id,
            similarity_type="sequence_similarity",
            score=score,
            identity=score,
            method=SimilarityMethod(
                engine="pandora.basic.ungapped_entity_identity.v1"
            ),
            provenance=SimilarityRelationshipProvenance(
                computed_at=computed_at
            ),
        )
    )
    print(f"  {source_id} vs {target_id}: identity={score:.3f}")

# 3. Connected-component clustering ----------------------------------------
clusters = cluster_similar_items(ENTRY_IDS, relationships, IDENTITY_THRESHOLD)
print(f"\n{len(clusters)} cluster(s) at threshold={IDENTITY_THRESHOLD}:")
for cluster in clusters:
    print(f"  {cluster.components}")

# 4. Leakage-safe train/val/test split — whole clusters move together, so
#    near-duplicate structures never end up split across partitions.
splits = partition_dataset(clusters, pct_train=0.5, pct_val=0.25, pct_test=0.25)
print("\nleakage-safe split:")
for split_name, entry_ids in splits.items():
    print(f"  {split_name}: {entry_ids}")

# 5. Dataset records for the whole batch ------------------------------------
all_chain_records = [
    record
    for entry_id in ENTRY_IDS
    for record in extract_chain_records(structures[entry_id])
]
print(
    f"\n{len(all_chain_records)} chain records across {len(ENTRY_IDS)} entries"
)
