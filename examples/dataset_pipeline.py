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
from pandora.datasets import (
    curate_structure,
    deduplicate_structures,
    extract_chain_records,
)
from pandora.export import write_json
from pandora.parsing import mmcif_to_structure
from pandora.provenance import build_dataset_manifest, build_provenance_bundle
from pandora.schemas.canonicalisation import canonicalisationPolicy
from pandora.schemas.dataset import DatasetCurationPolicy, DeduplicationRules
from pandora.schemas.similarity import (
    SimilarityMethod,
    SimilarityRelationship,
    SimilarityRelationshipProvenance,
)
from pandora.similarity import cluster_similar_items, partition_dataset

from os import listdir

MMCIF_DIR = Path("./datasets/dev/mmcif")
OUTPUT_DIR = Path("./datasets/output/dev/")
FILE_IDS = [file[0:-4] for file in listdir(MMCIF_DIR)]
# FILE_IDS = ["104m", "112m", "118l", "138l", "134l", "1ayi", "1a08", "1a1b"]
IDENTITY_THRESHOLD = 0.9

# No identifier remapping needed for this workflow — the default
# "preserve" policy is enough to resolve altlocs consistently before
# extracting per-chain sequences.
policy = canonicalisationPolicy(
    policy_id="dataset-default",
    policy_name="Default",
    policy_version="1.0.0",
)
curation_policy = DatasetCurationPolicy(
    policy_id="curation-default",
    policy_name="Default",
    policy_version="1.0.0",
)
dedup_rules = DeduplicationRules(enabled=True)

# 1. Parse + canonicalise every entry --------------------------------------

structures = {}
canonicalisation_provenance = {}
for file_id in FILE_IDS:
    structure, _, status = mmcif_to_structure(str(MMCIF_DIR / f"{file_id}.cif"))
    canonical, _, canon_prov = canonicalise_structure(structure, policy)
    structures[canonical.entry_id] = canonical
    canonicalisation_provenance[canonical.entry_id] = canon_prov
    print(
        f"[{structure.entry_id}] parsed+canonicalised: status={status} "
        f"chains={len(canonical.asym_units)}"
    )
ENTRY_IDS = list(structures)

# 1b. Curate + deduplicate ---------------------------------------------------
exclusions = []
for entry_id in list(ENTRY_IDS):
    curated, exclusion, _ = curate_structure(
        structures[entry_id], None, curation_policy
    )
    if curated is None:
        exclusions.append(exclusion)
        del structures[entry_id]
    else:
        structures[entry_id] = curated
ENTRY_IDS = list(structures)

retained, removed, dedup_prov = deduplicate_structures(
    list(structures.values()), dedup_rules
)
exclusions.extend(removed)
structures = {s.entry_id: s for s in retained}
ENTRY_IDS = list(structures)
print(
    f"\ncuration+dedup: {len(ENTRY_IDS)} retained, {len(exclusions)} excluded"
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
clusters, cluster_prov = cluster_similar_items(
    ENTRY_IDS, relationships, IDENTITY_THRESHOLD
)
print(f"\n{len(clusters)} cluster(s) at threshold={IDENTITY_THRESHOLD}:")
for cluster in clusters:
    print(f"  {cluster.components}")

# 4. Leakage-safe train/val/test split — whole clusters move together, so
#    near-duplicate structures never end up split across partitions.
splits, partition_prov = partition_dataset(
    clusters, pct_train=0.6, pct_val=0.2, pct_test=0.2
)
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

# 6. Dataset manifest — every step's provenance in one shareable file -------
manifest = build_dataset_manifest(
    dataset_id="dev-dataset",
    dataset_name="Development fixture dataset",
    dataset_version="1.0.0",
    curation_policy=curation_policy,
    canonicalisation_policy=policy,
    excluded=exclusions,
    deduplication=dedup_prov,
    clustering=cluster_prov,
    partition=partition_prov,
    splits=splits,
    structures=[
        build_provenance_bundle(
            structures[entry_id],
            canonicalisation=canonicalisation_provenance[entry_id],
        )
        for entry_id in ENTRY_IDS
    ],
)
write_json(manifest, OUTPUT_DIR / "dataset_manifest.json")
print(f"\ndataset manifest written to {OUTPUT_DIR / 'dataset_manifest.json'}")
