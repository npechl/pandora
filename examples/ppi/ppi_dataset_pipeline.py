"""Protein-protein interaction (PPI) dataset pipeline, built entirely from
functionality Pandora already provides — no new logic, just chaining
existing stages. It approximates ATOM3D's PIP task
(https://www.atom3d.ai/pip.html) rather than reproducing it: the unit of
data is a chain-chain `InterfaceRecord` (a pair of interacting chains plus
the residues on each side that contact the other), not ATOM3D's exact
residue-pair binary-classification framing — Pandora doesn't have a
residue-pair-contact or negative-sampling primitive (see the gap notes at
the bottom of this file).

Pipeline: parse+canonicalise -> curate/dedup -> keep only structures with
>=1 chain-chain interface -> extract_interface_records() per structure ->
entry-level sequence-identity clustering (MMseqs2) -> leakage-safe
30%-identity split -> dataset manifest.
"""

from __future__ import annotations

from pathlib import Path

from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import (
    curate_structure,
    deduplicate_structures,
    entry_sequences,
    extract_interface_records,
)
from pandora.export import write_json, write_records
from pandora.parsing import mmcif_to_structure
from pandora.provenance import build_dataset_manifest, build_provenance_bundle
from pandora.schemas.canonicalisation import canonicalisationPolicy
from pandora.schemas.dataset import DatasetCurationPolicy, DeduplicationRules
from pandora.schemas.structure import Structure
from pandora.similarity import (
    cluster_similar_items,
    compute_sequence_similarity,
    partition_dataset,
)

MMCIF_DIR = Path("./datasets/dev/mmcif")
OUTPUT_DIR = Path("./datasets/output/ppi/")

DISTANCE_CUTOFF = 4.0  # extract_interface_records's own default
IDENTITY_THRESHOLD = 0.3  # ATOM3D PIP's split cutoff

policy = canonicalisationPolicy(
    policy_id="ppi-default", policy_name="Default", policy_version="1.0.0"
)
curation_policy = DatasetCurationPolicy(
    policy_id="ppi-curation-default",
    policy_name="Default",
    policy_version="1.0.0",
)
dedup_rules = DeduplicationRules(enabled=True)


# 1. Parse + canonicalise every entry, skipping ones that fail
#    validation instead of aborting the whole batch (canonicalise_structure
#    raises ValueError on a failed validation status).
structures: dict[str, Structure] = {}
canonicalisation_provenance = {}
for path in sorted(MMCIF_DIR.glob("*.cif")):
    structure, _, status = mmcif_to_structure(str(path))
    try:
        canonical, _, canon_prov = canonicalise_structure(structure, policy)
    except ValueError as exc:
        print(f"[{structure.entry_id}] skipped: {exc}")
        continue
    structures[canonical.entry_id] = canonical
    canonicalisation_provenance[canonical.entry_id] = canon_prov
print(f"parsed+canonicalised: {len(structures)} structures")

# 1b. Curate + deduplicate.
exclusions = []
for entry_id in list(structures):
    curated, exclusion, _ = curate_structure(
        structures[entry_id], None, curation_policy
    )
    if curated is None:
        exclusions.append(exclusion)
        del structures[entry_id]
    else:
        structures[entry_id] = curated
retained, removed, dedup_prov = deduplicate_structures(
    list(structures.values()), dedup_rules
)
exclusions.extend(removed)
structures = {s.entry_id: s for s in retained}
print(f"curation+dedup: {len(structures)} retained, {len(exclusions)} excluded")

# 2. Keep only structures with >=1 chain-chain interface, and collect
#    Pandora's own InterfaceRecord for each one found. This is the PPI
#    unit of data: a pair of chains + the residues on each side that
#    contact the other — not a residue-pair binary label, since Pandora
#    has no primitive for that (see gap notes below).
interfaces = []
heteromeric: dict[str, Structure] = {}
for entry_id, structure in structures.items():
    records = extract_interface_records(structure, DISTANCE_CUTOFF)
    if not records:
        continue
    heteromeric[entry_id] = structure
    interfaces.extend(records)
structures = heteromeric
print(
    f"interfaces: {len(structures)} heteromeric structures, "
    f"{len(interfaces)} chain-chain interfaces"
)

# 3. Entry-level sequence identity -> similarity network -> clusters.
sequences = entry_sequences(structures)
relationships = compute_sequence_similarity(sequences)
item_ids = sorted(sequences)  # same ids compute_sequence_similarity used
clusters, cluster_prov = cluster_similar_items(
    item_ids, relationships, IDENTITY_THRESHOLD
)
print(
    f"{len(clusters)} sequence-identity cluster(s) at {IDENTITY_THRESHOLD:.0%}"
)

# 4. Leakage-safe split: whole clusters move together, so no two
#    complexes above the identity threshold land in different splits.
splits, partition_prov = partition_dataset(
    clusters, pct_train=0.7, pct_val=0.15, pct_test=0.15
)
entry_split = {
    entry_id: split_name
    for split_name, entry_ids in splits.items()
    for entry_id in entry_ids
}
for split_name, entry_ids in splits.items():
    print(f"  {split_name}: {len(entry_ids)} structures")

# 5. Write outputs — one interfaces file per split (InterfaceRecord has
#    no split field of its own, so this is the split without touching
#    the schema).
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
for split_name in splits:
    split_interfaces = [
        record
        for record in interfaces
        if entry_split.get(record.entry_id) == split_name
    ]
    if split_interfaces:
        write_records(
            split_interfaces, OUTPUT_DIR / f"interfaces_{split_name}.json"
        )

manifest = build_dataset_manifest(
    dataset_id="ppi-dev-dataset",
    dataset_name="PPI interface dataset (dev fixtures)",
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
        for entry_id in structures
    ],
)
write_json(manifest, OUTPUT_DIR / "dataset_manifest.json")
print(f"\nwrote {len(interfaces)} interfaces + manifest -> {OUTPUT_DIR}")

# Gaps vs. an ATOM3D-PIP-style dataset, using only what exists today:
#
# 1. No "requires an interface" / minimum-chain-count curation rule.
#    QualityRules has min_chain_length (longest chain's residue count) but
#    nothing about polymer chain count or interaction presence — filtering
#    to heteromeric structures had to happen as ad hoc script logic (step 2
#    above) instead of a declarative, reusable curation policy rule.
#
# 2. annotate_chain_interfaces/extract_interface_records report each
#    side's contact residues as two separate marginal sets, not the actual
#    residue-to-residue pairing. ATOM3D's task is a per-*pair* binary
#    label ("does residue i on chain A contact residue j on chain B?"),
#    which isn't recoverable from an InterfaceRecord alone.
#
# 3. No non-contact ("negative") sampling primitive anywhere in the
#    library. Every existing stage is deterministic filtering/transform;
#    there's no stochastic-sampling step in the pipeline at all.
