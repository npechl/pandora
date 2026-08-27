# PPI interface dataset

A narrated walkthrough of `examples/ppi_dataset_pipeline.py` — a protein-protein interaction (PPI) dataset built entirely from stages Pandora already provides, no new logic, just chaining existing ones together. It approximates [ATOM3D's PIP task](https://www.atom3d.ai/pip.html) rather than reproducing it exactly: the unit of data here is a chain-chain `InterfaceRecord` (a pair of interacting chains, plus the residues on each side that contact the other), not ATOM3D's residue-pair binary-classification framing.

## Prerequisites

Runs against the local fixtures in `datasets/dev/mmcif/`, so no network access is needed — but [`compute_sequence_similarity()`](../usage/similarity.md#sequence-similarity-mmseqs2) shells out to `mmseqs2`, so that binary needs to be on `PATH` (`conda install -c bioconda mmseqs2`; see [Installation](../getting-started/installation.md#external-binaries)).

## Run it

??? example "Full script"
    ```python linenums="1"
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
    ```

```text
parsed+canonicalised: 100 structures
curation+dedup: 100 retained, 0 excluded
interfaces: 45 heteromeric structures, 147 chain-chain interfaces
35 sequence-identity cluster(s) at 30%
  train: 31 structures
  val: 7 structures
  test: 7 structures

wrote 147 interfaces + manifest -> datasets/output/ppi
```

## Walking through the steps

### 1. Parse + canonicalise, skipping failures

`canonicalise_structure()` raises `ValueError` on a failed validation status — the script catches it per entry so one bad structure doesn't abort the whole batch:

```python linenums="1"
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
```

### 2. Curate + deduplicate

Both default policies here — an empty `DatasetCurationPolicy` and `DeduplicationRules(enabled=True)` — see [Curate one structure](../usage/datasets.md#curate-one-structure) and [Deduplicate a batch](../usage/datasets.md#deduplicate-a-batch) for what each rule group does and how to tighten them.

```python linenums="1"
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
```

### 3. Keep only heteromeric structures, collect interfaces

[`extract_interface_records()`](../usage/datasets.md#reshape-into-flat-records) reshapes [`annotate_chain_interfaces()`](../usage/annotation.md#chain-chain-interfaces)'s output; a structure with no chain-chain contacts (a monomer, or chains that never touch) contributes nothing and is dropped from the dataset entirely — this is the PPI-specific filter, and it's ad hoc script logic rather than a declarative curation rule (see gap 1 below):

```python linenums="1"
interfaces = []
heteromeric: dict[str, Structure] = {}
for entry_id, structure in structures.items():
    records = extract_interface_records(structure, DISTANCE_CUTOFF)
    if not records:
        continue
    heteromeric[entry_id] = structure
    interfaces.extend(records)
structures = heteromeric
```

### 4. Cluster and split, leakage-safe

Entry-level sequence identity (via `entry_sequences()` — one representative sequence per entry) feeds [`compute_sequence_similarity()`](../usage/similarity.md#sequence-similarity-mmseqs2), whose relationships feed [`cluster_similar_items()`](../usage/similarity.md#clustering) at ATOM3D PIP's own 30% identity threshold. [`partition_dataset()`](../usage/similarity.md#leakage-safe-partitioning) then moves whole clusters together, so no two complexes above that threshold end up split across train/val/test:

```python linenums="1"
sequences = entry_sequences(structures)
relationships = compute_sequence_similarity(sequences)
clusters, cluster_prov = cluster_similar_items(
    sorted(sequences), relationships, IDENTITY_THRESHOLD
)
splits, partition_prov = partition_dataset(
    clusters, pct_train=0.7, pct_val=0.15, pct_test=0.15
)
```

### 5. Build manifest & write outputs

`InterfaceRecord` has no split field of its own, so the script splits by filtering on `entry_split` instead of touching the schema, then writes one file per non-empty split plus a [dataset manifest](../usage/provenance.md#assemble-a-dataset-manifest) recording every policy, exclusion, and provenance record used along the way:

```python linenums="1"
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
```

```text
datasets/output/ppi/
├── interfaces_train.json
├── interfaces_val.json
├── interfaces_test.json
└── dataset_manifest.json
```

??? example "Manifest file"
    ```json linenums="1"
    {
    "dataset_id": "ppi-dev-dataset",
    "dataset_name": "PPI interface dataset (dev fixtures)",
    "dataset_version": "1.0.0",
    "pandora_version": "0.1.0",
    "generated_at": "2026-08-27T13:40:30.147529+00:00",
    "curation_policy": {
        "policy_id": "ppi-curation-default",
        "policy_name": "Default",
        "policy_version": "1.0.0",
        "description": "",
        "quality_rules": {
        "max_resolution": null,
        "null_resolution_behavior": "exclude",
        "min_chain_length": null,
        "include_experimental_methods": [],
        "exclude_experimental_methods": []
        },
        "organism_rules": {
        "include_taxa": [],
        "exclude_taxa": []
        },
        "content_rules": {
        "keep_ligands": true,
        "keep_waters": true,
        "keep_ions": true
        }
    },
    "canonicalisation_policy": {
        "policy_id": "ppi-default",
        "policy_name": "Default",
        "policy_version": "1.0.0",
        "description": "",
        "identifier_rules": {
        "chain_id": {
            "strategy": "preserve"
        },
        "residue_numbering": {
            "strategy": "preserve",
            "preserve_insertion_codes": true
        },
        "assembly_id": {
            "strategy": "preserve"
        }
        },
        "missing_data_rules": {
        "missing_atoms": {
            "strategy": "annotate",
            "allow_imputation": false,
            "record_missingness": true
        },
        "missing_residues": {
            "strategy": "annotate",
            "record_gaps": true
        },
        "incomplete_chains": {
            "strategy": "preserve"
        }
        },
        "altloc_rules": {
        "strategy": "select_best_occupancy",
        "tie_breaker": "alphabetical_first",
        "user_defined_altloc": null,
        "record_selection": true
        },
        "assembly_rules": {
        "strategy": "preserve_as_reported",
        "preferred_assembly_source": "author",
        "record_original_assembly_mapping": true
        },
        "entity_rules": {
        "strategy": "preserve",
        "preserve_original_entity_ids": true
        },
        "ligand_rules": {
        "strategy": "preserve",
        "keep_waters": true,
        "keep_ions": true,
        "keep_nonpolymer_ligands": true
        },
        "validation_rules": {
        "strictness": "moderate",
        "fail_on_unresolved_issues": true,
        "warnings_as_errors": false
        },
        "provenance_rules": {
        "record_original_mappings": true,
        "record_transforms": true,
        "record_policy_application": true,
        "emit_canonicalisation_report": false
        }
    },
    "excluded": [],
    "deduplication": {
        "deduplicated_at": "2026-08-27T13:40:28.348044+00:00",
        "enabled": true,
        "duplicates_found": 0
    },
    "clustering": {
        "clustered_at": "2026-08-27T13:40:30.146426+00:00",
        "threshold": 0.3,
        "n_relationships": 19,
        "n_clusters": 35,
        "similarity_method": {
        "engine": "MMseqs2",
        "version": "18.8cc5c",
        "parameters": {
            "sensitivity": 5.7,
            "mmseqs_bin": "mmseqs"
        }
        }
    },
    "partition": {
        "partitioned_at": "2026-08-27T13:40:30.146469+00:00",
        "pct_train": 0.7,
        "pct_val": 0.15,
        "pct_test": 0.15,
        "keep_similar_items": true,
        "split_sizes": {
        "train": 31,
        "val": 7,
        "test": 7
        }
    },
    "splits": {
        "train": [
        "13SD",
        "13SJ",
        "13SK",
        "13SR",
        "13SV",
        "1A2C",
        "1A3E",
        "1ABI",
        "1AVG",
        "1A08",
        "1A1B",
        "1A3N",
        "1ABW",
        "1AL4",
        "1AV2",
        "10TM",
        "11KB",
        "12HH",
        "12UM",
        "13DG",
        "179D",
        "1A02",
        "1A03",
        "1A7F",
        "1AC9",
        "1AUI",
        "1AX7",
        "1BU1",
        "1C9O",
        "1D1K",
        "1P58"
        ],
        "val": [
        "1ALL",
        "1AUU",
        "1BKV",
        "1BVI",
        "1CA7",
        "1I0T",
        "1TRJ"
        ],
        "test": [
        "1AO7",
        "1AUY",
        "1BMS",
        "1BYZ",
        "1CE0",
        "1LS2",
        "22JY"
        ]
    },
    "structures": [
        {
        "entry_id": "10TM",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147366+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.229025+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "11KB",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147378+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.303726+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "12HH",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147383+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.322423+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "12UM",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147386+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.333476+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "13DG",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147390+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.371535+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "13SD",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147393+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.699101+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "13SJ",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147396+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.725444+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "13SK",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147399+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.749490+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "13SR",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147402+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.774415+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "13SV",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147408+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.800050+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "179D",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147412+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.814653+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1A02",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147415+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.839286+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1A03",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147417+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.925100+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1A08",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147420+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.980795+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1A1B",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147423+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:26.994746+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1A2C",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147429+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.013854+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1A3E",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147432+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.031997+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1A3N",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147437+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.067539+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1A7F",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147440+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.199496+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1ABI",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147445+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.220438+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1ABW",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147448+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.255889+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1AC9",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147452+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.262685+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1AL4",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147454+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.344104+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1ALL",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147457+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.362344+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1AO7",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147460+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.398772+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1AUI",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147463+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.495855+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1AUU",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147466+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.529242+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1AUY",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147469+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.559004+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1AV2",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147472+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.633076+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1AVG",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147475+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.656718+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1AX7",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147478+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.758544+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1BKV",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147480+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.865407+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1BMS",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147483+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.894716+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1BU1",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147487+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.937810+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1BVI",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147490+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.963131+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1BYZ",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147495+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:27.984738+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1C9O",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147498+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:28.001886+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1CA7",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147501+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:28.019708+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1CE0",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147504+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:28.027855+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1D1K",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147507+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:28.082688+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1I0T",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147510+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:28.121683+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1LS2",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147513+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:28.291597+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1P58",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147519+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:28.320349+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "1TRJ",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147523+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:28.333389+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        },
        {
        "entry_id": "22JY",
        "pandora_version": "0.1.0",
        "generated_at": "2026-08-27T13:40:30.147524+00:00",
        "ingestion": null,
        "canonicalisation": {
            "canonicalised_at": "2026-08-27T13:40:28.347841+00:00",
            "policy_id": "ppi-default",
            "policy_name": "Default",
            "policy_version": "1.0.0",
            "transforms": [
            "missing_atoms:annotate",
            "missing_residues:annotate",
            "altloc:select_best_occupancy"
            ],
            "report": {}
        },
        "metadata_sources": [],
        "annotations": []
        }
    ]
    }
    ```