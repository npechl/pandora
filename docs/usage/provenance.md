# Provenance

`pandora.provenance` aggregates whatever provenance you already have in
hand — it never fetches, re-derives, or validates anything itself — into
one per-structure `ProvenanceBundle` or one dataset-wide
`DatasetManifest`. `reproduce_dataset()` is the one function that
actually re-runs the pipeline, replaying a `DatasetManifest` from
scratch. See [Architecture](../architecture.md#provenance-and-the-dataset-manifest)
for how these models relate, and
[Functions](../reference/functions.md#pandora.provenance) for full
signatures.

```python
from pathlib import Path
from pandora.ingestion import fetch_mmcif
from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.metadata import collect_metadata
from pandora.annotations import annotate_structure_counts
from pandora.schemas.canonicalisation import canonicalisationPolicy
from pandora.schemas.dataset import DatasetCurationPolicy

output_dir = Path("./datasets/output/provenance/")
policy = canonicalisationPolicy(
    policy_id="p", policy_name="p", policy_version="1.0.0"
)
curation_policy = DatasetCurationPolicy(
    policy_id="c", policy_name="Default", policy_version="1.0.0"
)

# Fetch (not skipped here, unlike the other usage pages — reproduce_dataset()
# below needs real ingestion provenance to replay from).
structures, ingestion_prov, canon_prov = {}, {}, {}
for entry_id in ["104m", "112m"]:
    ing = fetch_mmcif(entry_id, "pdbe", None, output_dir / "raw")
    structure, _, _ = mmcif_to_structure(
        str(output_dir / "raw" / f"{entry_id}.cif")
    )
    canonical, _, prov = canonicalise_structure(structure, policy)
    structures[canonical.entry_id] = canonical
    ingestion_prov[canonical.entry_id] = ing
    canon_prov[canonical.entry_id] = prov
```

## Flatten a MetadataRecord's provenance

`collect_metadata_provenance()` pulls the per-field `MetadataProvenance`
stamps out of a `MetadataRecord` into a flat list — what
`build_provenance_bundle()` stores as `metadata_sources`.

```python
from pandora.provenance import collect_metadata_provenance

metadata = collect_metadata(structures["104M"])
stamps = collect_metadata_provenance(metadata)
print(len(stamps), stamps[0])
# 13 source='mmcif' source_category='_struct' source_record_id=None
```

## Bundle one structure's provenance

`build_provenance_bundle()` assembles a `ProvenanceBundle` from
whichever stage provenance you pass — every argument is optional.

```python
from pandora.provenance import build_provenance_bundle

counts = annotate_structure_counts(structures["104M"])
bundle = build_provenance_bundle(
    structures["104M"],
    ingestion=ingestion_prov["104M"],
    canonicalisation=canon_prov["104M"],
    metadata=metadata,
    annotations=[counts],
)
print(bundle.entry_id, bundle.ingestion is not None, len(bundle.annotations))
# 104M True 1
```

## Assemble a dataset manifest

`build_dataset_manifest()` aggregates dataset-wide provenance —
policies (stored by value, not just id/version, so a rebuild has the
actual rules), exclusions, dedup/clustering/partition provenance, the
split assignment — plus every retained structure's `ProvenanceBundle`
into one `DatasetManifest`.

```python
from pandora.provenance import build_dataset_manifest
from pandora.export import write_json

bundles = [
    build_provenance_bundle(
        structures[entry_id],
        ingestion=ingestion_prov[entry_id],
        canonicalisation=canon_prov[entry_id],
    )
    for entry_id in structures
]
manifest = build_dataset_manifest(
    dataset_id="demo",
    dataset_name="Demo",
    dataset_version="1.0.0",
    canonicalisation_policy=policy,
    curation_policy=curation_policy,
    structures=bundles,
)
manifest_path = write_json(manifest, output_dir / "manifest.json")
print(manifest.dataset_id, len(manifest.structures), "structures")
# demo 2 structures
```

This is exactly what the [`pandora manifest`](cli.md#manifest) CLI
subcommand builds — reach for the CLI when you're assembling a manifest
from already-written stage outputs on disk, and this function directly
when you're scripting the whole pipeline in one process, as above.

## Replay a manifest from scratch

`reproduce_dataset()` re-fetches every structure via its bundle's
`ingestion` provenance, then replays canonicalisation, curation, dedup,
similarity/clustering, partitioning, and annotation exactly as the
manifest recorded them.

```python
from pandora.provenance import reproduce_dataset

reproduced, new_manifest = reproduce_dataset(
    manifest, output_dir / "reproduced"
)
print(len(reproduced), "structures reproduced")
# 2 structures reproduced
```

It's a best-effort re-run, not a guaranteed byte-identical rebuild —
diff `new_manifest` against `manifest` to see what changed. Two hard
requirements, both raising `ValueError` if unmet: every structure's
bundle needs `ingestion` provenance (this is why the walkthrough above
fetches for real, unlike the other usage pages), and reproducing
`clustering` needs `ClusteringProvenance.similarity_method` set. See
[`pandora reproduce`](cli.md#reproduce) for the CLI equivalent.
