# Architecture

Pandora (see [Home](../index.md)) is a toolkit of plain functions and typed
models that you combine as your dataset needs. The *data* those functions
pass around is what actually defines the library's shape. This page has two
parts: a table mapping each component's functions to the models they take and return, then entity-relationship
diagrams showing how those models reference each other. The diagrams are
generated straight from the models' type hints via
[erdantic](https://erdantic.drivendata.org/), so they can never drift from
the code. Regenerate them after changing a model:

```sh
uv run --extra docs python docs/scripts/generate_erd.py
```

For per-field descriptions, see [Schemas](../reference/schemas.md); for full
signatures and docstrings, see [Functions](../reference/functions.md).

## Functions and the models they move

Every function below is pure: it takes a `Structure` (or another typed
model) in and returns a new one, never mutating its input. Types are
narrowed for readability — optional/keyword arguments and provenance-only
parameters are omitted; follow a function's link in
[Functions](../reference/functions.md) for the exact signature. Every
`pandora <subcommand>` in the CLI wraps exactly one function below,
reading/writing the same types as files instead of in-memory objects — see
[Library vs. CLI](../usage/overview.md#library-vs-cli).

### Ingestion

| Function | Takes | Returns |
|---|---|---|
| `fetch_mmcif()` | `entry_id`, `provider`, `output_dir` | `IngestionProvenance` — writes the mmCIF file as a side effect |
| `fetch_list_mmcif()` | a list of `entry_id`s | `list[IngestionProvenance]` — calls `fetch_mmcif()` once per id |
| `ingest_local_mmcif()` | `path` to a file already on disk | `IngestionProvenance` — describes it without fetching/copying |
| `load_policy()` | `path` to a policy YAML | `canonicalisationPolicy` |

### Parsing

| Function | Takes | Returns |
|---|---|---|
| `mmcif_to_structure()` | `path_to_mmcif` | `(Structure | None, DiagnosticBundle, ResultStatus)` — the framework's central data model |

### Canonicalisation

| Function | Takes | Returns |
|---|---|---|
| `canonicalise_structure()` | `Structure`, [`canonicalisationPolicy`](#canonicalisation-policy) | `(Structure, CanonicalMappings, canonicalisationProvenance)` |
| `filter_ligands()` | `atoms`, `asym_units`, `entities`, `LigandRules` | `(atoms, asym_units)` filtered — one of `canonicalise_structure()`'s nine internal rule steps, also callable on its own |

### Metadata

| Function | Takes | Returns |
|---|---|---|
| `collect_metadata()` | `Structure` | `MetadataRecord` — entry/quality/taxonomy/entity/ligand/UniProt records, all source-backed |

`collect_metadata()` is a thin wrapper over the `extract_*()` functions in
`pandora.metadata.mmcif` (one per record type); call those directly if you
only need one piece. See [Metadata](../usage/metadata.md).

### Annotations

| Function | Takes | Returns |
|---|---|---|
| `annotate_structure_counts()` | `Structure` | `AnnotationLayer` (`type="structure_counts"`) |
| `annotate_ligand_contacts()` | `Structure` | `AnnotationLayer` (`type="ligand_contacts"`) |
| `annotate_chain_interfaces()` | `Structure` | `AnnotationLayer` (`type="chain_interfaces"`) |
| `annotate_pairwise_sequence_identity()` | two `Structure`s | `AnnotationLayer` (`type="pairwise_sequence_identity"`) |

Every `annotate_*()` function returns the same `AnnotationLayer` envelope
regardless of what it computed — `layer.data` is where the shape actually
differs; see [Annotation](../usage/annotation.md) for each one's `data` schema.

### Datasets

| Function | Takes | Returns |
|---|---|---|
| `curate_structure()` | `Structure`, `MetadataRecord | None`, `DatasetCurationPolicy` | `(Structure | None, ExclusionRecord | None, CurationProvenance)` — `None` structure means excluded |
| `deduplicate_structures()` | `list[Structure]`, `DeduplicationRules` | `(list[Structure], list[ExclusionRecord], DeduplicationProvenance)` |
| `extract_chain_records()` | `Structure` | `list[ChainRecord]` |
| `extract_residue_records()` | `Structure` | `list[ResidueRecord]` |
| `extract_interface_records()` | `Structure` | `list[InterfaceRecord]` — reshapes `annotate_chain_interfaces()`'s output |
| `entry_sequences()` | `dict[str, Structure]` | `dict[str, str]` — one representative sequence per entry, feeds `compute_sequence_similarity()` |

### Similarity

| Function | Takes | Returns |
|---|---|---|
| `compute_sequence_similarity()` | sequences (`dict[str, str]`, `list[ChainRecord]`, or a FASTA directory) | `list[SimilarityRelationship]` — via MMseqs2 |
| `compute_structure_similarity()` | structure file paths | `list[SimilarityRelationship]` — via Foldseek |
| `cluster_similar_items()` | `item_ids`, `list[SimilarityRelationship]`, `threshold` | `(list[SimilarityCluster], ClusteringProvenance)` |
| `partition_dataset()` | `list[SimilarityCluster]` | `(dict[str, list[str]], PartitionProvenance)` — the `train`/`val`/`test` split, whole clusters per side |

### Provenance

| Function | Takes | Returns |
|---|---|---|
| `build_provenance_bundle()` | one `Structure` + whatever ingestion/canonicalisation/metadata/annotation provenance you already have | `ProvenanceBundle` |
| `build_dataset_manifest()` | policies + every stage's provenance record above, plus a `ProvenanceBundle` per retained structure | `DatasetManifest` |
| `reproduce_dataset()` | a `DatasetManifest` | `(dict[str, Structure], DatasetManifest)` — replays every recorded step from just the manifest |

### Export

| Function | Takes | Returns |
|---|---|---|
| `structure_to_mmcif()` | `Structure` | an mmCIF file |
| `export_chain_mmcif()` | `Structure`, `chain_id` | an mmCIF file with only that chain's atoms/connections |
| `write_json()` | any pydantic `BaseModel` (e.g. a `DatasetManifest`) | a JSON file |
| `write_records()` | a `Sequence[BaseModel]` (e.g. `list[ChainRecord]`) | a JSON or Parquet file |

## The parsed molecule

`mmcif_to_structure()` (`pandora.parsing`) turns a raw mmCIF file into a
`Structure` — the framework's central data model. Every later stage
(`canonicalise_structure()`, `collect_metadata()`, the `annotate_*()`
functions, `structure_to_mmcif()`) takes a `Structure` in and, where it
transforms one, returns a new `Structure` rather than mutating in place.
`pandora.datasets.extract_*_records()` reshapes a canonical `Structure`
into flatter, ML-friendlier views (`ChainRecord`, `ResidueRecord`,
`InterfaceRecord`) for export.

<img src="../assets/diagrams/structure.svg" alt="Structure entity-relationship diagram" style="max-width: 100%;">

## Canonicalisation policy

`canonicalise_structure()` runs nine rule groups in a fixed order — chain
IDs, residue numbering, assemblies, entities, missing data, altlocs,
ligands, then validation — each configured by one rules sub-model bundled
into a single `canonicalisationPolicy`. This tree is exactly what a
policy YAML file (see [Policies](../reference/policies.md)) deserializes into.

<img src="../assets/diagrams/canonicalisation-policy.svg" alt="canonicalisationPolicy entity-relationship diagram" style="max-width: 100%;">

## Dataset curation policy

`curate_structure()` filters a canonical `Structure` by quality
(resolution, experimental method, chain length), source organism, and
non-polymer content — governed by a much smaller `DatasetCurationPolicy`.

<img src="../assets/diagrams/curation-policy.svg" alt="DatasetCurationPolicy entity-relationship diagram" style="max-width: 100%;">

## Provenance and the dataset manifest

Every stage produces its own provenance record; `build_provenance_bundle()`
collects one structure's worth (`ProvenanceBundle`) and
`build_dataset_manifest()` aggregates dataset-wide provenance —
curation, deduplication, clustering, partitioning — plus every retained
structure's bundle into one `DatasetManifest`, the single JSON file that
`reproduce_dataset()` can later replay from. The `canonicalisationPolicy`
and `DatasetCurationPolicy` boxes below are the same models diagrammed
above, shown collapsed here to keep this diagram readable.

<img src="../assets/diagrams/provenance-manifest.svg" alt="DatasetManifest entity-relationship diagram" style="max-width: 100%;">

## Standalone models

A few models aren't pictured because they don't reference — or get
referenced by — any other schema:

- **`AnnotationLayer`** (`pandora.schemas.annotation`) is the live result
  returned by each `annotate_*()` function; only its lighter-weight
  sibling `AnnotationProvenanceRecord` (pictured above, under
  `ProvenanceBundle`) gets persisted into a manifest.
- **`Diagnostic`/`DiagnosticBundle`** (`pandora.schemas.common`) carry
  parsing/canonicalisation warnings and errors; they're threaded through
  as function arguments, not stored on any other model.
- **`FetchOptions`** (`pandora.schemas.ingestion`) configures one
  `fetch_mmcif()` call and isn't retained afterward — only its result,
  `IngestionProvenance`, is.
