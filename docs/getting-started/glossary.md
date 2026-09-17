# Glossary

Pandora defines several terms of its own on top of raw mmCIF vocabulary. Each entry below links to where the concept actually lives in code, and to a usage guide or recipe that puts it to work — follow those if you want to see the term in action rather than just read about it. 

For the underlying pydantic models, see [Schemas](../reference/schemas.md); for every callable, see [Functions](../reference/functions.md).

## Data model

### Structure {: #structure }

Pandora's central typed model of one parsed mmCIF entry: atoms, entities, assemblies, plus the raw passthrough for everything else. Never mutated in place — every pipeline stage returns a new `Structure` via `.model_copy(update=...)`. See [Ingestion & parsing](../usage/ingestion.md).

### Entry {: #entry }

One PDB/PDBe deposition, identified by its 4-character `entry_id` (e.g. `104M`). One entry's mmCIF file maps to one [`Structure`](#structure).

### Entity {: #entity }

An mmCIF `entity` — one distinct polymer or non-polymer chemical component in the structure (e.g. "the protein", "a bound heme", "water") — as opposed to a [chain](#chain), which is one physical copy of an entity in the asymmetric unit.

### Chain {: #chain }

One physical copy of an [entity](#entity) in the asymmetric unit, identified by its `label_asym_id`. Two chains can share the same entity (e.g. a homodimer's two copies of the same protein) — [canonicalisation](#canonicalisation) can remap chain IDs; see `identifier_rules` in [Policies](../reference/policies.md).

### Raw category / raw passthrough {: #raw-passthrough }

Any mmCIF category `mmcif_to_structure()` doesn't promote to a typed field, kept verbatim in `Structure.raw`, readable via `raw_rows()`/`first_row()`. Values here are the literal CIF token as gemmi returns it — quoted strings keep their `'...'` delimiters.

## Canonicalisation {: #canonicalisation-section }

### Canonicalisation {: #canonicalisation }

Policy-driven normalization of a parsed [`Structure`](#structure): remapping chain IDs, resolving altlocs, handling missing atoms/residues, filtering ligands, and validating the result. The single entry point is `canonicalise_structure()`; see [Canonicalisation](../usage/canonicalisation.md).

### Canonicalisation policy {: #canonicalisation-policy }

The typed, versioned [policy](#policy) object (`canonicalisationPolicy`) that controls every canonicalisation rule group (`identifier_rules`, `missing_data_rules`, `altloc_rules`, `assembly_rules`, `entity_rules`, `ligand_rules`, `validation_rules`, `provenance_rules`). Fully documented field-by-field in [Policies](../reference/policies.md).

### CanonicalMappings {: #canonicalmappings }

The record of every id/coordinate change [canonicalisation](#canonicalisation) made (chain remaps, dropped atoms, altloc selections, ...), returned alongside the transformed `Structure` so the change is auditable rather than silent.

### Diagnostic / DiagnosticBundle {: #diagnostic }

A structured warning or error (`Diagnostic`) raised during canonicalisation or [validation](#validation), collected into a `DiagnosticBundle`. Only aggregate warning/error *counts* surface in provenance today, gated behind `provenance_rules.emit_canonicalisation_report`.

### Validation {: #validation }

The final canonicalisation step: computing a `"success"`/`"warning"`/`"failed"` status from the [diagnostics](#diagnostic) collected during the run, per `validation_rules.strictness`. A `"failed"` status raises `ValueError` when `fail_on_unresolved_issues=True`.

## Annotations

### Annotation / annotation layer {: #annotation }

A derived, per-entry or pairwise value computed from an already-[canonicalised](#canonicalisation) `Structure` — e.g. atom/chain counts, chain-chain interfaces, ligand contacts, or pairwise sequence identity. Produced by the `annotate_*()` functions in [`pandora.annotations`](../usage/annotation.md) and recorded as an `AnnotationLayer`.

## Curation, similarity & datasets

### Curation {: #curation }

Policy-driven filtering of a dataset — dropping structures that fail quality/size/resolution rules — via `curate_structure()`. See [Curate one structure](../usage/datasets.md#curate-one-structure).

### Deduplication {: #deduplication }

Removing structures that share the same `entry_id` from a batch, via `deduplicate_structures()`. Distinct from [clustering](#clustering), which groups by *similarity* rather than identical id.

### Similarity relationship {: #similarity-relationship }

One scored pairwise comparison (`SimilarityRelationship`) between two items, produced by `compute_sequence_similarity()` (MMseqs2) or `compute_structure_similarity()` (Foldseek). The input [clustering](#clustering) and [partitioning](#partition) work from.

### Cluster / clustering {: #clustering }

Grouping items (entries, chains, ...) that a [similarity relationship](#similarity-relationship) network connects into connected components, so that near-duplicates land in the same dataset split. Done by `cluster_similar_items()`; see [Clustering](../usage/similarity.md#clustering).

### Partition / split {: #partition }

Assigning clustered items to `train`/`val`/`test` sets such that no [cluster](#clustering) straddles a split boundary ("leakage-safe"), via `partition_dataset()`. See [Leakage-safe partitioning](../usage/similarity.md#leakage-safe-partitioning) and the worked example in the [PPI recipe](../recipes/ppi-01.md).

## Provenance

### Provenance bundle {: #provenance-bundle }

The per-structure provenance report (`ProvenanceBundle`) assembled by `build_provenance_bundle()`, aggregating whatever ingestion/canonicalisation/metadata/annotation provenance the caller already has for one structure. The building block a [dataset manifest](#dataset-manifest) is made of.

### Dataset manifest {: #dataset-manifest }

The dataset-level provenance report (`DatasetManifest`) produced by `build_dataset_manifest()`: the [curation](#curation) and [canonicalisation policies](#canonicalisation-policy) (stored by value), exclusion records, dedup/clustering/partition provenance, the split assignment, and every retained structure's [provenance bundle](#provenance-bundle), all in one JSON-serializable object. See [Assemble a dataset manifest](../usage/provenance.md#assemble-a-dataset-manifest).

## Documentation conventions

### Policy {: #policy }

Pandora's convention for making every stage's rules explicit and inspectable: a typed, versioned pydantic settings object (`canonicalisationPolicy`, `DatasetCurationPolicy`, ...) passed into the stage's function rather than hard-coded. Provenance records reference the policy that was applied via the compact `AppliedPolicyRef` (id, name, version).

### Recipe {: #recipe }

A narrated, end-to-end walkthrough doc under `docs/recipes/` (e.g. the [PPI interface dataset](../recipes/ppi-01.md) recipe) that chains multiple pipeline stages into one realistic workflow, as opposed to the single-stage guides under `docs/usage/`.
