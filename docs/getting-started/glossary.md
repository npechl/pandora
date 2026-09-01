# Glossary

Terms as Pandora actually uses them today. For the underlying pydantic
models, see [Schemas](../reference/schemas.md); for every callable, see
[Functions](../reference/functions.md).

**Altloc** — An alternate location for an atom (mmCIF `label_alt_id`),
recording two or more conformations of the same residue observed in the
crystal. `altloc_rules` in a [canonicalisation policy](../reference/policies.md#altloc_rules)
decides how these are resolved into a single kept conformer.

**Annotation / annotation layer** — A derived, per-entry or pairwise
value computed from an already-canonicalised `Structure` — e.g. atom/chain
counts, chain-chain interfaces, ligand contacts, or pairwise sequence
identity. Produced by the `annotate_*()` functions in
[`pandora.annotations`](../usage/annotation.md) and recorded as an
`AnnotationLayer`.

**Canonicalisation** — Policy-driven normalization of a parsed `Structure`:
remapping chain IDs, resolving altlocs, handling missing atoms/residues,
filtering ligands, and validating the result. The single entry point is
`canonicalise_structure()`; see [Canonicalisation](../usage/canonicalisation.md).

**Canonicalisation policy** — The typed, versioned settings object
(`canonicalisationPolicy`) that controls every canonicalisation rule group
(`identifier_rules`, `missing_data_rules`, `altloc_rules`, `assembly_rules`,
`entity_rules`, `ligand_rules`, `validation_rules`, `provenance_rules`).
Fully documented field-by-field in [Policies](../reference/policies.md).

**CanonicalMappings** — The record of every id/coordinate change
`canonicalise_structure()` made (chain remaps, dropped atoms, altloc
selections, ...), returned alongside the transformed `Structure` so the
change is auditable rather than silent.

**Cluster / clustering** — Grouping items (entries, chains, ...) that a
`SimilarityRelationship` network connects into connected components, so
that near-duplicates land in the same dataset split. Done by
`cluster_similar_items()`; see [Clustering](../usage/similarity.md#clustering).

**Curation** — Policy-driven filtering of a dataset — dropping structures
that fail quality/size/resolution rules — via `curate_structure()`. See
[Curate one structure](../usage/datasets.md#curate-one-structure).

**Dataset manifest** — The dataset-level provenance report
(`DatasetManifest`) produced by `build_dataset_manifest()`: the curation
and canonicalisation policies (stored by value), exclusion records,
dedup/clustering/partition provenance, the split assignment, and every
retained structure's `ProvenanceBundle`, all in one JSON-serializable
object. See [Assemble a dataset manifest](../usage/provenance.md#assemble-a-dataset-manifest).

**Deduplication** — Removing structures that share the same `entry_id`
from a batch, via `deduplicate_structures()`. Distinct from clustering,
which groups by *similarity* rather than identical id.

**Diagnostic / DiagnosticBundle** — A structured warning or error
(`Diagnostic`) raised during canonicalisation or validation, collected
into a `DiagnosticBundle`. Only aggregate warning/error *counts* surface
in provenance today, gated behind `provenance_rules.emit_canonicalisation_report`.

**Entity** — An mmCIF `entity` — one distinct polymer or non-polymer
chemical component in the structure (e.g. "the protein", "a bound heme",
"water"), as opposed to a **chain**, which is one physical copy of an
entity in the asymmetric unit.

**Entry** — One PDB/PDBe deposition, identified by its 4-character
`entry_id` (e.g. `104M`). One entry's mmCIF file maps to one `Structure`.

**Ligand** — A non-polymer entity (mmCIF `HETATM` records) other than
water — a cofactor, inhibitor, or other bound small molecule. Filtered or
kept per `ligand_rules`; see [Filter ligands directly](../usage/canonicalisation.md#filter-ligands-directly).

**Partition / split** — Assigning clustered items to `train`/`val`/`test`
sets such that no cluster straddles a split boundary ("leakage-safe"),
via `partition_dataset()`. See [Leakage-safe partitioning](../usage/similarity.md#leakage-safe-partitioning).

**Provenance bundle** — The per-structure provenance report
(`ProvenanceBundle`) assembled by `build_provenance_bundle()`, aggregating
whatever ingestion/canonicalisation/metadata/annotation provenance the
caller already has for one structure. The building block a dataset
manifest is made of.

**Raw category / raw passthrough** — Any mmCIF category `mmcif_to_structure()`
doesn't promote to a typed field is kept verbatim in `Structure.raw`,
readable via `raw_rows()`/`first_row()`. Values here are the literal CIF
token as gemmi returns it — quoted strings keep their `'...'` delimiters.

**Reproduce** — Best-effort replay of the whole
fetch → canonicalise → curate → dedup → similarity → cluster → partition →
annotate pipeline from a `DatasetManifest` alone, via `reproduce_dataset()`.
Not a guaranteed byte-identical rebuild — source data and external tools
can drift. See [Replay a manifest from scratch](../usage/provenance.md#replay-a-manifest-from-scratch).

**Similarity relationship** — One scored pairwise comparison
(`SimilarityRelationship`) between two items, produced by
`compute_sequence_similarity()` (MMseqs2) or `compute_structure_similarity()`
(Foldseek). The input clustering and partitioning work from.

**Structure** — Pandora's central typed model of one parsed mmCIF entry:
atoms, entities, assemblies, plus the raw passthrough for everything else.
Never mutated in place — every pipeline stage returns a new `Structure`
via `.model_copy(update=...)`.

**Validation** — The final canonicalisation step: computing a
`"success"`/`"warning"`/`"failed"` status from the diagnostics collected
during the run, per `validation_rules.strictness`. A `"failed"` status
raises `ValueError` when `fail_on_unresolved_issues=True`.
