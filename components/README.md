# Pandora — Architecture Overview

Pandora is a **Python library** for producing leakage-safe, ML-ready protein
structure datasets from PDB/PDBe mmCIF files.

The library is built from six independent, composable components (C01–C06).
Each component exposes a set of typed, policy-driven functions that can be
called individually or composed into a full pipeline. The "pipeline" is a
convenience — any function can be called in isolation, and external data can
enter the library at any stage via adapter functions.

---

## Architecture Diagram

```mermaid
flowchart TD
    %% ── Styles ─────────────────────────────────────────────────────────────
    classDef comp  fill:#f0f7ee,stroke:#4caf50,color:#1a1a1a
    classDef iface fill:#fff9e6,stroke:#d4a017,color:#1a1a1a,font-weight:bold
    classDef ext   fill:#e8f4f8,stroke:#4a9eca,color:#1a1a1a
    classDef artif fill:#fdebd0,stroke:#e67e22,color:#1a1a1a,font-weight:bold
    classDef store fill:#e8f8f5,stroke:#1abc9c,color:#1a1a1a,font-weight:bold
    classDef adapt fill:#fdedec,stroke:#e74c3c,color:#1a1a1a

    %% ── External inputs & entry-point adapters ──────────────────────────────
    SRC_MMCIF["mmCIF source\n(PDBe / PDB / local / URI / bytes)"]
    AD_PARSED["from_parsed_structure()\nBioPython · MDAnalysis · custom\n→ produces MmCIFIngestionResult"]
    SRC_META["UniProt · SIFTS · Taxonomy"]
    AD_CANON["from_canonical_structures()\nlist[CanonicalStructureResult]\n→ enters C04 directly"]
    AD_ANN["from_annotated_structures()\nlist[AnnotatedStructureWithPlugins]\n→ enters C04 directly"]
    AD_STORE["External DatasetStoreRef\n(pre-built Parquet store)\n→ enters C05 directly"]

    %% ── C01 ─────────────────────────────────────────────────────────────────
    C01["📥 C01 · Ingestion\nfetch_mmCIF()\nparse_mmCIF()\nvalidate_mmCIF()\ningest_mmCIF()  ★"]

    I01["MmCIFIngestionResult\n───────────────────────\nstatus · diagnostics · provenance\natoms · residues · chains\nentities · assemblies · ligands"]

    %% ── C02 ─────────────────────────────────────────────────────────────────
    C02["🔧 C02 · Canonicalization\ncanonicalize_structure()  ★\nvalidate_canonical_structure()"]

    I02["CanonicalStructureResult\n───────────────────────\ncanonical_structure\ncanonical_mappings\napplied_policy · provenance"]

    %% ── C03 ─────────────────────────────────────────────────────────────────
    C03["🏷 C03 · Metadata & Annotation\nretrieve_metadata()\nattach_metadata()\napply_annotation_plugins()  ★"]

    I03["AnnotatedStructureWithPlugins\n───────────────────────\ncanonical_structure_result\nmetadata_annotations\nderived_annotations · applied_plugins"]

    %% ── C04 ─────────────────────────────────────────────────────────────────
    C04["🗃 C04 · Curation\nbuild_dataset()  ★\nfilter_dataset()\ndeduplicate_dataset()\nextract_chain / interface / residue_records()\nmaterialize_dataset()"]

    STORE["📦 DatasetStore\nParquet files on disk\n(materialized mode only)"]

    I04["PandoraDataset\n───────────────────────\nDataset | ChainDataset\nInterfaceDataset | ResidueDataset\ngranularity\nmode: in_memory | materialized\nstore: DatasetStoreRef | null"]

    %% ── C05 ─────────────────────────────────────────────────────────────────
    C05["✂️ C05 · Splitting\ncompute_similarity_relationships()\n  sequences/structures → files → MMseqs2 / Foldseek\ncluster_similar_items()\npartition_dataset()  (greedy bin-packing)\nbuild_leakage_safe_dataset()  ★"]

    I05["LeakageSafeDataset\n───────────────────────\nsource_dataset_id\nsource_dataset: PandoraDataset | null\nsource_dataset_ref: DatasetStoreRef | null\npartitions: train / val / test\nsimilarity_network · clusters"]

    %% ── C06 ─────────────────────────────────────────────────────────────────
    C06["📋 C06 · Provenance\nbuild_provenance_bundle()\ngenerate_manifest()\ncompute_checksums()\nbuild_pandora_artifact()  ★"]

    I06["🏆 PandoraArtifact\n───────────────────────\nartifact_mode: embedded | by_reference\nmanifest · checksums · provenance_bundle\nembedded → LeakageSafeDataset in memory\nby_reference → ArtifactStore on disk"]

    %% ── Primary data flow (solid arrows = typed interface handoffs) ──────────
    SRC_MMCIF  -->  C01
    AD_PARSED  -->|"→ MmCIFIngestionResult"| I01
    C01        -->  I01

    I01        -->  C02
    C02        -->  I02

    SRC_META   -->  C03
    I02        -->  C03
    C03        -->  I03

    I03        -->  C04
    C04        -->|"in_memory"| I04
    C04       -.->|"materialized"| STORE
    STORE     -.->|"DatasetStoreRef"| I04

    I04        -->  C05
    C05        -->  I05

    I05        -->  C06
    C06        -->  I06

    %% ── Entry-point adapters (dashed = bypass upstream components) ───────────
    AD_CANON  -.->|"skip C01–C03"| C04
    AD_ANN    -.->|"skip C01–C03"| C04
    AD_STORE  -.->|"skip C01–C04"| I04

    %% ── Styles ──────────────────────────────────────────────────────────────
    class C01,C02,C03,C04,C05,C06 comp
    class I01,I02,I03,I04,I05 iface
    class I06 artif
    class STORE store
    class SRC_MMCIF,SRC_META ext
    class AD_PARSED,AD_CANON,AD_ANN,AD_STORE adapt
```

---

## Components at a glance

| Component | Key functions | Input | Output |
|-----------|--------------|-------|--------|
| **01 · Ingestion** | `fetch_mmCIF`, `parse_mmCIF`, `validate_mmCIF`, `ingest_mmCIF` | Raw mmCIF (PDBe/PDB/local/bytes) | `MmCIFIngestionResult` |
| **02 · Canonicalization** | `canonicalize_structure`, `validate_canonical_structure` | `MmCIFIngestionResult` + CanonicalizationPolicy | `CanonicalStructureResult` |
| **03 · Metadata & Annotation** | `attach_metadata`, `apply_annotation_plugins` | `CanonicalStructureResult` + Metadata/Plugin Policies | `AnnotatedStructureWithPlugins` |
| **04 · Curation** | `build_dataset`, `extract_chain_records`, `materialize_dataset` | `[AnnotatedStructureWithPlugins]` + CurationPolicy | `PandoraDataset` (in-memory or materialized) |
| **05 · Splitting** | `compute_similarity_relationships`, `cluster_similar_items`, `partition_dataset`, `build_leakage_safe_dataset` | `PandoraDataset` or `DatasetStoreRef` + LeakagePolicy | `LeakageSafeDataset` |
| **06 · Provenance** | `build_provenance_bundle`, `generate_manifest`, `build_pandora_artifact` | `LeakageSafeDataset` + ProvenancePolicy | `PandoraArtifact` (embedded or by-reference) |

---

## Library design principles

### 1. Library-first: every function is independently callable

Each component exposes typed, independently callable functions. Callers are
not required to run all preceding components. A researcher with existing
parsed structures can enter at C02; a team that already has canonical
structures can enter at C03; a team with a curated structure collection can
enter at C04 without using Pandora's ingestion layer at all.

The full pipeline is a convenience — it is implemented by composing the
same library functions that users can call directly.

```python
# Call any function independently
result = canonicalize_structure(my_ingestion_result, policy)
dataset = build_dataset(my_annotated_structures, curation_policy)
artifact = build_pandora_artifact(leakage_safe_dataset, prov_policy, export_policy)
```

### 2. Entry-point adapters for external data

Each component provides adapter functions that allow external data to enter
without requiring upstream Pandora types. Adapters wrap external data into
the expected Pandora schemas with sensible provenance defaults.

| Entry stage | Adapter | External input |
|-------------|---------|----------------|
| C01 (parse) | `from_raw_bytes(bytes)` | Pre-loaded mmCIF bytes |
| C02 (canonicalize) | `from_parsed_structure(structure, entry_id)` | BioPython / MDAnalysis / custom parsed structure |
| C04 (curate) | `from_canonical_structures(results, policy)` | list of CanonicalStructureResult without metadata |
| C04 (curate) | `from_annotated_structures(structures, policy)` | External annotated structure list |

### 3. Two execution modes for C04–C06

**In-memory mode** (default) — suitable for datasets up to ~10K structures:

```
C01 → C02 → C03 → C04 (list in RAM) → C05 → C06 (embedded artifact)
```

All intermediate results are held in memory. The final `PandoraArtifact`
embeds the full `LeakageSafeDataset` as a nested Python object.

**Materialized (streaming) mode** — required for datasets above ~10K structures:

```
C01 → C02 → C03 → C04 (write to DatasetStore) → C05 (reads IDs + files) → C06 (by-reference artifact)
```

C01-C04 processes one structure at a time and immediately writes each record
to a `DatasetStore` (Parquet files on disk). C05 reads only item IDs and
sequence/structure files from the store — it never loads full atom coordinate
objects into memory. C06 produces a manifest + split Parquet files rather than
a monolithic nested object.

### 4. C01–C04 is per-structure; C05–C06 is collection-level

Components 1–4 are **embarrassingly parallel** — each structure is processed
independently. They can be run as a stream: one structure in, one record out,
written to the store. No global state is needed.

Component 5 is the **global barrier** — it requires all item IDs and pairwise
similarity scores before it can cluster and partition. This is inherent to
leakage-safe splitting, not a design limitation. However, C05 only needs item
identifiers and sequence/structure files, not full atom-coordinate objects.

Component 6 is **lightweight metadata assembly** — it reads provenance fields
and partition lists, then generates manifests and checksums.

### 5. Policy-driven and reproducible

Every component accepts a typed policy object that makes all decisions
explicit, versionable, and reproducible. The identity of any generated
`PandoraArtifact` is fully determined by: source archive release, all policy
objects applied, and Pandora software version.

---

## Key design invariants

- **Components 1–3** are structure-centric. **Component 4** introduces
  `PandoraDataset` as the first-class collection object, at four granularity
  levels: structure, chain, interface, and residue.

- **Extraction is always downstream of curation** — `ChainDataset`,
  `InterfaceDataset`, and `ResidueDataset` are derived from a curated
  structure-level `Dataset`. Extraction is optional; the library proceeds at
  whichever granularity the user selects.

- **`PandoraDataset` is a discriminated union** —
  `Dataset | ChainDataset | InterfaceDataset | ResidueDataset`,
  discriminated by the `granularity` field. Components 5 and 6 accept any of
  these types. Each type carries a `mode` field (`in_memory | materialized`)
  that indicates whether records are embedded or held in a `DatasetStore`.

- **Metadata is attached, never embedded** — the canonical structure is never
  mutated after Component 02.

- **Component 05** wraps external engines (MMseqs2, Foldseek). In both modes,
  sequences and structures are extracted to files before invoking the external
  tool. The engine outputs are then read back and normalized to
  `SimilarityRelationship` records.

- **Component 06** supports two artifact modes. In `embedded` mode (default
  for small datasets), the `PandoraArtifact` holds the full
  `LeakageSafeDataset` in memory. In `by_reference` mode (large datasets),
  the artifact is a manifest + checksums pointing to split Parquet files in an
  `ArtifactStore` directory.

- **Component 06 provenance depth depends on granularity** — structure-level
  datasets yield full upstream provenance (ingestion → splitting);
  chain/interface/residue datasets yield only curation and splitting
  provenance (`UPSTREAM_PROVENANCE_NOT_EMBEDDED` — upstream
  `AnnotatedStructureWithPlugins` objects are not embedded at sub-structure
  granularities).

- **Component 06** only *records* — it never performs ingestion,
  transformation, or splitting.

- Every component accepts a **policy object** that makes all decisions
  explicit and reproducible.

- All batch variants (`ingest_list_mmCIF`, `canonicalize_many_structures`, …)
  are thin orchestration wrappers over the single-entry functions.
