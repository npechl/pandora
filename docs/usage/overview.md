# Overview

[Getting Started](../getting-started/overview.md) is a single walkthrough
of the pipeline end to end. This section is the per-stage reference
underneath it — one page per stage, covering every function's
signature, edge cases, and (where one exists) the equivalent `pandora`
CLI subcommand.

Most pages open with a short setup block (parse a fixture, canonicalise
it) that the rest of the page's examples build on — skim that first if
you're jumping straight to a page rather than reading in order.

## Library vs. CLI

Every operation that a `pandora` subcommand wraps is shown twice:

=== "`library`"

    Call the Python function directly — for scripting the pipeline in
    one process, or when you need something a subcommand doesn't
    expose (a keyword argument, an in-memory object).

=== "`cli`"

    The equivalent `pandora <subcommand>` invocation — for running a
    stage from the shell over a directory of files, chaining it with
    the next stage's `--input-dir`.

Both call the same underlying function; the CLI just reads/writes
files where the library form passes objects in memory. A handful of
steps (loading a policy, filtering ligands standalone, reshaping
records) have no CLI equivalent — those pages say so inline instead of
showing an empty tab.

## Pages

| Page | Covers |
|------|--------|
| [Ingestion](ingestion.md) | Fetch mmCIF files from PDBe/PDB (with on-disk caching), load a canonicalisation policy YAML |
| [Canonicalisation](canonicalisation.md) | Normalize a parsed `Structure` against a policy — chain IDs, altlocs, entities, ligands, and more |
| [Metadata](metadata.md) | Extract source-backed entry/quality/taxonomy/entity/ligand/UniProt records — library-only, no CLI subcommand |
| [Annotation](annotation.md) | Derived per-entry/pairwise layers: structure counts, ligand contacts, chain interfaces, sequence identity |
| [Datasets](datasets.md) | Curate and deduplicate a batch, reshape a canonical `Structure` into flat Chain/Residue/Interface records |
| [Similarity](similarity.md) | All-vs-all sequence/structure similarity, connected-component clustering, leakage-safe train/val/test partitioning |
| [Export](export.md) | Write a `Structure` (or the records from Datasets) back to mmCIF, JSON, or Parquet |
| [Provenance](provenance.md) | Assemble a per-structure `ProvenanceBundle` or dataset-wide `DatasetManifest`, and replay one from scratch |

See [Policies](../reference/policies.md) for every canonicalisation/curation
policy field, and [Functions](../reference/functions.md) /
[Schemas](../reference/schemas.md) for generated API reference.
