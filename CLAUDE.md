# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Conventions

This is primarily a Python project (with Markdown docs and YAML config). Follow PEP 8 conventions and use type hints for new Python code.

## What this is

Pandora is a Python library that turns raw PDB/PDBe mmCIF files into typed, policy-driven, ML-ready protein structure datasets. Every stage is a plain function: pass a `Structure` (or typed record) in, get one out — nothing is hidden behind a framework object or global state.

**Status:** ingestion, parsing, canonicalisation, metadata, annotations, export (`pandora/export/`), and a per-structure provenance bundle (`pandora/provenance/`) are implemented. Dataset curation (`pandora/datasets/curation.py`) and the CLI (`pandora/cli/app.py`) are still stubs (`# TODO` / `raise NotImplementedError`) — part of the roadmap, not a bug.

## Commands

```sh
uv sync --all-extras        # install everything into .venv/ (locked via uv.lock)

uv run pytest                                    # run all tests
uv run pytest tests/test_clustering.py::test_transitive_merge_and_isolate  # single test
uv run ruff check .                              # lint
uv run ruff format .                             # format (format --check . for CI-style check)
```

- The `similarity` extra in `pyproject.toml` is deliberately empty — `pandora.similarity.sequence`/`.structure` shell out to the `mmseqs`/`foldseek` binaries, which must be installed separately and be on `PATH`. Nothing in `tests/` currently exercises them.
- `pandora/ingestion` needs network access (fetches from PDBe/RCSB); `datasets/dev/mmcif/` has local fixture files for offline work — see `examples/overview.py` for an end-to-end run against them.
- Ruff's rule set is deliberately pinned in `pyproject.toml` (`select = ["E4", "E7", "E9", "F", "BLE001", "S110"]`) rather than relying on Ruff's version-dependent defaults; line length is 80 but `E501` is ignored (formatter's job).

## Architecture

### Layout: schemas vs. logic

- `pandora/schemas/` — pydantic models only, no logic. One module per pipeline stage: `structure.py` (the mmCIF data model — `Structure`, `AtomSiteRecord`, etc.), `canonicalisation.py` (policy + rule + mapping + provenance models), `metadata.py`, `annotation.py`, `similarity.py`, `dataset.py`, `ingestion.py`, `provenance.py` (`ProvenanceBundle`, `Checksums`, `AnnotationProvenanceRecord`), and `common.py` (shared `Diagnostic`/`DiagnosticBundle`).
- `pandora/{ingestion,parsing,canonicalisation,metadata,annotations,similarity,datasets,provenance}/` — logic, one package per stage, each built on its schema counterpart. Structures are never mutated in place; every transform returns a new object via pydantic's `.model_copy(update=...)`.

### The implemented pipeline

```
fetch_mmcif()          pandora.ingestion    — HTTP fetch from PDBe/PDB, with on-disk caching
mmcif_to_structure()   pandora.parsing      — gemmi-backed mmCIF -> Structure
canonicalise_structure() pandora.canonicalisation — policy-driven normalization
collect_metadata()     pandora.metadata     — source-backed entry/entity/quality/taxonomy records
annotate_*()           pandora.annotations  — derived per-entry/pairwise layers (counts, contacts, interfaces)
extract_*_records()    pandora.datasets     — reshape a canonical Structure into Chain/Residue/Interface records
compute_*_similarity() pandora.similarity   — MMseqs2/Foldseek wrappers -> SimilarityRelationship
cluster_similar_items() / partition_dataset() pandora.similarity — leakage-safe train/val/test splitting
build_provenance_bundle() pandora.provenance — aggregates ingestion/canonicalisation/metadata/annotation provenance + a structure checksum
structure_to_mmcif() / write_json() / write_records() pandora.export — serialize a Structure or records back to mmCIF/JSON
```

Each stage is independently callable — `examples/overview.py` shows the intended chaining, but nothing requires running the whole thing.

### Canonicalisation: one orchestrator over nine rule modules

`canonicalisation/canonicalise.py::canonicalise_structure(structure, policy)` is the single entry point. It runs, in a fixed order, one function per rule group — `chain_ids`, `residues`, `assemblies`, `entities`, `missing_data` (atoms/residues/incomplete-chains), `altlocs`, `ligands`, then `validation` — each living in its own sibling module. Every step returns `(transformed_data, mapping)`; `canonicalise_structure` assembles the mappings into `CanonicalMappings` and appends a transform label (e.g. `"chain_id:remap"`) whenever a rule deviates from "preserve".

`docs/policies.md` is the authoritative, implementation-accurate reference for every policy field — including explicit callouts for pieces accepted by the schema but not yet implemented (`missing_atoms.strategy: impute`, `assembly_rules.strategy: standardize_biological_assembly`, `validation_rules.strictness: permissive`). Check there before assuming a policy field does something.

Known sharp edge: `_validate()` (`canonicalisation/validation.py`) computes a `"failed"/"warning"/"success"` status from `validation_rules`, but `canonicalise_structure` discards that return value — a policy with `fail_on_unresolved_issues=True` cannot currently signal failure to the caller. The collected `DiagnosticBundle` itself is also not part of the function's return tuple; only aggregate warning/error *counts* surface, and only when `provenance_rules.emit_canonicalisation_report=True` (default `False`).

### Raw category passthrough

`mmcif_to_structure` promotes well-known mmCIF categories to typed records but keeps every other category verbatim in `Structure.raw: dict[str, list[dict[str, str | None]]]`. `metadata/utils.py::raw_rows()`/`first_row()` are the only read path into `.raw`; `metadata/mmcif.py`'s `extract_*` functions are its only consumers. String values read via gemmi (`_cs()` in `parsing/mmcif.py`) are the raw CIF token — quoted values keep their literal `'...'` delimiters unless unquoted, since gemmi's `find_value()`/loop access don't do it for you.

### Provenance is per-structure, not per-dataset

`pandora/provenance/manifest.py::build_provenance_bundle(structure, ...)` assembles a `ProvenanceBundle` from whatever provenance the caller already has in hand — `ingestion` (`IngestionProvenance`), `canonicalisation` (`canonicalisationProvenance`), `metadata` (`MetadataRecord`, flattened via `collect_metadata_provenance()`), and `annotations` (a list of `AnnotationLayer`s) — plus a SHA-256 checksum of the structure (`checksums.compute_checksum()`, over the model's sorted-key JSON dump). All arguments are optional; it does not fetch, re-derive, or validate anything itself. There is no dataset-level manifest, artifact export, or checksum-of-checksums yet — see the design-doc caveat below.

### External-tool wrappers (`pandora/similarity/`)

`sequence.py` (MMseqs2) and `structure.py` (Foldseek) both follow the same shape: validate the binary is on `PATH`, materialize inputs to a temp directory (FASTA for sequences, symlinked structure files for structures — or point directly at a directory the caller already has), shell out to the tool's `easy-search`, parse the tab-separated result, dedupe to the best hit per `(source_id, target_id)` pair (`source_id < target_id`), and return `SimilarityRelationship` objects. Neither tool is a Python dependency — both are external binaries the caller must install.

### Design-doc vs. implemented API — don't conflate them

`docs/components/*.md` documents an aspirational six-component architecture (C01–C06: `PandoraDataset`, `LeakageSafeDataset`, `DatasetStore`, `PandoraArtifact`, function names like `ingest_mmCIF`, `attach_metadata`, `build_leakage_safe_dataset`) that is **not** what exists in `pandora/` today — none of those types/functions are implemented. This includes `docs/components/06-provenance.md`: the implemented `pandora/provenance/` is a much smaller, per-structure `build_provenance_bundle()` (see above), not the `PandoraArtifact`/embedded-vs-by-reference/manifest-export system the doc describes, because that system depends on `PandoraDataset`/`LeakageSafeDataset`/`DatasetStore`, none of which exist. Treat `docs/components/` as a target design spec, not API documentation. `docs/policies.md`, by contrast, documents the real, implemented canonicalisation policy and stays accurate to the code (including "not implemented yet" notes). When in doubt about what's actually callable, check `pandora/*/__init__.py`'s `__all__` or `examples/overview.py`, not `docs/components/`.

## Workflow

Prefer running commands with Bash tool for verification, but batch related shell operations and prefer non-destructive read commands first.