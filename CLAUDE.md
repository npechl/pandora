# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What Pandora is

Pandora is a Python library of building blocks for preparing ML-ready
datasets from structural biology data (PDB/PDBe mmCIF today). It provides
typed data models, plain functions, and policies that users combine for
their own use case: parse a structure, normalise it, attach metadata and
annotations, filter and deduplicate, compute similarity, split without
leakage, record provenance, and export.

It is a toolkit, not a pipeline. The set of components and how they fit
together are still evolving. Do not assume a fixed stage order. Do not
add orchestrators, runners, registries, or framework base classes unless
asked. `examples/` shows some ways to combine the pieces; none of them
is the canonical flow.

## Architecture decisions

- **Schemas are separate from logic.** `pandora/schemas/` holds pydantic
  models only, no behaviour. Each logic package (`pandora/<package>/`)
  builds on its schema counterpart.
- **Functions, not objects with state.** Public API is plain functions
  that take a `Structure` or typed record and return a new one. No global
  state, no hidden configuration.
- **Never mutate inputs.** Return new objects via
  `.model_copy(update=...)`.
- **Behaviour is driven by policies.** Choices a user might want to vary
  (canonicalisation rules, curation filters, split settings) live in
  pydantic policy models, loadable from YAML via `load_policy()` (see
  `datasets/canonicalisation.yaml`; `docs/policies/*.yaml` only document
  each policy's field types).
  Do not hard-code them.
- **Report problems, don't swallow them.** Use `Diagnostic` /
  `DiagnosticBundle` (`schemas/common.py`) for recoverable issues. Catch
  named exception types only; Ruff enforces no blind `except Exception`
  (BLE001) and no silent `try/except/pass` (S110).
- **Provenance is a recipe, not a checksum.** Record source, policies,
  and tool versions so a dataset can be rebuilt. Pandora deliberately
  computes no content checksums.
- **Unknown mmCIF data is kept.** The parser promotes known categories to
  typed records and keeps every other category verbatim in
  `Structure.raw`.
- **Heavy or external tools stay optional.** Extra Python dependencies go
  behind an extra in `pyproject.toml` and are imported inside the function
  that needs them (see `export/records.py`). MMseqs2 and Foldseek are
  external binaries on `PATH`, not Python dependencies.

## Where to look

| Area | Logic | Models |
|---|---|---|
| Fetching files and policy loading | `pandora/ingestion/` | `schemas/ingestion.py` |
| mmCIF to `Structure` | `pandora/parsing/` | `schemas/structure.py` |
| Normalisation (chain IDs, residues, assemblies, altlocs, ligands, ...) | `pandora/canonicalisation/` (one module per rule group, `canonicalise.py` runs them) | `schemas/canonicalisation.py` |
| Entry, entity, quality, taxonomy, UniProt metadata | `pandora/metadata/` | `schemas/metadata.py` |
| Derived annotation layers | `pandora/annotations/` | `schemas/annotation.py` |
| Chain/residue/interface records, curation, dedup | `pandora/datasets/` | `schemas/dataset.py` |
| Similarity, clustering, leakage-safe splits | `pandora/similarity/` | `schemas/similarity.py` |
| Provenance bundles, dataset manifests, rebuilds | `pandora/provenance/` | `schemas/provenance.py` |
| mmCIF / JSON / JSONL / Parquet output | `pandora/export/` | — |
| Command-line wrappers | `pandora/cli/app.py` (argparse; I/O convention in `pandora/cli/README.md`) | — |

Other places:

- `docs/reference/policies.md` is the authoritative description of what
  every policy field does, including fields accepted but not implemented.
  Check it before assuming a policy field has an effect.
- `docs/contributing/architecture.md` maps functions to the models they
  take and return.
- `datasets/dev/mmcif/` holds offline mmCIF fixtures for tests and
  examples.

## Preferred libraries

- pydantic v2 for all data models.
- gemmi for reading and writing mmCIF.
- httpx for HTTP (the `ingestion` extra).
- pandas and pyarrow only for Parquet export (the `export` extra).
- PyYAML for policy files.
- argparse (stdlib) for the CLI.
- pytest for tests, Ruff for lint and format, uv for environments,
  Zensical for docs.

Reach for the stdlib or an existing dependency first. Ask before adding a
new dependency (for example BioPython, Biotite, or a CLI framework).

## Coding standards

- PEP 8, formatted by `ruff format`. Line length 80.
- Type hints on all new code. Start modules with
  `from __future__ import annotations`.
- Public functions get a Google-style docstring with `Args:`,
  `Returns:`, and `Raises:` where relevant. The docs site generates its
  reference pages from these. Internal helpers get one short line.
- Public names are exported through the package `__init__.py` and its
  `__all__`. Internal helpers start with `_`.
- Keep functions small and composable. Prefer a new function over a new
  flag on an existing one.

## Tests

- One test file per area in `tests/`, plain pytest functions.
- Use the fixtures in `datasets/dev/mmcif/`. Tests must not need network
  access or real MMseqs2/Foldseek binaries; fake the binary with a small
  shell script (see `tests/test_sequence_similarity.py`).
- Pin a known bug with `pytest.mark.xfail(strict=True, reason=...)` so
  the test starts failing once the bug is fixed.

## Commands

```sh
uv sync --all-extras                 # install into .venv/ from uv.lock
uv run pytest                        # all tests
uv run pytest tests/test_parsing.py::test_missing_model_falls_back_with_warning
uv run ruff format .                 # format
uv run ruff check .                  # lint
uv run zensical build --clean        # build the docs site
```

If you move the project directory, delete `.venv/` and run
`uv sync --all-extras` again. The venv scripts keep absolute paths, so
`uv run pytest` fails with `No module named 'pandora'`.

## Known issues

- `parsing/mmcif.py::_cs()` keeps CIF quote delimiters in string values
  (`"'X-ray diffraction'"`) and `;` markers in multi-line values,
  including sequences. `export/mmcif.py::_unwrap()` works around it.
  Tests in `tests/test_parsing.py` pin it as `xfail`.
- `canonicalise_structure` does not return its `DiagnosticBundle`; only
  warning/error counts surface, and only when
  `provenance_rules.emit_canonicalisation_report=True`.
- `reproduce_dataset` is a best-effort rebuild, not byte-identical.
  Upstream data and external tool versions can drift.

## Review checklist

Before calling a change done:

- [ ] `uv run pytest`, and
      `uv run ruff format --check .` pass.
- [ ] New behaviour has a test that uses local fixtures only.
- [ ] No input object is mutated; new objects come from `model_copy`.
- [ ] New models live in `pandora/schemas/`, new logic in the matching
      package.
- [ ] User-facing choices are policy fields, not hard-coded values.
- [ ] Public functions have type hints, a full docstring, and an
      `__all__` entry.
- [ ] No new dependency without agreement; optional ones sit behind an
      extra and are imported lazily.
- [ ] If a policy field's behaviour changed, `docs/reference/policies.md`
      is updated. If docs changed, the docs site builds.
- [ ] If schema relationships changed, the diagrams are regenerated:
      `uv run --extra docs python docs/scripts/generate_erd.py`.
