# CLI `canonicalise` / `curate` step outputs — design

Follows `2026-11-04-modular-canonicalisation-curation-design.md`, which
gave the library standalone steps and rebuilt PPI v2.1 part A as two
scripts (`examples/ppi/v2.1/canonicalise.py`, `curate.py`). This brings
the CLI in line, so running the two commands gives the same step
folders as the scripts.

## Intent

The CLI is a thin wrapper over the library (`pandora/cli/README.md`).
Today `pandora canonicalise` and `pandora curate` write an older shape:
a flat folder, one big JSON dict per step, every structure held in
memory, and a single canonicalisation error stops the run. After this
change a step's folder looks the same whether it was written by the CLI
or by a user's own loop.

## Decisions (user)

- The old output files **disappear**: `canonicalisation_provenance.json`,
  `canonicalisation_mappings.json`, `curation_provenance.json`,
  `curation_exclusions.json`. No compatibility period (Pandora 0.5.6).

## Goals

- `pandora canonicalise` writes, under `--output-dir`:

  | File | One per | Content |
  |---|---|---|
  | `structures/<mid>/<stem>.cif` | canonical entry | canonical structure |
  | `entries.<fmt>` | canonical entry | `CanonicalisationEntryRecord` |
  | `chain_copies.<fmt>` | symmetry-copy chain | `ChainCopyRecord` |
  | `mappings.jsonl` | canonical entry | `{"entry_id": ..., **CanonicalMappings}`, one line per entry |
  | `failures.jsonl` | input that raised | `path`, `error_type`, `message` |
  | `canonicalisation.json` | run | `StepSummary` |

- `pandora curate` writes, under `--output-dir`:

  | File | One per | Content |
  |---|---|---|
  | `structures/<mid>/<stem>.cif` | kept entry | curated structure |
  | `entries.<fmt>` | kept entry | `CurationEntryRecord` |
  | `chains.<fmt>` | kept chain | `DatasetChainRecord` (with `copy_of` from the input's `chain_copies`, if present) |
  | `exclusions.<fmt>` | dropped entry or chain | `ExclusionRecord` |
  | `failures.jsonl` | input that raised | as above |
  | `curation.json` | run | `StepSummary`; `previous_summary` = input's `canonicalisation.json` if present |

- `<fmt>` is `parquet` by default; `--table-format jsonl` for installs
  without the `export` extra. A table with no rows is **not written**
  (as in the scripts); the summary's counts still say 0.
- Both process one entry at a time and record a failing entry in
  `failures.jsonl` instead of stopping. Exception types caught:
  `ValueError, OSError, RuntimeError, KeyError, IndexError` (same as
  the scripts).
- `mappings.jsonl` keeps everything `canonicalisation_mappings.json`
  had (chain-id remaps, residue renumbering, modified residues, altloc
  choices, assembly copies), one JSON object per line so it streams.
- File stems keep the input file name's case (`10kz`), as in the
  scripts.

## Inputs

- `canonicalise --input-dir DIR`: every `*.cif` in `DIR` **or** in its
  `*/` mirror subfolders. New `--entries FILE` (a path list read with
  `read_entry_paths`) as an alternative; exactly one of the two.
- `curate --input-dir DIR`: if `DIR/structures/` exists (a canonicalise
  folder), read structures from there and `DIR/chain_copies.*`; else
  read `DIR` like `canonicalise` does.

## Other commands

- The shared reader `read_structures_dir` (`pandora/cli/_io.py`) accepts
  a flat folder, a mirror layout, or a step folder (`DIR/structures/`),
  so `similarity`, `annotate`, `dedup` and `manifest` keep working on
  `curate`'s output.
- `write_structures_dir` writes the mirror layout too (my call), so
  `dedup` and `reproduce` produce the same layout as the steps. One
  convention for every command that writes structures.
- `manifest`:
  - `--canonicalisation-provenance` (the old dict JSON) is replaced by
    `--canonicalisation-dir`, a canonicalise output folder. Its
    `entries` table gives each entry's `canonicalisationProvenance`
    (rebuilt from the flat row; `report` from `warnings`/`errors`).
  - `--exclusions` accepts the new `exclusions.<fmt>` tables.
  - New `--curation-dir`, filling `ProvenanceBundle.curation` per entry
    from its `entries` table.

## What Pandora gains

- `read_records(path, model) -> list[model]` in `pandora/export/records.py`
  (exported), the inverse of `write_records` for `.jsonl`, `.json`,
  `.parquet`; Parquet nulls come back as `None`, not NaN (the bug part B
  hit). The CLI and the PPI scripts use it to read tables back.
- Rebuilding provenance from the flat rows (for `manifest`): two small
  functions next to `canonicalisation_entry_record` /
  `curation_entry_record` that go the other way.

## Non-goals

- No FASTA or per-chain mmCIF in `curate`. Those are a PPI tool-input
  choice; they stay in `examples/ppi/v2.1/curate.py`.
- No change to the library step functions themselves.
- No parallelism.

## Docs and tests

- `tests/test_cli.py`: the chain test and any test reading old files
  move to the new files; new tests for `--entries`, a failing entry,
  `--table-format jsonl`, mirror-layout input to `similarity`.
- Every usage page and recipe that runs `pandora canonicalise` /
  `curate` / `manifest` or reads the old files is updated with real
  outputs (`docs/usage/canonicalisation.md`, `datasets.md`, `cli.md`,
  `provenance.md`, `ingestion.md`, `metadata.md`,
  `docs/recipes/benchmark-quality-filters.md`), plus
  `pandora/cli/README.md`. `jq` on a dict JSON becomes `jq` on JSONL
  (`--table-format jsonl` in docs examples so `jq` works without
  Parquet tooling), e.g. the transforms of one entry:
  `jq -c 'select(.entry_id=="104M") | .transforms' canonical/entries.jsonl`.
- The PPI scripts switch their Parquet reads to `read_records`.

## Calls made by Claude, to check

- `write_structures_dir` moves to the mirror layout for every command.
- Docs examples use `--table-format jsonl`, so readers can `jq` them.
- `manifest` takes step folders (`--canonicalisation-dir`,
  `--curation-dir`) rather than individual files.
