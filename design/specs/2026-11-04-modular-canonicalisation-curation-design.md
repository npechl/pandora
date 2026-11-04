# Modular canonicalisation and curation steps — design

Supersedes the output layout of
`2026-10-19-ppi-v21-prepare-entries-design.md` (part A). Part A's
filters and policies are unchanged.

## Intent

Someone preparing a structural-biology dataset should be able to write
the canonicalisation and curation stages with as little code as
possible, by calling Pandora functions in a loop they write themselves.
Today `examples/ppi/v2.1/a_prepare_entries.py` hand-builds that glue:
it combines parse, canonicalise and curate in one function, builds
chain rows as loose dicts, writes its own mmCIF tree, tallies counts,
and discards the provenance that Pandora already returns.

Each stage should be **a standalone step**: run it, and it leaves
behind its structures plus a provenance record that says how they were
made. A user can stop after any step, give its output to another
tool, or start the next step from it.

## Goals

- Canonicalisation and curation are two separate steps, each with its
  own outputs and its own provenance. Both steps can run in one loop or
  in two scripts with the same result.
- Canonical structures are written as **mmCIF** (user decision), so
  other tools can read them directly.
- Each step writes its provenance as files: a small step summary and
  per-entry rows.
- Pandora supplies the functions a step script needs. The user only
  writes the loop: iterating paths, catching failures, collecting rows.
- Tables are typed Pandora models, joinable on `(entry_id, chain_id)`.
- PPI v2.1 part A is split into two scripts, one per step, built only
  from these functions.

## Non-goals

- No pipeline, runner, registry, or `run_step()`-style orchestration.
  The loop stays in user code (CLAUDE.md, user confirmed).
- No ML-ready training export. That is the job of the final PPI part,
  after G. Part A's job is clean, traceable tables and tool inputs.
- No parallelisation, sharding, restarts, or merge (memory:
  no-parallelisation-yet).
- No change to the canonicalisation or curation rules themselves.

## Anti-goals

- A step whose output needs the previous step's in-memory objects to
  make sense. Everything the next step needs must be on disk.
- Provenance the user has to assemble by hand from parts.
- One large JSON that embeds every exclusion and every entry (today's
  `DatasetManifest` shape). It does not scale to the full PDB.

## Evidence the mmCIF hand-off works

Checked on 2026-11-04 against the 100 fixtures in
`examples/ppi/v2.1/entries.txt`: for each entry, canonicalise →
`structure_to_mmcif` → `mmcif_to_structure` → `curate_structure` gives
the same kept chains, exclusion records, and atom counts as curating
the in-memory canonical structure (99 entries, 0 differences). Curation
reads its metadata (method, resolution, R-factors) from raw categories,
and those survive the round trip.

Not carried by mmCIF: which chains are symmetry copies
(`CanonicalMappings.assembly_mapping`, `AssemblyChainCopy`). That goes
in the canonicalisation provenance rows (below).

Side finding: `22AG` now raises `ValueError: Validation of
canonicalised structure failed`. The 2026-10-28 run of part A reported
0 failures. Not part of this work; worth looking at separately.

## The two steps

### Step 1: canonicalise

Input: mmCIF paths and a `canonicalisationPolicy`.

Per entry: `mmcif_to_structure` → `canonicalise_structure` → write the
canonical mmCIF.

Output folder:

| File | One per | Content |
|---|---|---|
| `structures/<mid>/<entry>.cif` | canonical entry | canonical structure (mirror layout) |
| `canonicalisation.json` | step run | step summary (below) |
| `entries.parquet` | canonical entry | per-entry provenance row (below) |
| `chain_copies.parquet` | symmetry-copy chain | `entry_id`, `chain_id`, `source_chain_id`, `operator_id` |
| `failures.jsonl` | input that raised | `path`, `error_type`, `message` |

### Step 2: curate

Input: canonical structures (read from step 1's folder, or passed in
memory when both steps run in one loop) and a `DatasetCurationPolicy`.

Per entry: `collect_metadata` → `curate_structure` → write the curated
structure and rows.

Output folder:

| File | One per | Content |
|---|---|---|
| `structures/<mid>/<entry>.cif` | kept entry | curated structure (today's `biounits/`) |
| `curation.json` | step run | step summary |
| `entries.parquet` | kept entry | entry row: structure path, method, resolution, chain count, per-entry provenance |
| `chains.parquet` | kept chain | `ChainRecord` + completeness + method/resolution + structure path + `copy_of` |
| `exclusions.parquet` | dropped entry or chain | `ExclusionRecord` |
| `failures.jsonl` | input that raised | as above |

`copy_of` is the source chain id for symmetry copies (from step 1's
`chain_copies.parquet` or in-memory mappings), null otherwise.

### Step summary file (both steps)

A small JSON: the step name, the full policy as it ran, Pandora and
gemmi versions, the input it read (path list file or previous step's
folder plus that step's summary), the output folder, start/end
timestamps, and counts (in, kept, failed, excluded per reason code).
Step 2's summary refers to step 1's summary by path rather than
embedding it (my call), so the chain of recipes stays a chain of small
files. When both steps run in one loop, step 2 refers to the step 1
summary written in the same run.
It replaces `counts.tsv`, which can be dropped or derived from it.

The summary is the recipe: rerunning each step with its summary's
policy and input rebuilds the folder (best-effort, like
`reproduce_dataset`).

### Per-entry provenance rows

Step 1: the `canonicalisationProvenance` of each entry (policy
id/version, transforms, report). Step 2: the `CurationProvenance` of
each entry (policy id/version, `dropped_raw_categories`). Each row is
flattened into a table row keyed by `entry_id`.

## What Pandora gains

These are the user-facing pieces the step scripts need. The names are my
calls (the builder may rename them); the split into small public
functions follows CLAUDE.md.

1. **Mirror-layout paths**: one function mapping `(root, entry_id)` to
   `root/<mid>/<entry>.cif`, and one listing the entry files under such a
   root. Today this is hand-written in parts A and C.
2. *(Dropped in planning: `list_mirror_files` plus `mmcif_to_structure`
   covers it; see the plan.)* **Read a step's structures back**: iterate the canonical mmCIF files
   of a step folder as `Structure`s (a thin wrapper over
   `mmcif_to_structure` and item 1), with parse failures reported, not
   raised. Without it, step 2 run on its own needs extra glue.
3. **Step summary model + writer**: a pydantic model in
   `schemas/provenance.py` for the step summary, a function that builds
   it (fills versions and timestamps), and a writer. The policy is
   embedded in full.
4. **Per-entry provenance rows**: functions that turn a
   `canonicalisationProvenance` / `CurationProvenance` plus `entry_id`
   into a flat row model, so `write_records` can write them.
5. **Chain copy rows**: a function returning `AssemblyChainCopy` rows
   (with `entry_id`) from `CanonicalMappings`.
6. **Dataset chain row**: a typed model (in `schemas/dataset.py`) and a
   function that builds it from a curated structure and its metadata:
   `ChainRecord` fields + `ChainCompleteness` + method/resolution +
   `copy_of`. This replaces the loose dict in today's script.
7. **Curation provenance in `ProvenanceBundle`**: the bundle has a
   `canonicalisation` slot but none for curation. Add one.

Out of scope for Pandora (stays in the script): the loop, the
`try/except` around each entry, and the choice of which files to write.

Each new public function gets a docstring, an `__all__` entry, and a
test that uses fixtures only. `docs/contributing/architecture.md` and
the ERD are updated (CLAUDE.md review checklist).

## PPI v2.1 scripts

`a_prepare_entries.py` becomes two scripts:

- `canonicalise.py`: step 1, writes `output/canonicalisation/`.
- `curate.py`: step 2, reads `output/canonicalisation/`, writes `output/curation/`, plus
  the tool inputs part B and C need:
  - `chains.fasta`: kept chains, **excluding symmetry copies** (same
    sequence as their source, so they only add MMseqs2 work);
  - `chains/<mid>/<entry>_<chain>.cif`: per-chain files for Foldseek,
    also excluding copies.

Parts B and C switch their inputs from `output/a/` to `output/curation/`:
`chains.parquet` instead of `chains.jsonl`, and part B skips rows with
`copy_of` set. Part B's containment check uses cluster sets, so dropping
copies does not change its result.

Script names and output folders are user decisions (2026-11-04).

Each script clears its own output folder at the start, so stale
structure files from earlier runs cannot survive (review finding 1).

## Constraints

- Parquet needs the `export` extra (pandas, pyarrow). The step scripts
  require it; Pandora functions import it lazily as today.
- No new dependencies.
- Inputs are never mutated; new objects come from `model_copy`.
- Tests need no network and no external binaries.
- `examples/` is untracked and has no pytest (memory:
  tests-only-pandora). Only the new Pandora functions get tests.

## Calls made by Claude, to check

- Tables as Parquet rather than JSONL; `failures` stay JSONL (small,
  written while the loop runs, readable after a crash).
- Symmetry copies stay in `chains.parquet` (flagged by `copy_of`) but
  are left out of `chains.fasta` and the per-chain mmCIF files.
- `counts.tsv` is replaced by counts in the step summary.
- Output file and folder names use the entry id as the input file
  spells it, not uppercased (user decision, 2026-11-04). Pandora does
  not uppercase anything: `Structure.entry_id` comes from the file's
  `_entry.id`, which is `10KZ` in `0k/10kz.cif`. My reading of the
  decision: paths follow the input file name
  (`structures/0k/10kz.cif`, `chains/0k/10kz_A.cif`, FASTA id
  `10kz_A`); the `entry_id` column keeps the parsed value. Parts B
  and C build chain paths from the same rule.
- The function names in "What Pandora gains".

## Left to the builder

- Exact model and function names, and which module each lands in,
  following the "Where to look" table in CLAUDE.md.
- Whether item 2 yields `(path, Structure | None, diagnostics)` or
  similar.
- The column order of the tables.
