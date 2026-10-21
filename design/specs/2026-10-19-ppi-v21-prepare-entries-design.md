# PPI v2.1, part A: prepare entries — design

> **Update 2026-10-21 (user):** no parallelisation for now. Sharding,
> restart files (`started.txt`, `done.txt`), `manifest.json`, `--merge`
> and the Slurm script were removed. The script runs one plain loop and
> writes every output directly under `--out` (overwritten each run).
> Sections below that mention shards, restarts, merge or Codon are
> superseded.

## Context

The PPI v2.1 benchmark (`examples/ppi/benchmarks/ppi-v2.1-benchmark.md`)
is built in parts that hand files to each other, so each part can be
rerun alone and run as its own Codon job:

| Part | v2.1 steps | Summary |
|---|---|---|
| **A. Prepare entries** | Step 1 filters 3–9 | This spec |
| B. Sequence clusters | Step 2, Step 1 filter 1 | MMseqs2 90/90, representatives, biounit containment |
| C. Structure similarity | Step 3 | Foldseek on representatives, TM ≥ 0.5, cov ≥ 0.8 |
| D. Dimers | Step 4 | Import COCαDA contacts, any-contact and enough-contacts lists, size cap, interface gaps, dedup |
| E. Dimer clusters | Steps 5–6 | Dimers linked through chain TM hits, connected components |
| F. Positives | Steps 7, 9 | Pick positives per cluster, attach contacts |
| G. Negatives | Step 8 | Sample negatives, excluding any-contact pairs |

Scope decided with the user (2026-10-19): **no** PDB-REDO, PDBFixer,
DSSP or pH; the Step 1 quality filters **are** included. Interface
similarity in part E waits until v2.1 picks a method. Contacts come
from COCαDA, not Pandora (user, 2026-10-19): part A writes the files
COCαDA runs on, and part D imports its output.

## Goal

One script that turns a list of PDB entries into filtered, protein-only
biological-assembly chains, and reports how many entries and chains
each filter removed. It must run on a laptop
over the dev fixtures and on Codon over the whole PDB as a Slurm array.

## Non-goals

- Anything from parts B–G (clustering, dimers, negatives).
- Contacts. No Pandora contact computation and no COCαDA run; the
  user's colleagues run COCαDA on part A's biounit files.
- Biounit containment (Step 1 filter 1): it needs cross-entry sequence
  clusters, so it lives in part B.
- Poly-X tail trimming, PDB-REDO models, PDBFixer, DSSP, pH, FSC.
- A batch runner inside Pandora. Parallelism is the Slurm array; the
  script handles one shard.

## Anti-goals

- A crash on one bad entry kills a 10k-entry shard.
- Filter counts that can't be traced to a reason code.
- Hard-coded thresholds in the script: every v2.1 value lives in a
  policy YAML.

## Layout

```
examples/ppi/v2.1/
  README.md                       how to run A locally and on Codon
  policies/canonicalisation.yaml
  policies/curation.yaml
  a_prepare_entries.py
  slurm/a_prepare_entries.sh      example array job
```

Parts B–G add their own `b_*.py` … scripts here later.

## Inputs

```
python examples/ppi/v2.1/a_prepare_entries.py \
  --entries entries.txt            # one PDB id per line
  --mmcif-dir /path/to/mirror      # local files; fetched if absent (see below)
  --out out/a                      # output root
  --shard 3 --n-shards 200         # optional; default 0 of 1
```

- Shard `i` of `n` takes entries `i, i+n, i+2n, …` of the sorted list,
  so shards are balanced and deterministic.
- **Source of files:** if `<mmcif-dir>/<id>.cif` or `<id>.cif.gz`
  exists, use `ingest_local_mmcif`; else `fetch_mmcif` from RCSB into
  `--mmcif-dir`. On Codon, point `--mmcif-dir` at a local PDB mirror if
  EBI has one; otherwise pre-fetch (fetching from 200 array tasks at
  once is rude to RCSB).
- The input is the asymmetric-unit mmCIF. Assembly 1 is built by
  canonicalisation, not downloaded as `-assembly1.cif`, so assembly
  expansion and chain ids stay under Pandora's control and recorded in
  provenance.

## Per-entry flow

For each entry, in order:

1. **Ingest** (local or fetch) → path + `IngestionProvenance`.
2. **Parse** `mmcif_to_structure`.
3. **Canonicalise** with `policies/canonicalisation.yaml`:
   - `assembly_rules`: `standardize_biological_assembly`, source
     `first` (see Library change).
   - `identifier_rules`: `preserve` chain ids and `preserve` residue
     numbering, so `label_seq_id` still lines up with SEQRES
     (`chain_completeness` needs this).
   - `entity_rules`: `preserve` (not `merge_equivalent_entities`; see the
     deferred bug where completeness looks up SEQRES by the merged
     entity id).
   - `modified_residue_rules`: `map_to_parent` (MSE → MET).
   - `ligand_rules`: `filter`, drop waters, ions and ligands.
   - `altloc_rules`: `select_best_occupancy`.
4. **Metadata**: `collect_metadata` → `MetadataRecord` (quality,
   method, resolution) for curation.
5. **Curate** with `policies/curation.yaml`:

   | v2.1 filter | Field | Value |
   |---|---|---|
   | 3. Resolution | `max_resolution_by_method` | X-ray 2.5, EM 2.0 |
   | (implied) no NMR | `include_experimental_methods` | X-ray, EM |
   | 4. \|Rfree − Rwork\| | `max_r_free_gap` | 0.07 |
   | 5. Rfree, Rsym | `max_r_free`, `max_r_sym` | 0.35, 0.10 |
   | 6. Missing residues | `missing_residue_definition`, `max_missing_tail_fraction`, `max_missing_middle_fraction` | `incomplete_backbone`, 0.30, 0.10 |
   | 7. Non-standard residues | `exclude_nonstandard_residues` | true |
   | 8. No/mismatching SEQRES | covered by 6: a chain with no SEQRES or `SEQRES_MISMATCH` gets `NO_SEQRES` | — |
   | 9. ≥ 1 chain left | built in (`NO_CHAINS_LEFT`) | — |
   | Size guard | `max_atoms` | value chosen from the spike; see Codon |

6. **Protein only**: from the curated structure keep chains whose
   entity's `EntityPolyRecord.type` is `polypeptide(L)`. Non-protein polymers
   (DNA/RNA) are dropped from the outputs, not from the structure;
   the drop is counted as `NOT_PROTEIN`. An entry left with no protein
   chain is counted as `NO_PROTEIN_CHAINS` and writes nothing else.
7. **Write** for the protein chains (see Outputs):
   - `extract_chain_records` → chain records (sequence, lengths);
   - `chain_completeness` → missing N/C/middle counts (parts B and D
     need them: representative choice, interface gaps);
   - `export_chain_mmcif` → one file per chain (Foldseek input, part C);
   - `structure_to_mmcif` on the curated biounit restricted to its
     protein chains → one file per entry (COCαDA input). The script
     filters atoms and connections with `model_copy` and drops `raw`,
     the same way `export_chain_mmcif` does for one chain. Its chain ids
     are the ones in `chains.jsonl`, so part D can join COCαDA's
     contacts to Pandora's chains without remapping.

**Errors.** Each entry runs inside a `try` that catches named
exceptions only (gemmi/parse errors, `OSError`, `httpx.HTTPError`,
`ValueError`) and writes one line to `failures.jsonl`
(`entry_id`, stage, exception type, message). Unknown exceptions still
crash: they are bugs.

## Outputs

Per shard, under `<out>/shards/<i>-of-<n>/`:

| File | One row per | Content |
|---|---|---|
| `chains.jsonl` | kept protein chain | `ChainRecord` + `ChainCompleteness` fields + method and resolution |
| `chains.fasta` | kept protein chain | `>{entry}_{chain}` + sequence (MMseqs2 input, part B) |
| `exclusions.jsonl` | dropped entry or chain | `ExclusionRecord` fields; the script's own `NOT_PROTEIN` / `NO_PROTEIN_CHAINS` rows use the same fields but are plain dicts, since they aren't `ExclusionReason` values |
| `failures.jsonl` | entry that raised | see Errors |
| `done.txt` | finished entry | entry id |
| `manifest.json` | shard | policies, Pandora/gemmi versions, args, timestamps, counts |

Structure files are sharded by the entry id's middle two characters,
like the PDB archive, so no directory holds 1M files:

- chains: `<out>/chains/<entry[1:3]>/<entry>_<chain>.cif`;
- biounits: `<out>/biounits/<entry[1:3]>/<entry>.cif`.

**Restarts.** On start the script reads `done.txt` and `failures.jsonl`
and skips those entries, so a killed array task is resubmitted as is.
Rows for one entry are appended only after the whole entry finished,
so a kill mid-entry leaves no partial rows.

**Merge.** `a_prepare_entries.py --merge --out out/a` concatenates the
shards into `<out>/chains.jsonl`, `chains.fasta`, `exclusions.jsonl`, `failures.jsonl`, and writes `<out>/counts.tsv`:
entries in, entries kept, chains kept, and one row per reason code
(entry- and chain-level separately). This table answers v2.1 risk 2
("count what survives each filter").

## Codon

- Example `slurm/a_prepare_entries.sh`: `--array=0-199`, one shard per
  task, `uv run` from the repo checkout, `--mem` sized from the spike
  (6.5 GB to parse 3j3q, 2.44M atoms; ask 16 GB and let `max_atoms`
  stop the rest of the flow for the giants).
- Entries that exceed memory show up as killed tasks; the restart rule
  plus `--entries` of the leftovers handles them on a high-memory
  partition.
- The README says how to point `--mmcif-dir` at a mirror and how to
  check progress (`wc -l shards/*/done.txt`).

## Library change

`preferred_assembly_source: first` currently behaves like `author`
(author-determined, then software-determined, then first). Change it to
return the first listed assembly, which is assembly 1 in PDB files.
Update `docs/reference/policies.md` and add a test on a fixture whose
first assembly is not author-determined (or whose first and
author-determined assemblies differ).

No other library changes. Anything else part A needs exists.

## Tests

- `tests/test_ppi_v21.py`:
  - run part A on 3–5 dev fixtures (in-process, `--mmcif-dir
    datasets/dev/mmcif`, no network), check files exist, chain ids in
    `chains.fasta` match `chains.jsonl`, and each biounit file holds
    exactly that entry's chains from `chains.jsonl`;
  - a fixture that fails a known filter shows up in `exclusions.jsonl`
    with the expected reason code;
  - restart: run twice, second run processes nothing and outputs are
    unchanged;
  - two shards merged equal one unsharded run.
- The `first` assembly test above.

## Left to the builder

- Which dev fixtures to use for each test.
- The exact `max_atoms` value (from the spike numbers).
- Function boundaries inside the script.

## Decisions

Confirmed by the user (2026-10-19):

- No contacts in A; COCαDA provides them, imported in part D.
- Assembly 1 built from the asymmetric unit, not `-assembly1.cif`.
- Every chain mmCIF is exported in A, not only part B's representatives:
  it costs disk but avoids parsing every entry a second time.
- No poly-X tail trimming. Unobserved His-tags count as missing tails,
  so they push chains toward the 30% tail limit.
- `entity_rules: preserve`, avoiding the SEQRES lookup bug.

My calls:

- Protein-only is applied to outputs, not via a new canonicalisation
  rule.
- A writes one protein-only biounit mmCIF per entry as COCαDA's input,
  so COCαDA sees Pandora's chain ids.
