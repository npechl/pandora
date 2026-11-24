# CLI namespace and command interface

Date: 2026-11-19, revised 2026-11-24 after a command-by-command review.
Related: `pandora/cli/README.md` (current I/O convention),
`examples/ppi/v2.1/src/` (the main existing user of the CLI).

Decisions marked **(user)** came from the user; **(my call)** were
proposed by Claude and accepted, or left to Claude.

## Why

The `pandora` command grew one command at a time. It now has 16
top-level commands that mix verbs (`curate`, `cluster`), nouns
(`similarity`) and method names (`rebalance`). Flags differ from command
to command for the same kind of input, and some commands can only read
the output of one specific other command. New users can't guess where a
command lives or what it accepts (`rebalance` was hard to understand even
for the maintainer).

The goal is a namespace that is **clear, easy to extend, and made of
commands that each work as a separate tool (user)**. The model is MMseqs2:
every module reads data in a known format and writes data in a known
format, and does not care which module (or which other tool) produced its
input **(user)**.

## Goals

1. Group commands by the kind of data they write, under 7 noun groups.
2. Name commands after what they produce or the method they use;
   rename the ones whose names don't say that.
3. One input/output convention for every command: positional main input
   and output, named flags for other inputs, a self-describing output
   folder.
4. Every input accepts either a Pandora output folder or the bare file
   **(user)**, with one stated limit for structures (see "Two kinds of
   structures input").
5. Remove hidden dependencies between commands: no command finds an
   input by following `previous_summary`; every input is passed
   explicitly.
6. Shared flags (hit filters, `--table-format`, `--overwrite`) defined
   once and behaving the same everywhere.
7. Remove commands, flags and outputs that nothing uses or that work
   against the toolkit's promises (plain `dedup`, `cluster --pairs`,
   `partition --no-keep-similar-items`) **(user)**.

## Non-goals

- `pandora inspect DIR` (describe any output folder). Wanted later, not
  now **(user)**.
- Moving library code to match the new groups (for example the `pairs`
  logic out of `pandora/similarity/rebalance.py`). Later **(user)**.
- A dataset manifest on the CLI. The `manifest` group is removed for now
  and revisited later; the library code in `pandora/provenance/` stays
  **(user)**.
- New functionality, beyond the small input changes listed under each
  group below. The commands marked "later" in the tree are examples of
  where future work fits, not part of this design.

## Anti-goals

- **No use-case groups** (`pandora ppi`, `pandora antibody`). Recipes
  stay in `examples/` **(user)**.
- **No orchestrators or composite workflows** (MMseqs2's `easy-*`).
  CLAUDE.md rules them out.
- **No command that only accepts the output of one specific other
  command.** If a command needs data, it takes it as an input in a
  documented format.

## Constraints

- argparse only (CLAUDE.md), one CLI module per group with a
  `register()` function, as today.
- Breaking every existing invocation is accepted **(user)**. No
  deprecated aliases for the old names **(my call: the user accepted the
  break, and aliases add code for a project with few external users)**.
  This is a minor version bump: 0.7.x to 0.8.0.
- Tests use local fixtures only (CLAUDE.md).

## Approaches considered

1. **Use-case groups** (`pandora ppi ...`). Rejected: it conflicts with
   "toolkit, not a pipeline", hides reusable commands
   (`interface-gaps`, `filter`) from non-PPI users, and invites
   orchestrators.
2. **Object groups** (`pandora interfaces ...`). Rejected: interfaces
   are handled in many places (`clusters build --items interface`,
   `similarity interface-overlap`, annotation layers), so the group
   would split `annotations` and `similarity` and still not hold
   everything about interfaces.
3. **Noun groups named after the data written (chosen, user).** Each
   group is a data type; a new command goes in the group of the data it
   writes.

How far the change reaches:

- A. Command names only. Rejected: flags stay inconsistent.
- **B. Names plus one input/output convention (chosen, user).**
- C. B plus `pandora inspect`. Deferred (non-goal).

What counts as a valid input:

- A. Only a Pandora output folder.
- **B. A folder or the bare file (chosen, user).** Each command must be
  usable as a separate tool.
- C. A folder with `summary.json`; outside data only via `import`.
  Rejected: least independent.

Argument style:

- A. Named flags for everything.
- B. Fully positional, like MMseqs2. Rejected: `pairs assemble` would
  need six positional arguments in a fixed order.
- **C. Hybrid (chosen, user).** The main input(s) and the output are
  positional; every other input is a named flag.

## Design

### Rules for placing and naming commands (user)

1. A **group** is named after the kind of data its commands write, not
   what they read and not a use case. (`drop-contained` reads clusters
   but writes structures, so it is in `structures`.) A group also holds
   data picked directly out of its own type: representatives are picked
   out of clusters, so `representatives` is in `clusters` **(user)**.
2. A **command** is named after what it produces or the method it uses,
   preferably a noun. A few general verbs mean the same thing in every
   group: `import` (bring outside data in), `filter`, `convert`.
3. Variants that are **separate methods with different flags** are
   separate commands; variants that take the same flags share one
   command with a choice flag. This removes the "flag X only applies
   with Y" checks in `annotate compute --layer`, `similarity search
   --engine`, `dedup --contained`, `ingest --ids`, `representatives
   --items` and `partition --stratify`. Accepted exceptions, where the
   level of the data decides which input is needed **(my call)**:
   `similarity import` (`--structures` or `--contacts`, by the item
   level in `--description`) and `clusters build` (`--contacts` and
   `--chain-representatives` only for interface items).
4. Every group that writes layers has an `import` command, so outside
   tools are first-class inputs.
5. One CLI module per group: `pandora/cli/<group>.py` with a
   `register()`.
6. No use-case groups.

### The tree

```
pandora entries      fetch | scan                       (later: search)
pandora structures   canonicalise | curate | drop-contained | convert
pandora annotations  contacts | ligand-contacts | contact-counts |
                     interface-contacts | dimer-length | interface-gaps |
                     filter | import
pandora similarity   mmseqs | foldseek | interface-overlap | import
pandora clusters     build | representatives
pandora splits       assign                             (later: check)
pandora pairs        positives | negatives | assemble
```

`canonicalise` and `curate` keep their names: the same words name the
library packages, the policies (`canonicalisationPolicy`,
`DatasetCurationPolicy`) and the provenance step names, and `curate` does
more than drop rows, so the general verb `filter` would mislead
**(user)**. Their help text says what they do.

### Old to new

| Today | New |
|---|---|
| `fetch` | `entries fetch` |
| `ingest` | `entries scan` |
| `canonicalise` | `structures canonicalise` |
| `curate` | `structures curate` |
| `dedup` | removed: `canonicalise` and `curate` already drop repeated entry IDs **(user)** |
| `dedup --contained` | `structures drop-contained` |
| `export` | `structures convert` |
| `annotate compute --layer chain_interfaces` | `annotations contacts` |
| `annotate compute --layer ligand_contacts` | `annotations ligand-contacts` |
| `annotate count-contacts` | `annotations contact-counts` |
| `annotate label-interfaces --min-residues-per-chain` | `annotations interface-contacts` **(user)** |
| `annotate label-interfaces --max-total-residues` | `annotations dimer-length` **(user)** |
| `annotate interface-gaps` | `annotations interface-gaps` |
| `annotate filter-contacts` | `annotations filter` |
| `annotate import` | `annotations import` |
| `similarity search --engine mmseqs2` | `similarity mmseqs` |
| `similarity search --engine foldseek` | `similarity foldseek` |
| `similarity interfaces` | `similarity interface-overlap` |
| `similarity import` | `similarity import` |
| `cluster` | `clusters build` |
| `cluster --pairs` | removed: nothing reads `cluster_pairs.json` **(my call, not objected to)** |
| `representatives` | `clusters representatives` |
| `partition` | `splits assign` |
| `partition --no-keep-similar-items` | removed: `splits assign` is always leakage-safe **(user)** |
| `rebalance` | `pairs positives` (builds the dimer groups itself) **(user)** |
| `negatives` | `pairs negatives` |
| `assemble` | `pairs assemble` |
| `manifest` | removed for now **(user)** |
| `reproduce` | removed for now **(user)** |

`cluster`, `representatives` and `partition` wrote one file named by
`--output`. Their new commands write an `OUT` folder like every other
command (see Outputs).

### Argument convention (user)

Every command reads as:

```
pandora <group> <command> MAIN_INPUT... OUT [options]
```

- **Positional:** what the command processes, then where the result
  goes.
- **Named flags:** every other input, and every parameter.
- **Input flags are named after the data type** they take:
  `--structures`, `--contacts`, `--labels`, `--similarity`,
  `--clusters`, `--representatives`, `--splits`, `--negatives`. When a
  command takes two inputs of the same type, the second gets a prefix
  that says its role (`--chain-clusters`) **(my call)**.

"STRUCTURES (tables)" below means a structures folder with its tables
(see "Two kinds of structures input").

| Command | Positional | Named inputs and main options |
|---|---|---|
| `entries fetch` | `ID... OUT` | `--ids FILE`, `--provider`, `--no-cache`, `--allow-partial` |
| `entries scan` | `ROOT OUT` | `--glob` or `--path-template`, `--ids`, `--source-label`, `--release-date` |
| `structures canonicalise` | `STRUCTURES OUT` | `--policy` |
| `structures curate` | `STRUCTURES OUT` | `--policy` |
| `structures drop-contained` | `STRUCTURES (tables) OUT` | `--clusters`, `--policy`, `--disable` |
| `structures convert` | `FILE FILE` | — |
| `annotations contacts` | `STRUCTURES OUT` | `--distance-cutoff`, `--atom-set`, `--polymer-types` |
| `annotations ligand-contacts` | `STRUCTURES OUT` | `--distance-cutoff`, `--include-waters` |
| `annotations contact-counts` | `CONTACTS OUT` | `--level` |
| `annotations interface-contacts` | `CONTACTS OUT` | `--min-residues-per-chain` |
| `annotations dimer-length` | `CONTACTS OUT` | `--structures` (tables), `--max-total-residues` |
| `annotations interface-gaps` | `CONTACTS OUT` | `--structures` |
| `annotations filter` | `CONTACTS OUT` | `--labels` (repeatable), `--drop` (repeatable) |
| `annotations import` | `TABLE OUT` | `--description`, `--structures`, `--source-files` |
| `similarity mmseqs` | `STRUCTURES (tables) OUT` | `--chains`, `--sensitivity`, `--max-seqs`, `--bin` |
| `similarity foldseek` | `STRUCTURES (tables) OUT` | `--chains`, `--chain-files`, `--sensitivity`, `--alignment-type`, `--max-seqs`, `--exhaustive-search`, `--bin` |
| `similarity interface-overlap` | `CONTACTS OUT` | `--chain-clusters`, `--structures` (tables), `--min-identity` |
| `similarity import` | `TABLE OUT` | `--description`, `--structures` (tables) or `--contacts`, `--source-files` |
| `clusters build` | `SIMILARITY... OUT` | hit filters, `--similarity-items`, `--items`, `--structures` (tables), `--contacts`, `--chain-representatives` |
| `clusters representatives` | `CLUSTERS OUT` | `--structures` (tables), `--contacts`, `--similarity` + hit filters, `--policy` |
| `splits assign` | `CLUSTERS OUT` | `--pct-train`, `--pct-val`, `--pct-test`, `--labels` or `--stratify` |
| `pairs positives` | `SIMILARITY... OUT` | hit filters, `--similarity-items`, `--contacts`, `--structures` (tables), `--chain-clusters`, `--chain-representatives`, `--policy` |
| `pairs negatives` | `POSITIVES OUT` | `--representatives`, `--similarity` + hit filters, `--contacts`, `--seed` |
| `pairs assemble` | `POSITIVES OUT` | `--negatives`, `--splits`, `--representatives`, `--similarity` + hit filters |

### Per-command changes

**entries**

- `entries fetch` into a folder that already holds a fetch **adds to
  it**: it downloads only the missing files, then merges `entries.txt`
  and `ingestion_provenance.json` with what is there **(user)**. It is
  the one command exempt from the `--overwrite` rule, because its folder
  is also the download cache. It writes `summary.json` like every other
  command; the counts cover this run and the folder total **(my call)**.
- `entries fetch` takes IDs positionally, from `--ids FILE`, or both.
  With neither it exits 1 **(my call)**.
- `entries scan --ids` works with `--glob` as well as
  `--path-template`. With `--glob` it keeps the files whose name, before
  its first `.`, is one of the IDs (case-insensitive) **(my call)**.

**structures**

- Plain `dedup` is removed **(user)**.
- `chain_copies` is part of the structures folder format, next to the
  `entries` and `chains` tables. `structures curate` reads it from its
  input folder when present; bare files have no copies **(my call)**.
  This is the input's own data, not a `previous_summary` lookup, so goal
  5 holds.
- `structures drop-contained` needs a structures folder with tables.

**annotations**

- `label-interfaces` is split in two **(user)**:
  - `interface-contacts` labels `few_contacts` when either chain has
    fewer than `--min-residues-per-chain` residues in contact with the
    other, else `ok`. It needs only the contact layer.
  - `dimer-length` labels `too_large` when the two chains' sequence
    lengths add up to more than `--max-total-residues`, else `ok`. It
    reads `seqres_length` from the `chains` table of `--structures`.

  Both keep today's flag names and label names **(my call)**.
- `annotations filter` takes `--labels` repeatedly and drops an
  interface when any of the layers gives it a `--drop` label, so the PPI
  example filters once instead of twice **(user)**. The "no interface
  has this label" typo check runs over all the given layers together.

**similarity**

- `--chains` on `mmseqs` and `foldseek` accepts a plain ID list (one per
  line, `#` comments) or a representatives folder/file, whose
  representatives become the list. This removes the `jq` step in
  `examples/ppi/v2.1/src/03_structure_similarity.sh` **(my call)**.

**clusters**

- `clusters build` loses `--pairs` and `cluster_pairs.json`.
- `--items` and `--similarity-items` on `clusters build` mean different
  things: `--similarity-items` is the level of the input tables,
  `--items` the level clustered at (a chain layer can be rolled up to
  entries or to interfaces) **(my call)**.
- `clusters representatives` loses `--items`. The level comes from the
  flags: with `--contacts` it picks interfaces, without it chains. Both
  levels need `--structures` (tables), for chain values and for
  symmetry copies **(my call)**.

**splits**

- `--no-keep-similar-items` is removed **(user)**.
- `--labels L` places clusters so every label reaches every split.
  `--stratify L` splits each label's clusters separately (today's
  `--labels L --stratify`). The two are mutually exclusive in argparse,
  so the custom check goes **(my call)**.

**pairs**

- `pairs positives` builds the dimer groups itself as the connected
  components of the dimer graph from its similarity inputs, instead of
  taking `--clusters` plus the same layers again **(user)**. It writes
  the groups (limited to the positives) as `clusters.json`, so
  `splits assign` can take the positives folder directly, and the
  homomer/heteromer labels as the `kinds/` layer. The separate
  `clusters build` step for dimer groups (steps 5-6 of the PPI example)
  goes away.
- **Check before relying on it:** a test runs `clusters build --items
  interface` and `pairs positives` on the same fixture inputs and
  asserts the groups are identical.
- `pairs assemble` keeps warning when `--representatives` or
  `--similarity` differ from what `pairs negatives` used.

### Inputs: folder or bare file (user)

Each data type has one file format and one fixed file name. An input
accepts the folder that holds that file, or the file itself. The table
also defines what each command writes **(my call)**:

| Data type | File inside the folder | Bare-file forms accepted |
|---|---|---|
| entries | `entries.txt` | the `.txt` list |
| structures | `structures/<mid>/<entry>.cif` + `entries`, `chains`, `chain_copies` tables | at the entry points only: a flat or mirror folder of `.cif`/`.cif.gz` files, an `entries.txt` list, a single `.cif`/`.cif.gz` file |
| layer (contacts, labels, values) | `layer.<fmt>` | a `.parquet` or `.jsonl` table |
| similarity | `layer.parquet` | a `.parquet` table |
| clusters | `clusters.json` | the `.json` file |
| representatives | `representatives.json` | the `.json` file |
| splits | `splits.json` | the `.json` file |
| positives | `positives.<fmt>` (+ `clusters.json`, `kinds/`) | the table |
| negatives | `negatives.<fmt>` | the table |
| pairs | `pairs.<fmt>` | the table |

- A bare file is checked by its columns or schema; a missing column
  exits with one line naming it, as `import` does today.
- `summary.json`, when present, is read only for consistency warnings
  and for a similarity folder's item level. It is never required.
- A command that needs the structures behind a layer takes
  `--structures`; nothing follows `previous_summary` to find it
  (goal 5). The argument table lists every such command.
- A similarity table doesn't record whether its IDs are entries,
  chains or interfaces, and the IDs can't tell (entry IDs may contain
  `_`). Every command that takes similarity input (`clusters build`,
  `clusters representatives`, `pairs positives`, `pairs negatives`,
  `pairs assemble`) gets `--similarity-items {entry,chain,interface}`
  **(user)**. It is repeatable like the hit-filter flags: once for all
  similarity inputs, or once per input in order **(my call)**. It is
  required for a bare table and for a folder without `summary.json`.
  For a folder with one, the level comes from it, and a
  `--similarity-items` that disagrees exits **(my call)**.
- `annotations filter` today exits when the contacts and labels trace to
  different structures folders. With no tracing, it compares the
  `input` recorded in the `summary.json` files when all exist and warns
  if they differ; otherwise it skips the check **(my call)**.
- Structures folders are discovered with `.cif.gz` files as well as
  `*.cif` **(my call: fixes the inconsistency where the CLI README says
  gzipped inputs are read but folder discovery ignores them.)**

#### Two kinds of structures input (user)

Most commands after `curate` need the `entries` and `chains` tables
(sequences, lengths, resolution, copies), not only `.cif` files.

- **Entry points** accept any structures form, including bare `.cif`
  files and folders: `structures canonicalise`, `structures curate`,
  `annotations contacts`, `annotations ligand-contacts` and the
  `--structures` of `annotations interface-gaps`.
- **Every other structures input** needs a structures folder with
  tables. Without them the command exits 1 with
  `<group> <command>: <path> has no chains table; run
  \`pandora structures curate\` first`.

Goal 4 holds for every other data type.

### Outputs (user)

- `OUT` is always a folder holding the data under its fixed name plus a
  `summary.json` (a `StepSummary`, as today).
- Today's per-step summary names (`canonicalisation.json`,
  `curation.json`, `dedup.json`, `layer.json`, `ingestion.json`,
  `rebalance.json`, `negatives.json`, `assemble.json`) and side files
  (`cluster_provenance.json`, `representatives_provenance.json`,
  `partition_provenance.json`) all become `summary.json`. The step name
  is already recorded inside **(my call)**. `entries fetch`, which writes
  no summary today, gains one.
- Extra outputs keep their current names inside the folder:
  `failures.jsonl`, `mappings.jsonl`, `exclusions`, `contained`,
  `chain_copies`, `unmatched.<fmt>`, `source/`, `dropped.<fmt>`,
  `kinds/`, `ingestion_provenance.json`. `pairs positives` writes
  `clusters.json` instead of `groups.json`.
- **Exception:** `structures convert` writes exactly the file named, and
  takes exactly one file, because it is an edge converter like
  `mmseqs convertalis`. It is outside goal 4 and the folder rule.

### Shared behaviour (my call)

- `--overwrite` on every command except `entries fetch` (see above);
  without it, a command refuses an `OUT` that already holds its
  outputs, and with it deletes only those (for `structures/`, the whole
  tree).
- `--table-format {parquet,jsonl}` on every command that writes a
  table; similarity layers stay Parquet only.
- The hit-filter flags (`--hit-filter`, `--min-score`,
  `--min-coverage`) and `--similarity-items` are defined once and
  attached to every command that takes similarity input.
- One result line on stdout; warnings on stderr.
- Exit status 1 with a one-line `<group> <command>: ...` message for bad
  input. argparse's own usage errors stay at 2. `entries fetch` without
  `--allow-partial` exits with that one-line message instead of a
  traceback.
- Every group and command has `help=` text, and every flag has `help=`
  text.

### What changes besides `pandora/cli/`

- `pandora/cli/README.md`: rewritten around the new tree and the
  input/output convention.
- `docs/` pages that show CLI calls (including removing the `manifest`
  and `reproduce` CLI sections from `docs/usage/cli.md` and the CLI
  mention in `docs/usage/provenance.md`), `examples/API_SURFACE.md`,
  and the docs site builds.
- `examples/ppi/v2.1/src/*.sh`: updated to the new commands; the `jq`
  step, the second `filter` pass and the dimer-group `clusters build`
  step go away.
- `tests/test_cli.py`: updated, plus tests for the new behaviour (bare
  file inputs, explicit `--structures`, `summary.json`, additive fetch,
  repeatable `--labels`, the `pairs positives` group check).
- `pandora/cli/_io.layer_origin` is deleted once no command calls it
  **(my call)**.

## Later

- `entries search`, `splits check`, `pandora inspect`.
- A general `manifest` group: record the `summary.json` of any step
  folder, so a PPI dataset can be described too. Needs a `steps` list in
  `DatasetManifest`.
- Moving library code to match the groups.

## Left to the builder

- How the CLI modules share the "folder or bare file" lookup and the
  shared flag blocks (a helper in `_io.py` is the obvious place).
- The order of work and how to split it into commits.
- Wording of help texts.
