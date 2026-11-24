# CLI

The `pandora` command groups its subcommands by the kind of data they write. Run `pandora <group> <command> -h` for a command's arguments, or `pandora -h` for the groups. `docs/usage/cli.md` shows one example per command.

```
pandora entries      fetch | scan
pandora structures   canonicalise | curate | drop-contained | convert
pandora annotations  contacts | ligand-contacts | contact-counts | interface-contacts |
                     dimer-length | interface-gaps | filter | import
pandora similarity   mmseqs | foldseek | interface-overlap | import
pandora clusters     build | representatives
pandora splits       assign
pandora pairs        positives | negatives | assemble
```

A group is named after the data its commands write, not what they read: `structures drop-contained` reads clusters but writes structures. Each group is one module (`entries.py`, `structures.py`, ...) with a `register()`; `app.py` assembles them and `_io.py` holds the shared file helpers. A new command goes in the group of the data it writes.

## Argument convention

Every command reads as:

```
pandora <group> <command> MAIN_INPUT... OUT [options]
```

- Positional: what the command processes, then where the result goes.
- Named flags: every other input and every parameter. Input flags are named after the data type they take (`--structures`, `--contacts`, `--labels`, `--similarity`, `--clusters`, `--representatives`, `--splits`, `--negatives`). A second input of the same type gets a prefix for its role (`--chain-clusters`, `--chain-representatives`).
- No command finds an input by following a `previous_summary`; every input is passed explicitly.

## Outputs

`OUT` is always a folder. It holds the data under a fixed name plus `summary.json`, a `StepSummary` recording the step name, policy and settings, Pandora and gemmi versions, input, counts and (for information only) `previous_summary`. Extra outputs keep fixed names inside the folder: `failures.jsonl`, `mappings.jsonl`, `exclusions`, `contained`, `chain_copies`, `unmatched.<fmt>`, `source/`, `dropped.<fmt>`, `kinds/`, `ingestion_provenance.json`.

`structures convert` is the one exception: it takes exactly one file and writes exactly the file you name, like an edge converter.

## Inputs: a folder or the bare file

Each data type has one file format and one fixed file name. An input accepts the folder that holds that file, or the file itself.

| Data type | File inside the folder | Bare file accepted |
|---|---|---|
| entries | `entries.txt` | the `.txt` list |
| structures | `structures/<mid>/<entry>.cif` + `entries`, `chains`, `chain_copies` tables | entry points only: a flat or PDB-layout folder of `.cif`/`.cif.gz`, an `entries.txt` list, one `.cif`/`.cif.gz` |
| layer (contacts, labels, values) | `layer.<fmt>` | a `.parquet` or `.jsonl` table |
| similarity | `layer.parquet` | a `.parquet` table |
| clusters | `clusters.json` | the `.json` file |
| representatives | `representatives.json` | the `.json` file |
| splits | `splits.json` | the `.json` file |
| positives | `positives.<fmt>` (+ `clusters.json`, `kinds/`) | the table |
| negatives | `negatives.<fmt>` | the table |
| pairs | `pairs.<fmt>` | the table |

A bare file is checked by its columns or schema; a missing column exits with one line naming it. `summary.json` is read only for consistency warnings and for a similarity folder's item level; it is never required.

A similarity table doesn't record whether its IDs are entries, chains or interfaces. Every command that takes similarity input takes `--similarity-items {entry,chain,interface}`, once for all inputs or once per input in order. It is required for a bare table and for a folder without `summary.json`; a folder with one gets the level from it, and a flag that disagrees exits.

### Two kinds of structures input

- Entry points accept any structures form: `structures canonicalise`, `structures curate`, `annotations contacts`, `annotations ligand-contacts` and the `--structures` of `annotations interface-gaps`. Structures folders are discovered with `.cif.gz` files as well as `.cif`.
- Every other structures input needs a structures folder with its `entries` and `chains` tables (the output of `structures curate`, or of `structures drop-contained`). Without them the command exits 1 with `<group> <command>: <path> has no chains table; run \`pandora structures curate\` first`.

`entries scan OUT` and `entries fetch OUT` folders hold an `entries.txt`; pass that file (not the folder) as the structures input of an entry point.

## Shared flags and behaviour

- `--overwrite` on every command except `entries fetch`: without it a command refuses an `OUT` that already holds its outputs, with it deletes only those (for `structures/`, the whole tree). `entries fetch` adds to a folder that already holds a fetch, because the folder is also the download cache.
- `--table-format {parquet,jsonl}` on every command that writes a table (default Parquet; `jsonl` needs no `export` extra). Similarity layers are Parquet only. A table with no rows is not written.
- `--hit-filter`, `--min-score`, `--min-coverage` and `--similarity-items` are defined once and attached to every command that takes similarity input. Each can be given once for all inputs or once per input, in order.
- One result line on stdout; warnings go to stderr as `<group> <command>: warning: ...`.

## Exit codes

- 0: success.
- 1: bad input, with a one-line `<group> <command>: ...` message (missing column or table, an `OUT` that already holds outputs, an unknown label, ...). `entries fetch` without `--allow-partial` exits this way too.
- 2: argparse usage errors (unknown flag, missing positional argument).

## Scale

Commands handle one entry at a time so they scale to the whole PDB. Structures are parsed one at a time (`read_structures_dir`, or `structure_index` + `read_structure` to look an entry up by file name). Contact tables are streamed: `contact_rows` yields one `ResidueContactRow` at a time for the library functions that only keep per-interface totals, and `rows_by_entry` yields one entry's rows at a time, sorted by entry, for steps that need a whole entry. Tables as large as their input (`annotations contacts`, `import`, `filter`) are written row by row with `TableWriter`, which turns a column of mixed kinds into text; their rows come out grouped by entry, in sorted entry order. Don't load a whole contact table with `read_records`; at the full PDB it doesn't fit in memory.
