# Annotation layers from external tools — design

Roadmap sub-project 4 ("external tool outputs") from
`2026-10-07-contacts-design.md`, now scoped as a project of its own.
This document covers the whole project and its order. **Part 1 is
specified in detail.** Parts 2–6 are each brainstormed and specified
separately when their turn comes.

## Intent

Someone building a dataset uses several tools, and each one extracts a
different kind of information: COCaDA for interatomic contacts, MMseqs2
for sequence similarity, later FreeSASA, DSSP, PRODIGY and others. Each
writes its own file formats with its own chain and residue naming. Today
there is no way to bring that output into Pandora and use it alongside
the structures. The user does not know how to bridge COCaDA's output
with Pandora.

Pandora does not wrap these tools (exception: tools that are well known
and already wrapped, MMseqs2 and Foldseek). The user runs the tool on
Pandora's output structures. Pandora takes the result in, checks that
it refers to residues that exist, and records where it came from. The
user can then use it in later steps: curation, clustering, splitting,
export.

The key idea is **templates**: a small set of table shapes, each with a
meaning (residue contacts, per-residue values, chain labels, ...).
The user converts their tool's output into the template that fits.
Downstream steps work on templates, never on specific tools, so any
layer that fits a template can be used everywhere that template is
accepted.

## Goals (whole project)

- Every downstream use works on **any** layer of a fitting template, not
  on specific tools.
- Built-in and imported layers are interchangeable: Pandora's own
  annotations produce the same templates.
- A dataset built from several tools can be combined by joining layers
  on Pandora's chain and residue keys.

## Non-goals

- Wrapping or running external tools (beyond the existing MMseqs2 and
  Foldseek wrappers).
- Pandora-side parsers for each tool's native format. The user converts
  into a template. A tool-specific reader is added only when a format
  is common enough to be worth it, and then it emits a template table.
- Translating IDs for tools that were run on the **original** PDB files
  (through `mappings.jsonl`). Not in part 1. The supported path is
  running the tool on a Pandora step folder.

## Anti-goals

- A row silently attached to the wrong residue. Anything that doesn't
  match exactly, or matches more than one residue, is reported, never
  guessed.
- Pandora turning into a pipeline that drives tools.

## Order of work

1. **Foundation** (this document): templates, the layer folder, and
   `pandora` import of a user's table checked against a step folder.
2. **Built-in annotations produce templates.** Pandora's chain
   interfaces become a `residue_contacts` layer, and so on.
3. **Export joins layers** onto chain/residue/interface records, so the
   values reach the final dataset.
4. **Curation rules on layers**, e.g. drop chains with fewer than N
   contacts, or keep residues with a given label.
5. **Clustering from any contact layer.** Interface residues come from
   any `residue_contacts` layer, not only Pandora's own interfaces.
6. **Split balancing on a label** from any label layer (issue #95).

Export comes before curation because it is the quickest way to get a
combined COCaDA + MMseqs2 dataset out the other end.

The existing similarity import (`pandora similarity --precomputed-hits`)
is unchanged in part 1. It becomes the `similarity` template in a later
part.

## Part 1: foundation

### Templates

A template is a table shape: fixed **key columns** that say what a row
is about, plus required **value columns**. Any other columns in the
user's table are kept as-is and not interpreted.

| Template | One row per | Key columns | Required values |
|---|---|---|---|
| `residue_contacts` | residue pair | entry, side-1 residue, side-2 residue | none (`contact_type`, `distance` recognised if present) |
| `residue_values` | residue | entry, residue | `value` (number) |
| `residue_labels` | residue | entry, residue | `label` (string) |
| `chain_values` | chain | entry, chain | `value` (number) |
| `chain_labels` | chain | entry, chain | `label` (string) |
| `interface_labels` | chain pair | entry, chain 1, chain 2 | `label` (string) |
| `free_form` | any | at least `entry_id`; any chain/residue key columns present are checked | none |

`similarity` (item pairs) is listed for completeness, but it comes in a
later part, when the existing similarity import is moved over.

Builder's call: whether `residue_values` / `residue_labels` allow
several value columns in one table (e.g. `value_sasa`,
`value_plddt`), or one layer per value. Pick whichever keeps the
downstream parts simpler, and record it here.

**Free-form layers** are accepted. Their keys are checked like any
other layer, but nothing interprets their values: they are carried
along and, from part 3, joined in at export. If people keep importing
the same free-form shape, it should become a template.

### How rows identify chains and residues

Pandora names chains by `label_asym_id` and residues by
`label_seq_id`. These are mmCIF's own IDs, and `ChainRecord`,
`ResidueRecord` and `InterfaceRecord` already use them. Most tools use
the author IDs printed in PDB files: author chain, author residue
number, insertion code. Ligands and waters have no `label_seq_id` at
all.

The user's table may use **either** naming, per row:

- Pandora's: `chain_id`, `seq_id`
- Author's: `auth_chain_id`, `auth_seq_id`, `ins_code` (`ins_code`
  optional, empty means none)

For pairs, the same with suffixes `_1` / `_2` on each side
(`chain_id_1`, `auth_seq_id_2`, ...). Chain-level templates use only
the chain columns.

On import Pandora looks each row up in the structures of the step
folder and writes **both** namings in the stored layer table. From then
on, downstream steps use only Pandora's keys, and people reading the
table still see the familiar author numbers.

A row is **unmatched**, with a reason, when:
- its entry is not in the step folder;
- its chain or residue does not exist in that entry;
- it gives both namings and they disagree;
- it gives author IDs that point at more than one residue (this happens
  today when insertion codes are stripped, issue #125);
- a required column is missing or has the wrong type for the template.

**Which residues exist** (my call): a residue exists when it has
atoms in the step folder's structure. This is the same lookup
`datasets/records.py` uses for `ResidueRecord`. A row naming a residue
that is in the chain's sequence but has no coordinates is unmatched
with its own reason (`no_coordinates`), so sequence-based tools see
clearly what was left out instead of getting a generic "does not
exist". Downstream steps are structure-based, so they could not use
such rows anyway.

**Repeated keys** (my call). Parsing keeps one model per entry, so
repeats come only from the user's table:
- `residue_contacts`: several rows for the same residue pair are kept.
  One pair can legitimately have several contact types (e.g. a hydrogen
  bond and a hydrophobic contact).
- `free_form`: repeats are kept too, since its shape is unknown.
- every other template: the second and later rows for a key are
  unmatched with reason `duplicate_key`. The first row is kept. This
  makes the stored table one row per key, so part 3's export join can't
  multiply records.

**Pair order** (my call): sides are stored in the order the user gave.
Import does not reorder them or treat (A,B) and (B,A) as duplicates.
Downstream parts treat `residue_contacts` and `interface_labels` pairs
as unordered when they read them.

A table where a **required column is absent altogether**, or which
names an unknown template, is an error: the command stops before
writing anything.

How to identify a non-polymer residue (ligand, water) that has no
`seq_id`: author naming covers it. Builder's call whether the Pandora
naming also needs a way to address it; record the choice here.

### The command

```sh
pandora annotate import \
  --structures-dir curated/ \
  --table cocada_contacts.parquet \
  --description cocada.yaml \
  --source-files cocada_out/ \
  --output-dir annotations/cocada
```

- `--structures-dir`: the step folder the tool ran on. Any folder the
  other commands accept as structures works. If it holds a step summary
  (`curation.json`, `canonicalisation.json`), the layer's summary points
  to it.
- `--table`: the user's table, Parquet or JSON Lines (by extension),
  like Pandora's other tables.
- `--description`: a YAML the user writes once per layer:

  ```yaml
  name: cocada
  template: residue_contacts
  tool: COCaDA
  tool_version: "1.0"
  parameters:          # free key/value pairs, recorded as given
    cutoff: 5.0
  description: ""      # optional
  ```

  This is a pydantic model in `pandora/schemas/annotation.py`, loaded
  with `load_policy(path, model=...)`.
- `--source-files`: optional. A file or folder of the tool's original
  output, copied in untouched.
- `--output-dir`: where the layer folder is written. Refuses a folder
  that already holds a layer, unless `--overwrite` is given (same rule as
  `canonicalise` / `curate`).
- `--table-format`: as on the other step commands (Parquet default,
  `jsonl` without the `export` extra).

**`pandora annotate` becomes a group** (user's decision).
`pandora annotate` currently exists as a flat command (`--input-dir`,
`--layers`, `--pairwise`). It becomes `pandora annotate compute` with
the same options and behaviour, next to `pandora annotate import`. The
old form `pandora annotate --input-dir ...` stops working, with no
transition period (Pandora is 0.x, as in the CLI step-outputs change).
Update `pandora/cli/README.md`, `docs/usage/cli.md`, and any docs,
examples or tests that call it.

### The layer folder

```
annotations/cocada/
  layer.<fmt>        matched rows, both namings filled in, all user columns kept
  unmatched.<fmt>    rows that didn't match, with a `reason` column
  source/            the tool's original files (only with --source-files)
  layer.json         summary
```

- A table with no rows is not written (as in the other step folders).
- `layer.json` is a `StepSummary` (`step: "annotation_import"`), with:
  the description (name, template, tool, version, parameters) as its
  `policy`; `input` = the table path; `previous_summary` = the step
  folder's summary if present; `counts` = rows in, matched, unmatched
  (and unmatched per reason).
- Unmatched rows do not stop the import. The import carries on, and the
  counts are in the summary. No strict mode for now (my call). Add one
  if a real case needs it.

### Library side

The command is a thin wrapper, as with the other commands. The library
gets:

- the template and description models in `pandora/schemas/annotation.py`;
- a plain function in `pandora/annotations/` that takes the user's rows,
  the template, and the structures, and returns the matched rows with
  both namings plus the unmatched rows with reasons. It does no file
  I/O.

How this relates to the existing `AnnotationLayer` model (free-form
`data` dict per entry): part 1 does not change `AnnotationLayer`.
Whether imported layers are also exposed as `AnnotationLayer`s, or
`AnnotationLayer` moves to templates, is settled in part 2, when the
built-in annotations are reworked.

### Tests

Offline, using `datasets/dev/mmcif/` fixtures written into a step folder
in a temp dir:
- author-naming and Pandora-naming rows both resolve, and the stored
  table carries both;
- each unmatched reason is produced (including `no_coordinates` and
  `duplicate_key`);
- repeated pairs are kept in `residue_contacts`;
- a missing required column / unknown template stops the command;
- a pair template resolves both sides;
- `free_form` with only `entry_id`;
- `layer.json` counts and `previous_summary`.

## Decisions

The user decided:
- Templates, with users converting their tool output into one. Pandora
  does not parse each tool's own format.
- A user with a custom workflow provides the linking table. Without it,
  their data can't be imported.
- Data that fits no template is accepted as free-form.
- One folder per layer, like the canonicalise/curate step folders.
- Every downstream use (curation, clustering, splitting, export) must
  work on any layer of a fitting template.
- The order of parts 1–6 above.
- `pandora annotate` becomes a group: `annotate import` (new) and
  `annotate compute` (the current command, renamed, old form dropped).

My calls (Claude), which the user can overturn:
- Running the tool on a Pandora step folder is the supported path.
  Translating from original PDB files through `mappings.jsonl` is out
  of part 1.
- Either naming is accepted, and both are stored.
- Unmatched rows are reported, not fatal, with no strict mode.
- Free-form values are only carried and exported, not interpreted.
- The starter template list above.
- A residue exists only if it has atoms. Residues without coordinates
  are unmatched as `no_coordinates`.
- Repeated keys: kept for `residue_contacts` and `free_form`,
  `duplicate_key` for every other template.
- Entry IDs match case-insensitively; stored rows carry Pandora's.
- A chain given by author ID resolves among polymer chains only (an
  author chain is usually shared with its ligands and waters).
- Pair sides are stored as given, and read as unordered downstream.
