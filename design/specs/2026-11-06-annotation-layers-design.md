# Annotation layers from external tools — design

Roadmap sub-project 4 ("external tool outputs") from
`2026-10-07-contacts-design.md`, now scoped as a project of its own.
This document covers the whole project and its order. **Parts 1 and 2
are specified in detail.** Parts 3–6 are each brainstormed and
specified separately when their turn comes.

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

1. **Foundation**: templates, the layer folder, and `pandora` import of
   a user's table checked against a step folder. (Done.)
2. **Contacts are one layer.** Pandora's chain interfaces and ligand
   contacts produce `residue_contacts` layers, and interface records
   and interface clustering read any `residue_contacts` layer. (This
   absorbs the former "clustering from any contact layer" part.)
3. **Similarity is a layer.** A `similarity` template of item pairs,
   produced by the MMseqs2/Foldseek search, by `pandora similarity
   import` (replacing `--precomputed-hits`), and by the built-in
   pairwise sequence identity; `pandora cluster` reads any similarity
   layer.
4. **Export joins layers** onto chain/residue/interface records, so the
   values reach the final dataset. Entry-level structure counts become
   columns of an entry record here.
5. **Curation rules on layers**, e.g. drop chains with fewer than N
   contacts, or keep residues with a given label.
6. **Split balancing on a label** from any label layer (issue #95).

The guiding idea, confirmed in part 2's brainstorm: there is **one
layer per kind of information, not one per source**. Built-in
functions produce a layer directly; importing is the bridge that turns
another tool's output into the same layer. A layer is a typed table:
on disk a layer folder, in Python a list of rows of that template's
row model. Nothing downstream can tell a built-in layer from an
imported one apart from its recorded provenance. Extending Pandora for
a new kind of information means adding a template, not subclassing.

The existing similarity import (`pandora similarity --precomputed-hits`)
is unchanged in parts 1 and 2.

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

Recorded (part 1): **one layer per value.** A table has one `value`
(or `label`) column. Other columns are kept but not checked.

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

Recorded (part 1): **author naming only.** Non-polymer residues have
no `seq_id`, so a ligand or water row must use `auth_chain_id` +
`auth_seq_id` (+ `ins_code`). `no_coordinates` is reported for
Pandora naming from `_entity_poly_seq`. For author naming it is
reported only while `_pdbx_poly_seq_scheme` is present, which
canonicalisation may drop. Without it, the reason is
`residue_not_found`.

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
Settled in part 2: built-ins move to templates, and `AnnotationLayer`
is retired once nothing needs it.

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

## Part 2: contacts are one layer

### Intent

Today Pandora has two worlds. `annotate_chain_interfaces` and
`annotate_ligand_contacts` return an `AnnotationLayer` whose results sit
in a free-form `data` dict, and interface records
(`extract_interface_records`) and interface clustering
(`interface_residues_from_annotation`) read that dict. A COCaDA layer
imported through part 1 is a `residue_contacts` table that nothing
reads yet. Part 2 makes them one thing: both built-ins produce
`residue_contacts` rows, and both consumers take `residue_contacts`
rows, whoever produced them. COCaDA contacts and Pandora's own contacts
then give the same interface records and the same clustering input.

### Scope

In:
- `annotate_chain_interfaces` and `annotate_ligand_contacts` return
  `residue_contacts` rows.
- `extract_interface_records` and `interface_residues_from_annotation`
  take `residue_contacts` rows.
- `pandora annotate compute` writes a layer folder, one layer per run.
- `pandora cluster` takes a contact layer folder for interface coverage.
- `annotate_structure_counts` is removed.

Out:
- Pairwise sequence identity: unchanged on the old `AnnotationLayer`
  until part 3 (similarity).
- `reproduce` / `manifest`: left as they are, still on `AnnotationLayer`
  (they are likely to be removed; that is a separate decision). They
  lose `structure_counts`, and their replay of `chain_interfaces` /
  `ligand_contacts` follows whatever the new functions return (see
  "Builder's calls").
- The `AnnotationLayer` model itself stays until pairwise identity and
  reproduce/manifest no longer need it.
- Export of layers (part 4), curation on layers (part 5).

### The row model

Part 1 stores layer rows as `LayerRow` (any columns, checked by
`check_layer_columns`). Part 2 adds a typed row model for the one
template it touches, `ResidueContactRow`, in
`pandora/schemas/annotation.py` (my call; other templates get a row
model when a part needs one):

- `entry_id`
- per side `s` in `1`, `2`: `chain_id_s` (`label_asym_id`), `seq_id_s`
  (`label_seq_id`, `None` for non-polymer residues), `auth_chain_id_s`,
  `auth_seq_id_s`, `ins_code_s` (`None` for none, as part 1 stores
  it), `comp_id_s`
- `contact_type: str | None`, `distance: float | None`
- `extra="allow"`, so an imported table's own columns survive.

Column names are part 1's. `comp_id_s` is new and optional (my call):
both built-ins know it and it makes the table readable. Import fills
it from the structure when it resolves a row; part 1's
`resolve_layer_rows` does not do this yet, so its residue index gains
the comp ID. A `layer.<fmt>` written by part 1 reads into this model
unchanged except for `distance`: part 1 never type-checked it, so
import now marks a row whose `distance` is present but not a number as
`bad_value` (my call), and reading a stored layer coerces numeric
strings and fails loudly on anything else.

### Built-in producers

**`annotate_chain_interfaces(structure, distance_cutoff, atom_set,
polymer_types)`** keeps its arguments and contact detection, and
returns `list[ResidueContactRow]`: one row per contacting residue pair
across two polymer chains, side 1 the chain that sorts first (as
today). Both namings are filled from the atoms. `distance` is the
closest atom–atom distance for the pair (my call: cheap, already
computed while scanning, and makes the built-in layer at least as rich
as a typical COCaDA conversion). `contact_type` is `None`.

**`annotate_ligand_contacts(structure, distance_cutoff,
include_waters)`** returns `list[ResidueContactRow]`: side 1 the
polymer residue, side 2 the ligand residue (`seq_id_2` is `None`;
`auth_*` identifies it, as part 1 requires), `distance` the nearest
distance it reports today. `contact_type` is `None`.

A ligand is now a residue of a **non-polymer chain** (waters only with
`include_waters`), not any non-water `HETATM` (user's decision). Today
modified residues inside polymer chains (TYS in 1a2c, SME in 1ba6, HYP
in 1bkv, ...) count as ligands; as `residue_contacts` rows they would
look like polymer–polymer contacts and create false interfaces. They
are dropped from the ligand side, which changes today's ligand output
for those entries. This also makes `seq_id_2 = None` hold for every
ligand row.

Both stay plain functions with no file I/O. They do not go through
`resolve_layer_rows`: they build rows from the structure's own atoms,
so the rows already match.

### Consumers

**`extract_interface_records(rows, distance_cutoff=None,
atom_set=None)`** takes `residue_contacts` rows for any number of
entries. It groups them by entry and **unordered** chain pair, keeps
pairs where both sides are polymer residues on different chains
(`seq_id` present on both sides), and returns one `InterfaceRecord` per
pair. Side 1 of the record is the chain that sorts first. Rows that do
not qualify (ligand contacts, same-chain contacts) are skipped, not
errors: a ligand layer gives no interfaces rather than wrong ones.
Repeated rows for one residue pair (several contact types) count once.

`InterfaceRecord.distance_cutoff` becomes `float | None` (default
`None`) and `atom_set` becomes `ContactAtomSet | None`, keeping its
`"all"` default so records written before the field existed still
load as all-atom (pinned by
`test_interface_record_without_atom_set_loads_as_all`); callers pass
`None` explicitly when the setting is unknown. A
caller using the built-in passes the settings it ran with; for an
imported layer they are unknown unless the caller supplies them (e.g.
from the layer's `parameters`).

The interface residue IDs and `residue_pairs` keep today's
`label_asym_id:label_seq_id` string form, so existing
`InterfaceRecord` readers are unaffected.

Callers that today write `extract_interface_records(structure, 4.0)`
become `extract_interface_records(annotate_chain_interfaces(structure,
4.0), distance_cutoff=4.0)`. No convenience wrapper (my call; it is one
line).

**`interface_residues_from_annotation(structures, rows)`** takes
`residue_contacts` rows instead of `entry_id -> AnnotationLayer`, with
the same filter as above (polymer–polymer, different chains), and maps
each side's `seq_id` to Foldseek positions as today. Return value
unchanged. Builder's call: whether it reuses the interface-record
grouping or filters rows directly.

### The commands

**`pandora annotate compute`** writes one layer folder per run, the
same shape `annotate import` writes:

```sh
pandora annotate compute --layer chain_interfaces \
  --structures-dir curated/ \
  --distance-cutoff 5.0 --atom-set heavy \
  --output-dir annotations/interfaces
```

- `--layer`: `chain_interfaces` or `ligand_contacts`.
- `--structures-dir` replaces `--input-dir`, matching `annotate import`.
- Settings as flags: `--distance-cutoff`, `--atom-set`,
  `--polymer-types` (repeatable) for `chain_interfaces`;
  `--distance-cutoff`, `--include-waters` for `ligand_contacts`. Flags
  that don't apply to the chosen `--layer` are an error. Defaults are
  the functions' defaults.
- `--output-dir`, `--overwrite`, `--table-format` as on `import`.
- Writes `layer.<fmt>` and `layer.json`. No `unmatched` table (nothing
  can fail to match) and no `source/`.
- `layer.json` is a `StepSummary` with `step: "annotation_compute"`;
  its `policy` is a `LayerDescription` with `name` = the layer,
  `template: residue_contacts`, `tool: "pandora"`, `tool_version` =
  Pandora's version, `parameters` = the settings used. `input` and
  `previous_summary` as in `import`. `counts`: entries, rows.
- `--pairwise` is dropped from the CLI until part 3 (my call). The
  library function stays. `pandora manifest --annotations` is left as
  it is: it still reads an `AnnotationLayer` JSON, now only an old or
  hand-written one. Not worth reworking for commands likely to be
  removed; the docs example in `docs/usage/cli.md` loses that flag.
- `structure_counts` is no longer a choice.
- The old `annotations.json` output goes away (no transition, Pandora
  is 0.x).

**`pandora cluster --contacts <layer folder> --structures-dir <step
folder>`** replaces `--interface-residues <json>` (my call, no
transition). `--input-dir` stays the folder Foldseek ran on (the
per-chain files from `export_chain_mmcif`, stems `{entry}_{chain}`).
`--structures-dir` is the step folder the layer refers to; its entry
structures give the `seq_id` → Foldseek position mapping, and the
resulting `{entry}_{chain}` keys match the per-chain file stems.
Reading the per-chain files themselves would not work:
`read_structures_dir` keys by entry ID, so chains of one entry
overwrite each other. `--structures-dir` is required with `--contacts`. A layer whose `layer.json`
template is not `residue_contacts` is an error. The hand-written
positions JSON and `load_interface_residues` go away.

No command builds interface records today; part 2 does not add one
(export, part 4).

### Removed

`annotate_structure_counts`, from `pandora.annotations`,
`pandora.metadata`'s re-exports, `reproduce`'s replay table, and the
docs and examples that call it. Per-chain counts already live on chain
records; entry-level counts return as entry-record columns in part 4.

### Docs and examples to update

`docs/usage/annotation.md`, `docs/usage/datasets.md`,
`docs/usage/similarity.md`, `docs/usage/cli.md`,
`docs/usage/provenance.md`, `docs/getting-started/quickstart.md`,
`docs/index.md`, `docs/recipes/ppi-01.md`,
`docs/contributing/architecture.md`, `README.md`,
`pandora/cli/README.md`, and the scripts in `examples/` that call the
changed functions. Every snippet that reads `.data["interfaces"]` or
`.data["ligands"]` changes to rows. Regenerate the ERD
(`InterfaceRecord` and the new row model change).

### Tests

Offline, with `datasets/dev/mmcif/` fixtures:
- `annotate_chain_interfaces` rows: both namings filled, `distance` is
  the closest pair distance, same residue pairs as today's layer for the
  fixtures the current tests use;
- `annotate_ligand_contacts` rows: ligand side by author naming,
  `seq_id_2` empty; a modified residue in a polymer chain (e.g. 1a2c's
  TYS) is not a ligand;
- `extract_interface_records` gives the same records as today for the
  built-in layer, treats (A,B) and (B,A) rows as one pair, counts a
  repeated residue pair once, skips ligand and same-chain rows, and
  works on rows read back from an imported `layer.<fmt>`;
- `interface_residues_from_annotation` on built-in and imported rows
  gives the same positions;
- `annotate compute` writes `layer.<fmt>` + `layer.json` with the
  settings, rejects a flag that doesn't apply, refuses an existing layer
  without `--overwrite`;
- a layer from `annotate compute` and the same contacts brought in with
  `annotate import` give identical interface records;
- `cluster --contacts` with a fake Foldseek, and a non-contact layer is
  rejected.

### Builder's calls

- How reproduce/manifest keep working for `chain_interfaces` /
  `ligand_contacts` with the new return type (wrap the rows in an
  `AnnotationLayer` there, or drop those two from its replay table).
  Pick the smaller change; they are likely to be removed.
- Module layout inside `pandora/annotations/` and
  `pandora/datasets/records.py`.

## Decisions

The user decided (part 2):
- One layer per kind of information, not one per source. Built-ins
  produce the layer directly; import is the bridge for other tools.
  The `pandora` command stays enough to bring in another tool's output.
- Interface records and interface clustering read any
  `residue_contacts` layer, in part 2.
- Reproduce/manifest are left alone (likely to be removed later).
- Structure counts are dropped as an annotation; entry-level counts
  belong on an entry record at export.
- Ligand contacts count only residues of non-polymer chains as
  ligands; modified residues inside polymer chains are not ligands.
- Similarity becomes a layer (`pandora similarity import`), as its own
  part, next after part 2. Pairwise identity moves there.

My calls (part 2), which the user can overturn:
- A layer is a typed table (pydantic row model), not a class with
  behaviour; extending means a new template, not subclassing.
- `ResidueContactRow` added now; other templates' row models later.
- `comp_id_s` columns; closest-atom `distance` on chain-interface rows.
- `InterfaceRecord.distance_cutoff` / `atom_set` become optional.
- `annotate compute`: one layer per run, settings as flags,
  `--structures-dir`; `annotations.json` dropped.
- `cluster --contacts` replaces `--interface-residues`.
- No convenience wrapper for structure → interface records.

The user decided (part 1):
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
