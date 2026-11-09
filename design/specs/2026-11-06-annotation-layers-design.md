# Annotation layers from external tools — design

Roadmap sub-project 4 ("external tool outputs") from
`2026-10-07-contacts-design.md`, now scoped as a project of its own.
This document covers the whole project and its order. **Parts 1, 2
and 3 are specified in detail.** Parts 4–6 are each brainstormed and
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
4. **Layers are joined by the user.** Every layer can be joined onto
   chain/residue/interface records on shared keys; the docs show how,
   and Pandora adds no join step. Contact counts return as
   `residue_values` / `chain_values` layers. (Re-scoped 2026-11-09:
   was "export joins layers", with an entry record.)
5. **Interfaces are dropped by label.** A contact layer is filtered
   by an `interface_labels` layer into a new contact layer. (Re-scoped
   2026-11-09: was "curation rules on layers", e.g. drop chains with
   fewer than N contacts.)
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
is unchanged in parts 1 and 2, and replaced in part 3.

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
records; part 4 brings contact counts back as layers (it was planned
as entry-record columns; part 4 dropped the entry record).

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

## Part 3: similarity is a layer

### Intent

Similarity today is tied to two tools. `pandora similarity` runs
MMseqs2 or Foldseek and keeps the raw hit TSV; `--precomputed-hits`
accepts only a TSV with Pandora's exact MMseqs2/Foldseek columns; and
`pandora cluster` reads only that TSV. A user who computed similarity
another way (the motivating case: an R script that built an igraph
graph) cannot cluster it with Pandora. The items being compared also
vary: entries, chains, or interfaces.

Part 3 makes similarity a layer like the others: a `similarity`
template of item pairs. The built-in searches produce it, `pandora
similarity import` brings in anyone else's edge list, and `pandora
cluster` reads any similarity layer.

Two constraints shape everything below:

- **Whole-archive scale.** This must run on the whole PDBe archive:
  about 240k entries and around 1M polymer chains (my estimates). An
  all-vs-all search keeping up to 1,000 hits per query gives up to
  about 10⁹ rows. Nothing may hold all rows in memory, build a Python
  object per row, or parse every structure of the archive.
- **Modularity (user's decision).** Each step reads what an earlier
  step wrote and never recomputes it. Similarity import and clustering
  get their list of valid items from an earlier step's table or layer,
  never by parsing structures or detecting contacts again.

### Scope

In:
- the `similarity` template and its description field `items`;
- `pandora similarity search` and `pandora similarity import`, both
  writing a layer folder;
- `pandora cluster --similarity`, at any item level;
- a tool-neutral `HitFilter`;
- removing `annotate_pairwise_sequence_identity` (user's decision),
  `--precomputed-hits`, `SimilaritySearch`, `load_similarity_search`,
  the raw-TSV `iter_edges`, and `SimilarityRelationship` (unused).

Out:
- Running MMseqs2 or Foldseek at archive scale: that is the tool's job.
  Pandora reads what it writes.
- Parallel or sharded runs. One machine, one process, like the other
  steps.
- `partition` and `--pairs` / `pair_cluster_keys`: unchanged.
- Export of layers (part 4), curation on layers (part 5).

### The template

`similarity`: one row per hit between two items.

| Column | Type | Required |
|---|---|---|
| `item_id_1`, `item_id_2` | string | yes (keys) |
| `score` | number | yes |
| `identity` | number | recognised if present |
| `coverage_1`, `coverage_2` | number | recognised if present (coverage of item 1 / item 2) |
| `start_1`, `end_1`, `start_2`, `end_2` | integer | recognised if present (1-based alignment range on each item, used for interface coverage) |

Any other column is kept as given, as in the other templates.

**Item level.** Each layer holds one level, declared in the description
YAML with a new `items` field. `LayerDescription` gains
`items: Literal["entry", "chain", "interface"] | None`; it is required
when `template: similarity` and rejected for every other template.

```yaml
name: my_r_graph
template: similarity
items: chain
tool: R script
tool_version: "4.4"
parameters: {}
```

**Item IDs** are fixed-format strings, one format per level:

| `items` | ID form | Example |
|---|---|---|
| `entry` | `{entry}` | `1ABC` |
| `chain` | `{entry}_{chain}` | `1ABC_A` |
| `interface` | `{entry}_{chain1}_{chain2}` | `1ABC_A_B` |

- `chain` is Pandora's `label_asym_id` only. An ID string has nowhere to
  say which naming it uses, and Pandora's per-chain files
  (`export_chain_mmcif`, Foldseek's input) already use these (my call).
- Split from the right: chain IDs never contain `_`, but extended PDB
  entry IDs (`pdb_00001abc`) do.
- The entry part matches case-insensitively, and stored IDs carry
  Pandora's entry ID (as in part 1). Chain IDs are case-sensitive.
- Interface IDs are stored with the two chains sorted, so `1ABC_B_A`
  becomes `1ABC_A_B` and one interface has one ID (my call).
- "Structure" in the user's sense is the same as entry: a step folder
  holds one structure per entry (my assumption; if assemblies or
  models become separate items, that's a new level).

**Rows are not deduplicated or reordered.** A→B and B→A may both be
present, and self-hits are kept. Clustering is unaffected by either,
and deduplicating would mean holding every pair in memory (my call).
This differs from part 1's `duplicate_key` rule for single-key
templates, on purpose.

`SimilarityRow` (in `pandora/schemas/annotation.py`, `extra="allow"`)
describes the row and is used to check small tables and in tests. Bulk
reading never builds one per row.

### Where the item list comes from

| `items` | Item list read from | Written by |
|---|---|---|
| `entry` | `--structures-dir <step folder>`: its `entries.<fmt>`, rows that have a structure | `canonicalise` / `curate` |
| `chain` | `--structures-dir <step folder>`: its `chains.<fmt>` | `curate` |
| `interface` | `--contacts <residue_contacts layer folder>`: its distinct unordered pairs of polymer chains | `annotate compute` / `annotate import` |

- Only the tables are read, never the mmCIF files. A step folder
  without the needed table (e.g. chain items on a `canonicalise` folder,
  which writes no `chains` table) is an error that says which step
  writes it.
- Interface pairs use the same filter as `extract_interface_records`:
  both sides polymer residues (`seq_id` present) on different chains.
  The layer is read in one streaming pass over its entry and chain
  columns; the resulting set (at most a few million pairs) fits in
  memory. An interface is valid only if a contacts layer identified it
  (user's decision): there is no interface table of its own, since the
  contacts layer already lists the pairs (my call).
- The flag not needed for the level is an error, as in `annotate
  compute`.
- `layer.json` records which item list was used: its path and that
  folder's summary, alongside `previous_summary`.

### The commands

`pandora similarity` becomes a group, like `annotate`. The old flat
form stops working, with no transition (Pandora is 0.x).

**`pandora similarity search`** runs the tool and writes a layer
folder:

```sh
pandora similarity search --engine mmseqs2 \
  --structures-dir curated/ --output-dir similarity/seq
pandora similarity search --engine foldseek \
  --input-dir chains/ --structures-dir curated/ \
  --output-dir similarity/struct
```

- MMseqs2 searches **chains** (`items: chain`), with sequences read
  from the step folder's `chains.<fmt>`. Today's entry-level search
  (each entry's longest chain, `entry_sequences`) goes away;
  entry-level similarity can be imported (my call).
- Foldseek searches the per-chain files in `--input-dir`, as today
  (`items: chain`). `--structures-dir` gives the chain list the stems
  are checked against.
- Engine options as today: `--sensitivity`, `--max-seqs`,
  `--alignment-type`, `--exhaustive-search`, `--mmseqs-bin`,
  `--foldseek-bin`.
- The tool's TSV is converted to `layer.parquet` in one streaming pass
  and then deleted; the layer keeps every column. Column mapping:
  `query`→`item_id_1`, `target`→`item_id_2`, `fident`→`identity`,
  `qcov`→`coverage_1`, `tcov`→`coverage_2`, `qstart`/`qend`/`tstart`/
  `tend`→`start_1`/`end_1`/`start_2`/`end_2`; `score` = `fident` for
  MMseqs2 (today's behaviour) and `alntmscore` for Foldseek; `alnlen`,
  `qtmscore`, `ttmscore` kept as extra columns.
- Hits naming an item not in the item list go to `unmatched.parquet`,
  as in `import`.
- `layer.json`: `step: "similarity_search"`, `policy` a
  `LayerDescription` with `template: similarity`, `items: chain`,
  `tool` the engine, `tool_version` the tool's version, `parameters`
  the settings. `counts`: rows in, matched, unmatched (per reason).
- `--output-dir`, `--overwrite` as on `annotate`.

**`pandora similarity import`** brings in any edge list:

```sh
pandora similarity import \
  --table r_graph_edges.parquet --description r_graph.yaml \
  --structures-dir curated/ --output-dir similarity/r_graph
```

- `--table`: Parquet, or TSV with a header row. A user's own MMseqs2 or
  Foldseek TSV goes through here with its columns renamed into the
  template; this replaces `--precomputed-hits`.
- `--description`: the YAML above, with `items`.
- `--structures-dir` or `--contacts` for the item list (see above).
- `--source-files`, `--output-dir`, `--overwrite` as on `annotate
  import`.
- Writes `layer.parquet`, `unmatched.parquet` (with `reason`),
  `source/`, and `layer.json` (`step: "similarity_import"`).
- Unmatched reasons: `entry_not_found`, `chain_not_found`,
  `interface_not_found`, `bad_item_id` (doesn't fit the level's form),
  `bad_value` (a recognised column present but not a number). A missing
  required column, or `items` missing from the description, stops the
  command before anything is written, as in part 1.

**Parquet only** (my call). Similarity layers are always Parquet and
need pyarrow (the `export` extra). JSON Lines at 10⁹ rows isn't usable,
and one code path is simpler. There is no `--table-format` on these
commands; without the extra they stop with a message naming it.

**`pandora cluster`**:

```sh
pandora cluster --similarity similarity/struct \
  --structures-dir curated/ --min-score 0.5 --output clusters.json
```

- `--similarity <layer folder>` replaces `--search`. A layer whose
  `layer.json` template is not `similarity` is an error.
- The full item list comes from `--structures-dir` / `--contacts` as
  above, so items with no hits are single-item clusters. `--input-dir`
  and its `*.cif` listing go away.
- Clusters are at the layer's item level, except for one roll-up
  (user's decision): `--items entry` on a **chain** layer gives entry
  clusters. Each chain is linked to its entry in the union-find, so two
  entries share a cluster when any of their chains are similar. This is
  stricter against leakage than today's longest-chain-per-entry search,
  and needs only the one chain search. The entry list comes from
  `--structures-dir` (its `entries` and `chains` tables). Any other
  `--items` that differs from the layer's level is an error. The PPI
  recipe (`docs/recipes/ppi-01.md`, `examples/ppi_dataset_pipeline.py`)
  moves to a chain search plus this roll-up; its documented cluster
  count is re-measured, not kept.
- `cluster` may use a different item list from the one the layer was
  matched against (my call: search the archive once, cluster many
  curated subsets). Edges naming an item outside the list are dropped
  and counted (`n_edges_unknown_ids`, as today). Links that ran through
  a dropped item disappear.
- **Symmetry copies** (`chains.<fmt>` rows with `copy_of` set) are
  chain items, but `similarity search` leaves them out of the search:
  identical copies would fill the per-query `--max-seqs` slots and push
  out real homologues. `cluster`, with chain items, joins each copy to
  its source chain's cluster using `copy_of` from the same table (my
  call). An imported layer may name copies; they resolve like any chain.
- For interface items, `--contacts` gives the item list. For chain
  items, `--contacts` (with `--structures-dir`) still gives interface
  residues for `min_interface_coverage`, as in part 2.
- `--items entry|chain|interface`: optional, defaults to the layer's
  level; only `entry` on a chain layer differs (see the roll-up above).
- `--hit-filter`, `--min-score`, `--min-coverage`, `--pairs`,
  `--output` unchanged.

### Scaling

- **Reading:** layers are read with pyarrow in record batches, only the
  columns needed. The filter is applied to whole columns at once
  (pyarrow compute), so only rows that pass reach Python.
- **Clustering:** item IDs are mapped to integer indices and union-find
  runs on integer arrays, so memory grows with the number of items
  (about 1M), not edges. Expected minutes, not hours, for the archive;
  this is a guess to measure on a real run, not a requirement the tests
  check. `cluster_edges` keeps its contract (any iterable of id pairs).
- **Import and search conversion:** one streaming pass, batch by batch
  (pyarrow's streaming CSV reader for TSV); the item-ID check is a
  column-wise lookup per batch against the item set. Matched and
  unmatched rows are written batch by batch.
- **Known ceiling:** `min_interface_coverage` still needs the step
  folder's structures parsed to map `seq_id` to Foldseek positions
  (part 2's `interface_residues_from_annotation`). That doesn't reach
  archive scale; it is left as is until a real run needs it.

### HitFilter becomes tool-neutral

`HitFilter` applies to the template's columns, not to a tool's:

- `min_score`, `min_identity`, `min_coverage`: unchanged meaning, on
  `score`, `identity`, and coverage.
- `coverage_of`: `both` / `side_1` / `side_2` / `either` (was `query` /
  `target`).
- `tm_normalisation` is replaced by `score_column: str = "score"`, the
  column `min_score` tests. The built-in Foldseek layer keeps
  `qtmscore`/`ttmscore`, so `score_column: qtmscore` gives what
  `tm_normalisation: query` gave. The old `max` option is gone; add a
  column before import if needed (my call).
- A threshold on a column the layer doesn't have is an error before any
  row is read, as `iter_edges` does today for engine mismatches.
- `min_interface_coverage`: unchanged, needs `start_*`/`end_*` and
  interface residues.
- A NaN value fails every threshold, as today.

Update `docs/policies/*.yaml` and `docs/reference/policies.md` for the
changed fields.

### Library side

The commands are thin wrappers, as elsewhere. File reading and writing
stay out of the plain functions where possible, but streaming means the
similarity functions work on pyarrow record batches, not lists of rows.
Builder's call on exact names and signatures, within:

- a function that takes record batches of a user's table, the item
  level and the item set, and yields (matched, unmatched) batches;
- a function that takes record batches of a layer and a `HitFilter`
  (plus interface residues) and yields edges for `cluster_edges`;
- the MMseqs2/Foldseek search functions write the tool's TSV and
  report its path and settings; converting that TSV into a layer goes
  through the same column mapping and matching as import;
- functions that build the item set from an entries table, a chains
  table, or contact rows;
- pyarrow imported inside the functions that need it.

### Removed

- `annotate_pairwise_sequence_identity` (`pandora/annotations/
  pairwise.py`), its exports and docs (user's decision: an ungapped
  all-against-all in Python can't reach archive scale, and MMseqs2 is
  the built-in sequence route).
- `SimilaritySearch`, `load_similarity_search`, `HIT_COLUMNS`-based
  `iter_edges`, `cluster_similar_items` in its current form, and the
  `--precomputed-hits` flag.
- `SimilarityRelationship` and its provenance model (nothing uses
  them).
- `entry_sequences`, if nothing else needs it after the MMseqs2 search
  moves to chains.
- `reproduce` / `manifest`: they replay similarity through
  `SimilaritySearch` and `entry_sequences`. Builder's call, smallest
  change: drop the similarity step from their replay. They are likely
  to be removed (part 2 decision). `ClusteringProvenance.search`
  becomes the similarity layer's `layer.json` path or summary.

### Docs and examples to update

`docs/usage/similarity.md`, `docs/usage/cli.md`,
`docs/usage/annotation.md`, `docs/usage/export.md`,
`docs/usage/provenance.md`, `docs/getting-started/glossary.md`,
`docs/contributing/architecture.md`, `docs/contributing/contribute.md`,
`docs/reference/policies.md`, `docs/policies/*.yaml`,
`pandora/cli/README.md`, `README.md`, and the scripts in `examples/`
that call the changed functions. Regenerate the ERD
(`docs/scripts/generate_erd.py`, which names the removed models).

### Tests

Offline, with `datasets/dev/mmcif/` fixtures and fake MMseqs2/Foldseek
scripts as today:
- import at each level: entry, chain, interface IDs resolve; each
  unmatched reason is produced; interface IDs are sorted; the entry
  part matches case-insensitively; an extended entry ID with `_`
  parses;
- interface items come from a contacts layer only: a chain pair present
  in the step folder but not in the layer is `interface_not_found`;
- `items` required for `similarity`, rejected for other templates; a
  missing required column stops before writing;
- `similarity search` with fake MMseqs2 (chain sequences from
  `chains.<fmt>`) and Foldseek writes a layer with the column mapping
  above, and the raw TSV is gone;
- the same hits through `search` and through `import` (columns renamed)
  give identical clusters;
- `--items entry` on a chain layer: two entries linked only through a
  non-longest chain end up in one cluster;
- symmetry copies are left out of the fake search's input and land in
  their source's cluster; edges to items outside a smaller item list
  are dropped and counted;
- `cluster --similarity` at each level, singletons included; a
  non-similarity layer is rejected; the wrong item-list flag for the
  level is rejected;
- `HitFilter`: `score_column`, `coverage_of` sides, a threshold on an
  absent column is an error, NaN fails, `min_interface_coverage` with
  `--contacts`;
- batching: a layer written in several small record batches clusters
  the same as one batch (proves nothing depends on whole-table reads).

## Part 4: layers are joined by the user

### Intent

The roadmap line was "export joins layers onto records". Thinking as
the people who use the output (ML engineers training structure models,
structural biologists) changed that. They load one structure or one
table at a time and pick the features they want. A wide table per
level, whose columns change with every imported layer and which copies
chain values onto every residue, helps neither. Layers are already
Parquet tables keyed by the same IDs as the records, so joining is one
line in pandas, polars or DuckDB, and the user chooses which layers
to bring in.

So part 4 makes sure every layer **can** be joined onto the records,
and leaves the joining to the user. It also brings back the contact
counts that part 2 dropped, as layers.

### Scope

In:
- residue records gain the author keys they lack;
- a built-in step that turns a `residue_contacts` layer into a count
  layer;
- docs that show the joins.

Out:
- no join command and no join function (user decision);
- no entry record and no entry-level template (see "Entry level");
- similarity: unchanged. Cluster IDs come from `cluster`, and hit
  pairs stay their own layer.

### Joinable keys

Every template's key columns already match the records'
(`entry_id`, `chain_id`, `seq_id`; `chain_id_1`/`chain_id_2` for
interfaces), because import stores Pandora's naming. One gap:
`ResidueRecord` has `auth_seq_id` but no author chain or insertion
code, so a user who prefers author numbering can't join on it
unambiguously. Add to `ResidueRecord`:

- `auth_chain_id: str | None` (`auth_asym_id`)
- `ins_code: str | None` (`pdbx_PDB_ins_code`, `None` for none, as
  layers store it)

Both come from the residue's atoms. `extract_residue_records` groups by
`label_asym_id` + `label_seq_id` as now, and the new fields are read
from the first atom (my call). Records written before this change load
with both fields `None`.

### Contact counts as layers

A library function per level in `pandora/annotations/`, taking the
rows of one `residue_contacts` layer and returning rows of a value
layer. No file I/O, inputs not mutated.

- `count_residue_contacts(rows) -> list[LayerRow]`: a
  `residue_values` row per residue that appears on either side, with
  `value` = the number of **distinct partner residues**. A pair given
  several times (e.g. several contact types) counts once. The pair is
  unordered, as everywhere downstream.
- `count_chain_contacts(rows) -> list[LayerRow]`: a `chain_values`
  row per chain, with `value` = the number of **distinct residues of
  that chain** in contact.

Rows carry both namings, copied from the contact rows, so a ligand
residue (no `seq_id`) gets a row keyed by author naming, as part 1
allows. Residues and chains with no contacts get **no row** (my call):
the layer stays the size of the contact layer, and the docs show
`fillna(0)` after the join.

Command: `pandora annotate count-contacts --contacts <layer dir>
--level residue|chain --output-dir <dir> [--overwrite]
[--table-format]`. One layer per run, as `annotate compute`. It
writes an ordinary layer folder (`layer.<fmt>`, `layer.json`):
`StepSummary` with `step: "annotation_count"`, a `LayerDescription`
(`name`: `contact_count`, `template`: `residue_values` or
`chain_values`, `tool`: `pandora`, `parameters`: `level` and the
source layer's name), `input` = the contacts layer folder,
`previous_summary` = that folder's `layer.json`, `counts` = rows in
and rows out. It reads only the layer, never structures. A layer whose
template is not `residue_contacts` is an error before anything is
written.

### Entry level

No entry record and no `entry_values` template for now (user decision,
on my recommendation):
- per-entry tool output can already be imported as a `free_form` layer
  (only `entry_id` is required), which joins on `entry_id`;
- entry totals are a one-line sum of the chain count layer, shown in
  the docs.

If the same per-entry shape keeps arriving as `free_form`, add an
`entry_values` template then, by the rule part 1 set for free-form
layers. Entry metadata (method, resolution) stays on
`DatasetChainRecord`.

### Docs

A "Using layers in a dataset" section in `docs/usage/annotation.md`,
with full code inline and tested output, as the other usage examples
(library and CLI tabs where both exist):
- residue values onto residue records, on Pandora keys and on author
  keys;
- chain counts onto `chains.<fmt>`, with `fillna(0)`;
- interface labels onto interface records, matching the chain pair in
  either order (labels keep the user's order, records don't promise
  one);
- an entry total by summing the chain count layer.

`docs/usage/cli.md`, `pandora/cli/README.md`, `docs/usage/datasets.md`
(new residue record fields), `docs/contributing/architecture.md`.
Regenerate the ERD (`ResidueRecord` changes).

### Tests

Offline, with `datasets/dev/mmcif/` fixtures:
- residue records carry `auth_chain_id` and `ins_code`; a fixture with
  an insertion code shows it; old records without the fields still
  load;
- `count_residue_contacts`: a pair repeated with two contact types
  counts once; (A,B) and (B,A) count once; both sides get a row; a
  ligand side is keyed by author naming;
- `count_chain_contacts`: distinct residues per chain;
- the count layers pass `check_layer_columns` for their template;
- `annotate count-contacts` on a layer from `annotate compute` writes
  a layer folder whose summary links back; a non-contact layer is
  rejected before writing;
- inputs are not mutated.

## Part 5: interfaces are dropped by label

### Intent

The real case: a user imports COCaDA contacts, judges each interface
on its physicochemical properties (e.g. too few hydrogen bonds, too
hydrophobic) and wants the failing interfaces gone from everything
downstream: interface records, `cluster --contacts`, the split.

Following part 4, the user makes the judgement in pandas; Pandora
neither computes "H-bonds per interface" nor offers a rule language.
Pandora's job is to apply the verdict and record it. Every downstream
consumer already reads any `residue_contacts` layer, so the verdict is
applied there: the output is a **new contact layer** without the
dropped interfaces' rows. Interface records built from it, and
clustering, follow for free. The original layer is untouched.

The verdict itself is a layer, not a loose table (user decision): an
`interface_labels` layer (part 1 template) the user writes in pandas and
brings in with `annotate import`. So both inputs carry provenance.

**One label per interface** (user decision): import already rejects a
second label row for the same interface (`interface_labels` is unique
per key), and part 5 does not change that. A labels layer can still
use several values across interfaces (`ok`, `low_hbond`,
`too_hydrophobic`) to pick from per run. Several criteria for the same
interface are either folded into one label in pandas, or written as
one labels layer per criterion and applied by running
`filter-contacts` once per layer, each output feeding the next; each
run records its parent, so the chain stays readable.

### Scope

In:
- a library function that drops labelled interfaces from contact rows;
- `pandora annotate filter-contacts`, writing a new contact layer;
- docs and tests.

Out (my call; add each the same way when a real need appears):
- dropping chains or residues by label (`chain_labels`,
  `residue_labels`), or whole entries;
- thresholds on value layers (`chain_values` above N);
- filtering interface records directly, or modifying structures.

### Matching

- An interface is `(entry_id, unordered chain pair)` on Pandora chain
  IDs (`chain_id_1`, `chain_id_2`), which import stores for both
  layers. A label on 1abc B–A drops contact rows written A–B.
- Matching is exactly that, with no special cases: ligand contact rows
  match a label on (polymer chain, ligand chain); intra-chain rows
  match only a label on (A, A). (My call: simpler than excluding them,
  and a label on (A, A) is the user asking for it.)
- Import keeps one label row per interface, but its duplicate check
  is on the chain pair in the order given, so A–B and B–A can both
  slip through. If they do, the interface is dropped if **any** of its
  label rows carries a `--drop` label (my call; no error).
- **Unlabelled interfaces are kept** (user decision): the filter
  drops only what was named. The count of unlabelled interfaces is
  reported so a typo or a labels layer from another run is visible.
- A `--drop` value found on no row of the labels layer is an error
  before anything is written (my call: almost always a typo).

### Library side

In `pandora/annotations/`, exported:

`drop_labelled_interfaces(contacts, labels, drop) -> list[ResidueContactRow]`
- `contacts`: rows of a `residue_contacts` layer;
- `labels`: rows of an `interface_labels` layer (`LayerRow`, no typed
  row model needed);
- `drop`: the label values that mean drop;
- returns the kept contact rows in input order. No file I/O, inputs not
  mutated. Raises `ValueError` for a `drop` value absent from `labels`.

The counts for the summary are the builder's choice of shape (a second
small function or computed in the CLI).

### The command

`pandora annotate filter-contacts --contacts <layer dir> --labels
<layer dir> --drop <label> [--drop <label> ...] --output-dir <dir>
[--overwrite] [--table-format]`

- `--contacts` must be a `residue_contacts` layer
  (`check_contact_layer`), `--labels` an `interface_labels` layer;
  otherwise an error before anything is written. `--output-dir` may
  not be either input folder (as `count-contacts`).
- Reads only the two layers, never structures.
- **Both layers must come from the same structures folder** (user
  decision): chain IDs can change between step folders (e.g. raw vs
  canonicalised), and a label on A–B would then silently drop a
  different interface. Each layer's origin is found by following
  `previous_summary` back through `layer.json` files until it reaches
  something that is not a layer summary: for an imported layer, the
  structures folder's step summary (`curation.json` /
  `canonicalisation.json`); for a layer from `annotate compute`, its
  `input` (the structures folder). Compare resolved folders. Different
  origins are an error before anything is written. If an origin cannot
  be found (an import against a folder with no step summary records
  none today), print a warning and carry on (my call: refusing would
  block raw-folder layers outright). The library function does no such
  check; it sees only rows.
- Writes an ordinary layer folder (`layer.<fmt>`, `layer.json`), as
  `count-contacts` does: `StepSummary` with `step:
  "annotation_filter"`, a `LayerDescription` (`name`:
  `filtered_contacts`, `template`: `residue_contacts`, `tool`:
  `pandora`, `parameters`: the contacts layer's name, the labels
  layer's name and folder, the `drop` labels), `input` = the contacts
  folder, `previous_summary` = the contacts folder's `layer.json`
  (which chains back to the original import).
- `counts`: `rows_in`, `rows_out`, `interfaces_in`,
  `interfaces_dropped`, `interfaces_unlabelled`, and
  `labels_unmatched` (labelled interfaces not in the contact layer).
- Prints one line: interfaces dropped, kept, unlabelled.

### Docs

In `docs/usage/annotation.md`, after "Using layers in a dataset": the
COCaDA-style flow end to end, full code inline with tested output,
library and CLI tabs: count a contact type per interface in pandas,
write the labels table, `annotate import` it as `interface_labels`,
`annotate filter-contacts`, then build interface records / `cluster
--contacts` from the new layer. Contact types in the fixture rows are
written by hand (the built-in layers leave `contact_type` empty). Also
`docs/usage/cli.md`, `pandora/cli/README.md`,
`docs/contributing/architecture.md`.

### Tests

Offline, with `datasets/dev/mmcif/` fixtures:
- a labelled interface's rows go, in both chain orders; other rows stay
  in order;
- unlabelled interfaces are kept; label rows written as both A–B and
  B–A, one dropping and one keeping, drop the interface;
- a ligand interface and an (A, A) label match as described;
- an unknown `--drop` label raises before writing;
- CLI: layers traced to different structures folders are rejected
  before writing; an untraceable origin warns and runs;
- inputs are not mutated;
- CLI: on a layer from `annotate compute` plus an imported labels
  layer, writes a layer folder whose summary links back and whose
  counts are right; wrong templates and `--output-dir` equal to an
  input are rejected before writing; `cluster --contacts` accepts the
  output.

## Decisions

The user decided (part 5):
- The real need is dropping interfaces on their physicochemical
  properties (e.g. from COCaDA contact types), not chains or entries.
- The drop applies to the contact layer, producing a new layer, so
  every downstream consumer follows (my recommendation, accepted over
  filtering interface records).
- A filter step that records its parent layers, rather than "filter in
  pandas and re-import".
- The verdict comes in as an `interface_labels` layer, not a plain
  table.
- Interfaces without a label are kept, and their count reported (my
  recommendation, accepted).
- One label per interface; import is unchanged. Several criteria are
  folded in pandas or applied by chained runs (my recommendation,
  accepted).
- The two layers must trace back to the same structures folder, or
  the command refuses (my recommendation, accepted).

My calls (part 5), which the user can overturn:
- Command name `annotate filter-contacts`, `--drop` repeatable.
- No special case for intra-chain or ligand rows.
- A–B and B–A label rows that both slip past import: any `--drop`
  label drops the interface.
- An untraceable layer origin warns rather than refuses.
- An unknown `--drop` label is an error.
- Output described as `tool: pandora`, `name: filtered_contacts`, as
  `count-contacts` does.
- Chains, residues, entries and value thresholds are out of scope.

The user decided (part 4):
- Export is general purpose, not shaped for one benchmark.
- Layers stay separate files; users join them. No wide table per
  level, no join command (and, following from that, no join function:
  the docs show the pandas join).
- Both per-residue values and contact counts are worth having.
- No entry template or entry record for now (my recommendation,
  accepted).

My calls (part 4), which the user can overturn:
- Counts are layers made by `annotate count-contacts`, one level per
  run, not record columns.
- Residue count = distinct partner residues; chain count = distinct
  residues in contact; no rows for zero.
- `ResidueRecord` gains `auth_chain_id` and `ins_code`.

The user decided (part 3):
- Similarity input can come from any source (e.g. an R script's edge
  list), so the template is generic: two items and a score.
- Items can be entries, chains, or interfaces.
- It must run on the whole PDBe archive.
- Modularity: each step reads earlier outputs and never recomputes.
  An interface item is valid only if a contacts layer identified it.
- `annotate_pairwise_sequence_identity` is removed.
- Entry clusters come from rolling a chain layer up to entries in
  `cluster`; the longest-chain entry search goes.

My calls (part 3), which the user can overturn:
- One item level per layer, set by `items` in the description.
- Fixed-format string IDs; Pandora chain IDs only; interface IDs
  sorted at import.
- No deduplication or reordering of rows; thresholds applied at
  clustering, not import.
- Item lists from step-folder tables and contact layers, never from
  parsing structures; no interface table of its own.
- `similarity` becomes a group (`search`, `import`); the raw TSV is
  converted and deleted; `--precomputed-hits` and `SimilaritySearch`
  go.
- MMseqs2 searches chains, without symmetry copies; `cluster` puts
  copies in their source's cluster.
- `cluster` accepts a smaller item list than the layer was matched
  against; edges outside it are dropped and counted.
- Similarity layers are Parquet only (the `export` extra).
- `HitFilter`: `score_column` replaces `tm_normalisation` (no `max`);
  `coverage_of` uses sides.
- `min_interface_coverage` keeps parsing structures (known ceiling).

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
