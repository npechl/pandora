# CLI `ingest` finds input files — design

Follows `2026-11-04-cli-step-outputs-design.md`. `pandora canonicalise`
reads a folder (flat or PDB mirror layout) or a list of paths
(`--entries`). Real archives use other layouts, e.g. the PDBe release
tree on NFS:

```
/nfs/production/gerard/pdbe/release/data/entry/00/100d/clean_mmcif/100d_updated.cif.gz
```

## Intent

Finding the input files becomes its own step: `pandora ingest`. It
writes `entries.txt`, which `pandora canonicalise --entries` already
reads. Every source (a local mirror, the PDBe release tree, `fetch`)
then reaches `canonicalise` in the same way, and the steps after
`ingest` don't need to know any folder layout.

Speed matters. The PDBe tree holds about 230k entry folders on NFS,
where every directory read or file stat costs a network round trip.
Opening each entry folder (`find`, `ls -R`, a `*/*/clean_mmcif/*`
glob) is the slow part. Listing only the top two levels gives every
entry ID in about 1,300 directory reads, and each file path follows
from its ID.

## Decisions (user)

- `ingest` records provenance once for the whole run (the input root,
  how files were found, an optional label and release date). It does
  **not** open, parse or stat individual files, and writes no
  per-entry record. A file that is missing or broken shows up later as
  a `canonicalise` failure.
- `ingest` no longer writes `ingestion_provenance.json`. `manifest` and
  `reproduce` may be dropped later, so nothing here is designed to feed
  them. `manifest --ingestion-provenance` is optional and still accepts
  `fetch`'s file. Neither command changes.
- The `--glob` option added to `canonicalise` and `curate` earlier in
  this branch (uncommitted) is **removed**. Finding files is `ingest`'s
  job.
- Kept from that same uncommitted change: `canonicalise` and `curate`
  name output files after the parsed entry ID (lowercase), not the
  input file name, so `100d_updated.cif.gz` becomes
  `structures/00/100d.cif`. Duplicates are detected by entry ID and
  recorded as `DuplicateEntry`. This replaces the 2026-11-04 rule
  "file stems keep the input file name's case".

## The command

Exactly one of `--glob` or `--path-template` is required:

```sh
# any layout; reads every folder the pattern passes through
pandora ingest --input-dir mirror/ --glob '*/*.cif' --output-dir ingest/

# ID-derived paths; lists only the folders above the ID
pandora ingest \
  --input-dir /nfs/production/gerard/pdbe/release/data/entry \
  --path-template '{mid}/{id}/clean_mmcif/{id}_updated.cif.gz' \
  --source-label pdbe_release_2026-11-05 --release-date 2026-11-05 \
  --output-dir ingest/

# known IDs; no filesystem access under --input-dir at all
pandora ingest --input-dir ... --path-template '...' --ids ids.txt \
  --output-dir ingest/

pandora canonicalise --entries ingest/entries.txt --policy ... --output-dir canon/
```

Options:

| Option | Meaning |
|---|---|
| `--input-dir` | Root of the archive. Required. |
| `--glob PATTERN` | `Path.glob` pattern relative to `--input-dir`. Matches are listed as found, sorted. |
| `--path-template T` | Path relative to `--input-dir` with placeholders `{id}` (lowercase entry ID) and `{mid}` (`id[1:3]`, the PDB mirror folder). |
| `--ids FILE` | Only with `--path-template`. One entry ID per line; blank lines and `#` lines skipped; lowercased. The paths are built from the IDs without listing anything. |
| `--source-label`, `--release-date` | Free text stored in the summary as given. Optional. |
| `--output-dir` | Required. |
| `--overwrite` | As in the other step commands. |

### `--path-template` without `--ids`

The first path segment that contains `{id}` is the **ID segment**.

1. The segments before it are listed level by level, with `{mid}`
   matching any folder name and literal segments matching only
   themselves.
2. The names in each listed folder that match the ID segment (with
   `{id}` as the capture) give the IDs, e.g. folder `100d` for segment
   `{id}`, or file `104m.cif` for segment `{id}.cif`.
3. The full path is built from the template for each ID. Nothing below
   the ID segment is read.

So `{mid}/{id}/clean_mmcif/{id}_updated.cif.gz` reads `entry/` and
each `entry/<mid>/` (about 1,300 reads for the full archive), and
`{mid}/{id}.cif` works on a plain PDB mirror (including
`datasets/dev/mmcif/`).

A template without `{id}` is an error (`SystemExit` with a message).

Calls made in the spec (not discussed with the user):

- Without `--ids`, `{id}` keeps the case found on disk, so a mirror
  with uppercase folder names still works. With `--ids`, IDs are
  lowercased, which is what the PDB and PDBe archives use.
- A pattern or template that finds nothing still writes both files
  (`entries_found: 0`) and says so on stdout. That is not an error.
- `--input-dir` is resolved to an absolute path before the paths are
  written, so `entries.txt` works from any folder.
- One of `--glob` / `--path-template` is required, with no default.
  The old `pandora ingest --input-dir X` becomes
  `--glob '*.cif'`. No compatibility period, as in the 2026-11-04
  change.
- An ID segment with another placeholder in it (`{mid}{id}`) is an
  error.

## Output folder

| File | Content |
|---|---|
| `entries.txt` | One absolute path per line, sorted, no header. The format `read_entry_paths` reads. |
| `ingestion.json` | `StepSummary`, with `step="ingestion"`; `policy` = the discovery options below; `input` = `--input-dir`; `counts = {"entries_found": n}`. |

No `failures.jsonl`: nothing is checked per entry, so `ingest` has no
per-entry failures. (Earlier in the chat this file was planned for
`--ids` misses; checking them would mean one stat per ID, which goes
against the decision above.)

`--overwrite` and the refusal to write into a folder that holds step
outputs work as in `canonicalise`, but for `entries.txt` and
`ingestion.json`. `canonicalise --overwrite` must **not** start deleting
an `entries.txt` it finds in its own output folder.

## Linking the steps

When `canonicalise` reads `--entries FILE` and `FILE`'s folder holds
`ingestion.json`, its summary's `previous_summary` points there, just
as `curate` points to `canonicalisation.json`. The chain becomes
ingestion → canonicalisation → curation.

## `fetch`

`pandora fetch` also writes `entries.txt` in its `--output-dir` (the
paths of the files it saved, sorted). Its `ingestion_provenance.json`
is unchanged: per-entry records make sense there, because each file
comes from its own download.

## Library side

Per CLAUDE.md the CLI stays a thin wrapper. In `pandora/ingestion/paths.py`,
next to `read_entry_paths` and `mirror_path`:

- `template_entry_paths(root, template, ids=None) -> list[Path]`: the
  algorithm above.
- `write_entry_paths(paths, path) -> None`: the counterpart of
  `read_entry_paths`.

Both are public (docstrings, `__all__`). `--glob` is
`sorted(root.glob(pattern))` and needs no function of its own.

The discovery options go in a pydantic model in
`pandora/schemas/ingestion.py` (next to `FetchOptions`), so
`build_step_summary` can embed them as `policy`. Fields: `input_dir`,
`glob`, `path_template`, `ids_file`, `source_label`, `release_date`.
They come from CLI flags, not a YAML policy file, because they describe
one run, not a reusable rule set.

`ingest_local_mmcif()` stays as it is; the CLI stops calling it.

## Non-goals

- Checking that listed files exist, or reading anything inside them.
- Parallel listing, sharding or restart (see the "no parallelisation
  yet" decision for PPI v2.1).
- Presets for named archives (`--layout pdbe-release`). The template
  covers them.
- Remote listings (FTP/HTTP indexes, holdings files).
- Changes to `manifest`, `reproduce` or `IngestionProvenance`.
- `--entries` on `curate`; it reads the `canonicalise` step folder.

## Tests and docs

- `tests/test_ingestion.py`: `template_entry_paths` with and without
  `ids` on a small tree under `tmp_path` that mirrors the PDBe layout
  (gzipped `clean_mmcif/*_updated.cif.gz` plus a stray file the
  template must not pick up), and `{mid}/{id}.cif` on
  `datasets/dev/mmcif/`; `write_entry_paths` round-trips through
  `read_entry_paths`.
- `tests/test_cli.py`: `ingest` with `--glob` and with
  `--path-template`, then `canonicalise --entries` on its output,
  checking `previous_summary` and the entry-ID file names; `fetch`
  writing `entries.txt` with downloads faked as in the existing fetch
  tests. Revert the `--glob` test added to `canonicalise`/`curate`.
- `docs/usage/ingestion.md`: rewrite "Ingest already-downloaded files"
  with a library tab and a cli tab, real outputs on
  `datasets/dev/mmcif/`, tested by `tests/test_docs_examples.py`. Add
  the PDBe template as a non-executed example. Update
  `docs/usage/cli.md`, `pandora/cli/README.md` and
  `docs/contributing/architecture.md`; build the site.
