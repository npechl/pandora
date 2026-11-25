# CLI Namespace Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the 16 flat `pandora` commands with 7 noun groups
(`entries`, `structures`, `annotations`, `similarity`, `clusters`,
`splits`, `pairs`) that share one input/output convention.

**Architecture:** One CLI module per group (`pandora/cli/<group>.py`
with `register()`), assembled by `pandora/cli/app.py`. Shared input
resolution, flag blocks and output-folder handling live in
`pandora/cli/_io.py`. Library code outside `pandora/cli/` does not move;
the only library-facing change is `pairs positives` building its groups
with the existing `cluster_edges`.

**Tech Stack:** Python 3, argparse, pydantic v2, pytest, uv, Ruff,
Zensical.

**Spec:** `.superpowers/specs/2026-11-19-cli-namespace-design.md`.
Read it in full first: the argument table, the per-command changes and
the inputs table are normative and not repeated here.

## Global Constraints

- Work on branch `cli-namespace`, never on `main` **(user)**.
- argparse only. One module per group: `entries.py`, `structures.py`,
  `annotations.py`, `similarity.py`, `clusters.py`, `splits.py`,
  `pairs.py`, each with `register(subparsers)`.
- `pandora <group> <command> MAIN_INPUT... OUT [options]`: main inputs
  and `OUT` are positional, every other input is a named flag named after
  its data type.
- Every `OUT` is a folder holding its data under the fixed name from the
  spec's inputs table, plus `summary.json` (a `StepSummary`). The only
  exception is `structures convert`.
- Bad input exits with status 1 and one line starting
  `<group> <command>: `. argparse usage errors stay at 2.
- One result line on stdout, warnings on stderr.
- Every group, command and flag has `help=` text.
- No deprecated aliases for old command names.
- `from __future__ import annotations`, type hints, line length 80,
  `ruff format`, ruff's complexity limits (extract helpers). Internal
  helpers start with `_` and get a one-line docstring.
- No new dependencies. Tests use local fixtures only (`datasets/dev/mmcif/`,
  fake binaries as in `tests/test_sequence_similarity.py`).
- `examples/` is git-ignored: edit `examples/ppi/v2.1/src/*.sh`, never
  commit them.
- Version becomes `0.8.0`.

## Review Focus

1. **Folder vs bare file on every input.** `splits assign positives/`
   must find `positives/clusters.json`; `clusters representatives
   seq90/clusters.json` must work like `seq90/`. Tests in Tasks 1, 6, 7.
2. **Similarity input with no `summary.json` and no
   `--similarity-items`** exits 1 naming `--similarity-items`, and does
   not guess from the IDs. Test in Task 1.
3. **Re-running `entries fetch` with overlapping IDs and one failing ID
   under `--allow-partial`.** `entries.txt` is the union of old and new
   entries, without the failed one, and earlier provenance is kept. Test
   in Task 2.
4. **`--overwrite` on a structures `OUT`** removes the whole
   `structures/` tree and the step's tables, and leaves unrelated files
   in `OUT`. Test in Task 3.
5. **Running from another working directory.** Layers no longer trace
   their origin, so a command given only absolute or relative explicit
   paths works from any cwd (replaces
   `test_annotate_filter_contacts_from_another_cwd`). Test in Task 4.

---

### Task 0: Branch and commit the spec and plan

**Files:**
- Add: `.superpowers/specs/2026-11-19-cli-namespace-design.md`,
  `.superpowers/plans/2026-11-24-cli-namespace.md`

- [ ] **Step 1:** `git switch -c cli-namespace`
- [ ] **Step 2: Commit**

```bash
git add .superpowers/specs/2026-11-19-cli-namespace-design.md .superpowers/plans/2026-11-24-cli-namespace.md
git commit -m "docs(superpowers): CLI namespace spec and plan"
```

---

### Task 1: Shared input and output helpers

Everything later tasks share. Old commands keep working through this
task: the summary reader accepts the legacy names until Task 9 removes
that fallback.

**Files:**
- Modify: `pandora/cli/_io.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces (all in `pandora/cli/_io.py`):
  - `SUMMARY = "summary.json"`
  - `read_summary(folder: Path) -> dict[str, Any] | None`: parsed
    `folder/summary.json`. Until Task 9 it falls back to `layer.json`
    and the names in `STEP_SUMMARIES`; returns None when none exists.
  - `write_summary(out: Path, summary: StepSummary) -> None`: writes
    `out/summary.json`.
  - `input_summary(folder: Path) -> Path | None`: the summary file of an
    input folder, for `previous_summary` (a record only, never followed).
  - `resolve_file(command: str, path: Path, name: str) -> Path`: `path`
    when it is a file; `path / name` when `path` is a folder holding it;
    otherwise exits `"{command}: {path} has no {name}"`.
  - `resolve_table(command: str, path: Path, name: str) -> Path`: same,
    for a table: a `.parquet`/`.jsonl` file, or `find_table(path, name)`.
  - `require_tables(command: str, folder: Path) -> tuple[Path, Path]`:
    the `entries` and `chains` tables of a structures folder; otherwise
    exits `"{command}: {folder} has no chains table; run \`pandora
    structures curate\` first"`.
  - `add_out(p: ArgumentParser, *, tables: bool = True) -> None`: adds
    positional `out`, `--overwrite`, and (when `tables`)
    `--table-format {parquet,jsonl}` (default `parquet`).
  - `add_similarity_inputs(p: ArgumentParser, *, positional: bool) ->
    None`: positional `similarity` (`nargs="+"`) or repeatable
    `--similarity`, plus repeatable `--similarity-items
    {entry,chain,interface}`, `--hit-filter`, `--min-score`,
    `--min-coverage`.
  - `SimilarityInput` (a `NamedTuple`): `path: Path`, `table: Path |
    None` (None: an empty layer), `items: SimilarityItems`,
    `hit_filter: HitFilter`.
  - `similarity_inputs(command: str, args: Namespace) ->
    list[SimilarityInput]`: one per similarity input. `items` comes from
    the folder's summary, or from `--similarity-items` (once for all, or
    once per input in order); a disagreement exits; a bare table or a
    folder without a summary with no `--similarity-items` exits naming
    the flag. Hit filters as today's `_hit_filters`.
  - `prepare_output_dir(command: str, out: Path, outputs: Sequence[str],
    overwrite: bool) -> None`: today's function with the command prefix
    in its message.
- Moves here from `pandora/cli/similarity.py` (renamed public, same
  behaviour): `per_layer`, `hit_filters`, `layer_edges`.
- `structure_files` also finds `*.cif.gz` in flat folders.

- [ ] **Step 1: Write the failing tests**

```python
def test_resolve_file_accepts_folder_or_file(tmp_path):
    (tmp_path / "clusters.json").write_text("[]")
    assert _io.resolve_file("x y", tmp_path, "clusters.json") == tmp_path / "clusters.json"
    assert (
        _io.resolve_file("x y", tmp_path / "clusters.json", "clusters.json")
        == tmp_path / "clusters.json"
    )
    with pytest.raises(SystemExit, match="x y: .* has no clusters.json"):
        _io.resolve_file("x y", tmp_path / "empty", "clusters.json")


def test_require_tables_names_curate(tmp_path):
    with pytest.raises(SystemExit, match="run `pandora structures curate` first"):
        _io.require_tables("similarity mmseqs", tmp_path)


def test_similarity_inputs_bare_table_needs_items(tmp_path):
    # a layer.parquet copied out of its folder, no --similarity-items
    with pytest.raises(SystemExit, match="--similarity-items"):
        _io.similarity_inputs(
            "clusters build",
            _ns(similarity=[bare], similarity_items=None, min_score=[0.5]),
        )


def test_similarity_inputs_items_disagreeing_with_summary_exits(tmp_path): ...
def test_similarity_inputs_items_once_or_per_input(tmp_path): ...
def test_structure_files_finds_gzipped_flat_files(tmp_path): ...
```

(`_ns` builds an `argparse.Namespace` with the flag defaults; `bare` is
an entry-level `layer.parquet` written by `_similarity_layer` and moved
out of its folder.)

- [ ] **Step 2:** `uv run pytest tests/test_cli.py -k "resolve_file or require_tables or similarity_inputs or gzipped" -v`. Expected: FAIL (names not defined).
- [ ] **Step 3:** Implement the interfaces above in `pandora/cli/_io.py`.
  Update today's callers of `prepare_output_dir`, `_hit_filters`,
  `_per_layer` and `_layer_edges` to the moved names so the suite stays
  green.
- [ ] **Step 4:** `uv run pytest`. Expected: all pass.
- [ ] **Step 5: Commit** `refactor(cli): shared input, output and similarity-flag helpers`

---

### Task 2: `entries` group

**Files:**
- Create: `pandora/cli/entries.py` (from `pandora/cli/ingestion.py`)
- Delete: `pandora/cli/ingestion.py`
- Modify: `pandora/cli/app.py` (register `entries`, drop `ingestion`)
- Test: `tests/test_cli.py` (`_ingest`, `test_ingest_*`,
  `test_fetch_writes_entries_txt` move to `entries scan` / `entries fetch`)

**Interfaces:**
- Produces: `entries fetch ID... OUT`, `entries scan ROOT OUT` per the
  spec's argument table. Both write `entries.txt` and `summary.json`
  (step names `fetch` and `ingestion`).

Behaviour, from the spec's "entries" section:
- `fetch`: IDs from positionals (`nargs="*"`) and/or `--ids FILE` (one
  per line, `#` comments). Neither → exit 1. Into a folder that already
  holds `entries.txt`: merge it and `ingestion_provenance.json` (this
  run wins on the same ID); no `--overwrite` flag. A fetch error without
  `--allow-partial` exits 1 with `entries fetch: <ID>: <error>` (no
  traceback). Summary counts: `fetched` (this run), `failed`,
  `entries_total`.
- `scan --ids` works with `--glob`: keep the files whose name before its
  first `.` equals an ID (case-insensitive). Remove today's
  `--ids needs --path-template` check.

- [ ] **Step 1: Write the failing tests**
  - `test_entries_fetch_adds_to_existing_folder`: fetch `104m`, then
    `112m` into the same `OUT` (monkeypatch `fetch_mmcif` as today's
    `test_fetch_writes_entries_txt` does); `entries.txt` lists both;
    `summary.json` counts `{"fetched": 1, "entries_total": 2, ...}`.
  - `test_entries_fetch_partial_rerun_keeps_previous` (Review Focus 3):
    second run `112m bad1 --allow-partial`; `entries.txt` has 104m and
    112m only; `ingestion_provenance.json` still has 104M.
  - `test_entries_fetch_without_ids_exits`: `main(["entries", "fetch", str(out)])`
    raises `SystemExit` matching `entries fetch: `.
  - `test_entries_fetch_failure_is_one_line`: failing fetch without
    `--allow-partial` → `SystemExit` whose message starts `entries fetch: BAD1`.
  - `test_entries_scan_glob_with_ids`: `--glob '*/*.cif' --ids ids.txt`
    keeps only the listed entries.
  - Migrate the existing `test_ingest_*` tests; delete
    `test_ingest_rejects_ids_without_template`.
- [ ] **Step 2:** run them; expected FAIL (`invalid choice: 'entries'`).
- [ ] **Step 3:** implement `pandora/cli/entries.py`; register it in `app.py`.
- [ ] **Step 4:** `uv run pytest`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(cli)!: entries group (fetch adds to its folder, scan --ids with --glob)`

---

### Task 3: `structures` group

**Files:**
- Create: `pandora/cli/structures.py` (from `canonicalisation.py`,
  `datasets.py`, `export.py`)
- Delete: `pandora/cli/canonicalisation.py`, `pandora/cli/datasets.py`,
  `pandora/cli/export.py`
- Modify: `pandora/cli/app.py`, `pandora/cli/_io.py` (`prepare_step_dir`
  clears `summary.json`)
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: `structures canonicalise STRUCTURES OUT --policy`,
  `structures curate STRUCTURES OUT --policy`,
  `structures drop-contained STRUCTURES OUT --clusters [--policy]
  [--disable]`, `structures convert FILE FILE`. Step folders as today,
  with `summary.json` in place of `canonicalisation.json`,
  `curation.json` and `dedup.json`.
- Consumes: `resolve_file`, `require_tables`, `add_out`,
  `write_summary`, `input_summary` (Task 1).

Behaviour:
- `STRUCTURES` for `canonicalise`/`curate`: a structures folder, a flat
  or mirror folder of `.cif`/`.cif.gz`, an `entries.txt` list, or one
  `.cif`/`.cif.gz` file (replaces `--input-dir` / `--entries`).
- `curate` reads `chain_copies` from its input folder when present (as
  today).
- `drop-contained` takes `--clusters` as a folder or `clusters.json`,
  and needs `require_tables` on its input.
- Plain `dedup` is deleted, with `deduplicate_structures` CLI wiring and
  `dedup_provenance.json` / `dedup_exclusions.json` outputs. (The
  library function stays.)

- [ ] **Step 1: Write the failing tests**
  - `test_structures_canonicalise_single_file_input`: `STRUCTURES` is
    `FIXTURES_DIR / "04/104m.cif"`; `OUT/structures/04/104m.cif` and
    `OUT/summary.json` exist.
  - `test_structures_drop_contained_needs_tables`: a bare `.cif` folder
    exits with `structures drop-contained: ... run \`pandora structures curate\` first`.
  - `test_structures_drop_contained_accepts_clusters_folder`: `--clusters`
    given as the folder holding `clusters.json`.
  - `test_structures_overwrite_keeps_unrelated_files` (Review Focus 4):
    curate into `OUT` that also holds `notes.txt`; rerun with
    `--overwrite`; `structures/` is rebuilt, `notes.txt` is untouched.
  - `test_plain_dedup_is_gone`: `main(["dedup", ...])` exits 2.
  - Migrate the canonicalise/curate/dedup-contained/export tests and the
    `_canonicalise`, `_canonicalise_and_curate` helpers. Delete the
    plain-dedup assertions in `test_canonicalise_curate_dedup_export_chain`.
- [ ] **Step 2:** run; expected FAIL.
- [ ] **Step 3:** implement `pandora/cli/structures.py`; register it.
- [ ] **Step 4:** `uv run pytest`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(cli)!: structures group; drop plain dedup`

---

### Task 4: `annotations` group

**Files:**
- Create: `pandora/cli/annotations.py` (from `annotate.py`)
- Delete: `pandora/cli/annotate.py`
- Modify: `pandora/cli/app.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces the eight commands of the spec's argument table. Every layer
  folder holds `layer.<fmt>` and `summary.json` (no `layer.json`).
- Consumes: Task 1 helpers. `interface-gaps --structures` accepts any
  structures form (entry point); `dimer-length --structures` needs
  `require_tables`.

Behaviour:
- `contacts` / `ligand-contacts` call `annotate_chain_interfaces` /
  `annotate_ligand_contacts` with their own flags only (no
  "does not apply to --layer" check).
- `interface-contacts`: `label_interfaces(rows, {}, min_residues, None)`;
  `parameters.labels = ["ok", "few_contacts"]`.
- `dimer-length`: `label_interfaces(rows, lengths, None, max_total)`
  with lengths from the chains table; `parameters.labels = ["ok",
  "too_large"]`; keep the `unknown_length` count.
- `filter CONTACTS OUT --labels L [--labels L2 ...] --drop X [...]`: each
  `--labels` must be an `interface_labels` layer; labels from all layers
  are combined per entry before `drop_labelled_interfaces`; the typo
  check uses the union of present and known labels. Origin check: when
  every input has a summary, compare their recorded `input`; differing
  → warning on stderr; otherwise no check. `layer_origin` is no longer
  called.

- [ ] **Step 1: Write the failing tests**
  - `test_annotations_interface_contacts_labels_few_contacts` and
    `test_annotations_dimer_length_labels_too_large`: split from
    `test_annotate_label_interfaces_writes_layer`, same fixture, same
    expected labels per half.
  - `test_annotations_dimer_length_needs_tables`.
  - `test_annotations_filter_two_label_layers_one_pass`: contacts,
    `interface-contacts` labels and `interface-gaps` labels; one
    `filter` with both `--labels` and `--drop few_contacts --drop
    interface_gap` gives the same rows as the two-pass run.
  - `test_annotations_filter_works_from_another_cwd` (Review Focus 5):
    replaces `test_annotate_filter_contacts_from_another_cwd`.
  - Delete `test_annotate_compute_rejects_flag_for_other_layer`,
    `test_annotate_filter_contacts_rejects_different_origins`,
    `test_annotate_label_interfaces_untraceable_origin_exits`,
    `test_annotate_old_flat_form_is_gone`; migrate the rest.
- [ ] **Step 2:** run; expected FAIL.
- [ ] **Step 3:** implement `pandora/cli/annotations.py`; register it.
- [ ] **Step 4:** `uv run pytest`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(cli)!: annotations group; split label-interfaces; repeatable filter --labels`

---

### Task 5: `similarity` group

**Files:**
- Modify: `pandora/cli/similarity.py` (keep only `mmseqs`, `foldseek`,
  `interface-overlap`, `import`; cluster/representatives/partition code
  moves in Tasks 6–7)
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: similarity layer folders `layer.parquet` + `summary.json`,
  whose summary policy records `items` (read by `similarity_inputs`).

Behaviour:
- `mmseqs` and `foldseek` are separate parsers; `--bin` defaults
  `mmseqs` / `foldseek`; `--sensitivity`/`--max-seqs` defaults 5.7/300
  and 9.5/1000; `--alignment-type`, `--exhaustive-search` and
  `--chain-files DIR` (today's `--input-dir`) are foldseek-only.
  `--chains` with `--chain-files` still exits (today's check).
- `--chains PATH`: a representatives folder or `.json` (use its
  `representative` IDs) or a text list.
- `interface-overlap CONTACTS OUT --chain-clusters --structures
  --min-identity`: sequences from `require_tables(--structures)`;
  `--chain-clusters` via `resolve_file(..., "clusters.json")`.
- `import TABLE OUT --description ... (--structures | --contacts)`.

- [ ] **Step 1: Write the failing tests**
  - `test_similarity_chains_accepts_representatives_folder`: a
    `representatives.json` with two representatives limits the MMseqs2
    search to those chains (fake binary as in today's mmseqs tests).
  - `test_similarity_mmseqs_has_no_alignment_type`: `--alignment-type`
    on `similarity mmseqs` exits 2.
  - `test_similarity_interface_overlap_needs_structures_tables`.
  - Delete `test_similarity_search_mmseqs_rejects_input_dir`,
    `test_similarity_interfaces_untraceable_origin_exits`; migrate the rest.
- [ ] **Step 2:** run; expected FAIL.
- [ ] **Step 3:** implement.
- [ ] **Step 4:** `uv run pytest`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(cli)!: similarity group (mmseqs, foldseek, interface-overlap, import)`

---

### Task 6: `clusters` group

**Files:**
- Create: `pandora/cli/clusters.py` (`_cmd_cluster`,
  `_cmd_representatives` and their helpers from `similarity.py`)
- Modify: `pandora/cli/similarity.py`, `pandora/cli/app.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: `clusters build SIMILARITY... OUT` writing `clusters.json` +
  `summary.json` (policy: today's `ClusteringProvenance`);
  `clusters representatives CLUSTERS OUT` writing
  `representatives.json` + `summary.json` (policy: today's
  `RepresentativesProvenance`).
- Consumes: `similarity_inputs`, `layer_edges`, `resolve_file`,
  `require_tables` (Task 1).

Behaviour:
- `build`: no `--pairs`; `--items` keeps today's roll-up rules, compared
  against each input's level from `similarity_inputs`.
- `representatives`: no `--items`; interface level when `--contacts` is
  given, else chain; `--structures` required for both levels (resolution
  and copies come from its chains table, replacing `layer_origin`).
  `--similarity` inputs must be at that level.

- [ ] **Step 1: Write the failing tests**
  - `test_clusters_build_writes_folder`: `OUT/clusters.json`,
    `OUT/summary.json` with `step == "cluster"`.
  - `test_clusters_representatives_accepts_bare_clusters_file` (Review
    Focus 1).
  - `test_clusters_representatives_interface_level_from_contacts`:
    `--contacts` + `--structures`, no `--items`.
  - Delete `test_cluster_with_pairs`,
    `test_representatives_interface_items_need_contacts`,
    `test_representatives_interfaces_reject_structures_dir`; migrate the rest.
- [ ] **Step 2:** run; expected FAIL.
- [ ] **Step 3:** implement `pandora/cli/clusters.py`; register it.
- [ ] **Step 4:** `uv run pytest`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(cli)!: clusters group; drop --pairs and representatives --items`

---

### Task 7: `splits` group

**Files:**
- Create: `pandora/cli/splits.py` (`_cmd_partition`, `_partition_labels`)
- Modify: `pandora/cli/similarity.py`, `pandora/cli/app.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: `splits assign CLUSTERS OUT` writing `splits.json` +
  `summary.json` (policy: `PartitionProvenance`).

Behaviour: `--labels L` and `--stratify L` in an argparse mutually
exclusive group; `--stratify L` calls `partition_dataset(...,
labels=..., stratify=True)`. No `--no-keep-similar-items`
(`keep_similar_items` stays at its default `True`). `CLUSTERS` via
`resolve_file(..., "clusters.json")`; `L` via `resolve_table(..., "layer")`.

- [ ] **Step 1: Write the failing tests**
  - `test_splits_assign_takes_positives_folder` (Review Focus 1): a
    folder holding `clusters.json` and `kinds/`;
    `splits assign POS OUT --stratify POS/kinds` writes `OUT/splits.json`.
  - `test_splits_labels_and_stratify_exclusive`: both flags → exit 2.
  - Delete `test_partition_labels_rejects_no_keep`,
    `test_partition_stratify_needs_labels`; migrate the rest.
- [ ] **Step 2:** run; expected FAIL.
- [ ] **Step 3:** implement.
- [ ] **Step 4:** `uv run pytest`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(cli)!: splits group; always leakage-safe`

---

### Task 8: `pairs` group

**Files:**
- Create: `pandora/cli/pairs.py` (from `sampling.py`)
- Delete: `pandora/cli/sampling.py`
- Modify: `pandora/cli/app.py`
- Test: `tests/test_cli.py`

**Interfaces:**
- Produces: `pairs positives SIMILARITY... OUT` writing `positives.<fmt>`,
  `dropped.<fmt>`, `clusters.json`, `kinds/` (`layer.<fmt>` +
  `summary.json`), `summary.json`; `pairs negatives POSITIVES OUT`;
  `pairs assemble POSITIVES OUT`.
- Consumes: `similarity_inputs`, `layer_edges`, `resolve_file`,
  `resolve_table`, `require_tables` (Task 1).

Behaviour:
- `positives`: dimers are all interface items of `--contacts`; build the
  adjacency with `dimer_adjacency` as today, then the groups with
  `cluster_edges(dimers, ((a, b) for a, ns in adjacency.items() for b in ns))`,
  then `rebalance(groups, ...)`. `--structures` (tables) replaces
  `layer_origin` for interface values. No `--clusters`.
- `negatives` / `assemble`: `POSITIVES` via `resolve_table(...,
  "positives")`; `--negatives`, `--splits`, `--representatives` via
  `resolve_file`/`resolve_table`. The "differs from what negatives used"
  warning reads the negatives `summary.json` when present.

- [ ] **Step 1: Write the groups-equivalence test first**

```python
def test_pairs_positives_groups_match_clusters_build(tmp_path):
    # same contacts, seq90 representatives, interface and chain layers
    # and hit filters as examples/ppi/v2.1 steps 5-7, on the 1a2c fixture
    main(
        [
            "clusters",
            "build",
            iface,
            chain,
            str(grp),
            "--items",
            "interface",
            ...,
        ]
    )
    main(["pairs", "positives", iface, chain, str(pos), ...])
    built = json.loads((grp / "clusters.json").read_text())
    kept = {p["dimer"] for p in _jsonl(pos / "positives.jsonl")}
    restricted = [sorted(set(c["components"]) & kept) for c in built]
    assert sorted(c for c in restricted if c) == sorted(
        sorted(c["components"]) for c in json.loads((pos / "clusters.json").read_text())
    )
```

  If it fails on the fixture, stop and report to the user before going
  on: the `pairs` redesign rests on it.
- [ ] **Step 2: More failing tests**
  - `test_pairs_positives_writes_clusters_json_and_kinds` (replaces
    `test_rebalance_writes_groups_and_kinds`).
  - `test_pairs_positives_needs_structures_tables`.
  - Delete `test_rebalance_inputs_mismatch_exits`,
    `test_rebalance_dimer_missing_from_contacts_exits`,
    `test_negatives_needs_rebalance_folder`; rewrite `_rebalance_inputs`
    / `_rebalance` without `groups`; migrate the rest.
- [ ] **Step 3:** run; expected FAIL.
- [ ] **Step 4:** implement `pandora/cli/pairs.py`; register it.
- [ ] **Step 5:** `uv run pytest`. Expected: all pass.
- [ ] **Step 6: Commit** `feat(cli)!: pairs group; positives builds its own groups`

---

### Task 9: Remove `manifest`, finish the cutover

**Files:**
- Delete: `pandora/cli/provenance.py`
- Modify: `pandora/cli/app.py` (module order: entries, structures,
  annotations, similarity, clusters, splits, pairs; `app()` docstring
  lists the groups), `pandora/cli/_io.py` (delete `layer_origin`,
  `load_pairs` if unused, `STEP_SUMMARIES` and the legacy-name fallback
  in `read_summary`), `pyproject.toml` (`version = "0.8.0"`), `uv.lock`
  (`uv lock`)
- Test: `tests/test_cli.py`

- [ ] **Step 1: Write the failing test**

```python
_GROUPS = ["entries", "structures", "annotations", "similarity",
           "clusters", "splits", "pairs"]

def test_parser_registers_the_groups():
    help_text = _build_parser().format_help()
    for name in _GROUPS:
        assert name in help_text
    for gone in ("manifest", "reproduce", "rebalance", "partition", "ingest"):
        with pytest.raises(SystemExit) as exc:
            main([gone])
        assert exc.value.code == 2

def test_every_flag_has_help():
    # walk every subparser's actions; assert action.help for each
    # non-positional, non-help action
```

  Delete `_SUBCOMMANDS`, `test_parser_registers_all_subcommands`,
  `test_manifest_*`, `test_reproduce_*`.
- [ ] **Step 2:** run; expected FAIL.
- [ ] **Step 3:** implement the deletions and `app.py` changes. `grep -rn
  "layer.json\|canonicalisation.json\|curation.json\|rebalance.json\|negatives.json\|assemble.json\|_provenance.json" pandora tests`
  finds nothing that names an old summary (`ingestion_provenance.json`
  stays).
- [ ] **Step 4:** `uv run pytest && uv run ruff check . && uv run ruff format --check .`. Expected: all pass.
- [ ] **Step 5: Commit** `feat(cli)!: remove manifest group; 0.8.0`

---

### Task 10: Docs and examples

**Files:**
- Modify: `pandora/cli/README.md` (rewrite around the tree and the I/O
  convention: positional inputs/OUT, folder-or-file table, two kinds of
  structures input, `summary.json`, shared flags), `docs/usage/cli.md`
  (new commands; remove the `manifest`/`reproduce` sections),
  `docs/usage/*.md` pages that show CLI calls (`grep -rn "pandora "
  docs`), `docs/usage/provenance.md` (drop the CLI mention),
  `examples/API_SURFACE.md`
- Modify, do not commit: `examples/ppi/v2.1/src/*.sh` and
  `examples/ppi/v2.1/README.md` (drop the `jq` step, the second
  `filter`, the dimer-group `clusters build`; `splits assign` takes the
  positives folder)

- [ ] **Step 1:** update the docs and examples.
- [ ] **Step 2:** `grep -rnE "pandora (fetch|ingest|canonicalise|curate|dedup|export|annotate|cluster|representatives|partition|rebalance|negatives|assemble|manifest|reproduce)\b" docs pandora/cli/README.md examples`. Expected: no matches.
- [ ] **Step 3:** run the docs-build skill (`uv run zensical build --clean`). Expected: builds with no broken links.
- [ ] **Step 4:** `bash -n examples/ppi/v2.1/src/*.sh`. Expected: no syntax errors.
- [ ] **Step 5: Commit** `docs(cli): document the group namespace` (docs and `pandora/cli/README.md` only).
