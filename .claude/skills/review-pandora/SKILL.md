---
name: review-pandora
description: Use when the user asks to review Pandora for bugs, find gaps in its documentation or in the project, audit the codebase before a release, look for issues to file, or asks "what's broken", "what's missing", "find bugs", "review the docs", or "/review-pandora".
---

# Reviewing Pandora

Pandora's bugs hide in two places: where the code disagrees with what it
claims (docstrings, `docs/reference/policies.md`, policy fields), and where
tests only agree with themselves (synthetic round-trips that never touch a
real fixture). A review finds leads, confirms each one by reproducing it,
removes duplicates of existing issues, and reports. **It does not fix
anything.** Fixing or filing issues happens only when the user asks.

## 1. Run the scan

```bash
uv run python .claude/skills/review-pandora/scripts/scan.py           # ~3 min, all fixtures
uv run python .claude/skills/review-pandora/scripts/scan.py --limit 10 # quick pass
uv run python .claude/skills/review-pandora/scripts/scan.py --skip-sweep
```

If `uv run` fails with `No module named 'pandora'`, the venv is stale:
delete `.venv/` and run `uv sync --all-extras`.

Every line is a **lead**, not a finding. How to read each section:

| Section | Usually means | False positives |
|---|---|---|
| Policy fields no logic reads | Field is accepted and ignored | Already marked "unused" in `policies.md`: check before reporting |
| Policy fields missing from policies.md | Docs gap | None known |
| docs/policies YAML drift | Docs describe fields that don't exist | None known |
| Packages missing from functions.md | Public API not on the reference page | None known |
| Typed categories parsed as empty | Parser reads the wrong tags | A fixture where the category truly has no usable rows |
| Canonicalisation invariant failures | Stale references, empty entries, exports that don't re-parse | `raises ValueError` under `validation_rules.strictness=strict` or `warnings_as_errors` is the policy doing its job, unless the entry is clean |
| Policy values that never change the output | Ignored value | The value only matters together with another field (`keep_waters` needs `ligand_rules.strategy: filter`), only changes diagnostics (not returned, #69), or no fixture triggers it |

The sweep changes one canonicalisation field at a time from the defaults.
It doesn't cover curation, similarity, splits, provenance or the CLI;
lenses 4–6 below do.

## 2. Look through the lenses

Scan output covers what can be checked mechanically. Then read code for:

1. **Claims vs code.** Grep docstrings and comments for absolute words
   (`infinite`, `always`, `never`, `every`, `all`) and check the code keeps
   the promise at the edges (empty input, more than 26²
   ids, missing optional column).
2. **Tests that agree with themselves.** For each parser/exporter path,
   is there a test on a real fixture from `datasets/dev/mmcif/`? A test
   that writes then reads back with the same code proves nothing about
   real files.
3. **Docs vs behaviour.** For each "not implemented" / "currently unused"
   note in `docs/reference/policies.md`, check it's still true: stale
   notes are gaps too. Check the "Known issues" list in `CLAUDE.md` the same way.
4. **Anything that removes things.** Curation, dedup, ligand filtering,
   chain dropping: afterwards, does every reference (assemblies, bonds,
   secondary structure, `raw` categories, entities) still point at
   something that exists?
5. **Level mismatches.** Similarity, clustering, splits and
   `reproduce_dataset`: do the item ids (entry vs `chain_item_id`) match
   across the steps, and does the manifest record enough to rebuild?
6. **Docs that run.** Run the commands in `README.md` and
   `docs/getting-started/`, and the inline examples in `docs/usage/`.
   Check that every `pandora/cli` subcommand is documented.

## 3. Confirm each lead

- Write the smallest repro on a local fixture (a few lines of Python).
  A lead you can't reproduce gets dropped, not reported as "possible".
- To tell whether a lead is new, `git stash` the working tree and rerun.
- Search open and closed issues for it:
  `gh issue list --state all --search "<keyword>"`. Link the match
  instead of reporting it again. If `gh` fails, say so; don't skip this step.

## 4. Report

Ranked most severe first: wrong or lost data, then crashes, then ignored
settings, then docs.

| # | Kind | Finding | Evidence | Repro | Existing issue |
|---|---|---|---|---|---|
| 1 | bug / docs / gap | one sentence | `path:line` | command or snippet | `#N` or — |

After the table: how many leads you dropped and why (one line each), then
offer to open issues. New issues follow `.github/ISSUE_TEMPLATE/bug.md`
(Description, Expected, Actual, How to reproduce, Additional context) with
the `bug`, `enhancement` or `documentation` label.

## Common mistakes

- Reporting scan lines without a repro.
- Fixing bugs mid-review: the user decides what gets fixed.
- Reporting a value with "no effect" without checking whether it depends on another field.
- Filing a duplicate because the issue search was skipped.
- Counting a pre-existing problem as a regression from recent changes.
