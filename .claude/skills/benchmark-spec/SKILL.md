---
name: benchmark-spec
description: Write a benchmark dataset description (Markdown with a YAML header) from the project template in this skill's folder (assets/template.md). Use when the user asks to create, convert, fill in or document a benchmark/dataset spec, e.g. "convert ppi-v2.1.md using the template", "turn these notes on Pinder into a benchmark file", or "write a benchmark description for X".
---

# Benchmark spec

Produce one Markdown file per benchmark, following `assets/template.md` (in this skill's folder, `.claude/skills/benchmark-spec/assets/template.md`) exactly: same sections, same order, same field names.

## Steps

1. **Read the template** `assets/template.md` (in this skill's folder, `.claude/skills/benchmark-spec/assets/template.md`). Take the section list and field names from there, not from memory.
2. **Read the whole source**: a workflow doc, pasted notes, a paper summary or code. If the source is a file, re-read it now. The user may have edited it since you last saw it.
3. **Pick the output path.** Never overwrite the source. The default is `<source-dir>/<name>-benchmark.md`, or `examples/ppi/<name>-benchmark.md` for pasted notes. Use `<name>-v<version>-benchmark.md` when the source is a versioned workflow.
4. **Fill in every section** using the rules below.
5. **Check that the YAML header parses:**
   ```sh
   uv run python -c "import yaml,sys; print(yaml.safe_load(open(sys.argv[1]).read().split('---')[1]))" <file>
   ```
6. **Report back** in a few lines: the file path, the steps it contains, what came from outside the source (see "Honesty"), and the most important open decisions.

## Rules

**YAML header**
- `name`: kebab-case identifier.
- `task`: one line, input → output.
- `description`: 2–4 sentences.
- `tags`: 4–7 kebab-case tags covering task type, biological problem, input modality and data source.
- `version`: the workflow's own version. When describing someone else's dataset, use `"1.0"` with the comment `# version of this description, not of the upstream dataset`.

**Steps**
- Keep the source's step numbering, so readers can cross-reference the two documents.
- Every step has all the template's fields. Write `none` for a field that doesn't apply; never delete it.
- Rewrite the procedure as ordered, precise actions, keeping the source's parameter values and tool names.
- When a source rule was replaced (e.g. "(old rule) → new rule"), show the old rule struck through with `~~…~~`, and name the step that replaces it.
- Parameters table, `Status` column:
  - `fixed`: the source states the value.
  - `open`: the source marks it undecided ("?", "???", "check distribution") or leaves it out.
  - `choice`: a user-selectable option.
  - Add what's open inside the cell, e.g. `fixed; normalisation: open`.
- The **In plain words** box explains the step for a reader with no structural-biology background.

**Sections the source doesn't cover** (Splits, Output, Reproducibility, Data sources)
- Write `Not specified in <source>.`
- Then list what the steps imply, or what should be recorded, labelled as *implied* or *suggested*.

**Glossary**
- Define every domain term used in the steps.
- If a sibling benchmark file already defines a term, reuse its wording.

**Limitations**
- Use the template's item format: `L<n>. <title>`, then *Where / What goes wrong / Effect / Mitigation*.
- When the source already has a limitations list (e.g. an earlier review), keep its IDs and wording instead of renumbering.
- Cover the dataset-building workflow only, not how the dataset is later used in training or evaluation.

**Open decisions**
- A checklist of every `open` item from the steps and sections, each tagged `(Step N)`.

**Changelog**
- One row per version. Use the date the version was written (commit date if known). Use `—` when unknown.

## Honesty

- Procedures, parameters and tool names come **only from the source**. If the source is silent, write `not covered in notes` (or `not specified in <source>`), and add the item to Open decisions. Don't fill it from general knowledge.
- The goal lines, plain-words boxes, glossary and limitations may use your own knowledge. Tell the user which items are your interpretation, not the source's.
- If the source contradicts itself (e.g. two different cutoffs for the same thing), don't pick one. Record the conflict as an open question, and say so in the report.
