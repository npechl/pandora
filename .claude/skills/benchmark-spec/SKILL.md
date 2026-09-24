---
name: benchmark-spec
description: Write or enrich a benchmark dataset description (Markdown with a YAML header) from the project template in this skill's folder (assets/template.md). Use when the user asks to create, convert, fill in, update or document a benchmark/dataset spec, from notes, a workflow doc, a paper or a GitHub repo, e.g. "convert ppi-v2.1.md using the template", "turn these notes on Pinder into a benchmark file", "write a benchmark description for X from github.com/org/repo", or "fill the gaps in pinder-benchmark.md from the repo".
---

# Benchmark spec

Produce one Markdown file per benchmark, following `assets/template.md` in this skill's folder exactly: same sections, same order, same field names.

There are two modes:

- **New:** the source is notes, a workflow doc, a paper or a GitHub repo, and you write a new file. Follow "Steps".
- **Enrich:** the user gives an existing `*-benchmark.md` plus a new source, and you fill its gaps in place. Follow "Enriching an existing spec".

## Steps

1. **Read the template.** Take the section list and field names from it, not from memory.
2. **Read the whole source.** If it's a file, re-read it now; the user may have edited it since you last saw it. If it's a GitHub repo, follow "Reading a GitHub repo".
3. **Pick the output path.** Never overwrite the source.
   - Source is a file: `<source-dir>/<name>-benchmark.md`, or `<name>-v<version>-benchmark.md` when the source is a versioned workflow.
   - Source is pasted notes, a repo or a paper: ask the user which directory to use.
4. **Fill in every section** using the rules below.
5. **Check that the YAML header parses:**
   ```sh
   uv run python -c "import yaml,sys; print(yaml.safe_load(open(sys.argv[1]).read().split('---')[1]))" <file>
   ```
6. **Report back** in a few lines: the file path, the steps it contains, what came from outside the source (see "Honesty"), and the most important open decisions.

## Enriching an existing spec

1. **Read the spec and the new source.** For a repo, follow "Reading a GitHub repo".
2. **List every gap** in the spec: `not covered in …` / `not specified in …` text, `open` in a Status column, each **Open questions** entry, and every unchecked item under Open decisions.
3. **For each gap the new source answers:** replace the placeholder text with the answer plus a citation (see "Citations"). Set Status to `fixed` if the source states the value. Tick the matching Open decision (`- [x]`) and add `→ <answer> (<citation>)`.
4. **When the new source contradicts the spec,** keep the spec's text. Add the conflict as an open question citing both sources.
5. **Leave everything else alone.** That covers gaps the source doesn't answer, wording, step numbering and limitation IDs.
6. **Update the header.** Add the new source to the `*Source:*` line, bump `version` by a minor step (`1.0` → `1.1`), and add a Changelog row that names the source (with the commit SHA for a repo).
7. **Check the YAML** (Steps, item 5). Report: how many gaps were filled, how many are still open, and any conflicts.

## Reading a GitHub repo

1. **Shallow-clone** it into the scratchpad directory, then record the commit:
   ```sh
   gh repo clone <owner>/<repo> <scratchpad>/<repo> -- --depth 1
   git -C <scratchpad>/<repo> rev-parse --short HEAD
   ```
   If `gh` isn't authenticated, use `git clone --depth 1 https://github.com/<owner>/<repo>`.
2. **Read in this order,** and stop once every template section is covered or no file is left:
   1. README, and any docs/ folder or linked docs site.
   2. Config and parameter files (YAML/JSON/TOML), plus CLI and function defaults.
   3. Build or pipeline scripts. Cutoffs, filters and tie-breaks often appear only here.
   4. Split, index and manifest files. Their names give the Output files; row counts give Size.
   5. LICENSE, `CITATION.cff` and the paper link.
   If the repo is too large to read in this order, send an `Explore` agent to find where each open parameter is set. Ask it for `path:line` locations, and check each one yourself before citing it.
3. **Values found in code count as the source.** When the README and the code disagree, record the conflict; don't pick one.
4. Don't run the pipeline or download the dataset unless the user asks.

## Citations

- Notes or a workflow doc: no citation is needed; the file is the whole source.
- Repo: `(path:line@sha)`, e.g. `(src/pinder/core/utils.py:42@a1b2c3d)`.
- Paper: `(<first author> <year>, §3.2)`, with the DOI in Data sources or Reproducibility.

## Rules

**YAML header and source line**
- `name`: kebab-case identifier.
- `task`: one line, input → output.
- `description`: 2–4 sentences.
- `tags`: 4–7 kebab-case tags covering task type, biological problem, input modality and data source.
- `version`: the workflow's own version. When describing someone else's dataset, use `"1.0"` with the comment `# version of this description, not of the upstream dataset`.
- The `*Source:*` line under the title names what the description was built from: a notes file, `<owner>/<repo>@<sha>`, or a paper DOI.

**Steps**
- Keep the source's step numbering, so readers can cross-reference the two documents. When the source is a repo with no numbered steps, number the pipeline stages in execution order.
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
- Use the template's item format: `L<n>. <title>`, then *Where / What goes wrong / Effect / Mitigation*. All four are required; write `Mitigation: none.` rather than dropping it.
- When the source already has a limitations list (e.g. an earlier review), keep its IDs and wording instead of renumbering.
- Cover the dataset-building workflow only, not how the dataset is later used in training or evaluation.

**Open decisions**
- A checklist of every `open` item from the steps and sections, each tagged `(Step N)`.

**Changelog**
- One row per version. Use the date the version was written (commit date if known), or `—` when unknown. For a repo source, include the commit SHA.

## Honesty

- Procedures, parameters and tool names come **only from the source**: notes, repo files or paper. If the source is silent, write `not covered in <source>` and add the item to Open decisions. Don't fill it from general knowledge.
- The goal lines, plain-words boxes, glossary and limitations may use your own knowledge. Tell the user which items are your interpretation, not the source's.
- If the source contradicts itself (e.g. two different cutoffs for the same thing, or a README that disagrees with the code), don't pick one. Record the conflict as an open question, and say so in the report.
