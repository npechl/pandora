---
name: ""          # short identifier, e.g. ppi-benchmark
version: ""       # e.g. "2.1"
task: ""          # one line: what a model predicts, e.g. "Given two protein chains, predict whether they interact"
description: >
  2–4 sentences: what the dataset is for, what one example is,
  where the data comes from, and whether it's for training, testing or both.
tags: []          # e.g. [protein-protein-interaction, binary-classification, structure, pdb]
                  # suggested kinds: task type, biological problem, input modality, data source
---

# <Dataset name> (v<version>)

<!--
How to use this template:
- Replace every <placeholder> and delete the HTML comments.
- Delete a section only if it truly doesn't apply; otherwise write "None".
- Anything not decided yet goes in the Status column as "open", and in "Open decisions" at the end.
-->

## Overview

- **Task:** <input → output, e.g. "two protein chains → interact yes/no">
- **One example is:** <e.g. "a pair of chains plus their interface residues">
- **Labels:** <what positives/negatives or target values are, and how each is obtained>
- **Intended use:** <train + test / test-only; which kind of model it targets>
- **Size (after all steps):** <number of examples, class balance>

## Glossary

<!-- Domain terms a non-specialist needs to follow the steps. -->

| Term | Meaning |
|------|---------|
| <term> | <plain-language meaning> |

## Data sources

| Source | What is taken | Version / snapshot date | Access |
|--------|---------------|-------------------------|--------|
| <e.g. PDB> | <e.g. biological assembly 1, mmCIF> | <e.g. 2026-09-01> | <URL or API> |

## Steps

<!-- Copy this block once per step. -->

### Step 1: <short name>

**Goal:** <one sentence: why this step exists>

**Input:** <what it receives, and from which step> — **Output:** <what it produces>

**Procedure:**
1. <precise, ordered actions; say what is removed, kept or computed, and in what order>

**Parameters:**

| Parameter | Value | Status | Why this value |
|-----------|-------|--------|----------------|
| <e.g. max resolution (X-ray)> | <2.5 Å> | fixed / open | <reason or reference> |

**Tools:** <tool name + version + non-default options, or "None">

> **In plain words:** <the same step explained for a non-specialist>

**Open questions:** <or "None">

## Splits

- **Method:** <how examples are assigned to train/validation/test, or "test-only">
- **Grouping unit:** <what is kept together to avoid leakage, e.g. similarity clusters from Step N>
- **Ratios:** <e.g. 80/10/10>
- **Negatives:** <which split each negative goes to>

## Output

**Files:** <paths and formats>

| Field | Type | Meaning |
|-------|------|---------|
| <field> | <type> | <meaning> |

## Reproducibility

- **Data snapshot:** <date/release of every source>
- **Tool versions:** <tool: version, one per line>
- **Random seeds:** <value and where it's used>
- **Tie-breaking rules:** <how ties are resolved wherever "pick one" happens>
- **Code:** <repo, commit or script that builds the dataset>

## Limitations

<!-- One item per limitation. -->

- **L1. <title>.** *Where:* <step>. *What goes wrong:* <problem>. *Effect:* <consequence for the dataset>. *Mitigation:* <suggestion, or "none">.

## Open decisions

- [ ] <decision> (Step N)

## Changelog

| Version | Date | Changes |
|---------|------|---------|
| <x.y> | <YYYY-MM-DD> | <what changed and why> |
