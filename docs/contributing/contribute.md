# Contributing

For environment setup, branching, and how to open a PR, see
[`CONTRIBUTING.md`](https://github.com/npechl/pandora/blob/main/CONTRIBUTING.md)
in the repo root. This page covers what to do for each kind of
contribution.

Before diving in, [Architecture](architecture.md) maps every public
function to the data model it takes and returns — the fastest way to see
where a new piece belongs.

## Where to start

You do not need structural-biology experience for every contribution.
Choose an area that matches your background:

| Background | Starting points | Examples to follow |
|---|---|---|
| Python, no mmCIF knowledge needed | [Similarity and splits](../usage/similarity.md): clustering prepared edges and partitioning clusters | [Clustering tests](https://github.com/npechl/pandora/blob/main/tests/test_clustering.py), [partition tests](https://github.com/npechl/pandora/blob/main/tests/test_partition.py) |
| Python, no mmCIF knowledge needed | [Export](../usage/export.md): formatting typed records as JSON, JSONL or Parquet | [Export tests](https://github.com/npechl/pandora/blob/main/tests/test_export.py) |
| Command-line tools, no mmCIF knowledge needed | [CLI](../usage/cli.md): argument handling and file I/O around existing functions | [CLI conventions](https://github.com/npechl/pandora/blob/main/pandora/cli/README.md), [CLI tests](https://github.com/npechl/pandora/blob/main/tests/test_cli.py) |
| Comfortable with Pandora's typed models and policies | [Curation rules](../usage/datasets.md) and [annotation layers](../usage/annotation.md); follow an existing rule or `annotate_*()` function | [Curation tests](https://github.com/npechl/pandora/blob/main/tests/test_curation.py), [annotation tests](https://github.com/npechl/pandora/blob/main/tests/test_annotations.py) |
| Familiar with mmCIF and structure identifiers | [Parsing and canonicalisation](../usage/canonicalisation.md); read the [glossary](../getting-started/glossary.md) and [architecture](architecture.md) first | Local fixtures in `datasets/dev/mmcif/` and the matching tests |

Clustering and partition tests can work with prepared edges and typed
records; they do not require running MMseqs2 or Foldseek. Parquet export
needs the optional `export` dependencies. The full test suite uses local
fixtures and fakes external similarity tools rather than requiring live
services or those binaries.

Browse [open good first issues](https://github.com/npechl/pandora/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22)
for scoped tasks, and check the discussion and linked PRs before starting
so you do not duplicate work already in progress. Use the workflow below
once you have chosen a task.

## Ways to contribute

| I want to... | Go to |
|---|---|
| Report something broken | [Report a bug](#report-a-bug) |
| Propose something new | [Request a feature](#request-a-feature) |
| Fix an existing bug | [Fix a bug](#fix-a-bug) |
| Add a function to `pandora/ingestion`, `parsing`, `canonicalisation`, `metadata`, `datasets`, `similarity`, `provenance`, or `export` | [Add a function to an existing component](#add-a-function-to-an-existing-component) |
| Add a new `annotate_*()` function | [Add a new annotation](#add-a-new-annotation) |
| Write a runnable script or a narrated write-up | [Write an example or recipe](#write-an-example-or-recipe) |
| Add a bundled fixture dataset | [Build a new dataset](#build-a-new-dataset) |
| Fix or extend the docs site | [Update the docs](#update-the-docs) |

## Report a bug

[Open an issue](https://github.com/npechl/pandora/issues/new?template=bug.md)
with the bug template: expected vs. actual behavior, steps to reproduce,
and your Pandora version/commit. Attach the mmCIF file or policy YAML if
the bug depends on one.

## Request a feature

[Open an issue](https://github.com/npechl/pandora/issues/new?template=feature_request.md)
with the feature template: what problem it solves, and a proposed
approach if you have one.

## Fix a bug

1. Find the function via [Functions](../reference/functions.md) or
   `grep -rn` in `pandora/`.
2. Write a failing test in `tests/` that reproduces it (or extend an
   existing one).
3. Fix it, then verify:

   ```sh
   uv run pytest
   uv run ruff format .
   uv run ruff check .
   ```

4. Follow `CONTRIBUTING.md` to open the PR.

## Add a function to an existing component

Every stage follows the same split: **schema in `pandora/schemas/`,
logic in `pandora/<package>/`** (see [Architecture](architecture.md)).

1. Add/extend the Pydantic model in `pandora/schemas/<module>.py`, if
   the function needs a new typed input or output.
2. Add the function in `pandora/<package>/<module>.py`. Keep it a pure
   function: take a `Structure` (or typed record) in, return a new one
   — never mutate in place.
3. Export it from `pandora/<package>/__init__.py`'s `__all__` if it's
   meant to be public API.
4. Docstring it: exported functions get a full description with
   `Args:`/`Returns:` (and `Raises:` if applicable); internal helpers
   get one brief line. This is what powers
   [Functions](../reference/functions.md) and [Schemas](../reference/schemas.md)
   — nothing else to do to get a new function or model documented there.
5. Add a test in `tests/`.
6. If you changed which schemas reference each other, regenerate the
   architecture diagrams:

   ```sh
   uv run --extra docs python docs/scripts/generate_erd.py
   ```

## Add a new annotation

Entry-level annotations live in `pandora/annotations/entry.py`,
pairwise ones in `pandora/annotations/pairwise.py`. Both follow the
same shape as the existing `annotate_*()` functions: take a
`Structure` (or two), return an `AnnotationLayer`
(`pandora/schemas/annotation.py`). Register it in
`pandora/annotations/__init__.py`'s `__all__`.

There's no formal plugin-registration system yet. If you want a plugin protocol that lets
annotations be registered/discovered instead of hand-imported, that's
an open, unclaimed piece of the design.

## Write an example or recipe

Examples in `examples/` are self-contained scripts that run against
the local fixtures in `datasets/dev/mmcif/` — no network or external
binaries required. Follow the existing scripts' shape (see
`examples/README.md`), then run it to confirm it works:

```sh
uv run python examples/your_script.py
```

A [recipe](../recipes/ppi-01.md) is a narrated, task-oriented write-up
("how do I build a PPI dataset") built on top of an example script —
[PPI interface dataset](../recipes/ppi-01.md) is the only one so far,
and a second one covering a different task (e.g. clustering by
structural rather than sequence similarity, or a single-organism
dataset) is a welcome contribution. To add one:

1. Write the runnable script under `examples/` first, if one doesn't
   already exist for the workflow.
2. Add a narrated page under `docs/recipes/` that walks through it
   step by step — see `ppi-01.md`'s source for the shape (prerequisites,
   a collapsible full script, expected output, then a per-step
   breakdown with links to the [Usage](../usage/overview.md) pages
   each step comes from).
3. Register the new page in `zensical.toml`'s `nav`, under the
   `"Recipes"` list — a page not listed there builds but won't appear
   anywhere in the site.

## Build a new dataset

Follow the pattern in `datasets/scripts/build_dev_dataset.py`: fetch
via `pandora.ingestion.fetch_list_mmcif`, write mmCIF files plus a
`manifest.json` under `datasets/<name>/`. Only the script and the
manifest are committed — raw mmCIF files are gitignored, so anyone can
regenerate the dataset from the script.

```sh
uv run python datasets/scripts/build_dev_dataset.py
```

## Update the docs

Preview changes locally, then verify the build is clean before
pushing (CI deploys with `zensical build --clean` on merge to `main`):

```sh
uv run zensical serve      # live preview at http://localhost:8000
uv run zensical build --strict
```

You don't need to hand-write API reference content — `Functions` and
`Schemas` are pulled straight from docstrings on every build.
