# Contributing

For environment setup, branching, and how to open a PR, see
[`CONTRIBUTING.md`](https://github.com/npechl/pandora/blob/main/CONTRIBUTING.md)
in the repo root. This page covers what to do for each kind of
contribution.

Before diving in, [Architecture](getting-started/architecture.md) shows how the data
models fit together — the fastest way to see where a new piece belongs.

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

1. Find the function via [Functions](reference/functions.md) or
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
logic in `pandora/<package>/`** (see [Architecture](getting-started/architecture.md)).

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
   [Functions](reference/functions.md) and [Schemas](reference/schemas.md)
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

There's no formal plugin-registration system yet —
`pandora/annotations/plugins.py` and `pandora/annotations/base.py` are
placeholders for one. If you want a plugin protocol that lets
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

[Recipes](recipes/recipes.md) is a currently-empty stub for narrated,
task-oriented write-ups ("how do I do X") built on top of an example
script — the first one is a welcome contribution. Add a page under
`docs/recipes/` and link it from `docs/recipes/recipes.md`.

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
