# Overview

Pandora is developed in the open, and contributions are welcome —
bug reports, feature requests, code, examples, or docs.

New to the codebase? [Architecture](architecture.md)
is the fastest way to see how the schemas and stages fit together
before you go looking for where a change belongs.

## Where to start

- **Environment setup, branching, and opening a PR** — the mechanics
  live in [`CONTRIBUTING.md`](https://github.com/npechl/pandora/blob/main/CONTRIBUTING.md)
  in the repo root: cloning, `uv sync --all-extras`, running the test
  suite and linters, and the commit/push/PR flow.
- **What to do for each kind of contribution** —
  [How to contribute](contribute.md) covers reporting a bug, requesting
  a feature, fixing a bug, adding a function to an existing component,
  adding a new annotation, writing an example or recipe, building a
  dataset, and updating the docs themselves.

Every contribution ends the same way — lint, test, then a PR:

```sh
uv run ruff format .
uv run ruff check .
uv run pytest
```
