# Contributing to `pandora` development

## Prerequisites

1. Install [`git`](https://git-scm.com/install/).
2. Install [`uv`](https://docs.astral.sh/uv/) (installs and manages Python itself, no separate Python install needed). 
3. Install the [`gh` CLI](https://cli.github.com/) (optional, only needed for `gh pr create` in the last step).

## Get the code

Fork and clone the repo.

```sh
git clone https://github.com/<your-username>/pandora.git
cd pandora
```

## Install the project (locked, reproducible env)

```sh
uv sync --all-extras
```

This reads `uv.lock` + `pyproject.toml`, creates `.venv/`, and installs every extra (ingestion, similarity, annotations, cli, dev, docs) so nothing is missing later.

The similarity components need the `mmseqs2`/`foldseek` binaries on PATH separately if you touch that code. You can install these by following instructions here: [MMseqs2](https://github.com/soedinglab/mmseqs2) / [Foldseek](https://github.com/steineggerlab/foldseek).

### Verify the install works

```sh
uv run pytest
uv run ruff check .
uv run ruff format --check .
```

## Make a branch

```sh
git checkout -b your-feature-name
```

## Do the work, then re-check before committing

Make any changes in the code. Recheck after your changes and fix any potential issues.

```sh
uv run ruff format .
uv run ruff check .
uv run pytest
```

## Commit and push

```sh
git add <your-changes>
git commit -m "short description of the changes"
git push -u origin your-feature-name
```

## Open a PR

```sh
gh pr create --fill
```

(or push and open the PR from GitHub's UI). CI (lint, test, `pip-audit` security scan) runs automatically on the PR.

