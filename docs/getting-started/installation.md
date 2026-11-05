# Installation

Not yet on PyPI, so install from a checkout of the repo. You need
Python 3.11 or newer.

## With conda

`environment.yml` creates a `pandora` env with Python and Pandora (plus
the `ingestion` and `export` extras) in one step:

```sh
git clone https://github.com/npechl/pandora.git
cd pandora
conda env create -f environment.yml
conda activate pandora
pandora --help
```

To add MMseqs2 and Foldseek for similarity, uncomment them in
`environment.yml` before creating the env, or see
[External binaries](#external-binaries).

## With pip

Inside an activated virtual environment or conda env. Use
`python -m pip`, not `pip`, so the package goes into the active
environment's Python:

```sh
git clone https://github.com/npechl/pandora.git
cd pandora

# base install
python -m pip install -e .

# install full toolkit
pip install -e ".[full]"

# install multiple components
pip install -e ".[similarity,export]"

# install development
pip install -e ".[dev]"
```

The base install includes `pydantic`, `pyyaml`, and `gemmi` — enough to **parse mmCIF files you already have, canonicalise them, collect metadata, and run annotations**. Everything else below is an add-on for one specific capability.

## Check the install

```sh
python -m pandora --help    # always uses the active environment
which pandora               # should point inside your environment
```

If `pandora --help` fails with `ModuleNotFoundError: No module named
'pandora'` (or `'pydantic'`) but `python -m pandora --help` works, your
shell is running a `pandora` script left over from another Python
install. Run `rehash` (zsh) or `hash -r` (bash), or delete the stale
script that `which pandora` shows.

## Update

Both installs above are editable (`-e`), so the environment runs the code in your checkout. To move to a newer Pandora, update the checkout and refresh the dependencies:

```sh
cd pandora
git pull

# conda env: also picks up dependency changes in environment.yml
conda env update -f environment.yml --prune

# pip: re-run the install with the extras you used before
python -m pip install -e ".[full]"
```

`git pull` alone updates the code. The second step is needed when a release adds or changes a dependency or a command. Check the installed version with:

```sh
python -m pip show pandora
```

## Components

`pandora/` is one package per component (see the [Overview](../index.md) diagram), but most stages need nothing beyond the base install — only `ingestion`, `export`, and `similarity` pull in extra dependencies. Rows below with no `pip install` command aren't real `pyproject.toml` extras; running `pip install -e ".[parsing]"` (for example) would just error with "extra not found".

| Install | Description | Adds |
|---------|-------------|------|
| `pip install -e .` | **parsing** — turn a raw mmCIF file into Pandora's typed `Structure`, via `mmcif_to_structure()` | nothing — `gemmi` is a base dependency |
| `pip install -e ".[ingestion]"` | **ingestion** — fetch mmCIF files from PDBe/RCSB myself | `httpx`, for `fetch_mmcif()` / `fetch_list_mmcif()` |
| `pip install -e .` | **canonicalisation** — apply a policy to normalize chain IDs, altlocs, entities, ligands, etc. via `canonicalise_structure()` | nothing |
| `pip install -e .` | **metadata** — collect source-backed entry/quality/taxonomy/entity/ligand/UniProt-mapping records via `collect_metadata()` | nothing |
| `pip install -e .` | **annotations** — derived per-entry/pairwise layers: structure counts, ligand contacts, chain interfaces, sequence identity | nothing |
| `pip install -e .` | **datasets** — reshape a canonical `Structure` into Chain/Residue/Interface records via `extract_*_records()`, and curate/deduplicate a dataset via `curate_structure()` / `deduplicate_structures()` | nothing |
| `pip install -e ".[similarity]"` + `mmseqs2`/`foldseek` on `PATH` | **similarity** — compute sequence/structure similarity | nothing pip-managed — see below |
| `pip install -e .` | **provenance** — build a per-structure `ProvenanceBundle` (ingestion/canonicalisation/metadata/annotation provenance) via `build_provenance_bundle()` | nothing |
| `pip install -e .` | **schemas** — the typed Pydantic models (`Structure`, policies, records, ...) every other package builds on; not something you install on its own | nothing — always included |
| `pip install -e .` | **cli** — the `pandora` command, exposing each component as a subcommand (`fetch`, `ingest`, `canonicalise`, `curate`, `dedup`, `similarity`, `cluster`, `partition`, `annotate`, `manifest`, `reproduce`, `export`) | nothing to install the CLI itself — `fetch`/`reproduce` need the `ingestion` extra, `similarity` needs `mmseqs2`/`foldseek` on `PATH` |
| `pip install -e ".[export]"` | **export** — write datasets to Parquet | `pandas` + `pyarrow`, for `write_records(..., "*.parquet")`. Exporting to mmCIF/JSON needs no extra |
| `pip install -e ".[dev]"` | Run the test suite / lint | `pytest`, `pytest-cov`, `ruff` |
| `pip install -e ".[full]"` | All extras combined | `ingestion` + `export` + `dev` + `docs` |

If you're not fetching files over the network, computing similarity, or exporting to Parquet, the base install is all you need — parsing, canonicalisation, metadata, annotations, dataset records, and provenance bundles all work with no extras.

## External binaries

### Similarity

`pandora.similarity.sequence` (MMseqs2) and `pandora.similarity.structure` (Foldseek) shell out to command-line tools that aren't distributed on PyPI. Install them separately and make sure they're on `PATH`, e.g. via `conda`:

```bash
conda install -c bioconda mmseqs2 foldseek
```

Nothing in the test suite exercises these two modules, so a missing binary only breaks the specific call that needs it — the rest of Pandora is unaffected.
