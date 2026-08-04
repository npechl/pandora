# Installation

Not yet on PyPI — install from a checkout of the repo:

```sh
git clone https://github.com/npechl/pandora.git
cd pandora

# base install
pip install -e .

# install full toolkit
pip install -e ".[full]"

# install multiple components
pip install -e ".[similarity,export]"

# install development
pip install -e ".[dev]"
```

The base install includes `pydantic`, `pyyaml`, and `gemmi` — enough to
**parse mmCIF files you already have, canonicalise them, collect
metadata, and run annotations**. Everything else below is an add-on
for one specific capability.

## Components

`pandora/` is one package per pipeline stage (see the [Overview](index.md)
diagram), but most stages need nothing beyond the base install — only
`ingestion`, `export`, and `similarity` pull in extra dependencies.
Rows below with no `pip install` command aren't real `pyproject.toml`
extras; running `pip install -e ".[parsing]"` (for example) would just
error with "extra not found".

| Install | Description | Adds |
|---------|-------------|------|
| `pip install -e .` | **parsing** — turn a raw mmCIF file into Pandora's typed `Structure`, via `mmcif_to_structure()` | nothing — `gemmi` is a base dependency |
| `pip install -e ".[ingestion]"` | **ingestion** — fetch mmCIF files from PDBe/RCSB myself | `httpx`, for `fetch_mmcif()` / `fetch_list_mmcif()` |
| `pip install -e .` | **canonicalisation** — apply a policy to normalize chain IDs, altlocs, entities, ligands, etc. via `canonicalise_structure()` | nothing |
| `pip install -e .` | **metadata** — collect source-backed entry/quality/taxonomy/entity/ligand/UniProt-mapping records via `collect_metadata()` | nothing |
| *(base install*)* | **annotations** — derived per-entry/pairwise layers: structure counts, ligand contacts, chain interfaces, sequence identity | nothing yet — see note below |
| `pip install -e .` | **datasets** — reshape a canonical `Structure` into Chain/Residue/Interface records via `extract_*_records()`. Dataset curation (filtering, dedup, splitting) is **not implemented yet** | nothing |
| `pip install -e ".[similarity]"` + `mmseqs2`/`foldseek` on `PATH` | **similarity** — compute sequence/structure similarity | nothing pip-managed — see below |
| `pip install -e .` | **provenance** — build a per-structure `ProvenanceBundle` (ingestion/canonicalisation/metadata/annotation provenance + a checksum) via `build_provenance_bundle()` | nothing |
| `pip install -e .` | **schemas** — the typed Pydantic models (`Structure`, policies, records, ...) every other package builds on; not something you install on its own | nothing — always included |
| **not implemented yet** | **cli** — `pandora.cli.app()` currently just raises `NotImplementedError` | — |
| `pip install -e ".[export]"` | **export** — write datasets to Parquet | `pandas` + `pyarrow`, for `write_records(..., "*.parquet")`. Exporting to mmCIF/JSON needs no extra |
| `pip install -e ".[dev]"` | Run the test suite / lint | `pytest`, `pytest-cov`, `ruff` |
| `pip install -e ".[full]"` | All extras combined | `ingestion` + `annotations` + `export` + `dev` + `docs` |

If you're not fetching files over the network, computing similarity,
or exporting to Parquet, the base install is all you need — parsing,
canonicalisation, metadata, annotations, dataset records, and
provenance bundles all work with no extras.

!!! note "`annotations` extra"
    `pip install -e ".[annotations]"` adds `freesasa`, but no
    currently implemented annotation function uses it yet — it's
    reserved for a solvent-accessible-surface-area annotation that
    isn't built. You don't need it today.

## External binaries

### Similarity

`pandora.similarity.sequence` (MMseqs2) and `pandora.similarity.structure`
(Foldseek) shell out to command-line tools that aren't distributed on
PyPI. Install them separately and make sure they're on `PATH`, e.g.
via conda:

```bash
conda install -c bioconda mmseqs2 foldseek
```

Nothing in the test suite exercises these two modules, so a missing
binary only breaks the specific call that needs it — the rest of
Pandora is unaffected.
