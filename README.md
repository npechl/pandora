# Pandora

Pandora is a Python toolkit for building ML-ready datasets from structural
biology data. It gives you typed data models and plain functions to parse
PDB/PDBe mmCIF files, normalise them with explicit policies, attach metadata
and annotations, curate and split datasets without leakage, and record
enough provenance to rebuild the result.

> [!WARNING]
> Pandora is under active development. The API can change between
> versions, and the set of components is still growing.

**Documentation:** <https://npechl.github.io/pandora/>

## Why Pandora

Structural data arrives as loosely typed, inconsistently annotated mmCIF.
Every ML project then writes its own ad hoc cleanup, and small choices
(which altloc, which assembly, how to split similar chains) silently change
the dataset. Pandora makes those choices explicit and reusable:

- **Typed models.** Everything is a pydantic model, from `Structure` to
  dataset records and manifests.
- **Plain functions.** Pass a `Structure` or record in, get a new one out.
  No framework object, no global state, no mutation.
- **Policies, not hard-coded rules.** Canonicalisation and curation are
  driven by policy models you can write in YAML and share.
- **Provenance.** Every step can record its source, policy, and tool
  versions, so a dataset can be rebuilt from its manifest.

Use the components you need, in the order you need them.

## Components

| Package | What it does |
|---|---|
| `pandora.ingestion` | Fetch mmCIF files from PDBe/RCSB with on-disk caching; search RCSB/PDBe; load policies |
| `pandora.parsing` | Parse mmCIF into a typed `Structure` (via gemmi), keeping unknown categories |
| `pandora.canonicalisation` | Policy-driven normalisation: chain IDs, residue numbering, assemblies, altlocs, missing data, entities, ligands, validation |
| `pandora.metadata` | Entry, quality, taxonomy, entity, ligand, and UniProt-mapping records |
| `pandora.annotations` | Derived layers: structure counts, ligand contacts, chain interfaces, pairwise sequence identity |
| `pandora.datasets` | Chain/residue/interface records; policy-driven curation and deduplication |
| `pandora.similarity` | MMseqs2/Foldseek similarity, clustering, and leakage-safe train/val/test splits |
| `pandora.provenance` | Per-structure provenance bundles, dataset manifests, and `reproduce_dataset()` |
| `pandora.export` | Write structures to mmCIF and records to JSON, JSONL, or Parquet |
| `pandora.schemas` | The pydantic models every other package builds on |

Most components are also available as `pandora` CLI subcommands (see
[CLI](#command-line)).

## Install

Pandora is not on PyPI yet. Install it from source (Python 3.11 or newer):

```bash
git clone https://github.com/npechl/pandora.git
cd pandora
pip install -e .            # base install
pip install -e ".[full]"    # everything
```

The base install (`pydantic`, `pyyaml`, `gemmi`) covers parsing,
canonicalisation, metadata, annotations, dataset records, curation, and
provenance. Optional extras:

| Extra | Adds | Needed for |
|---|---|---|
| `ingestion` | `httpx` | Fetching and searching PDBe/RCSB |
| `export` | `pandas`, `pyarrow` | Writing Parquet |
| `similarity` | nothing | Marker only: install the `mmseqs` and `foldseek` binaries yourself and put them on `PATH` |
| `dev` | `pytest`, `pytest-cov`, `ruff` | Tests and linting |
| `docs` | `mkdocstrings`, `erdantic` | Building the docs and diagrams |
| `full` | all of the above | |

### With uv

[uv](https://docs.astral.sh/uv/) can install straight from GitHub, no
clone needed:

```bash
# add Pandora as a dependency of your own uv project
uv add "pandora[ingestion,export] @ git+https://github.com/npechl/pandora"

# or install only the `pandora` command-line tool
uv tool install "pandora[ingestion] @ git+https://github.com/npechl/pandora"
```

To work on Pandora itself, clone the repo and install from the lockfile:

```bash
uv sync --all-extras
```

## Example

This runs offline against a structure bundled in the repo (sperm whale
myoglobin, `104m`):

```python
from pandora.annotations import annotate_ligand_contacts
from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import extract_chain_records
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import canonicalisationPolicy

structure, diagnostics, status = mmcif_to_structure(
    "datasets/dev/mmcif/104m.cif"
)

policy = canonicalisationPolicy(
    policy_id="default", policy_name="Default", policy_version="1.0.0"
)
canonical, mappings, provenance = canonicalise_structure(structure, policy)

metadata = collect_metadata(canonical)
print(metadata.quality.resolution, metadata.taxonomies[0].ncbi_taxon_id)
# 1.71 9755

contacts = annotate_ligand_contacts(canonical)
print(
    [
        (l["ligand_comp_id"], l["contact_count"])
        for l in contacts.data["ligands"]
    ]
)
# [('SO4', 4), ('HEM', 15), ('NBN', 4)]

chains = extract_chain_records(canonical)
print([(c.chain_id, c.residue_count) for c in chains])
# [('A', 153)]
```

More complete, runnable scripts are in [`examples/`](examples/), including
multi-structure datasets with clustering and leakage-safe splits.

## Command line

The `pandora` command exposes the components as subcommands that read and
write directories and JSON files:

```text
fetch  ingest  canonicalise  curate  dedup  similarity  cluster
partition  annotate  manifest  reproduce  export
```

Run `pandora <subcommand> -h` for arguments. See the
[CLI guide](https://npechl.github.io/pandora/usage/cli/) for worked
examples.

## Documentation

The [documentation site](https://npechl.github.io/pandora/) has a
quickstart, a guide per component, the full policy reference, generated
function and schema references, and end-to-end recipes.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for setup and the PR workflow. In
short:

```bash
uv sync --all-extras
uv run pytest
uv run ruff format .
uv run ruff check .
```

## License

[CC0 1.0 Universal](LICENSE)
