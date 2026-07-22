# Pandora

Pandora is a Python library for turning raw PDB/PDBe data into typed, policy-driven, ML-ready protein structure datasets.

Each component, ingestion, parsing, canonicalisation, metadata, annotation, includes plain functions you can call on its own, or chain into a pipeline. Nothing is hidden behind a framework object: you pass a `Structure` in, you get a `Structure` (or a typed record) out.

> **Status:** under active development. Ingestion, parsing, canonicalisation,
> metadata, and annotations are implemented today; dataset curation,
> similarity/splitting, and provenance manifests are still on the roadmap.

## Install

Pandora is not yet published to PyPI or any other software distibution repository. Please install from source:

```bash
git clone https://github.com/npechl/pandora.git
cd pandora
pip install -e ".[full]"
```

Optional extras (see `pyproject.toml`): `ingestion`, `annotations`, `cli`,
`dev`, `docs`.

## Quick start

```python
from pathlib import Path

from pandora.ingestion import fetch_mmcif
from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.metadata import collect_metadata
from pandora.schemas.canonicalisation import canonicalisationPolicy

entry_id = "1cbs"
mmcif_dir = Path("./mmcif")

provenance = fetch_mmcif(
    entry_id=entry_id,
    provider="pdbe",
    source_uri=None,
    output_dir=mmcif_dir,
)

structure, diagnostics, status = mmcif_to_structure(
    str(mmcif_dir / f"{entry_id}.cif")
)

policy = canonicalisationPolicy(
    policy_id="quickstart",
    policy_name="Quickstart",
    policy_version="1.0.0",
)
canonical, mappings, canon_provenance = canonicalise_structure(structure, policy)

metadata = collect_metadata(canonical)
```

## Documentation

*TODO*

## Contributing

*TODO*

## License

[CC0 1.0 Universal](LICENSE).
