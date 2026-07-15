# Pandora

Pandora is a Python library for turning raw PDB/PDBe data into typed, policy-driven, ML-ready protein structure datasets.

Each component, ingestion, parsing, canonicalisation, metadata, annotation, is a plain function you can call on its own, or chain into a pipeline. Nothing is hidden behind a framework object: you pass a `Structure` in, you get a `Structure` (or a typed record) out.

!!! warning
    Pandora is under active development. Ingestion, parsing, canonicalisation, metadata, and annotations are implemented today; dataset curation, similarity/splitting, and provenance manifests are still on the roadmap.

## Install

```bash
# The following commands are not working. 

# pip install pandora[ingestion]      # fetch_mmcif() / fetch_list_mmcif()
# pip install pandora[annotations]    # freesasa-backed annotations
# pip install pandora[full]           # everything, including the dev/test extras
```

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

# 1. Ingestion — fetch the raw mmCIF file, with on-disk caching.
provenance = fetch_mmcif(
    entry_id=entry_id,
    provider="pdbe",
    source_uri=None,
    output_dir=mmcif_dir,
)

# 2. Parsing — raw mmCIF -> Pandora's typed Structure.
structure, diagnostics, status = mmcif_to_structure(
    str(mmcif_dir / f"{entry_id}.cif")
)

# 3. Canonicalisation — apply a policy (chain IDs, altlocs, ligands, ...).
policy = canonicalisationPolicy(
    policy_id="quickstart",
    policy_name="Quickstart",
    policy_version="1.0.0",
)
canonical, mappings, canon_provenance = canonicalise_structure(structure, policy)

# 4. Metadata — source-backed entry/entity/quality/taxonomy records.
metadata = collect_metadata(canonical)
```

## Policies
Use [Policies](policies.md) for how to configure canonicalisation

## Recipes
See [Recipes](recipes.md) for more end-to-end configurations of how to build datasets. 

## Reference
The full function and schema reference is under [Reference](reference.md).

## Comment