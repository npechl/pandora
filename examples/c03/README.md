# C03 — Metadata & Annotation

Retrieves external metadata (PDBe, SIFTS, UniProt, taxonomy) and attaches it to a
`CanonicalStructureResult` (from C02). Optionally executes annotation plugins that
produce derived `AnnotationLayer` records.

> **Status:** API wired; metadata retrieval and plugin execution are stubs returning empty records.

## Prerequisites

```bash
pip install -e ".[ingestion]"
```

## attach_metadata.py

Runs C01 → C02 → C03 for one entry.
Accepts an optional entry ID argument (default: `1cbs`).

```bash
python attach_metadata.py
python attach_metadata.py 4hhb
```

Expected output:

```
[C01] ingestion:        success
[C02] canonicalization: success
[C03] metadata attached
      entry_id:  1cbs
      retrieved: 2026-06-10T...

NOTE: C03 API calls (PDBe, SIFTS, UniProt) are stubs — metadata fields are empty until implemented.
```

## Batch mode

```python
from pandora.c03_metadata import attach_metadata_many
from pandora.schemas.c03_metadata import MetadataIntegrationPolicy

policy = MetadataIntegrationPolicy(
    policy_id="meta-v1",
    policy_name="Default Metadata",
    policy_version="1.0.0",
)

batch_result = attach_metadata_many(canonical_results, policy)
print(batch_result.summary)
```
