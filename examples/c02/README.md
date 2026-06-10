# C02 — Canonicalization

Converts a `ParsedStructure` (from C01) into a `CanonicalStructure` by applying a
`CanonicalizationPolicy`. Planned normalization steps: chain IDs, residue numbering,
assembly selection, altloc resolution, entity normalization, and ligand filtering.

> **Status:** API wired; normalization sub-steps are stubs returning empty mappings.

## Prerequisites

```bash
pip install -e ".[ingestion]"
```

## canonicalize.py

Runs C01 ingestion then C02 canonicalization for one entry.
Accepts an optional entry ID argument (default: `1cbs`).

```bash
python canonicalize.py
python canonicalize.py 4hhb
```

Expected output:

```
[C01] 1cbs ingestion status: success
[C02] canonicalization status: success
      chains:    1
      residues:  238
      atoms:     1213
      policy:    canon-v1 v1.0.0

NOTE: C02 normalization steps are stubs — mappings are empty until implemented.
```

## Custom policy

```python
from pandora.schemas.c02_canonicalization import CanonicalizationPolicy

policy = CanonicalizationPolicy(
    policy_id="my-policy",
    policy_name="My Canonicalization Policy",
    policy_version="2.0.0",
)
```
