# C04 — Dataset Curation

Applies selection, quality, content, and organism filters to a list of
`AnnotatedStructureWithPlugins` records (from C03), then extracts chain-,
interface-, or residue-level sub-datasets via a `DatasetCurationPolicy`.

> **Status:** API wired; filtering, deduplication, and extraction rules are stubs
> returning the full input unmodified and empty sub-datasets.

## Prerequisites

```bash
pip install -e ".[ingestion]"
```

## build_dataset.py

Runs C01 → C02 → C03 → C04 for a list of entries, then extracts chains and residues.
Accepts an optional list of entry IDs (default: `1cbs 4hhb 1tup`).

```bash
python build_dataset.py
python build_dataset.py 1cbs 4hhb 1tup 2hho
```

Expected output:

```
  prepared 1cbs: ingestion=success
  prepared 4hhb: ingestion=success
  prepared 1tup: ingestion=success

[C04] Dataset:    example-dataset
      structures: 3
      chains:     0 (stub — implement extraction)
      residues:   0 (stub — implement extraction)

NOTE: filtering, deduplication, and extraction rules are stubs.
```

## Granularity options

```python
from pandora.c04_curation import (
    build_dataset,
    extract_chains,
    extract_interfaces,
    extract_residues,
)

dataset        = build_dataset(annotated_structures, policy, ...)
chain_dataset  = extract_chains(dataset, policy)
iface_dataset  = extract_interfaces(dataset, policy)         # requires FreeSASA
residue_dataset = extract_residues(chain_dataset, policy)   # source: Dataset | ChainDataset
```
