# C05 — Similarity & Leakage-Safe Splitting

Computes pairwise similarity (MMseqs2 for sequence, Foldseek for structure),
builds a similarity network, clusters it, and assigns clusters to
train / validation / test partitions so that no similar items leak across splits.

> **Status:** API wired; MMseqs2/Foldseek subprocess calls are stubs returning no
> relationships. Clustering produces one singleton cluster per item. Partitioning
> is sequential (not cluster-aware yet).

## Prerequisites

```bash
pip install -e ".[ingestion]"
# For real similarity: MMseqs2 and Foldseek binaries must be in PATH
```

## create_splits.py

Runs the full C01–C05 chain for a list of entries.
Accepts an optional list of entry IDs (default: `1cbs 4hhb 1tup 2hho 1a4w`).

```bash
python create_splits.py
python create_splits.py 1cbs 4hhb 1tup 2hho 1a4w 2ci2
```

Expected output:

```
  prepared 1cbs
  prepared 4hhb
  prepared 1tup
  prepared 2hho
  prepared 1a4w

[C05] Split: split-example:split
      train:      4  (80%)
      validation: 0  (0%)
      test:       1  (20%)
      clusters:   5

NOTE: MMseqs2/Foldseek similarity and cluster-aware assignment are stubs.
```

## Controlling split fractions

```python
from pandora.schemas.c05_splitting import LeakagePolicy, PartitionRules

policy = LeakagePolicy(
    policy_id="split-v1",
    policy_name="70/15/15 Split",
    policy_version="1.0.0",
    partition_rules=PartitionRules(
        train_fraction=0.70,
        validation_fraction=0.15,
        test_fraction=0.15,
    ),
)
```
