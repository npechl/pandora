# C01 — mmCIF Ingestion

Fetches raw mmCIF files from PDBe or RCSB, parses them with **gemmi** into a
typed `ParsedStructure` (atoms → residues → chains → entities → assemblies),
and validates the result. Parsed files are cached under `~/.pandora/cache/mmcif/`.

## Prerequisites

```bash
pip install -e ".[ingestion]"
```

## fetch_single_entry.py

Fetch and parse one PDB entry. Accepts an optional entry ID argument (default: `1cbs`).

```bash
python fetch_single_entry.py
python fetch_single_entry.py 4hhb
```

Expected output:

```
Entry:      1cbs
Status:     success
From cache: False
Chains:     1
Residues:   238
Atoms:      1213
Entities:   3
Assemblies: 1
Ligands:    1

  chain 'A'   type=polymer       entity='1'

  entity '1'  type=polymer       seq_len=137
  entity '2'  type=non-polymer   seq_len=0
  entity '3'  type=water         seq_len=0
```

A second run of the same entry will show `From cache: True` — no network call is made.

## ingest_batch.py

Ingest a list of PDB entries sequentially. Accepts an optional list of entry IDs
(default: `1cbs 4hhb 1tup`).

```bash
python ingest_batch.py
python ingest_batch.py 1cbs 4hhb 1tup 2hho
```

Expected output:

```
Total:   3
Success: 3
Warning: 0
Failed:  0

  1cbs    success   chains=1  atoms=  1213  cached=True
  4hhb    success   chains=4  atoms=  4779  cached=False
  1tup    success   chains=5  atoms=  5828  cached=False
```

## Using a local file

```python
from pandora.c01_ingestion import ingest_mmcif
from pandora.schemas.c01_ingestion import MmCIFIngestionInput

result = ingest_mmcif(MmCIFIngestionInput(
    entry_id="1cbs",
    provider="local",
    source_uri="/path/to/1cbs.cif",
))
```

## Passing raw content directly

```python
raw = open("1cbs.cif").read()

result = ingest_mmcif(MmCIFIngestionInput(
    entry_id="1cbs",
    provider="local",
    raw_content=raw,   # skips fetch entirely
))
```
