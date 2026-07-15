# Examples

## `examples/overview.py`

Fetches PDB entry `1cbs` from PDBe, parses it, and runs it through three
different canonicalisation policies (preserve-as-reported, remap +
renumber, and author-centric identifiers). It's the same walkthrough
used in [Policies](policies.md#example-three-policies-one-structure).

```bash
python examples/overview.py
```

Requires network access on first run (the fetched mmCIF file is cached
under `examples/mmcif/` afterwards, so subsequent runs are offline).

## Dev dataset

`scripts/build_dev_dataset.py` builds the small, checked-in mmCIF sample
set under `datasets/dev/` used for local testing without hitting the
network. See that script for the entry list and how to regenerate it.
