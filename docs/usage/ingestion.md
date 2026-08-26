# Ingestion

`pandora.ingestion` fetches raw mmCIF files from PDBe/PDB (with on-disk
caching) and loads canonicalisation policy YAML files. See
[Functions](../reference/functions.md#pandora.ingestion) for full
signatures.

## Fetch one entry

`fetch_mmcif()` downloads (or reuses a cached copy of) one entry and
returns provenance describing where it came from.

```python
from pandora.ingestion import fetch_mmcif
from pathlib import Path

output_dir = Path("./datasets/output/fetched")
provenance = fetch_mmcif("1crn", "pdbe", None, output_dir)

print(provenance.source_uri)
# https://www.ebi.ac.uk/pdbe/entry-files/download/1crn_updated.cif
print(provenance.from_cache)
# False — running this again reuses the cached file on disk and prints True
```

`provider` is `"pdbe"` or `"pdb"`; pass an explicit `source_uri` instead
to fetch from a mirror or local file server.

## Fetch a batch, tolerating failures

`fetch_list_mmcif()` calls `fetch_mmcif()` once per id.
`FetchOptions(allow_partial=True)` skips entries that fail instead of
aborting the whole batch — useful when fetching a list from an external
search result that might include a typo or a withdrawn entry.

```python
from pandora.ingestion import fetch_list_mmcif
from pandora.schemas.ingestion import FetchOptions

provenance = fetch_list_mmcif(
    ["1crn", "not-a-real-id"],
    "pdbe",
    output_dir,
    fetch_options=FetchOptions(allow_partial=True),
)
print(len(provenance))
# 1 — the bad id was skipped, not raised
```

Without `allow_partial=True` (the default), a single bad id raises and
stops the batch.

### Cache control

`FetchOptions` also controls staleness: `max_age_seconds` plus
`stale_behavior` (`"use_stale"` / `"warn"` / `"fail"`) decide what
happens when a cached file is older than that. `use_cache=False` skips
the cache entirely and always re-downloads.

## Load a canonicalisation policy

`load_policy()` reads a policy YAML file into a validated
`canonicalisationPolicy` — the same object the `pandora canonicalise`
CLI subcommand builds internally, and what
[`canonicalise_structure()`](canonicalisation.md) expects.

```python
from pandora.ingestion import load_policy

policy = load_policy("datasets/canonicalisation.yaml")
print(policy.policy_id, policy.ligand_rules.strategy)
# overview-remap filter
```

See [Policies](../policies.md) for every field a policy YAML can set.
