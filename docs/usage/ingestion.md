# Ingestion

`pandora.ingestion` fetches raw mmCIF files from PDBe/PDB (with on-disk
caching) and loads canonicalisation policy YAML files. See
[Functions](../reference/functions.md#pandora.ingestion) for full
signatures.

## Fetch one entry

=== "`library`"

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

    `provenance.revision_date` carries the entry's own most recent revision
    date (from RCSB's mmCIF file directly, or PDBe's summary API — PDBe's
    mmCIF download doesn't include revision history itself), so a later
    re-fetch of the same entry can detect drift; it's `None` when a
    `source_uri` mirror doesn't expose one.

=== "`cli`"

    Download raw mmCIF files (with on-disk caching):

    ```bash
    pandora fetch 104m 112m 118l 138l 1ayi --output-dir raw/
    # fetched 5/5 entries -> raw/
    ```

    `raw/ingestion_provenance.json` is written alongside the files — keep
    it, `manifest` needs it later if you want `reproduce` to work.

## Find entry ids to fetch

No CLI subcommand for this yet (`search_rcsb`/`search_pdbe` output feeds
straight into `fetch`'s `entry_ids` argument, or `fetch_list_mmcif` in
code) — if you don't already have a list of ids, `search_rcsb()` and
`search_pdbe()` query each provider's own search API and return matching
entry ids. Each takes that provider's *native* query syntax verbatim —
Pandora doesn't translate between the two (they don't even share field
names for the same concept), so pick the one whose query language you
want to write.

```python
from pandora.ingestion import search_rcsb, search_pdbe

# RCSB: a query tree (https://search.rcsb.org/#search-api)
rcsb_ids = search_rcsb(
    {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "exptl.method",
                    "operator": "exact_match",
                    "value": "X-RAY DIFFRACTION",
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": 1.2,
                },
            },
        ],
    },
    rows=50,
)

# PDBe: a raw Solr/Lucene query string (https://www.ebi.ac.uk/pdbe/api/doc/search.html)
pdbe_ids = search_pdbe(
    'experimental_method:"X-ray diffraction" AND resolution:[0 TO 1.2]',
    rows=50,
)
```

Both are a single page (`rows`/`start`) of the full result set, not
auto-paginated — loop and increment `start` yourself for more than `rows`
matches. PDBe's Solr endpoint doesn't reject a malformed query (an
unparseable clause is silently dropped rather than erroring), so
sanity-check the returned count against what you expect.

## Ingest already-downloaded files

Pandora has no bulk-download step of its own — for a full PDB snapshot,
use an existing mirror (e.g. `rsync.rcsb.org`) rather than fetching entries
one at a time. `ingest_local_mmcif()` records provenance for files you
already have on disk, without fetching or copying anything.

=== "`library`"

    ```python
    from pandora.ingestion import ingest_local_mmcif
    from pathlib import Path

    provenance = ingest_local_mmcif(
        Path("./mirror/pdb_snapshot_2024-01-29/1crn.cif")
    )
    print(provenance.provider)
    # local
    ```

    Pass `source_uri` to record a more meaningful label than the raw path,
    e.g. `"pdb_snapshot_2024-01-29"`.

=== "`cli`"

    ```bash
    pandora ingest --input-dir mirror/pdb_snapshot_2024-01-29/ --output-dir raw/
    # ingested 5/5 entries -> raw/
    ```

    Same `ingestion_provenance.json` shape as `fetch`, so `manifest` reads
    it the same way regardless of which one produced it. `--source-uri`
    applies one label to every file in `--input-dir`.

## Fetch a batch, tolerating failures

=== "`library`"

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

=== "`cli`"

    `--allow-partial` skips entries that fail instead of aborting the batch:

    ```bash
    pandora fetch 1crn not-a-real-id --allow-partial --output-dir raw/
    # warning: failed to fetch not-a-real-id: ...
    # fetched 1/2 entries -> raw/
    ```

    Without `--allow-partial` (the default), a single bad id raises and
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
[`canonicalise_structure()`](canonicalisation.md) expects. There's no
standalone CLI subcommand for this step — `pandora canonicalise` calls
it internally from its `--policy` flag.

```python
from pandora.ingestion import load_policy

policy = load_policy("datasets/canonicalisation.yaml")
print(policy.policy_id, policy.ligand_rules.strategy)
# overview-remap filter
```

See [Policies](../reference/policies.md) for every field a policy YAML can set.
