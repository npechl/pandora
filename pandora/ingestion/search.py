from __future__ import annotations

from pandora._util import require

_RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
_PDBE_SEARCH_URL = "https://www.ebi.ac.uk/pdbe/search/pdb/select"


def search_rcsb(
    query: dict,
    return_type: str = "entry",
    rows: int = 100,
    start: int = 0,
) -> list[str]:
    """Search RCSB's Search API for matching entry ids.

    `query` is RCSB's own query tree, passed through verbatim — Pandora
    doesn't validate or translate it. See
    https://search.rcsb.org/#search-api for the query language
    (terminal/group nodes, services, attributes, operators) and
    https://search.rcsb.org/rcsbsearch/v2/metadata/schema for the full
    list of searchable attributes. The returned ids feed directly into
    `fetch_mmcif`/`fetch_list_mmcif` with `provider="pdb"`.

    Args:
        query: An RCSB query tree, e.g. `{"type": "terminal", "service":
            "text", "parameters": {"attribute": "exptl.method",
            "operator": "exact_match", "value": "X-RAY DIFFRACTION"}}`.
        return_type: What kind of id to return — "entry" for PDB entry
            ids (the default), or one of RCSB's other result types
            ("polymer_entity", "assembly", ...).
        rows: Maximum number of ids to return in this page.
        start: Offset into the full result set, for paging through more
            than `rows` matches.

    Returns:
        Matching identifiers, in RCSB's own casing (e.g. "104M"). Empty
        list if nothing matches.

    Raises:
        RuntimeError: If the request fails, including a malformed query
            — RCSB validates the query tree server-side and returns its
            own error message on a bad attribute/operator.
    """

    body = {
        "query": query,
        "return_type": return_type,
        "request_options": {"paginate": {"start": start, "rows": rows}},
    }
    httpx = require("httpx", "ingestion")
    try:
        resp = httpx.post(_RCSB_SEARCH_URL, json=body, timeout=60.0)
        resp.raise_for_status()
    except httpx.HTTPStatusError as exc:
        raise RuntimeError(
            f"HTTP {exc.response.status_code} querying RCSB search API: "
            f"{exc.response.text}"
        ) from exc
    except httpx.RequestError as exc:
        raise RuntimeError(
            f"Network error querying RCSB search API: {exc}"
        ) from exc

    if resp.status_code == 204:
        return []
    return [hit["identifier"] for hit in resp.json().get("result_set", [])]


def search_pdbe(
    query: str,
    rows: int = 100,
    start: int = 0,
) -> list[str]:
    """Search PDBe's Solr-based search API for matching entry ids.

    `query` is a raw Solr/Lucene query string, passed through verbatim as
    the `q` parameter — Pandora doesn't validate or translate it. See
    https://www.ebi.ac.uk/pdbe/api/doc/search.html for the query syntax
    and the full list of searchable fields. The returned ids feed
    directly into `fetch_mmcif`/`fetch_list_mmcif` with
    `provider="pdbe"`.

    Unlike RCSB's search API, PDBe's Solr endpoint doesn't reject a
    malformed query — an unparseable clause is silently dropped rather
    than erroring, which can make a typo'd query match far more than
    intended. Sanity-check the returned count against what you expect.

    Args:
        query: A Solr query string, e.g. `'experimental_method:"X-ray
            diffraction" AND resolution:[0 TO 1.5]'`.
        rows: Maximum number of ids to return in this page.
        start: Offset into the full result set, for paging through more
            than `rows` matches.

    Returns:
        Matching PDB ids, lowercase (PDBe's own convention), deduplicated
        — PDBe's Solr index can list the same id more than once (e.g. one
        row per assembly/entity). Empty list if nothing matches.

    Raises:
        RuntimeError: If the request fails.
    """

    params = {
        "q": query,
        "fl": "pdb_id",
        "rows": rows,
        "start": start,
        "wt": "json",
    }
    httpx = require("httpx", "ingestion")
    try:
        resp = httpx.get(_PDBE_SEARCH_URL, params=params, timeout=60.0)
        resp.raise_for_status()
    except httpx.HTTPStatusError as exc:
        raise RuntimeError(
            f"HTTP {exc.response.status_code} querying PDBe search API: "
            f"{exc.response.text}"
        ) from exc
    except httpx.RequestError as exc:
        raise RuntimeError(
            f"Network error querying PDBe search API: {exc}"
        ) from exc

    docs = resp.json().get("response", {}).get("docs", [])
    # PDBe's Solr index can return the same pdb_id more than once (e.g.
    # once per assembly/entity row); dedupe, preserving result order.
    return list(dict.fromkeys(doc["pdb_id"] for doc in docs if "pdb_id" in doc))
