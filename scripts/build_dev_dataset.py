"""Build a small, stratified mmCIF dataset for pandora development.

Queries the RCSB Search API (https://search.rcsb.org) for entries in a
handful of buckets chosen to exercise different pandora code paths (altlocs,
multimers, ligands/ions, NMR multi-model, cryo-EM assemblies, ...), then
fetches the mmCIF files via pandora.ingestion.fetch_list_mmcif.

Downloaded files are NOT committed to git (see .gitignore) -- only this
script and the resulting manifest.json are. Re-run this script to
regenerate the dataset from scratch, or a fresh sample of it.

Usage:
    python scripts/build_dev_dataset.py
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import json
import random
import urllib.request

from pandora.ingestion import fetch_list_mmcif
from pandora.schemas.ingestion import FetchOptions

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "dev-data" / "mmcif"
MANIFEST_PATH = OUTPUT_DIR.parent / "manifest.json"

RANDOM_SEED = 20260714
ROWS_PER_BUCKET = 200  # fetched from the API before sub-sampling

# atom-count ceilings keep files small for fast iterative dev, per bucket
# criteria chosen to hit distinct pandora canonicalisation/parsing paths.
BUCKETS: dict[str, tuple[list[dict], int]] = {
    "small_apo_monomer": (
        [
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
                    "attribute": "rcsb_entry_info.deposited_polymer_entity_instance_count",
                    "operator": "equals",
                    "value": 1,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                    "operator": "equals",
                    "value": 0,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_atom_count",
                    "operator": "less_or_equal",
                    "value": 2500,
                },
            },
        ],
        16,
    ),
    "small_ligand_monomer": (
        [
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
                    "attribute": "rcsb_entry_info.deposited_polymer_entity_instance_count",
                    "operator": "equals",
                    "value": 1,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                    "operator": "greater",
                    "value": 0,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_atom_count",
                    "operator": "less_or_equal",
                    "value": 3000,
                },
            },
        ],
        18,
    ),
    "homo_multimer": (
        [
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
                    "attribute": "rcsb_entry_info.polymer_entity_count",
                    "operator": "equals",
                    "value": 1,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_polymer_entity_instance_count",
                    "operator": "greater",
                    "value": 2,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_atom_count",
                    "operator": "less_or_equal",
                    "value": 6000,
                },
            },
        ],
        16,
    ),
    "hetero_multimer": (
        [
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
                    "attribute": "rcsb_entry_info.polymer_entity_count",
                    "operator": "greater",
                    "value": 1,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_atom_count",
                    "operator": "less_or_equal",
                    "value": 6000,
                },
            },
        ],
        16,
    ),
    "nmr_multimodel": (
        [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "exptl.method",
                    "operator": "exact_match",
                    "value": "SOLUTION NMR",
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_model_count",
                    "operator": "greater",
                    "value": 5,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_atom_count",
                    "operator": "less_or_equal",
                    "value": 2000,
                },
            },
        ],
        16,
    ),
    "cryoem_small": (
        [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "exptl.method",
                    "operator": "exact_match",
                    "value": "ELECTRON MICROSCOPY",
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_atom_count",
                    "operator": "less_or_equal",
                    "value": 8000,
                },
            },
        ],
        8,
    ),
    "high_res_altloc_prone": (
        [
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
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.deposited_atom_count",
                    "operator": "less_or_equal",
                    "value": 3000,
                },
            },
        ],
        10,
    ),
}


def _search(nodes: list[dict], rows: int) -> list[str]:
    body = {
        "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows}},
    }
    req = urllib.request.Request(
        SEARCH_URL,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        raw = resp.read()
    if not raw:
        return []
    return [r["identifier"] for r in json.loads(raw).get("result_set", [])]


def select_entries() -> dict[str, list[str]]:
    rng = random.Random(RANDOM_SEED)
    chosen: dict[str, list[str]] = {}
    seen: set[str] = set()
    for bucket, (nodes, count) in BUCKETS.items():
        candidates = [
            c for c in _search(nodes, ROWS_PER_BUCKET) if c not in seen
        ]
        rng.shuffle(candidates)
        picked = candidates[:count]
        seen.update(picked)
        chosen[bucket] = picked
    return chosen


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    buckets = select_entries()
    all_ids = sorted({eid for ids in buckets.values() for eid in ids})
    print(f"Selected {len(all_ids)} entries across {len(buckets)} buckets.")

    provenance = fetch_list_mmcif(
        entry_ids=all_ids,
        provider="pdbe",
        source_uri=None,
        output_dir=OUTPUT_DIR,
        fetch_options=FetchOptions(allow_partial=True),
    )
    print(f"Fetched {len(provenance)}/{len(all_ids)} entries.")

    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "seed": RANDOM_SEED,
        "provider": "pdbe",
        "buckets": buckets,
        "requested_count": len(all_ids),
        "fetched_count": len(provenance),
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))
    print(f"Manifest written to {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
