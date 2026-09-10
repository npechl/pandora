# CLI

The `pandora` command exposes the pipeline as subcommands: `fetch`, `ingest`, `canonicalise`, `curate`, `dedup`, `similarity`, `cluster`, `partition`, `annotate`, `manifest`, `reproduce`, `export`. Run `pandora <subcommand> -h` for its arguments, or `pandora -h` for the full list.

`fetch` downloads from PDBe/PDB; `ingest` instead records the same shape of `ingestion_provenance.json` for mmCIF files you already have on disk (e.g. a local bulk PDB mirror) — it doesn't download or copy anything, just describes what's already there.

Structure-transforming stages (`canonicalise`, `curate`, `dedup`) read and write directories of mmCIF files, one per entry_id. Provenance-producing stages write one JSON file per stage into their `--output`/`--output-dir` (e.g. `canonicalisation_provenance.json`), which downstream stages — `manifest` in articular — read back in.
