# CLI

The `pandora` command exposes the pipeline as subcommands: `fetch`,
`canonicalise`, `curate`, `dedup`, `similarity`, `cluster`, `partition`,
`annotate`, `manifest`, `reproduce`, `export`. Run `pandora <subcommand> -h`
for its arguments, or `pandora -h` for the full list.

Structure-transforming stages (`canonicalise`, `curate`, `dedup`) read and
write directories of mmCIF files, one per entry_id. Provenance-producing
stages write one JSON file per stage into their `--output`/`--output-dir`
(e.g. `canonicalisation_provenance.json`), which downstream stages —
`manifest` in particular — read back in.
