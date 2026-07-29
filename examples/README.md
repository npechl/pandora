# Examples

Every script below runs directly against local fixtures in
`datasets/dev/mmcif/` — no network access or external binaries needed.

- `overview.py` — minimal end-to-end walkthrough (parsing,
  canonicalisation, metadata, annotations) on two small monomers.
- `multichain_ligand_complex.py` — the fuller pipeline on a richer
  structure (1a3n, a hemoglobin tetramer with 4 heme groups):
  chain-chain interfaces, ligand contacts, dataset record extraction,
  every `pandora.export` format, and a `pandora.provenance` bundle.
- `dataset_pipeline.py` — multi-structure batch processing: an
  all-vs-all sequence-similarity network built from Pandora's own
  annotation (no MMseqs2/Foldseek required), connected-component
  clustering, and a leakage-safe train/val/test split.