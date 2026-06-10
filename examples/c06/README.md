# C06 — Provenance & Reproducibility

Assembles full pipeline provenance, generates a manifest with dataset summary and
checksums, and seals a `PandoraArtifact` that bundles the leakage-safe dataset with
its reproducibility record.

> **Status:** API wired; SHA-256 checksum computation, lineage graph export, and full
> provenance traversal are stubs returning empty/zeroed values.

## Prerequisites

```bash
pip install -e ".[ingestion]"
```

## finalize_artifact.py

Runs the full C01–C06 chain and prints the sealed artifact summary.
Accepts an optional list of entry IDs (default: `1cbs 4hhb 1tup 2hho 1a4w`).

```bash
python finalize_artifact.py
python finalize_artifact.py 1cbs 4hhb 1tup 2hho 1a4w 2ci2
```

Expected output:

```
  prepared 1cbs
  prepared 4hhb
  ...

[C06] Artifact:  pandora-artifact-001
      name:      Example Artifact
      manifest:  <uuid4>
      version:   0.1.0
      generated: 2026-06-10T...
      train:     4
      val:       0
      test:      1

NOTE: SHA-256 checksums, lineage export, and full provenance traversal are stubs.
```

## Serializing the artifact to JSON

```python
import json

artifact_json = artifact.model_dump_json(indent=2)
print(artifact_json)

# or write to disk
with open("pandora-artifact-001.json", "w") as f:
    f.write(artifact_json)
```
