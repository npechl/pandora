# What is Pandora?

This walks through the implemented pipeline end to end: parse an mmCIF
file, canonicalise it, collect its metadata, and compute a couple of
annotations. It uses a file already checked into the repo, so it runs
fully offline — no network access, no PDB account.

Pandora's base install (`pip install -e .`) and a checkout of this repo
are all you need — steps 1-4 below use no optional extras. The last,
optional section (fetching a file yourself over the network) needs the
`ingestion` extra instead; see [Installation](installation.md).

## The file we'll use

`datasets/dev/mmcif/104m.cif` is one of ~100 small mmCIF files bundled
under `datasets/dev/mmcif/` for local development and docs — sperm
whale myoglobin bound to a small-molecule ligand (heme). Any entry
under that directory works the same way; swap the entry ID if you want
to follow along with a different one.

## 1. Parse

`mmcif_to_structure()` reads the raw mmCIF file and returns Pandora's
typed `Structure`, plus a diagnostics bundle and a status
(`"success"`, `"warning"`, or `"failed"`):

```python
from pathlib import Path
from pandora.parsing import mmcif_to_structure

mmcif_dir = Path("datasets/dev/mmcif")
structure, diagnostics, status = mmcif_to_structure(str(mmcif_dir / "104m.cif"))

print(status, len(structure.atoms))
# success 1450
```

## 2. Canonicalise

`canonicalise_structure()` applies a policy — chain ID remapping,
altloc resolution, ligand filtering, and more — and returns the
canonical `Structure` alongside the mappings back to the original
identifiers and a provenance record of what actually changed:

```python
from pandora.canonicalisation import canonicalise_structure
from pandora.schemas.canonicalisation import canonicalisationPolicy

policy = canonicalisationPolicy(
    policy_id="getting-started",
    policy_name="Default",
    policy_version="1.0.0",
)
canonical, mappings, canon_prov = canonicalise_structure(structure, policy)

print(canon_prov.transforms)
# ['missing_atoms:annotate', 'missing_residues:annotate', 'altloc:select_best_occupancy']
```

An unconfigured policy still resolves alternate locations (the default
`altloc_rules` strategy) and annotates missing atoms/residues rather
than silently dropping them — nothing here changed chain IDs or
residue numbers, because this entry didn't need it. See
[Policies](../reference/policies.md) for every field you can configure.

## 3. Collect metadata

`collect_metadata()` pulls the source-backed entry, quality, taxonomy,
entity, ligand, and UniProt-mapping records straight from the mmCIF
categories Pandora kept:

```python
from pandora.metadata import collect_metadata

metadata = collect_metadata(canonical)

print(metadata.entry.title)
# SPERM WHALE MYOGLOBIN N-BUTYL ISOCYANIDE AT PH 7.0
print(metadata.quality.experimental_method, metadata.quality.resolution)
# X-ray diffraction 1.71
```

## 4. Annotate

Annotations are derived layers computed from the canonical structure.
`annotate_structure_counts()` gives you a quick per-entry summary;
`annotate_ligand_contacts()` finds polymer residues near each ligand:

```python
from pandora.annotations import (
    annotate_structure_counts,
    annotate_ligand_contacts,
)

counts = annotate_structure_counts(canonical)
print(counts.data["atom_count"], counts.data["entity_type_counts"])
# 1450 {'polymer': 1, 'non-polymer': 3, 'water': 1}

contacts = annotate_ligand_contacts(canonical)
for ligand in contacts.data["ligands"]:
    print(ligand["ligand_comp_id"], ligand["contact_count"])
# SO4 4
# HEM 15
# NBN 4
```

`HEM` (the heme group) has the most polymer residues in contact,
which is what you'd expect from myoglobin's oxygen-binding pocket.

## Putting it together

```python
from pathlib import Path

from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.metadata import collect_metadata
from pandora.annotations import (
    annotate_structure_counts,
    annotate_ligand_contacts,
)
from pandora.schemas.canonicalisation import canonicalisationPolicy

mmcif_dir = Path("datasets/dev/mmcif")
structure, diagnostics, status = mmcif_to_structure(str(mmcif_dir / "104m.cif"))

policy = canonicalisationPolicy(
    policy_id="getting-started",
    policy_name="Default",
    policy_version="1.0.0",
)
canonical, mappings, canon_prov = canonicalise_structure(structure, policy)

metadata = collect_metadata(canonical)
counts = annotate_structure_counts(canonical)
contacts = annotate_ligand_contacts(canonical)
```

## Fetching a file yourself

Instead of a bundled fixture, `fetch_mmcif()` downloads an entry from
PDBe or RCSB (needs network access, caches to disk so re-runs are
offline):

```python
from pandora.ingestion import fetch_mmcif

fetch_mmcif(
    entry_id="1cbs",
    provider="pdbe",
    source_uri=None,
    output_dir=mmcif_dir,
)
```

Its return value is provenance about the download, not the structure
itself — feed the file it just wrote to `mmcif_to_structure()` as
above.

## Next steps

- [Policies](../reference/policies.md) — every canonicalisation policy field, with
  worked examples.
- [Usage](../usage/ingestion.md) — a per-stage reference for ingestion,
  canonicalisation, metadata, annotation, similarity, and the CLI.
- [Recipes](../recipes/ppi-01.md) — end-to-end multi-structure workflows
  (clustering, leakage-safe splits).
- [Examples](../examples.md) — runnable scripts under `examples/`.

## Comment
