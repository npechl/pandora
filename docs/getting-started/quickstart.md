# Quickstart

This walks through one typical workflow: parse an mmCIF file, canonicalise it, collect its metadata, and compute a couple of annotations. It uses files already checked into the repo, so it runs fully offline — no network access, no PDB account. Steps 1-4 walk through a single structure so each stage's output is easy to read; step 5 then scales the same functions up to all ~100 bundled entries at once.

Pandora's base install (`pip install -e .`) and a checkout of this repo are all you need — steps 1-5 below use no optional extras. The last, optional section (fetching a file yourself over the network) needs the `ingestion` extra instead; see [Installation](installation.md).

## The fixtures we'll use

`datasets/dev/mmcif/` bundles ~100 small mmCIF files for local development and docs — real PDB entries, no network access needed. Steps 1-4 use `104m.cif`, sperm whale myoglobin bound to a small-molecule ligand (heme); any other entry under that directory works the same way if you want to follow along with a different one.

## 1. Parse

`mmcif_to_structure()` reads the raw mmCIF file and returns Pandora's typed `Structure`, plus a diagnostics bundle and a status (`"success"`, `"warning"`, or `"failed"`):

```python
from pathlib import Path
from pandora.parsing import mmcif_to_structure

mmcif_dir = Path("datasets/dev/mmcif")
structure, diagnostics, status = mmcif_to_structure(str(mmcif_dir / "104m.cif"))

print(status, len(structure.atoms))
# success 1450
```

## 2. Canonicalise

`canonicalise_structure()` applies a policy — chain ID remapping, altloc resolution, ligand filtering, and more — and returns the canonical `Structure` alongside the mappings back to the original identifiers and a provenance record of what actually changed:

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

An unconfigured policy still resolves alternate locations (the default `altloc_rules` strategy) and annotates missing atoms/residues rather than silently dropping them — nothing here changed chain IDs or residue numbers, because this entry didn't need it. See [Policies](../reference/policies.md) for every field you can configure.

## 3. Collect metadata

`collect_metadata()` pulls the source-backed entry, quality, taxonomy, entity, ligand, and UniProt-mapping records straight from the mmCIF categories Pandora kept:

```python
from pandora.metadata import collect_metadata

metadata = collect_metadata(canonical)

print(metadata.entry.title)
# SPERM WHALE MYOGLOBIN N-BUTYL ISOCYANIDE AT PH 7.0
print(metadata.quality.experimental_method, metadata.quality.resolution)
# X-ray diffraction 1.71
```

## 4. Annotate

Annotations are derived layers computed from the canonical structure. `annotate_structure_counts()` gives you a quick per-entry summary; `annotate_ligand_contacts()` finds polymer residues near each ligand:

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

`HEM` (the heme group) has the most polymer residues in contact, which is what you'd expect from myoglobin's oxygen-binding pocket.

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

## 5. Scale it up: the whole directory

Nothing about steps 1-4 is single-file-only — each function takes one `Structure` and returns one `Structure`. Loop over all ~100 fixtures under `datasets/dev/mmcif/` and feed each through the same functions to get dataset-level numbers instead of one entry's:

```python
from collections import Counter
from pathlib import Path

from pandora.annotations import annotate_ligand_contacts
from pandora.canonicalisation import canonicalise_structure
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import canonicalisationPolicy

mmcif_dir = Path("datasets/dev/mmcif")
policy = canonicalisationPolicy(
    policy_id="quickstart", policy_name="Default", policy_version="1.0.0"
)

# canonicalise_structure() raises ValueError on a "failed" validation
# status — catch it per entry so one bad structure doesn't abort the batch.
structures = {}
for path in sorted(mmcif_dir.glob("*.cif")):
    structure, _, status = mmcif_to_structure(str(path))
    try:
        canonical, _, canon_prov = canonicalise_structure(structure, policy)
    except ValueError as exc:
        print(f"[{structure.entry_id}] skipped: {exc}")
        continue
    structures[canonical.entry_id] = canonical

print(f"parsed+canonicalised: {len(structures)} structures")
# parsed+canonicalised: 100 structures

resolutions = []
for structure in structures.values():
    resolution = collect_metadata(structure).quality.resolution
    if resolution is not None:
        resolutions.append(resolution)

print(
    f"resolution range: {min(resolutions):.2f}-{max(resolutions):.2f} Å ({len(resolutions)} entries)"
)
# resolution range: 0.60-16.80 Å (79 entries)

ligand_counts = Counter()
for structure in structures.values():
    for ligand in annotate_ligand_contacts(structure).data["ligands"]:
        ligand_counts[ligand["ligand_comp_id"]] += 1

print(ligand_counts.most_common(5))
# [('ZN', 56), ('HEM', 11), ('CA', 11), ('POL', 10), ('CL', 9)]
```

Zinc alone shows up in 56 of the 100 entries — this fixture set leans heavily on metalloproteins and small ligand-bound structures.

`curate_structure()` and `deduplicate_structures()` (both pure Python — no external tools) filter and dedupe that same batch:

```python
from pandora.datasets import curate_structure, deduplicate_structures
from pandora.schemas.dataset import DatasetCurationPolicy, DeduplicationRules

curation_policy = DatasetCurationPolicy(
    policy_id="quickstart-curation",
    policy_name="Default",
    policy_version="1.0.0",
)
exclusions = []
for entry_id in list(structures):
    curated, exclusion, _ = curate_structure(
        structures[entry_id], None, curation_policy
    )
    if curated is None:
        exclusions.append(exclusion)
        del structures[entry_id]
    else:
        structures[entry_id] = curated

retained, removed, dedup_prov = deduplicate_structures(
    list(structures.values()), DeduplicationRules(enabled=True)
)
exclusions.extend(removed)
print(f"curation+dedup: {len(retained)} retained, {len(exclusions)} excluded")
# curation+dedup: 100 retained, 0 excluded
```

This dev fixture set is already clean, so nothing gets excluded at default policy settings — see [Policies](../reference/policies.md) for the quality/organism/content rules that would start filtering it. See [PPI interface dataset](../recipes/ppi-01.md) for what comes after this: reshaping into flat records, sequence-similarity clustering, and a leakage-safe train/val/test split.

## Fetching a file yourself

Instead of a bundled fixture, `fetch_mmcif()` downloads an entry from PDBe or RCSB (needs network access, caches to disk so re-runs are offline):

```python
from pandora.ingestion import fetch_mmcif

fetch_mmcif(
    entry_id="1cbs",
    provider="pdbe",
    source_uri=None,
    output_dir=mmcif_dir,
)
```

Its return value is provenance about the download, not the structure itself — feed the file it just wrote to `mmcif_to_structure()` as above.

## Next steps

- [Policies](../reference/policies.md) — every canonicalisation policy field, with worked examples.
- [Usage](../usage/ingestion.md) — a per-stage reference for ingestion, canonicalisation, metadata, annotation, similarity, and the CLI.
- [Recipes](../recipes/ppi-01.md) — end-to-end multi-structure workflows (clustering, leakage-safe splits).
