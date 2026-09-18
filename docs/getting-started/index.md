# What is `Pandora`?

Pandora turns raw PDBe mmCIF files into typed, policy-driven, ML-ready protein structure datasets. Structural biology data usually arrives as loosely-typed, inconsistently-annotated mmCIF — Pandora parses it into a typed `Structure`, applies a configurable policy to normalize it (chain IDs, altlocs, ligands, missing atoms), and derives the metadata, annotations, and provenance records you need for downstream ML work.

!!! warning "Status"
    Pandora is under active development. Ingestion, parsing,
    canonicalisation, metadata, annotations, export, dataset curation,
    provenance (per-structure bundles, dataset manifests, and
    reproducing a dataset from one), and the `pandora` CLI are all
    implemented today.

Every stage is a plain function: pass a `Structure` (or a typed record) in, get one out. Nothing is hidden behind a framework object or global state, so you can call one stage on its own or chain all of them into a pipeline.

## See it in action

This runs fully offline against a fixture already checked into the repo — sperm whale myoglobin bound to a heme ligand:

```python
from pathlib import Path
from pandora.parsing import mmcif_to_structure
from pandora.canonicalisation import canonicalise_structure
from pandora.metadata import collect_metadata
from pandora.annotations import annotate_ligand_contacts
from pandora.schemas.canonicalisation import canonicalisationPolicy

mmcif_dir = Path("datasets/dev/mmcif")
structure, diagnostics, status = mmcif_to_structure(str(mmcif_dir / "104m.cif"))

policy = canonicalisationPolicy(
    policy_id="quickstart", policy_name="Default", policy_version="1.0.0"
)
canonical, mappings, canon_prov = canonicalise_structure(structure, policy)

metadata = collect_metadata(canonical)
contacts = annotate_ligand_contacts(canonical)

print(metadata.entry.title)
# SPERM WHALE MYOGLOBIN N-BUTYL ISOCYANIDE AT PH 7.0
print([(l["ligand_comp_id"], l["contact_count"]) for l in contacts.data["ligands"]])
# [('SO4', 4), ('HEM', 15), ('NBN', 4)]
```

`HEM` (the heme group) has the most polymer residues in contact — exactly what you'd expect from myoglobin's oxygen-binding pocket. See [Usage](../usage/overview.md) for this same pipeline broken down stage by stage, library and CLI side by side.

## Install

Not yet on PyPI — install from source:

```bash
git clone https://github.com/npechl/pandora.git
cd pandora

# Base install
pip install -e .
```

See [Installation](installation.md) for what each optional extra (`ingestion`, `export`, `similarity`, ...) unlocks.

## Where to go next

- **New here?** [Quickstart](quickstart.md) — a five-minute, fully offline walkthrough using the bundled sample data.
- **Want the full pipeline, stage by stage?** [Usage](../usage/overview.md) has the library call and its CLI equivalent side by side.
- **Configuring canonicalisation?** [Policies](../reference/policies.md) covers every policy field.
- **Building a training dataset?** [Recipes](../recipes/ppi-01.md) has end-to-end multi-structure workflows — clustering, leakage-safe splits.
- **Looking for a specific function?** [Functions](../reference/functions.md) and [Data models](../reference/schemas.md) are generated straight from the docstrings.
- **Hit an unfamiliar term?** Check the [Glossary](glossary.md).
</content>
