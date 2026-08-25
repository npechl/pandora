# Metadata

`pandora.metadata` extracts source-backed metadata (never derived —
everything here comes straight from an mmCIF category) from a parsed
`Structure`. See
[Functions](../reference/functions.md#pandora.metadata) for full
signatures, and [Schemas](../reference/schemas.md) for every record's
fields.

## Collect everything at once

`collect_metadata()` is the orchestrator: it calls every `extract_*()`
function below and bundles the results into one `MetadataRecord`.

```python
from pandora.parsing import mmcif_to_structure
from pandora.metadata import collect_metadata

structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/104m.cif")
metadata = collect_metadata(structure)

print(metadata.entry.title)
# 'SPERM WHALE MYOGLOBIN N-BUTYL ISOCYANIDE AT PH 7.0'
print([ligand.comp_id for ligand in metadata.ligands])
# ['SO4', 'HEM', 'NBN']
print(len(metadata.taxonomies))
# 1
```

## Extract one category at a time

Each `extract_*()` function reads one or a few specific mmCIF
categories and is independently callable — reach for these when you
only need one piece, not the full `MetadataRecord`.

```python
from pandora.metadata.mmcif import (
    extract_entry_metadata,
    extract_entity_metadata,
    extract_ligand_metadata,
    extract_quality,
    extract_taxonomies,
    extract_taxonomy,
    extract_uniprot_mappings,
)

entry = extract_entry_metadata(structure)
print(entry.title, entry.doi)
# 'SPERM WHALE MYOGLOBIN N-BUTYL ISOCYANIDE AT PH 7.0' None

entities = extract_entity_metadata(structure)
print([(e.entity_id, e.entity_type) for e in entities])
# [('1', 'polymer'), ('2', 'non-polymer'), ('3', 'non-polymer'), ('4', 'non-polymer'), ('5', 'water')]

ligands = extract_ligand_metadata(structure)
print([(l.comp_id, l.name) for l in ligands])
# [('SO4', "'SULFATE ION'"), ('HEM', "'PROTOPORPHYRIN IX CONTAINING FE'"), ...]

quality = extract_quality(structure)
print(quality.experimental_method, quality.resolution)
# 'X-ray diffraction' 1.71

taxonomies = extract_taxonomies(structure)
print([t.organism_scientific for t in taxonomies])
# ["'Physeter catodon'"]

taxon = extract_taxonomy(structure)  # first taxonomy record, or None
print(taxon.organism_scientific)
# 'Physeter catodon'

mappings = extract_uniprot_mappings(structure)
print([m.accession for m in mappings])
# ['P02185', 'P02185']
```

!!! note "Quoted values"
    `name`/`organism_scientific` above keep their literal `'...'`
    delimiters (e.g. `"'SULFATE ION'"`) — this is the raw mmCIF token,
    not a formatting bug. See the note on `parsing/mmcif.py::_cs` in
    `CLAUDE.md` for why.

## Read a raw category directly

Not every mmCIF category has a typed extractor yet.
`extract_metadata_category()` is the escape hatch: it reads any
category straight out of `Structure.raw`.

```python
from pandora.metadata.mmcif import extract_metadata_category

rows = extract_metadata_category(structure, "_exptl")
print(rows)
# [{'entry_id': '104M', 'method': "'X-ray diffraction'", 'crystals_number': '1'}]
```

Pass `columns=[...]` to restrict which fields come back. `extract_metadata()`
is an older-name alias for the same function, kept for compatibility.
