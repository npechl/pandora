from pathlib import Path

from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure

MMCIF_PATH = (
    Path(__file__).parent.parent / "datasets" / "dev" / "mmcif" / "104m.cif"
)


def test_collect_metadata_from_fixture():
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))

    metadata = collect_metadata(structure)

    assert metadata.entry_id == "104M"
    assert metadata.quality.resolution == 1.71
    assert metadata.quality.r_work == 0.154
    assert metadata.quality.r_free == 0.225
    assert [t.ncbi_taxon_id for t in metadata.taxonomies] == [9755]
    assert [(e.entity_id, e.chain_ids) for e in metadata.entities] == [
        ("1", ["A"]),
        ("2", ["B"]),
        ("3", ["C"]),
        ("4", ["D"]),
        ("5", ["E"]),
    ]
    assert [lig.comp_id for lig in metadata.ligands] == ["SO4", "HEM", "NBN"]
    assert {m.accession for m in metadata.uniprot_mappings} == {"P02185"}
    assert metadata.raw_categories == sorted(structure.raw)
