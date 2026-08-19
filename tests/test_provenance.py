from pathlib import Path

from pandora.parsing import mmcif_to_structure
from pandora.provenance import build_provenance_bundle
from pandora.schemas.ingestion import IngestionProvenance

MMCIF_PATH = (
    Path(__file__).parent.parent / "datasets" / "dev" / "mmcif" / "104m.cif"
)


def test_build_provenance_bundle():
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    ingestion = IngestionProvenance(provider="pdbe", from_cache=True)

    bundle = build_provenance_bundle(structure, ingestion=ingestion)

    assert bundle.entry_id == structure.entry_id
    assert bundle.ingestion == ingestion
    assert bundle.canonicalisation is None
