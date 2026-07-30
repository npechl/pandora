from pathlib import Path

from pandora.parsing import mmcif_to_structure
from pandora.provenance import build_provenance_bundle, compute_checksum
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
    assert bundle.checksums.structure_checksum == compute_checksum(structure)


def test_compute_checksum_is_stable_and_content_sensitive():
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    other, _, _ = mmcif_to_structure(str(MMCIF_PATH))

    assert compute_checksum(structure) == compute_checksum(other)

    mutated = structure.model_copy(update={"entry_id": "different"})
    assert compute_checksum(mutated) != compute_checksum(structure)
