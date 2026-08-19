from pathlib import Path

from pandora.parsing import mmcif_to_structure
from pandora.provenance import build_dataset_manifest, build_provenance_bundle
from pandora.schemas.dataset import (
    DatasetCurationPolicy,
    DeduplicationProvenance,
)
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


def test_build_dataset_manifest():
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    bundle = build_provenance_bundle(structure)
    policy = DatasetCurationPolicy(
        policy_id="p", policy_name="p", policy_version="1.0.0"
    )
    dedup_prov = DeduplicationProvenance(
        deduplicated_at="2026-01-01T00:00:00+00:00", enabled=True
    )

    manifest = build_dataset_manifest(
        dataset_id="d1",
        dataset_name="Dataset One",
        dataset_version="1.0.0",
        curation_policy=policy,
        deduplication=dedup_prov,
        splits={"train": [structure.entry_id]},
        structures=[bundle],
    )

    assert manifest.dataset_id == "d1"
    assert manifest.curation_policy == policy
    assert manifest.deduplication == dedup_prov
    assert manifest.splits == {"train": [structure.entry_id]}
    assert manifest.structures == [bundle]
    assert manifest.excluded == []
