import shutil
from pathlib import Path

from pandora.annotations.entry import annotate_structure_counts
from pandora.parsing import mmcif_to_structure
from pandora.provenance import (
    build_dataset_manifest,
    build_provenance_bundle,
    reproduce_dataset,
)
from pandora.schemas.dataset import (
    DatasetCurationPolicy,
    DeduplicationProvenance,
)
from pandora.schemas.ingestion import IngestionProvenance
from pandora.schemas.provenance import (
    AnnotationProvenanceRecord,
    ProvenanceBundle,
)

MMCIF_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"
MMCIF_PATH = MMCIF_DIR / "104m.cif"


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


def test_reproduce_dataset(tmp_path, monkeypatch):
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    entry_id = structure.entry_id

    def fake_fetch_mmcif(
        entry_id, provider, source_uri, output_dir, fetch_options=None
    ):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(MMCIF_PATH, output_dir / f"{entry_id.lower()}.cif")
        return IngestionProvenance(provider=provider, source_uri=source_uri)

    monkeypatch.setattr(
        "pandora.provenance.reproduce.fetch_mmcif", fake_fetch_mmcif
    )

    policy = DatasetCurationPolicy(
        policy_id="p", policy_name="p", policy_version="1.0.0"
    )
    manifest = build_dataset_manifest(
        dataset_id="d1",
        dataset_name="Dataset One",
        dataset_version="1.0.0",
        curation_policy=policy,
        structures=[
            ProvenanceBundle(
                entry_id=entry_id,
                pandora_version="0.0.0",
                generated_at="2026-01-01T00:00:00+00:00",
                ingestion=IngestionProvenance(provider="pdbe"),
            )
        ],
    )

    structures, new_manifest = reproduce_dataset(manifest, tmp_path)

    assert set(structures) == {entry_id}
    assert new_manifest.dataset_id == "d1"
    assert new_manifest.curation_policy == policy
    assert new_manifest.excluded == []
    assert [b.entry_id for b in new_manifest.structures] == [entry_id]
    assert new_manifest.structures[0].ingestion.provider == "pdbe"


def test_reproduce_dataset_regenerates_annotations(tmp_path, monkeypatch):
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    entry_id = structure.entry_id
    original_layer = annotate_structure_counts(structure)

    def fake_fetch_mmcif(
        entry_id, provider, source_uri, output_dir, fetch_options=None
    ):
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(MMCIF_PATH, output_dir / f"{entry_id.lower()}.cif")
        return IngestionProvenance(provider=provider, source_uri=source_uri)

    monkeypatch.setattr(
        "pandora.provenance.reproduce.fetch_mmcif", fake_fetch_mmcif
    )

    manifest = build_dataset_manifest(
        dataset_id="d1",
        dataset_name="Dataset One",
        dataset_version="1.0.0",
        structures=[
            ProvenanceBundle(
                entry_id=entry_id,
                pandora_version="0.0.0",
                generated_at="2026-01-01T00:00:00+00:00",
                ingestion=IngestionProvenance(provider="pdbe"),
                annotations=[
                    AnnotationProvenanceRecord(
                        layer_name=original_layer.layer_name,
                        layer_type=original_layer.layer_type,
                        method=original_layer.method,
                        target_ids=original_layer.target_ids,
                        parameters=original_layer.parameters,
                    )
                ],
            )
        ],
    )

    _, new_manifest = reproduce_dataset(manifest, tmp_path)

    new_bundle = new_manifest.structures[0]
    assert len(new_bundle.annotations) == 1
    assert new_bundle.annotations[0].layer_type == "structure_counts"


def test_reproduce_dataset_requires_ingestion_provenance(tmp_path):
    manifest = build_dataset_manifest(
        dataset_id="d1",
        dataset_name="Dataset One",
        dataset_version="1.0.0",
        structures=[
            ProvenanceBundle(
                entry_id="1abc",
                pandora_version="0.0.0",
                generated_at="2026-01-01T00:00:00+00:00",
            )
        ],
    )

    try:
        reproduce_dataset(manifest, tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for missing ingestion")
