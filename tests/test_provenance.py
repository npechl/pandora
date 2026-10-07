import shutil
import stat
from pathlib import Path

import pytest

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

from pandora.schemas.similarity import (
    ClusteringProvenance,
    HitFilter,
    SimilaritySearch,
)
from pandora.similarity.hits import HIT_COLUMNS

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

    with pytest.raises(ValueError):
        reproduce_dataset(manifest, tmp_path)


FAKE_MMSEQS = """#!/bin/sh
if [ "$1" = version ]; then echo 15.6f452; exit 0; fi
printf '104M\\t104M\\t1.0\\t10\\t1.0\\t1.0\\n' > "$4"
"""


def _patch_fetch(monkeypatch):
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


def _manifest_with_clustering(clustering):
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    return build_dataset_manifest(
        dataset_id="d1",
        dataset_name="Dataset One",
        dataset_version="1.0.0",
        clustering=clustering,
        structures=[
            ProvenanceBundle(
                entry_id=structure.entry_id,
                pandora_version="0.0.0",
                generated_at="2026-01-01T00:00:00+00:00",
                ingestion=IngestionProvenance(provider="pdbe"),
            )
        ],
    )


def _clustering(search, hit_filter=None):
    return ClusteringProvenance(
        clustered_at="2026-01-01T00:00:00+00:00",
        hit_filter=hit_filter or HitFilter(min_score=0.9),
        search=search,
        n_clusters=1,
    )


def test_reproduce_reruns_computed_search(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    mmseqs = bin_dir / "mmseqs"
    mmseqs.write_text(FAKE_MMSEQS)
    mmseqs.chmod(mmseqs.stat().st_mode | stat.S_IEXEC)
    search = SimilaritySearch(
        engine="MMseqs2",
        hits_path="/somewhere/else/hits.tsv",
        columns=list(HIT_COLUMNS["MMseqs2"]),
        parameters={
            "mmseqs_bin": str(mmseqs),
            "sensitivity": 5.7,
            "max_seqs": 300,
            "mmseqs_options": [],
        },
    )
    out = tmp_path / "out"

    _, new_manifest = reproduce_dataset(
        _manifest_with_clustering(_clustering(search)), out
    )

    assert (out / "hits.tsv").exists()
    assert new_manifest.clustering.search.hits_path == str(out / "hits.tsv")
    assert new_manifest.clustering.search.parameters == search.parameters
    assert new_manifest.clustering.hit_filter == HitFilter(min_score=0.9)


def test_reproduce_reuses_existing_precomputed_hits(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    hits = tmp_path / "pre.tsv"
    hits.write_text("104M\t104M\t1.0\t10\t1.0\t1.0\n")
    search = SimilaritySearch(
        engine="MMseqs2",
        hits_path=str(hits),
        columns=list(HIT_COLUMNS["MMseqs2"]),
        origin="precomputed",
    )

    _, new_manifest = reproduce_dataset(
        _manifest_with_clustering(_clustering(search)), tmp_path / "out"
    )

    assert new_manifest.clustering.search == search


def test_reproduce_precomputed_hits_missing_raises(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    search = SimilaritySearch(
        engine="MMseqs2",
        hits_path=str(tmp_path / "gone.tsv"),
        columns=list(HIT_COLUMNS["MMseqs2"]),
        origin="precomputed",
    )

    with pytest.raises(ValueError, match="precomputed"):
        reproduce_dataset(
            _manifest_with_clustering(_clustering(search)), tmp_path / "out"
        )


def test_reproduce_old_manifest_clustering_raises(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    old = ClusteringProvenance.model_validate(
        {
            "clustered_at": "2026-01-01T00:00:00+00:00",
            "threshold": 0.9,
            "n_relationships": 1,
            "n_clusters": 1,
        }
    )

    with pytest.raises(ValueError, match="did not record"):
        reproduce_dataset(_manifest_with_clustering(old), tmp_path / "out")


def test_reproduce_interface_filter_raises(tmp_path, monkeypatch):
    _patch_fetch(monkeypatch)
    search = SimilaritySearch(
        engine="Foldseek",
        hits_path=str(tmp_path / "hits.tsv"),
        columns=list(HIT_COLUMNS["Foldseek"]),
    )
    clustering = _clustering(search, HitFilter(min_interface_coverage=0.5))

    with pytest.raises(ValueError, match="interface"):
        reproduce_dataset(
            _manifest_with_clustering(clustering), tmp_path / "out"
        )
