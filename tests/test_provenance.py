import shutil
import stat
from pathlib import Path

import gemmi
import pytest

import pandora
from pandora.annotations.entry import annotate_structure_counts
from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import curate_structure, extract_chain_records
from pandora.export import write_json, write_records
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.provenance import (
    build_dataset_manifest,
    build_provenance_bundle,
    build_step_summary,
    canonicalisation_entry_record,
    canonicalisation_provenance_from_record,
    curation_entry_record,
    curation_provenance_from_record,
    reproduce_dataset,
)
from pandora.schemas.canonicalisation import canonicalisationPolicy
from pandora.schemas.dataset import (
    DatasetCurationPolicy,
    DeduplicationProvenance,
)
from pandora.schemas.ingestion import IngestionProvenance
from pandora.schemas.provenance import (
    AnnotationProvenanceRecord,
    ProvenanceBundle,
    StepSummary,
)
from pandora.schemas.similarity import (
    ClusteringProvenance,
    HitFilter,
    SimilaritySearch,
)
from pandora.similarity.hits import HIT_COLUMNS

MMCIF_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"
MMCIF_PATH = MMCIF_DIR / "04" / "104m.cif"


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


def test_build_step_summary_embeds_policy_and_versions(tmp_path):
    policy = DatasetCurationPolicy(
        policy_id="p", policy_name="n", policy_version="1"
    )
    summary = build_step_summary(
        "curation",
        policy,
        input=tmp_path / "in",
        output=tmp_path / "out",
        started_at="2026-11-04T00:00:00Z",
        counts={"entries_in": 2},
        previous_summary=tmp_path / "canonicalisation.json",
    )
    assert summary.step == "curation"
    assert summary.policy["policy_id"] == "p"
    assert summary.pandora_version == pandora.__version__
    assert summary.gemmi_version == gemmi.__version__
    assert summary.input == str(tmp_path / "in")
    assert summary.output == str(tmp_path / "out")
    assert summary.previous_summary == str(tmp_path / "canonicalisation.json")
    assert summary.counts == {"entries_in": 2}
    path = write_json(summary, tmp_path / "curation.json")
    assert StepSummary.model_validate_json(path.read_text()) == summary


def _curation_policy() -> DatasetCurationPolicy:
    return DatasetCurationPolicy(
        policy_id="p", policy_name="n", policy_version="1"
    )


def test_entry_records_write_to_parquet(tmp_path):
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / "ay" / "1ayi.cif"))
    canon, _, cprov = canonicalise_structure(
        structure,
        canonicalisationPolicy(
            policy_id="c", policy_name="n", policy_version="1"
        ),
    )
    metadata = collect_metadata(canon)
    curated, _, kprov = curate_structure(canon, metadata, _curation_policy())
    assert curated is not None
    c_row = canonicalisation_entry_record(
        canon, cprov, tmp_path / "ay" / "1ayi.cif"
    )
    k_row = curation_entry_record(curated, metadata, kprov)
    assert c_row.entry_id == canon.entry_id
    assert c_row.structure_path == str(tmp_path / "ay" / "1ayi.cif")
    assert c_row.policy_id == cprov.policy_id
    assert c_row.transforms == cprov.transforms
    assert k_row.structure_path is None
    assert k_row.chain_count == len(extract_chain_records(curated))
    assert k_row.resolution == metadata.quality.resolution
    assert k_row.policy_id == "p"
    # cprov.report may be {}; the flat row must still write to Parquet.
    write_records([c_row], tmp_path / "c.parquet")
    write_records([k_row], tmp_path / "k.parquet")


def test_provenance_bundle_carries_curation():
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / "ay" / "1ayi.cif"))
    _, _, kprov = curate_structure(
        structure, collect_metadata(structure), _curation_policy()
    )
    bundle = build_provenance_bundle(structure, curation=kprov)
    assert bundle.curation == kprov


def test_provenance_from_record_round_trips():
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / "ay" / "1ayi.cif"))
    canon, _, cprov = canonicalise_structure(
        structure,
        canonicalisationPolicy(
            policy_id="c", policy_name="n", policy_version="1"
        ),
    )
    metadata = collect_metadata(canon)
    curated, _, kprov = curate_structure(canon, metadata, _curation_policy())
    assert (
        canonicalisation_provenance_from_record(
            canonicalisation_entry_record(canon, cprov)
        )
        == cprov
    )
    assert (
        curation_provenance_from_record(
            curation_entry_record(curated, metadata, kprov)
        )
        == kprov
    )
    reported = cprov.model_copy(update={"report": {"warnings": 2, "errors": 0}})
    assert (
        canonicalisation_provenance_from_record(
            canonicalisation_entry_record(canon, reported)
        )
        == reported
    )
