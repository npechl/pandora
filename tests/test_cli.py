import json
import shutil
import stat
from pathlib import Path

import pytest

from pandora.cli.app import _build_parser, main
from pandora.schemas.similarity import SimilaritySearch
from pandora.similarity.hits import HIT_COLUMNS

FIXTURES_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"
CANON_POLICY = (
    Path(__file__).parent.parent / "datasets" / "canonicalisation.yaml"
)

_SUBCOMMANDS = [
    "fetch",
    "ingest",
    "canonicalise",
    "curate",
    "dedup",
    "similarity",
    "cluster",
    "partition",
    "annotate",
    "manifest",
    "reproduce",
    "export",
]


def test_parser_registers_all_subcommands():
    help_text = _build_parser().format_help()
    for name in _SUBCOMMANDS:
        assert name in help_text


def test_ingest_local_directory(tmp_path):
    mirror_dir = tmp_path / "mirror"
    mirror_dir.mkdir()
    (mirror_dir / "104m.cif").write_text(
        (FIXTURES_DIR / "104m.cif").read_text()
    )

    output_dir = tmp_path / "raw"
    main(
        [
            "ingest",
            "--input-dir",
            str(mirror_dir),
            "--output-dir",
            str(output_dir),
            "--source-uri",
            "dev-fixture-snapshot",
        ]
    )

    prov = json.loads((output_dir / "ingestion_provenance.json").read_text())
    assert prov["104M"]["provider"] == "local"
    assert prov["104M"]["source_uri"] == "dev-fixture-snapshot"


def test_canonicalise_curate_dedup_export_chain(tmp_path):
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    (raw_dir / "104m.cif").write_text((FIXTURES_DIR / "104m.cif").read_text())

    canonical_dir = tmp_path / "canonical"
    main(
        [
            "canonicalise",
            "--input-dir",
            str(raw_dir),
            "--policy",
            str(CANON_POLICY),
            "--output-dir",
            str(canonical_dir),
        ]
    )
    assert (canonical_dir / "104m.cif").exists()
    canon_prov = json.loads(
        (canonical_dir / "canonicalisation_provenance.json").read_text()
    )
    assert "104M" in canon_prov
    canon_mappings = json.loads(
        (canonical_dir / "canonicalisation_mappings.json").read_text()
    )
    assert canon_mappings["104M"]["chain_id_mapping"]["items"]

    curation_policy = tmp_path / "curation.yaml"
    curation_policy.write_text(
        "policy_id: test\npolicy_name: Test\npolicy_version: 1.0.0\n"
    )
    curated_dir = tmp_path / "curated"
    main(
        [
            "curate",
            "--input-dir",
            str(canonical_dir),
            "--policy",
            str(curation_policy),
            "--output-dir",
            str(curated_dir),
        ]
    )
    assert (curated_dir / "104m.cif").exists()

    deduped_dir = tmp_path / "deduped"
    main(
        [
            "dedup",
            "--input-dir",
            str(curated_dir),
            "--output-dir",
            str(deduped_dir),
        ]
    )
    assert (deduped_dir / "104m.cif").exists()

    exported = tmp_path / "104m.json"
    main(
        [
            "export",
            "--input",
            str(deduped_dir / "104m.cif"),
            "--output",
            str(exported),
        ]
    )
    assert json.loads(exported.read_text())["entry_id"] == "104M"


def _cluster_inputs(tmp_path):
    input_dir = tmp_path / "clustered"
    input_dir.mkdir()
    for name in ("104m", "112m", "118l"):
        (input_dir / f"{name}.cif").touch()
    hits = tmp_path / "hits.tsv"
    hits.write_text("104M\t112M\t0.95\t100\t0.95\t0.95\n")
    search = tmp_path / "hits.tsv.search.json"
    search.write_text(
        SimilaritySearch(
            engine="MMseqs2",
            hits_path=str(hits),
            columns=list(HIT_COLUMNS["MMseqs2"]),
        ).model_dump_json()
    )
    return input_dir, search


def test_cluster_with_pairs(tmp_path):
    input_dir, search = _cluster_inputs(tmp_path)
    pairs = tmp_path / "pairs.json"
    pairs.write_text(json.dumps([["104M", "112M"], ["104M", "118L"]]))
    output = tmp_path / "clusters.json"

    main(
        [
            "cluster",
            "--input-dir",
            str(input_dir),
            "--search",
            str(search),
            "--min-score",
            "0.9",
            "--pairs",
            str(pairs),
            "--output",
            str(output),
        ]
    )

    cluster_pairs = json.loads(
        (output.parent / "cluster_pairs.json").read_text()
    )
    by_pair = {(p["item_id_1"], p["item_id_2"]): p for p in cluster_pairs}
    same = by_pair[("104M", "112M")]
    assert same["cluster_id_1"] == same["cluster_id_2"]
    different = by_pair[("104M", "118L")]
    assert different["cluster_id_1"] != different["cluster_id_2"]


def test_cluster_reads_hit_filter_yaml_and_records_it(tmp_path):
    input_dir, search = _cluster_inputs(tmp_path)
    hit_filter = tmp_path / "filter.yaml"
    hit_filter.write_text("min_score: 0.99\nmin_coverage: 0.5\n")
    output = tmp_path / "clusters.json"

    main(
        [
            "cluster",
            "--input-dir",
            str(input_dir),
            "--search",
            str(search),
            "--hit-filter",
            str(hit_filter),
            "--output",
            str(output),
        ]
    )

    clusters = json.loads(output.read_text())
    assert all(c["n_components"] == 1 for c in clusters)  # 0.95 < 0.99
    provenance = json.loads(
        (output.parent / "cluster_provenance.json").read_text()
    )
    assert provenance["hit_filter"]["min_score"] == 0.99
    assert provenance["hit_filter"]["min_coverage"] == 0.5
    assert provenance["search"]["engine"] == "MMseqs2"


def test_cluster_flag_overrides_hit_filter_yaml(tmp_path):
    input_dir, search = _cluster_inputs(tmp_path)
    hit_filter = tmp_path / "filter.yaml"
    hit_filter.write_text("min_score: 0.99\n")
    output = tmp_path / "clusters.json"

    main(
        [
            "cluster",
            "--input-dir",
            str(input_dir),
            "--search",
            str(search),
            "--hit-filter",
            str(hit_filter),
            "--min-score",
            "0.9",
            "--output",
            str(output),
        ]
    )

    clusters = json.loads(output.read_text())
    assert ["104M", "112M"] in [c["components"] for c in clusters]


FAKE_MMSEQS_CLI = """#!/bin/sh
if [ "$1" = version ]; then echo 15.6f452; exit 0; fi
echo "$@" > "$(dirname "$0")/args.txt"
printf '104M\\t112M\\t0.99\\t100\\t0.95\\t0.95\\n' > "$4"
"""


def test_similarity_writes_hits_and_search_record(tmp_path):
    mmseqs = tmp_path / "mmseqs"
    mmseqs.write_text(FAKE_MMSEQS_CLI)
    mmseqs.chmod(mmseqs.stat().st_mode | stat.S_IEXEC)
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    for name in ("104m", "112m"):
        shutil.copy(FIXTURES_DIR / f"{name}.cif", input_dir / f"{name}.cif")
    output = tmp_path / "hits.tsv"

    main(
        [
            "similarity",
            "--input-dir",
            str(input_dir),
            "--engine",
            "mmseqs2",
            "--mmseqs-bin",
            str(mmseqs),
            "--max-seqs",
            "500",
            "--output",
            str(output),
        ]
    )

    record = json.loads((tmp_path / "hits.tsv.search.json").read_text())
    assert output.read_text().startswith("104M\t112M")
    assert record["hits_path"] == str(output)
    assert record["parameters"]["max_seqs"] == 500
    assert record["origin"] == "computed"


def test_similarity_wraps_precomputed_hits(tmp_path):
    hits = tmp_path / "elsewhere.tsv"
    hits.write_text("104M\t112M\t0.99\t100\t0.95\t0.95\n")
    output = tmp_path / "hits.tsv"

    main(
        [
            "similarity",
            "--engine",
            "mmseqs2",
            "--precomputed-hits",
            str(hits),
            "--output",
            str(output),
        ]
    )

    record = json.loads((tmp_path / "hits.tsv.search.json").read_text())
    assert record["origin"] == "precomputed"
    assert record["hits_path"] == str(hits)


def test_similarity_needs_input_dir_without_precomputed_hits(tmp_path):
    with pytest.raises(SystemExit, match="--input-dir"):
        main(
            [
                "similarity",
                "--engine",
                "mmseqs2",
                "--output",
                str(tmp_path / "hits.tsv"),
            ]
        )


def test_cluster_requires_a_threshold(tmp_path):
    input_dir, search = _cluster_inputs(tmp_path)

    with pytest.raises(SystemExit, match="--min-score"):
        main(
            [
                "cluster",
                "--input-dir",
                str(input_dir),
                "--search",
                str(search),
                "--output",
                str(tmp_path / "clusters.json"),
            ]
        )


def test_similarity_records_absolute_hits_path(tmp_path, monkeypatch):
    (tmp_path / "elsewhere.tsv").write_text(
        "104M\t112M\t0.99\t100\t0.95\t0.95\n"
    )
    monkeypatch.chdir(tmp_path)

    main(
        [
            "similarity",
            "--engine",
            "mmseqs2",
            "--precomputed-hits",
            "elsewhere.tsv",
            "--output",
            "hits.tsv",
        ]
    )

    record = json.loads((tmp_path / "hits.tsv.search.json").read_text())
    assert record["hits_path"] == str((tmp_path / "elsewhere.tsv").resolve())
