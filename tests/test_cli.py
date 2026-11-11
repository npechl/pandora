import gzip
import json
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from pandora.annotations import (
    count_chain_contacts,
    count_residue_contacts,
    drop_labelled_interfaces,
    label_interfaces,
)
from pandora.cli._io import (
    layer_origin,
    read_structures_dir,
    write_structures_dir,
)
from pandora.cli.app import _build_parser, main
from pandora.datasets import extract_interface_records
from pandora.export import read_records
from pandora.ingestion import read_entry_paths
from pandora.parsing import mmcif_to_structure
from pandora.schemas.annotation import LayerRow, ResidueContactRow
from pandora.schemas.provenance import CanonicalisationEntryRecord
from pandora.schemas.similarity import SimilarityCluster
from pandora.similarity import interface_similarity

FIXTURES_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"
CANON_POLICY = (
    Path(__file__).parent.parent / "datasets" / "canonicalisation.yaml"
)

_SUBCOMMANDS = [
    "fetch",
    "representatives",
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


PDBE_TEMPLATE = "{mid}/{id}/clean_mmcif/{id}_updated.cif.gz"


def _pdbe_tree(root):
    """A PDBe-release-like tree with two gzipped entries and strays."""

    for mid, entry in [("04", "104m"), ("0k", "10kz")]:
        folder = root / mid / entry / "clean_mmcif"
        folder.mkdir(parents=True)
        data = (FIXTURES_DIR / mid / f"{entry}.cif").read_bytes()
        (folder / f"{entry}_updated.cif.gz").write_bytes(gzip.compress(data))
        (folder.parent / "other.cif").write_bytes(data)
    (root / "README").write_text("")


def _ingest(*args):
    main(["ingest", *args])


def test_ingest_path_template_feeds_canonicalise(tmp_path, monkeypatch):
    _pdbe_tree(tmp_path / "entry")
    monkeypatch.chdir(tmp_path)
    _ingest(
        "--input-dir",
        "entry",
        "--path-template",
        PDBE_TEMPLATE,
        "--source-label",
        "lbl",
        "--output-dir",
        "ing",
    )
    paths = read_entry_paths(tmp_path / "ing" / "entries.txt")
    assert len(paths) == 2 and all(p.is_absolute() for p in paths)
    summary = json.loads((tmp_path / "ing" / "ingestion.json").read_text())
    assert summary["step"] == "ingestion"
    assert summary["counts"] == {"entries_found": 2}
    assert summary["policy"]["path_template"] == PDBE_TEMPLATE
    assert summary["policy"]["source_label"] == "lbl"

    main(
        [
            "canonicalise",
            "--entries",
            "ing/entries.txt",
            "--policy",
            str(CANON_POLICY),
            "--output-dir",
            "canon",
            "--table-format",
            "jsonl",
        ]
    )
    assert (tmp_path / "canon" / "structures" / "04" / "104m.cif").exists()
    canon = json.loads(
        (tmp_path / "canon" / "canonicalisation.json").read_text()
    )
    assert canon["previous_summary"].endswith("ing/ingestion.json")


def test_ingest_glob(tmp_path):
    _ingest(
        "--input-dir",
        str(FIXTURES_DIR),
        "--glob",
        "*/*.cif",
        "--output-dir",
        str(tmp_path / "ing"),
    )
    paths = read_entry_paths(tmp_path / "ing" / "entries.txt")
    assert paths == sorted(FIXTURES_DIR.resolve().glob("*/*.cif"))


def test_ingest_ids_file(tmp_path):
    ids = tmp_path / "ids.txt"
    ids.write_text("# c\n\n104M\n")
    (tmp_path / "entry").mkdir()
    _ingest(
        "--input-dir",
        str(tmp_path / "entry"),
        "--path-template",
        PDBE_TEMPLATE,
        "--ids",
        str(ids),
        "--output-dir",
        str(tmp_path / "ing"),
    )
    paths = read_entry_paths(tmp_path / "ing" / "entries.txt")
    assert [str(p) for p in paths] == [
        str(tmp_path / "entry/04/104m/clean_mmcif/104m_updated.cif.gz")
    ]


def test_ingest_nothing_found(tmp_path):
    (tmp_path / "entry").mkdir()
    _ingest(
        "--input-dir",
        str(tmp_path / "entry"),
        "--path-template",
        PDBE_TEMPLATE,
        "--output-dir",
        str(tmp_path / "ing"),
    )
    assert (tmp_path / "ing" / "entries.txt").read_text() == ""
    summary = json.loads((tmp_path / "ing" / "ingestion.json").read_text())
    assert summary["counts"] == {"entries_found": 0}


def test_ingest_rejects_ids_without_template(tmp_path):
    ids = tmp_path / "ids.txt"
    ids.write_text("104m\n")
    with pytest.raises(SystemExit):
        _ingest(
            "--input-dir",
            str(tmp_path),
            "--glob",
            "*.cif",
            "--ids",
            str(ids),
            "--output-dir",
            str(tmp_path / "ing"),
        )


def test_canonicalise_overwrite_keeps_entries_txt(tmp_path):
    _pdbe_tree(tmp_path / "entry")
    out = tmp_path / "work"
    _ingest(
        "--input-dir",
        str(tmp_path / "entry"),
        "--path-template",
        PDBE_TEMPLATE,
        "--output-dir",
        str(out),
    )
    for extra in ([], ["--overwrite"]):
        main(
            [
                "canonicalise",
                "--entries",
                str(out / "entries.txt"),
                "--policy",
                str(CANON_POLICY),
                "--output-dir",
                str(out),
                "--table-format",
                "jsonl",
                *extra,
            ]
        )
    assert len(read_entry_paths(out / "entries.txt")) == 2


def _canonicalise_and_curate(tmp_path, fixture="04/104m.cif", canon=None):
    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copy(FIXTURES_DIR / fixture, raw / Path(fixture).name)
    curation_policy = tmp_path / "curation.yaml"
    curation_policy.write_text(
        "policy_id: test\npolicy_name: Test\npolicy_version: 1.0.0\n"
    )
    canonical_dir = tmp_path / "canonical"
    curated_dir = tmp_path / "curated"
    main(
        [
            "canonicalise",
            "--input-dir",
            str(raw),
            "--policy",
            str(canon or CANON_POLICY),
            "--output-dir",
            str(canonical_dir),
            "--table-format",
            "jsonl",
        ]
    )
    main(
        [
            "curate",
            "--input-dir",
            str(canonical_dir),
            "--policy",
            str(curation_policy),
            "--output-dir",
            str(curated_dir),
            "--table-format",
            "jsonl",
        ]
    )
    return canonical_dir, curated_dir


def test_canonicalise_curate_dedup_export_chain(tmp_path):
    canonical_dir, curated_dir = _canonicalise_and_curate(tmp_path)
    assert (curated_dir / "structures" / "04" / "104m.cif").exists()
    chains = _jsonl(curated_dir / "chains.jsonl")
    assert chains and all(c["entry_id"] == "104M" for c in chains)
    assert all("copy_of" in c for c in chains)
    summary = json.loads((curated_dir / "curation.json").read_text())
    assert summary["previous_summary"] == str(
        canonical_dir / "canonicalisation.json"
    )
    assert not (curated_dir / "curation_provenance.json").exists()
    assert not (curated_dir / "curation_exclusions.json").exists()

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
    assert (deduped_dir / "04" / "104m.cif").exists()

    exported = tmp_path / "104m.json"
    main(
        [
            "export",
            "--input",
            str(deduped_dir / "04" / "104m.cif"),
            "--output",
            str(exported),
        ]
    )
    assert json.loads(exported.read_text())["entry_id"] == "104M"


def _similarity_layer(folder, rows, items):
    """A similarity layer folder holding `rows` (None: no matched rows)."""

    import pyarrow as pa
    import pyarrow.parquet as pq

    folder.mkdir(parents=True)
    if rows:
        pq.write_table(pa.Table.from_pylist(rows), folder / "layer.parquet")
    (folder / "layer.json").write_text(
        json.dumps({"policy": {"template": "similarity", "items": items}})
    )
    return folder


def _cluster_inputs(tmp_path, rows=None, items="entry"):
    """A tables-only step folder (entries 104M, 112M, 118L; chains 104M_A,
    112M_A, 112M_B copy of A, 118L_A) and a similarity layer."""

    step = tmp_path / "step"
    step.mkdir()
    _write_jsonl(
        step / "entries.jsonl",
        [
            {"entry_id": e, "structure_path": f"{e}.cif"}
            for e in ("104M", "112M", "118L")
        ],
    )
    _write_jsonl(
        step / "chains.jsonl",
        [
            {"entry_id": "104M", "chain_id": "A"},
            {"entry_id": "112M", "chain_id": "A"},
            {"entry_id": "112M", "chain_id": "B", "copy_of": "A"},
            {"entry_id": "118L", "chain_id": "A"},
        ],
    )
    if rows is None:
        rows = [
            {
                "item_id_1": "104M",
                "item_id_2": "112M",
                "score": 0.95,
                "coverage_1": 0.95,
                "coverage_2": 0.95,
            }
        ]
    layer = _similarity_layer(tmp_path / "layer", rows, items)
    return step, layer


def _cluster(step, layer, output, *extra):
    main(
        ["cluster", "--similarity", str(layer), "--structures-dir", str(step)]
        + ["--output", str(output), *extra]
    )
    return [c["components"] for c in json.loads(output.read_text())]


def test_cluster_with_pairs(tmp_path):
    step, layer = _cluster_inputs(tmp_path)
    pairs = tmp_path / "pairs.json"
    pairs.write_text(json.dumps([["104M", "112M"], ["104M", "118L"]]))
    output = tmp_path / "clusters.json"

    _cluster(step, layer, output, "--min-score", "0.9", "--pairs", str(pairs))

    cluster_pairs = json.loads(
        (output.parent / "cluster_pairs.json").read_text()
    )
    by_pair = {(p["item_id_1"], p["item_id_2"]): p for p in cluster_pairs}
    same = by_pair[("104M", "112M")]
    assert same["cluster_id_1"] == same["cluster_id_2"]
    different = by_pair[("104M", "118L")]
    assert different["cluster_id_1"] != different["cluster_id_2"]


def test_cluster_reads_hit_filter_yaml_and_records_it(tmp_path):
    step, layer = _cluster_inputs(tmp_path)
    hit_filter = tmp_path / "filter.yaml"
    hit_filter.write_text("min_score: 0.99\nmin_coverage: 0.5\n")
    output = tmp_path / "clusters.json"

    clusters = _cluster(step, layer, output, "--hit-filter", str(hit_filter))

    assert clusters == [["104M"], ["112M"], ["118L"]]  # 0.95 < 0.99
    provenance = json.loads(
        (output.parent / "cluster_provenance.json").read_text()
    )
    [recorded] = provenance["layers"]
    assert recorded["hit_filter"]["min_score"] == 0.99
    assert recorded["hit_filter"]["min_coverage"] == 0.5
    assert recorded["similarity_layer"] == str(layer / "layer.json")
    assert recorded["items"] == "entry"
    assert "search" not in provenance


def test_cluster_flag_overrides_hit_filter_yaml(tmp_path):
    step, layer = _cluster_inputs(tmp_path)
    hit_filter = tmp_path / "filter.yaml"
    hit_filter.write_text("min_score: 0.99\n")

    clusters = _cluster(
        step,
        layer,
        tmp_path / "clusters.json",
        "--hit-filter",
        str(hit_filter),
        "--min-score",
        "0.9",
    )

    assert ["104M", "112M"] in clusters


def test_cluster_requires_a_threshold(tmp_path):
    step, layer = _cluster_inputs(tmp_path)

    with pytest.raises(SystemExit, match="--min-score"):
        _cluster(step, layer, tmp_path / "clusters.json")


def _chain_rows(*pairs):
    return [{"item_id_1": a, "item_id_2": b, "score": 0.9} for a, b in pairs]


def test_cluster_chain_layer_with_copies_and_singletons(tmp_path):
    step, layer = _cluster_inputs(
        tmp_path, _chain_rows(("104M_A", "112M_A")), items="chain"
    )

    clusters = _cluster(step, layer, tmp_path / "c.json", "--min-score", "0.5")

    assert clusters == [["104M_A", "112M_A", "112M_B"], ["118L_A"]]


def test_cluster_rolls_chain_layer_up_to_entries(tmp_path):
    step, layer = _cluster_inputs(
        tmp_path, _chain_rows(("112M_B", "118L_A")), items="chain"
    )

    clusters = _cluster(
        step,
        layer,
        tmp_path / "c.json",
        "--min-score",
        "0.5",
        "--items",
        "entry",
    )

    assert clusters == [["104M"], ["112M", "118L"]]


def test_cluster_layer_with_no_matched_rows_gives_singletons(tmp_path):
    step, layer = _cluster_inputs(tmp_path, rows=[], items="entry")

    clusters = _cluster(step, layer, tmp_path / "c.json", "--min-score", "0.5")

    assert clusters == [["104M"], ["112M"], ["118L"]]


def test_cluster_drops_edges_outside_smaller_item_list(tmp_path):
    step, layer = _cluster_inputs(
        tmp_path,
        [{"item_id_1": "104M", "item_id_2": "999X", "score": 0.9}],
    )
    output = tmp_path / "c.json"

    _cluster(step, layer, output, "--min-score", "0.5")

    provenance = json.loads((tmp_path / "cluster_provenance.json").read_text())
    assert provenance["n_edges_unknown_ids"] == 1


def test_cluster_rejects_non_similarity_layer(tmp_path):
    step, _ = _cluster_inputs(tmp_path)
    labels = tmp_path / "labels"
    labels.mkdir()
    (labels / "layer.json").write_text(
        json.dumps({"policy": {"template": "chain_labels"}})
    )
    with pytest.raises(SystemExit, match="similarity"):
        _cluster(step, labels, tmp_path / "c.json", "--min-score", "0.5")


@pytest.mark.parametrize("items", ["interface", "chain"])
def test_cluster_rejects_items_other_than_entry_rollup(tmp_path, items):
    step, layer = _cluster_inputs(tmp_path)  # an entry layer
    with pytest.raises(SystemExit, match="--items"):
        _cluster(
            step,
            layer,
            tmp_path / "c.json",
            "--min-score",
            "0.5",
            "--items",
            items,
        )


def test_python_m_pandora_runs_cli():
    result = subprocess.run(
        [sys.executable, "-m", "pandora", "--help"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert "canonicalise" in result.stdout


def test_read_structures_dir_reads_flat_mirror_and_step_folders(tmp_path):
    src = FIXTURES_DIR / "04" / "104m.cif"
    flat = tmp_path / "flat"
    flat.mkdir()
    shutil.copy(src, flat / "104m.cif")
    mirror = tmp_path / "mirror"
    (mirror / "04").mkdir(parents=True)
    shutil.copy(src, mirror / "04" / "104m.cif")
    step = tmp_path / "step"
    (step / "structures" / "04").mkdir(parents=True)
    shutil.copy(src, step / "structures" / "04" / "104m.cif")
    for folder in (flat, mirror, step):
        assert list(read_structures_dir(folder)) == ["104M"]


def test_write_structures_dir_uses_mirror_layout(tmp_path):
    structure, _, _ = mmcif_to_structure(str(FIXTURES_DIR / "04" / "104m.cif"))
    write_structures_dir({"104M": structure}, tmp_path)
    assert (tmp_path / "04" / "104m.cif").exists()


def _jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def test_canonicalise_writes_step_folder(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copy(FIXTURES_DIR / "04" / "104m.cif", raw / "104m.cif")
    (raw / "bad.cif").write_text("garbage\n")
    out = tmp_path / "canon"
    main(
        [
            "canonicalise",
            "--input-dir",
            str(raw),
            "--policy",
            str(CANON_POLICY),
            "--output-dir",
            str(out),
            "--table-format",
            "jsonl",
        ]
    )
    assert (out / "structures" / "04" / "104m.cif").exists()
    entries = _jsonl(out / "entries.jsonl")
    assert [e["entry_id"] for e in entries] == ["104M"]
    assert entries[0]["structure_path"] == "structures/04/104m.cif"
    mappings = _jsonl(out / "mappings.jsonl")
    assert mappings[0]["entry_id"] == "104M"
    assert mappings[0]["chain_id_mapping"]["items"]
    failures = _jsonl(out / "failures.jsonl")
    assert [Path(f["path"]).name for f in failures] == ["bad.cif"]
    summary = json.loads((out / "canonicalisation.json").read_text())
    assert summary["counts"] == {
        "entries_failed": 1,
        "entries_in": 2,
        "entries_kept": 1,
    }
    assert not (out / "canonicalisation_provenance.json").exists()
    assert not (out / "canonicalisation_mappings.json").exists()


def test_canonicalise_reads_entries_file_and_writes_parquet(tmp_path):
    entries = tmp_path / "entries.txt"
    entries.write_text(str(FIXTURES_DIR / "04" / "104m.cif") + "\n")
    out = tmp_path / "canon"
    main(
        [
            "canonicalise",
            "--entries",
            str(entries),
            "--policy",
            str(CANON_POLICY),
            "--output-dir",
            str(out),
        ]
    )
    rows = read_records(out / "entries.parquet", CanonicalisationEntryRecord)
    assert [r.entry_id for r in rows] == ["104M"]
    summary = json.loads((out / "canonicalisation.json").read_text())
    assert summary["input"] == str(entries)


def test_curate_flat_folder_and_everything_excluded(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copy(FIXTURES_DIR / "04" / "104m.cif", raw / "104m.cif")
    strict = tmp_path / "strict.yaml"
    strict.write_text(
        "policy_id: s\npolicy_name: s\npolicy_version: '1'\n"
        "quality_rules:\n  max_resolution: 1.0\n"
    )
    out = tmp_path / "strict"
    main(
        [
            "curate",
            "--input-dir",
            str(raw),
            "--policy",
            str(strict),
            "--output-dir",
            str(out),
            "--table-format",
            "jsonl",
        ]
    )
    assert not (out / "entries.jsonl").exists()
    assert not (out / "chains.jsonl").exists()
    exclusions = _jsonl(out / "exclusions.jsonl")
    assert exclusions[0]["reason_code"] == "RESOLUTION_THRESHOLD"
    summary = json.loads((out / "curation.json").read_text())
    assert summary["counts"]["entries_kept"] == 0
    assert summary["previous_summary"] is None


def test_curate_fills_copy_of_from_canonicalise_folder(tmp_path):
    # 1ayi's assembly 1 has a symmetry copy A_2 of chain A.
    canon = tmp_path / "canon.yaml"
    canon.write_text(
        "policy_id: c\npolicy_name: c\npolicy_version: '1'\n"
        "assembly_rules:\n  strategy: standardize_biological_assembly\n"
        "  preferred_assembly_source: first\n"
    )
    _, curated_dir = _canonicalise_and_curate(
        tmp_path, "ay/1ayi.cif", canon=canon
    )
    rows = {r["chain_id"]: r for r in _jsonl(curated_dir / "chains.jsonl")}
    assert rows["A-2"]["copy_of"] == "A"
    assert rows["A"]["copy_of"] is None


def test_manifest_reads_step_folders(tmp_path):
    canonical_dir, curated_dir = _canonicalise_and_curate(tmp_path)
    strict = tmp_path / "strict.yaml"
    strict.write_text(
        "policy_id: s\npolicy_name: s\npolicy_version: '1'\n"
        "quality_rules:\n  max_resolution: 1.0\n"
    )
    excluded_dir = tmp_path / "excluded"
    main(
        [
            "curate",
            "--input-dir",
            str(canonical_dir),
            "--policy",
            str(strict),
            "--output-dir",
            str(excluded_dir),
            "--table-format",
            "jsonl",
        ]
    )
    manifest = tmp_path / "manifest.json"
    main(
        [
            "manifest",
            "--dataset-id",
            "d",
            "--dataset-name",
            "d",
            "--dataset-version",
            "1",
            "--structures-dir",
            str(curated_dir),
            "--canonicalisation-dir",
            str(canonical_dir),
            "--curation-dir",
            str(curated_dir),
            "--exclusions",
            str(excluded_dir / "exclusions.jsonl"),
            "--output",
            str(manifest),
        ]
    )
    m = json.loads(manifest.read_text())
    bundle = m["structures"][0]
    assert bundle["entry_id"] == "104M"
    assert bundle["canonicalisation"]["policy_id"] == "overview-remap"
    assert bundle["curation"]["policy_id"] == "test"
    assert [e["reason_code"] for e in m["excluded"]] == ["RESOLUTION_THRESHOLD"]


def test_steps_accept_relative_output_dirs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    canonical_dir, curated_dir = _canonicalise_and_curate(Path("."))
    assert _jsonl(canonical_dir / "entries.jsonl")[0]["structure_path"] == (
        "structures/04/104m.cif"
    )
    assert _jsonl(curated_dir / "entries.jsonl")[0]["structure_path"] == (
        "structures/04/104m.cif"
    )


def _canonicalise(raw, out, *extra):
    main(
        [
            "canonicalise",
            "--input-dir",
            str(raw),
            "--policy",
            str(CANON_POLICY),
            "--output-dir",
            str(out),
            "--table-format",
            "jsonl",
            *extra,
        ]
    )


def test_step_refuses_non_empty_output_dir_unless_overwrite(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copy(FIXTURES_DIR / "04" / "104m.cif", raw / "104m.cif")
    out = tmp_path / "canon"
    _canonicalise(raw, out)
    (raw / "104m.cif").unlink()
    shutil.copy(FIXTURES_DIR / "0k" / "10kz.cif", raw / "10kz.cif")
    with pytest.raises(SystemExit):
        _canonicalise(raw, out)
    (out / "notes.txt").write_text("mine")
    _canonicalise(raw, out, "--overwrite")
    assert not (out / "structures" / "04" / "104m.cif").exists()
    assert [e["entry_id"] for e in _jsonl(out / "entries.jsonl")] == ["10KZ"]
    assert (out / "notes.txt").read_text() == "mine"


def test_step_checks_parquet_support_before_running(tmp_path, monkeypatch):
    raw = tmp_path / "raw"
    raw.mkdir()
    shutil.copy(FIXTURES_DIR / "04" / "104m.cif", raw / "104m.cif")
    monkeypatch.setitem(sys.modules, "pandas", None)
    out = tmp_path / "canon"
    with pytest.raises(SystemExit) as exc:
        main(
            [
                "canonicalise",
                "--input-dir",
                str(raw),
                "--policy",
                str(CANON_POLICY),
                "--output-dir",
                str(out),
            ]
        )
    assert "--table-format jsonl" in str(exc.value)
    assert not (out / "structures").exists()


def test_step_skips_duplicate_entries(tmp_path):
    raw = tmp_path / "raw"
    (raw / "04").mkdir(parents=True)
    shutil.copy(FIXTURES_DIR / "04" / "104m.cif", raw / "104m.cif")
    shutil.copy(FIXTURES_DIR / "04" / "104m.cif", raw / "04" / "104m.cif")
    out = tmp_path / "canon"
    _canonicalise(raw, out)
    assert [e["entry_id"] for e in _jsonl(out / "entries.jsonl")] == ["104M"]
    failures = _jsonl(out / "failures.jsonl")
    assert [f["error_type"] for f in failures] == ["DuplicateEntry"]


def test_step_names_outputs_by_entry_id(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    data = (FIXTURES_DIR / "04" / "104m.cif").read_bytes()
    (raw / "104m_updated.cif.gz").write_bytes(gzip.compress(data))
    (raw / "entries.txt").write_text("104m_updated.cif.gz\n")
    canon = tmp_path / "canon"
    main(
        [
            "canonicalise",
            "--entries",
            str(raw / "entries.txt"),
            "--policy",
            str(CANON_POLICY),
            "--output-dir",
            str(canon),
            "--table-format",
            "jsonl",
        ]
    )
    assert (canon / "structures" / "04" / "104m.cif").exists()
    assert [e["entry_id"] for e in _jsonl(canon / "entries.jsonl")] == ["104M"]


def test_manifest_rejects_folder_that_is_not_a_step_folder(tmp_path):
    canonical_dir, curated_dir = _canonicalise_and_curate(tmp_path)
    with pytest.raises(SystemExit):
        main(
            [
                "manifest",
                "--dataset-id",
                "d",
                "--dataset-name",
                "d",
                "--dataset-version",
                "1",
                "--structures-dir",
                str(curated_dir),
                "--canonicalisation-dir",
                str(tmp_path / "typo"),
                "--output",
                str(tmp_path / "m.json"),
            ]
        )


def test_fetch_writes_entries_txt(tmp_path, monkeypatch):
    data = (FIXTURES_DIR / "04" / "104m.cif").read_bytes()
    monkeypatch.setattr("pandora.ingestion.mmcif._download", lambda url: data)
    raw = tmp_path / "raw"
    main(
        ["fetch", "104m", "--provider", "pdb", "--no-cache"]
        + ["--output-dir", str(raw)]
    )
    assert read_entry_paths(raw / "entries.txt") == [
        (raw / "104m.cif").resolve()
    ]


def test_ingest_rejects_missing_input_dir(tmp_path):
    for find in (["--glob", "*.cif"], ["--path-template", "{id}.cif"]):
        with pytest.raises(SystemExit, match="not a folder"):
            _ingest(
                "--input-dir",
                str(tmp_path / "missing"),
                *find,
                "--output-dir",
                str(tmp_path / "ing"),
            )


def test_curate_duplicate_check_ignores_case(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    text = (FIXTURES_DIR / "04" / "104m.cif").read_text()
    (raw / "a.cif").write_text(text)
    (raw / "b.cif").write_text(text.replace("data_104M", "data_104m", 1))
    policy = tmp_path / "curate.yaml"
    policy.write_text("policy_id: c\npolicy_name: c\npolicy_version: '1'\n")
    out = tmp_path / "curated"
    main(
        ["curate", "--input-dir", str(raw), "--policy", str(policy)]
        + ["--output-dir", str(out), "--table-format", "jsonl"]
    )
    failures = _jsonl(out / "failures.jsonl")
    assert [f["error_type"] for f in failures] == ["DuplicateEntry"]


def _step_folder(tmp_path):
    """A step folder holding 1a2c; returns its path."""

    structures = tmp_path / "curated" / "structures" / "a2"
    structures.mkdir(parents=True)
    shutil.copy(FIXTURES_DIR / "a2" / "1a2c.cif", structures / "1a2c.cif")
    return tmp_path / "curated"


def _compute(folder, out, *extra):
    main(
        ["annotate", "compute", "--structures-dir", str(folder)]
        + ["--output-dir", str(out), "--table-format", "jsonl", *extra]
    )


def test_annotate_compute_writes_contact_layer(tmp_path):
    out = tmp_path / "annotations" / "interfaces"
    _compute(
        _step_folder(tmp_path),
        out,
        "--layer",
        "chain_interfaces",
        "--atom-set",
        "heavy",
    )

    rows = read_records(out / "layer.jsonl", ResidueContactRow)
    assert len(rows) == 74 + 21 + 25
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "annotation_compute"
    assert summary["policy"]["template"] == "residue_contacts"
    assert summary["policy"]["tool"] == "pandora"
    assert summary["policy"]["parameters"] == {
        "distance_cutoff": 4.0,
        "atom_set": "heavy",
        "polymer_types": None,
    }
    assert summary["counts"] == {"entries": 1, "rows": 120}


def test_annotate_compute_ligand_layer_settings(tmp_path):
    out = tmp_path / "ligands"
    _compute(
        _step_folder(tmp_path),
        out,
        "--layer",
        "ligand_contacts",
        "--distance-cutoff",
        "5.0",
        "--include-waters",
    )

    summary = json.loads((out / "layer.json").read_text())
    assert summary["policy"]["parameters"] == {
        "distance_cutoff": 5.0,
        "include_waters": True,
    }
    assert read_records(out / "layer.jsonl", ResidueContactRow)


def test_annotate_compute_rejects_flag_for_other_layer(tmp_path):
    out = tmp_path / "out"
    with pytest.raises(SystemExit, match="--include-waters"):
        _compute(
            _step_folder(tmp_path),
            out,
            "--layer",
            "chain_interfaces",
            "--include-waters",
        )
    assert not out.exists()


def test_annotate_compute_refuses_existing_layer(tmp_path):
    folder = _step_folder(tmp_path)
    out = tmp_path / "out"
    _compute(folder, out, "--layer", "ligand_contacts")
    with pytest.raises(SystemExit):
        _compute(folder, out, "--layer", "ligand_contacts")
    _compute(folder, out, "--layer", "ligand_contacts", "--overwrite")


def test_compute_and_import_give_same_interface_records(tmp_path):
    folder = _step_folder(tmp_path)
    computed = tmp_path / "computed"
    _compute(folder, computed, "--layer", "chain_interfaces")
    # The same contacts as a tool would give them: author naming only,
    # lowercase entry id, sides swapped.
    table = tmp_path / "tool.jsonl"
    table.write_text(
        "".join(
            json.dumps(
                {
                    "entry_id": r.entry_id.lower(),
                    "auth_chain_id_1": r.auth_chain_id_2,
                    "auth_seq_id_1": r.auth_seq_id_2,
                    "ins_code_1": r.ins_code_2,
                    "auth_chain_id_2": r.auth_chain_id_1,
                    "auth_seq_id_2": r.auth_seq_id_1,
                    "ins_code_2": r.ins_code_1,
                }
            )
            + "\n"
            for r in read_records(computed / "layer.jsonl", ResidueContactRow)
        )
    )
    (tmp_path / "d.yaml").write_text(
        "name: tool\ntemplate: residue_contacts\ntool: T\n"
    )
    imported = tmp_path / "imported"
    main(
        ["annotate", "import", "--structures-dir", str(folder)]
        + ["--table", str(table), "--description", str(tmp_path / "d.yaml")]
        + ["--output-dir", str(imported), "--table-format", "jsonl"]
    )

    a = extract_interface_records(
        read_records(computed / "layer.jsonl", ResidueContactRow)
    )
    b = extract_interface_records(
        read_records(imported / "layer.jsonl", ResidueContactRow)
    )
    assert a == b
    assert len(a) == 3


class _Captured(Exception):
    pass


def test_cluster_reads_contact_layer(tmp_path, monkeypatch):
    folder = _step_folder(tmp_path)
    _write_jsonl(
        folder / "chains.jsonl",
        [{"entry_id": "1A2C", "chain_id": c} for c in "ABCD"],
    )
    contacts = tmp_path / "interfaces"
    _compute(folder, contacts, "--layer", "chain_interfaces")
    row = {
        "item_id_1": "1A2C_A",
        "item_id_2": "1A2C_B",
        "score": 0.9,
        "start_1": 1,
        "end_1": 50,
        "start_2": 1,
        "end_2": 50,
    }
    layer = _similarity_layer(tmp_path / "struct", [row], "chain")
    captured = {}

    def fake_edges(batches, hit_filter, interface_residues=None):
        captured["residues"] = interface_residues
        raise _Captured

    monkeypatch.setattr("pandora.cli.similarity.similarity_edges", fake_edges)
    hit_filter = tmp_path / "filter.yaml"
    hit_filter.write_text("min_interface_coverage: 0.5\n")
    with pytest.raises(_Captured):
        main(
            ["cluster", "--similarity", str(layer)]
            + ["--hit-filter", str(hit_filter)]
            + ["--contacts", str(contacts), "--structures-dir", str(folder)]
            + ["--output", str(tmp_path / "c.json")]
        )

    residues = captured["residues"]
    assert set(residues) == {"1A2C_A", "1A2C_B", "1A2C_C", "1A2C_D"}
    assert all(residues.values())


def test_cluster_interface_layer_takes_items_from_contacts(tmp_path):
    contacts = tmp_path / "contacts"
    contacts.mkdir()
    _write_jsonl(
        contacts / "layer.jsonl",
        [
            {"entry_id": "1A2C", "chain_id_1": a, "chain_id_2": b}
            | {"seq_id_1": 1, "seq_id_2": 2}
            for a, b in (("A", "B"), ("B", "C"), ("B", "D"))
        ],
    )
    (contacts / "layer.json").write_text(
        json.dumps({"policy": {"template": "residue_contacts"}})
    )
    layer = _similarity_layer(
        tmp_path / "ifaces",
        [{"item_id_1": "1A2C_A_B", "item_id_2": "1A2C_B_C", "score": 0.9}],
        "interface",
    )
    output = tmp_path / "c.json"
    main(
        ["cluster", "--similarity", str(layer), "--min-score", "0.5"]
        + ["--contacts", str(contacts), "--output", str(output)]
    )

    clusters = [c["components"] for c in json.loads(output.read_text())]
    assert clusters == [["1A2C_A_B", "1A2C_B_C"], ["1A2C_B_D"]]


def _rollup_contacts(tmp_path):
    """A contact layer with interfaces 1AAA_A_B, 2BBB_C_D and 3CCC_E_F."""

    contacts = tmp_path / "contacts"
    contacts.mkdir()
    _write_jsonl(
        contacts / "layer.jsonl",
        [
            {"entry_id": e, "chain_id_1": a, "chain_id_2": b}
            | {"seq_id_1": 1, "seq_id_2": 2}
            for e, a, b in (("1AAA", "A", "B"), ("2BBB", "C", "D"))
            + (("3CCC", "E", "F"),)
        ],
    )
    (contacts / "layer.json").write_text(
        json.dumps({"policy": {"template": "residue_contacts"}})
    )
    return contacts


def _rollup_cluster(tmp_path, contacts, *extra):
    output = tmp_path / "c.json"
    main(
        ["cluster", "--items", "interface", "--contacts", str(contacts)]
        + ["--output", str(output), *extra]
    )
    provenance = json.loads((tmp_path / "cluster_provenance.json").read_text())
    return [c["components"] for c in json.loads(output.read_text())], provenance


def test_cluster_rolls_chain_layer_up_to_interfaces(tmp_path):
    contacts = _rollup_contacts(tmp_path)
    tm = _similarity_layer(
        tmp_path / "tm",
        _chain_rows(("1AAA_A", "2BBB_D"), ("9ZZZ_A", "3CCC_E"))
        + _chain_rows(("9ZZZ_A", "1AAA_A")),
        "chain",
    )

    clusters, provenance = _rollup_cluster(
        tmp_path, contacts, "--similarity", str(tm), "--min-score", "0.5"
    )

    # 9ZZZ_A is in no interface, so it does not join 3CCC_E_F.
    assert clusters == [["1AAA_A_B", "2BBB_C_D"], ["3CCC_E_F"]]
    assert provenance["layers"][0]["items"] == "chain"
    assert provenance["layers"][0]["n_edges"] == 3


def _chain_reps(path, mapping):
    """A representatives.json mapping each chain to its representative."""

    groups: dict[str, list[str]] = {}
    for chain, rep in mapping.items():
        groups.setdefault(rep, []).append(chain)
    path.write_text(
        json.dumps(
            [
                {"representative": r, "members": m, "n_members": len(m)}
                for r, m in groups.items()
            ]
        )
    )
    return path


_ROLLUP_CHAINS = ["1AAA_A", "1AAA_B", "2BBB_C", "2BBB_D", "3CCC_E", "3CCC_F"]


def test_cluster_chain_representatives_joins_members(tmp_path):
    contacts = _rollup_contacts(tmp_path)
    mapping = dict(zip(_ROLLUP_CHAINS, _ROLLUP_CHAINS, strict=True))
    mapping |= {"1AAA_A": "9R_X", "2BBB_D": "8R_Y"}
    reps = _chain_reps(tmp_path / "reps.json", mapping)
    tm = _similarity_layer(
        tmp_path / "tm", _chain_rows(("9R_X", "8R_Y")), "chain"
    )

    clusters, provenance = _rollup_cluster(
        tmp_path,
        contacts,
        *["--similarity", str(tm), "--min-score", "0.5"],
        *["--chain-representatives", str(reps)],
    )

    assert clusters == [["1AAA_A_B", "2BBB_C_D"], ["3CCC_E_F"]]
    assert provenance["chain_representatives"] == str(reps)


def test_cluster_chain_representatives_without_chain_layer(tmp_path):
    # Dimers sharing a representative link even with no chain layer.
    contacts = _rollup_contacts(tmp_path)
    mapping = dict(zip(_ROLLUP_CHAINS, _ROLLUP_CHAINS, strict=True))
    mapping["2BBB_C"] = "1AAA_A"
    reps = _chain_reps(tmp_path / "reps.json", mapping)
    iface = _similarity_layer(tmp_path / "iface", None, "interface")

    clusters, _ = _rollup_cluster(
        tmp_path,
        contacts,
        *["--similarity", str(iface), "--min-score", "0.5"],
        *["--chain-representatives", str(reps)],
    )

    assert clusters == [["1AAA_A_B", "2BBB_C_D"], ["3CCC_E_F"]]


def test_cluster_chain_representatives_missing_chain_exits(tmp_path):
    contacts = _rollup_contacts(tmp_path)
    mapping = {c: c for c in _ROLLUP_CHAINS if c != "3CCC_F"}
    reps = _chain_reps(tmp_path / "reps.json", mapping)
    tm = _similarity_layer(tmp_path / "tm", _chain_rows(), "chain")
    with pytest.raises(SystemExit, match="3CCC_F"):
        _rollup_cluster(
            tmp_path,
            contacts,
            *["--similarity", str(tm), "--min-score", "0.5"],
            *["--chain-representatives", str(reps)],
        )


def test_cluster_chain_representatives_needs_interface_items(tmp_path):
    step, layer = _cluster_inputs(
        tmp_path, _chain_rows(("104M_A", "112M_A")), items="chain"
    )
    reps = _chain_reps(tmp_path / "reps.json", {"104M_A": "104M_A"})
    with pytest.raises(SystemExit, match="--items interface"):
        _cluster(
            step,
            layer,
            tmp_path / "c.json",
            *["--min-score", "0.5", "--chain-representatives", str(reps)],
        )


def test_cluster_joins_several_layers(tmp_path):
    contacts = _rollup_contacts(tmp_path)
    tm = _similarity_layer(
        tmp_path / "tm", _chain_rows(("1AAA_A", "2BBB_D")), "chain"
    )
    rows = [
        {"item_id_1": "2BBB_C_D", "item_id_2": "3CCC_E_F", "score": 0.3},
    ]
    iface = _similarity_layer(tmp_path / "iface", rows, "interface")
    args = ["--similarity", str(tm), "--min-score", "0.5"]
    args += ["--similarity", str(iface), "--min-score", "0.2"]

    clusters, provenance = _rollup_cluster(tmp_path, contacts, *args)

    assert clusters == [["1AAA_A_B", "2BBB_C_D", "3CCC_E_F"]]
    assert [
        (layer["hit_filter"]["min_score"], layer["n_edges"])
        for layer in provenance["layers"]
    ] == [(0.5, 1), (0.2, 1)]

    # One --min-score applies to every layer: 0.3 < 0.5 drops the link.
    clusters, _ = _rollup_cluster(
        tmp_path,
        contacts,
        *["--similarity", str(tm), "--similarity", str(iface)],
        *["--min-score", "0.5"],
    )
    assert clusters == [["1AAA_A_B", "2BBB_C_D"], ["3CCC_E_F"]]


def test_cluster_threshold_count_must_match_layers(tmp_path):
    contacts = _rollup_contacts(tmp_path)
    tm = _similarity_layer(tmp_path / "tm", _chain_rows(), "chain")
    with pytest.raises(SystemExit, match="once per --similarity"):
        _rollup_cluster(
            tmp_path,
            contacts,
            *["--similarity", str(tm), "--similarity", str(tm)],
            *["--min-score", "0.5", "--min-score", "0.5"],
            *["--min-score", "0.5"],
        )


def test_cluster_rejects_non_contact_layer(tmp_path):
    step, layer = _cluster_inputs(tmp_path, _chain_rows(), items="chain")
    labels = tmp_path / "labels"
    labels.mkdir()
    (labels / "layer.json").write_text(
        json.dumps({"policy": {"template": "chain_labels"}})
    )
    with pytest.raises(SystemExit, match="residue_contacts"):
        _cluster(
            step,
            layer,
            tmp_path / "c.json",
            "--min-score",
            "0.5",
            "--contacts",
            str(labels),
        )


def test_cluster_entry_layer_needs_structures_dir(tmp_path):
    _, layer = _cluster_inputs(tmp_path)
    with pytest.raises(SystemExit, match="--structures-dir"):
        main(
            ["cluster", "--similarity", str(layer), "--min-score", "0.5"]
            + ["--output", str(tmp_path / "c.json")]
        )


def test_annotate_old_flat_form_is_gone(tmp_path):
    with pytest.raises(SystemExit):
        main(
            [
                "annotate",
                "--input-dir",
                str(tmp_path),
                "--output-dir",
                str(tmp_path),
            ]
        )


_LAYER_ROWS = [
    '{"entry_id":"1a2c","auth_chain_id":"L","auth_seq_id":"1",'
    '"ins_code":"H","value":1.0}',
    '{"entry_id":"1A2C","chain_id":"A","seq_id":2,"value":2.0}',
    '{"entry_id":"1A2C","auth_chain_id":"H","auth_seq_id":"147",'
    '"ins_code":"A","value":3.0}',
    '{"entry_id":"9XYZ","chain_id":"A","seq_id":1,"value":4.0}',
]


def _layer_inputs(tmp_path):
    """A curated step folder with 1a2c, a 4-row residue_values table, a
    description and a raw tool output; returns the import CLI args."""

    structures = tmp_path / "curated" / "structures" / "a2"
    structures.mkdir(parents=True)
    shutil.copy(FIXTURES_DIR / "a2" / "1a2c.cif", structures / "1a2c.cif")
    (tmp_path / "curated" / "curation.json").write_text("{}")
    (tmp_path / "cocada.yaml").write_text(
        "name: t\ntemplate: residue_values\ntool: T\n"
        'tool_version: "1"\nparameters: {cutoff: 5.0}\n'
    )
    (tmp_path / "raw").mkdir()
    (tmp_path / "raw" / "out.txt").write_text("raw tool output\n")
    (tmp_path / "t.jsonl").write_text("\n".join(_LAYER_ROWS) + "\n")
    return [
        "--structures-dir",
        str(tmp_path / "curated"),
        "--table",
        str(tmp_path / "t.jsonl"),
        "--description",
        str(tmp_path / "cocada.yaml"),
        "--source-files",
        str(tmp_path / "raw"),
        "--output-dir",
        str(tmp_path / "annotations" / "t"),
        "--table-format",
        "jsonl",
    ]


def test_annotate_import_writes_layer_folder(tmp_path):
    main(["annotate", "import", *_layer_inputs(tmp_path)])
    out = tmp_path / "annotations" / "t"
    matched = read_records(out / "layer.jsonl", LayerRow)
    assert [r.model_dump()["seq_id"] for r in matched] == [1, 2]
    unmatched = [
        r.model_dump()["reason"]
        for r in read_records(out / "unmatched.jsonl", LayerRow)
    ]
    assert unmatched == ["no_coordinates", "entry_not_found"]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "annotation_import"
    assert summary["policy"]["template"] == "residue_values"
    assert summary["counts"] == {
        "rows_in": 4,
        "rows_matched": 2,
        "rows_unmatched": 2,
        "unmatched:entry_not_found": 1,
        "unmatched:no_coordinates": 1,
    }
    assert summary["previous_summary"].endswith("curation.json")
    assert (out / "source" / "out.txt").exists()


def test_annotate_import_refuses_existing_layer_without_overwrite(tmp_path):
    args = _layer_inputs(tmp_path)
    main(["annotate", "import", *args])
    with pytest.raises(SystemExit):
        main(["annotate", "import", *args])
    main(["annotate", "import", *args, "--overwrite"])


def test_annotate_import_missing_value_column_writes_nothing(tmp_path):
    args = _layer_inputs(tmp_path)
    (tmp_path / "t.jsonl").write_text(
        '{"entry_id":"1A2C","chain_id":"A","seq_id":1}\n'
    )
    with pytest.raises(SystemExit, match="value"):
        main(["annotate", "import", *args])
    assert not (tmp_path / "annotations" / "t").exists()


def test_annotate_import_empty_table(tmp_path):
    args = _layer_inputs(tmp_path)
    (tmp_path / "t.jsonl").write_text("")
    main(["annotate", "import", *args])
    out = tmp_path / "annotations" / "t"
    summary = json.loads((out / "layer.json").read_text())
    assert summary["counts"]["rows_in"] == 0
    assert not (out / "layer.jsonl").exists()


def test_annotate_import_parquet_with_mixed_type_columns(tmp_path):
    args = _layer_inputs(tmp_path)[:-2]  # default --table-format parquet
    (tmp_path / "t.jsonl").write_text(
        '{"entry_id":"1A2C","chain_id":"A","seq_id":1,"value":1,"note":1}\n'
        '{"entry_id":"1A2C","chain_id":"A","seq_id":"x","value":1,'
        '"note":"a"}\n'
        '{"entry_id":"1A2C","chain_id":"A","seq_id":999,"value":"n/a"}\n'
    )
    main(["annotate", "import", *args])
    out = tmp_path / "annotations" / "t"
    reasons = [
        r.model_dump()["reason"]
        for r in read_records(out / "unmatched.parquet", LayerRow)
    ]
    assert reasons == ["bad_value", "residue_not_found"]
    assert (out / "layer.json").exists()


def test_annotate_import_refuses_similarity_layer(tmp_path):
    args = _layer_inputs(tmp_path)
    (tmp_path / "cocada.yaml").write_text(
        "name: s\ntemplate: similarity\nitems: chain\ntool: R\n"
    )
    (tmp_path / "t.jsonl").write_text("")
    with pytest.raises(SystemExit, match="pandora similarity import"):
        main(["annotate", "import", *args])
    assert not (tmp_path / "annotations" / "t").exists()


# Stands in for `mmseqs`: keeps the query FASTA, writes two hits.
FAKE_MMSEQS_CLI = """#!/bin/sh
if [ "$1" = version ]; then echo 15.6f452; exit 0; fi
cp "$2" "$(dirname "$0")/query.fasta"
printf '104M_A\\t112M_A\\t0.99\\t100\\t0.95\\t0.9\\n' > "$4"
printf '104M_A\\t999X_A\\t0.5\\t100\\t0.95\\t0.9\\n' >> "$4"
"""

# Stands in for `foldseek`: one hit between two chain files.
FAKE_FOLDSEEK_CLI = """#!/bin/sh
if [ "$1" = version ]; then echo 10.941cd33; exit 0; fi
printf '104m_A\\t112m_A\\t0.5\\t100\\t0.9\\t0.8\\t0.7\\t0.6\\t0.6\\t1\\t100\\t1\\t90\\n' > "$4"
"""


def _write_jsonl(path, rows):
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))


def _similarity_step(tmp_path):
    """A step folder holding only tables: entries 104M, 112M; chains
    104M_A, 112M_A and 112M_B (a symmetry copy of 112M_A)."""

    folder = tmp_path / "curated"
    folder.mkdir()
    _write_jsonl(
        folder / "entries.jsonl",
        [
            {"entry_id": "104M", "structure_path": "s/104m.cif"},
            {"entry_id": "112M", "structure_path": "s/112m.cif"},
            {"entry_id": "118L", "structure_path": None},
        ],
    )
    _write_jsonl(
        folder / "chains.jsonl",
        [
            {"entry_id": "104M", "chain_id": "A", "sequence": "MKV"},
            {"entry_id": "112M", "chain_id": "A", "sequence": "MKL"},
            {
                "entry_id": "112M",
                "chain_id": "B",
                "sequence": "MKL",
                "copy_of": "A",
            },
        ],
    )
    (folder / "curation.json").write_text("{}")
    return folder


def _executable(path, text):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def _search_mmseqs(tmp_path, out, *extra):
    mmseqs = _executable(tmp_path / "mmseqs", FAKE_MMSEQS_CLI)
    main(
        [
            "similarity",
            "search",
            "--engine",
            "mmseqs2",
            "--mmseqs-bin",
            str(mmseqs),
            "--structures-dir",
            str(_similarity_step(tmp_path)),
            "--output-dir",
            str(out),
            *extra,
        ]
    )


def _parquet_rows(path):
    import pyarrow.parquet as pq

    return pq.read_table(path).to_pylist()


def test_similarity_search_mmseqs_writes_layer_without_copies(tmp_path):
    out = tmp_path / "similarity" / "seq"
    _search_mmseqs(tmp_path, out, "--max-seqs", "500")

    query = (tmp_path / "query.fasta").read_text()
    assert ">104M_A" in query
    assert ">112M_B" not in query
    rows = _parquet_rows(out / "layer.parquet")
    assert [(r["item_id_1"], r["item_id_2"]) for r in rows] == [
        ("104M_A", "112M_A")
    ]
    assert rows[0]["score"] == rows[0]["identity"] == 0.99
    assert rows[0]["coverage_2"] == 0.9
    unmatched = _parquet_rows(out / "unmatched.parquet")
    assert [r["reason"] for r in unmatched] == ["entry_not_found"]
    assert not (out / "hits.tsv").exists()
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "similarity_search"
    assert summary["policy"]["items"] == "chain"
    assert summary["policy"]["tool"] == "MMseqs2"
    assert summary["policy"]["parameters"]["max_seqs"] == 500
    assert summary["previous_summary"].endswith("curation.json")
    assert summary["counts"] == {
        "chains_searched": 2,
        "rows_in": 2,
        "rows_matched": 1,
        "rows_unmatched": 1,
        "unmatched:entry_not_found": 1,
    }


def test_similarity_search_refuses_existing_layer(tmp_path):
    out = tmp_path / "similarity" / "seq"
    _search_mmseqs(tmp_path, out)
    shutil.rmtree(tmp_path / "curated")
    with pytest.raises(SystemExit, match="overwrite"):
        _search_mmseqs(tmp_path, out)


def test_similarity_search_foldseek_matches_file_stems(tmp_path):
    foldseek = _executable(tmp_path / "foldseek", FAKE_FOLDSEEK_CLI)
    chains = tmp_path / "chains"
    chains.mkdir()
    out = tmp_path / "similarity" / "struct"
    main(
        [
            "similarity",
            "search",
            "--engine",
            "foldseek",
            "--foldseek-bin",
            str(foldseek),
            "--input-dir",
            str(chains),
            "--structures-dir",
            str(_similarity_step(tmp_path)),
            "--output-dir",
            str(out),
        ]
    )

    rows = _parquet_rows(out / "layer.parquet")
    assert [(r["item_id_1"], r["item_id_2"]) for r in rows] == [
        ("104M_A", "112M_A")
    ]
    assert rows[0]["score"] == 0.7
    assert rows[0]["end_2"] == 90


def test_similarity_search_mmseqs_rejects_input_dir(tmp_path):
    with pytest.raises(SystemExit, match="--input-dir"):
        _search_mmseqs(tmp_path, tmp_path / "out", "--input-dir", str(tmp_path))


def test_similarity_search_mmseqs_limits_to_chains(tmp_path):
    chains = tmp_path / "reps.txt"
    chains.write_text("# reps\n104M_A \n\n")
    out = tmp_path / "seq"
    _search_mmseqs(tmp_path, out, "--chains", str(chains))

    query = (tmp_path / "query.fasta").read_text()
    assert ">104M_A" in query
    assert ">112M_A" not in query
    summary = json.loads((out / "layer.json").read_text())
    assert summary["counts"]["chains_searched"] == 1
    assert summary["input"] == str(chains)


def test_similarity_search_chains_unknown_id_exits(tmp_path):
    chains = tmp_path / "reps.txt"
    chains.write_text("104M_A\n9XXX_A\n")
    with pytest.raises(SystemExit, match="9XXX_A"):
        _search_mmseqs(tmp_path, tmp_path / "out", "--chains", str(chains))


def _search_foldseek(tmp_path, structures_dir, out, *extra):
    main(
        ["similarity", "search", "--engine", "foldseek"]
        + ["--foldseek-bin", str(tmp_path / "foldseek")]
        + ["--structures-dir", str(structures_dir)]
        + ["--output-dir", str(out), *extra]
    )


def test_similarity_search_chains_rejects_input_dir(tmp_path):
    _executable(tmp_path / "foldseek", FAKE_FOLDSEEK_CLI)
    chains = tmp_path / "reps.txt"
    chains.write_text("104M_A\n")
    with pytest.raises(SystemExit, match="--chains"):
        _search_foldseek(
            tmp_path,
            _similarity_step(tmp_path),
            tmp_path / "out",
            *["--chains", str(chains), "--input-dir", str(tmp_path)],
        )


# Stands in for `foldseek`: lists the files it was given, no hits.
FAKE_FOLDSEEK_LIST = """#!/bin/sh
if [ "$1" = version ]; then echo 10.941cd33; exit 0; fi
ls "$2" > "$(dirname "$0")/files.txt"
for f in "$2"/*; do grep -c '^ATOM\\|^HETATM' "$f"; done \
  > "$(dirname "$0")/atoms.txt"
grep -h '^ATOM\\|^HETATM' "$2"/* | awk '{print $7}' | sort -u \
  > "$(dirname "$0")/asyms.txt"
: > "$4"
"""


def test_similarity_search_foldseek_writes_chain_files(tmp_path):
    _, curated = _canonicalise_and_curate(tmp_path, "a2/1a2c.cif")
    chains = _jsonl(curated / "chains.jsonl")
    ids = sorted(
        f"{c['entry_id']}_{c['chain_id']}" for c in chains if not c["copy_of"]
    )
    assert len(ids) > 1
    _executable(tmp_path / "foldseek", FAKE_FOLDSEEK_LIST)
    _search_foldseek(tmp_path, curated, tmp_path / "tm")

    files = (tmp_path / "files.txt").read_text().split()
    assert files == [f"{i}.cif" for i in ids]
    assert all(int(n) > 0 for n in (tmp_path / "atoms.txt").read_text().split())

    # With --chains, only the listed chain is written.
    listed = tmp_path / "reps.txt"
    listed.write_text(ids[0] + "\n")
    _search_foldseek(
        tmp_path, curated, tmp_path / "tm1", "--chains", str(listed)
    )
    assert (tmp_path / "files.txt").read_text().split() == [f"{ids[0]}.cif"]
    one_chain = ids[0].split("_", 1)[1]
    assert (tmp_path / "asyms.txt").read_text().split() == [one_chain]


def _import(tmp_path, items, *flags, rows=None):
    (tmp_path / "graph.yaml").write_text(
        f"name: r_graph\ntemplate: similarity\n{items}tool: R\n"
    )
    table = tmp_path / "edges.tsv"
    table.write_text(
        "item_id_1\titem_id_2\tscore\tweight\n"
        + "".join("\t".join(r) + "\n" for r in rows or [])
    )
    out = tmp_path / "similarity" / "r_graph"
    main(
        [
            "similarity",
            "import",
            "--table",
            str(table),
            "--description",
            str(tmp_path / "graph.yaml"),
            "--output-dir",
            str(out),
            *flags,
        ]
    )
    return out


def test_similarity_import_tsv_entry_level(tmp_path):
    out = _import(
        tmp_path,
        "items: entry\n",
        "--structures-dir",
        str(_similarity_step(tmp_path)),
        rows=[("104m", "112M", "0.8", "w"), ("104M", "999X", "0.5", "w")],
    )

    rows = _parquet_rows(out / "layer.parquet")
    assert rows == [
        {"item_id_1": "104M", "item_id_2": "112M", "score": 0.8, "weight": "w"}
    ]
    assert [r["reason"] for r in _parquet_rows(out / "unmatched.parquet")] == [
        "entry_not_found"
    ]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "similarity_import"
    assert summary["input"].endswith("edges.tsv")


def test_similarity_import_118l_without_structure_is_not_an_item(tmp_path):
    out = _import(
        tmp_path,
        "items: entry\n",
        "--structures-dir",
        str(_similarity_step(tmp_path)),
        rows=[("104M", "118L", "0.8", "w")],
    )

    assert not (out / "layer.parquet").exists()
    assert [r["reason"] for r in _parquet_rows(out / "unmatched.parquet")] == [
        "entry_not_found"
    ]


def test_similarity_import_interface_level_from_contacts(tmp_path):
    contacts = tmp_path / "contacts"
    contacts.mkdir()
    _write_jsonl(
        contacts / "layer.jsonl",
        [
            {
                "entry_id": "1A2C",
                "chain_id_1": "B",
                "chain_id_2": "A",
                "seq_id_1": 3,
                "seq_id_2": 9,
            }
        ],
    )
    (contacts / "layer.json").write_text(
        json.dumps({"policy": {"template": "residue_contacts"}})
    )
    out = _import(
        tmp_path,
        "items: interface\n",
        "--contacts",
        str(contacts),
        rows=[
            ("1A2C_B_A", "1A2C_A_B", "1", "w"),
            ("1A2C_A_C", "1A2C_A_B", "1", "w"),
        ],
    )

    assert [r["item_id_1"] for r in _parquet_rows(out / "layer.parquet")] == [
        "1A2C_A_B"
    ]
    assert [r["reason"] for r in _parquet_rows(out / "unmatched.parquet")] == [
        "interface_not_found"
    ]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["previous_summary"] == str(contacts / "layer.json")


def test_similarity_import_interface_level_needs_contacts(tmp_path):
    with pytest.raises(SystemExit, match="--contacts"):
        _import(
            tmp_path,
            "items: interface\n",
            "--structures-dir",
            str(_similarity_step(tmp_path)),
        )


def test_similarity_import_chain_level_without_chains_table(tmp_path):
    folder = _similarity_step(tmp_path)
    (folder / "chains.jsonl").unlink()
    with pytest.raises(SystemExit, match="pandora curate"):
        _import(tmp_path, "items: chain\n", "--structures-dir", str(folder))


def test_similarity_import_missing_items_stops_before_writing(tmp_path):
    with pytest.raises(SystemExit, match="items"):
        _import(
            tmp_path, "", "--structures-dir", str(_similarity_step(tmp_path))
        )
    assert not (tmp_path / "similarity" / "r_graph").exists()


def test_similarity_import_missing_score_stops_before_writing(tmp_path):
    (tmp_path / "g.tsv").write_text("item_id_1\titem_id_2\n104M\t112M\n")
    (tmp_path / "g.yaml").write_text(
        "name: g\ntemplate: similarity\nitems: entry\ntool: R\n"
    )
    with pytest.raises(SystemExit, match="score"):
        main(
            [
                "similarity",
                "import",
                "--table",
                str(tmp_path / "g.tsv"),
                "--description",
                str(tmp_path / "g.yaml"),
                "--structures-dir",
                str(_similarity_step(tmp_path)),
                "--output-dir",
                str(tmp_path / "out"),
            ]
        )
    assert not (tmp_path / "out").exists()


def test_cluster_interface_filter_without_contacts_exits(tmp_path):
    step, layer = _cluster_inputs(
        tmp_path, _chain_rows(("104M_A", "112M_A")), items="chain"
    )
    hit_filter = tmp_path / "filter.yaml"
    hit_filter.write_text("min_interface_coverage: 0.5\n")
    with pytest.raises(SystemExit, match="interface_residues"):
        _cluster(
            step, layer, tmp_path / "c.json", "--hit-filter", str(hit_filter)
        )


def _count(contacts, out, level):
    main(
        ["annotate", "count-contacts", "--contacts", str(contacts)]
        + ["--level", level, "--output-dir", str(out)]
        + ["--table-format", "jsonl"]
    )


@pytest.mark.parametrize(
    "level, template, function",
    [
        ("residue", "residue_values", count_residue_contacts),
        ("chain", "chain_values", count_chain_contacts),
    ],
)
def test_annotate_count_contacts_writes_value_layer(
    tmp_path, level, template, function
):
    contacts = tmp_path / "interfaces"
    _compute(_step_folder(tmp_path), contacts, "--layer", "chain_interfaces")
    out = tmp_path / f"counts_{level}"
    _count(contacts, out, level)

    rows = read_records(contacts / "layer.jsonl", ResidueContactRow)
    written = read_records(out / "layer.jsonl", LayerRow)
    assert written
    assert [r.model_dump() for r in written] == [
        r.model_dump() for r in function(rows)
    ]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "annotation_count"
    assert summary["policy"]["template"] == template
    assert summary["policy"]["name"] == "contact_count"
    assert summary["policy"]["tool"] == "pandora"
    assert summary["policy"]["parameters"] == {
        "level": level,
        "contacts_layer": "chain_interfaces",
    }
    assert summary["counts"] == {
        "rows_in": len(rows),
        "rows_out": len(written),
    }
    assert summary["previous_summary"] == str(contacts / "layer.json")


def test_annotate_count_contacts_rejects_non_contact_layer(tmp_path):
    contacts = tmp_path / "interfaces"
    _compute(_step_folder(tmp_path), contacts, "--layer", "chain_interfaces")
    counts = tmp_path / "counts"
    _count(contacts, counts, "chain")
    out = tmp_path / "again"
    with pytest.raises(SystemExit, match="needs residue_contacts"):
        _count(counts, out, "chain")
    assert not out.exists()


def test_annotate_count_contacts_empty_layer(tmp_path):
    contacts = tmp_path / "empty"
    contacts.mkdir()
    (contacts / "layer.json").write_text(
        json.dumps({"policy": {"name": "x", "template": "residue_contacts"}})
    )
    out = tmp_path / "counts"
    _count(contacts, out, "residue")
    assert sorted(p.name for p in out.iterdir()) == ["layer.json"]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["counts"] == {"rows_in": 0, "rows_out": 0}


def test_annotate_count_contacts_refuses_its_input_as_output(tmp_path):
    contacts = tmp_path / "interfaces"
    _compute(_step_folder(tmp_path), contacts, "--layer", "chain_interfaces")
    before = (contacts / "layer.jsonl").read_text()
    with pytest.raises(SystemExit, match="is the --contacts folder"):
        main(
            ["annotate", "count-contacts", "--contacts", str(contacts)]
            + ["--level", "chain", "--output-dir", str(contacts)]
            + ["--overwrite", "--table-format", "jsonl"]
        )
    assert (contacts / "layer.jsonl").read_text() == before


def _import_labels(tmp_path, step, rows, out, template="interface_labels"):
    table = tmp_path / f"{out.name}.jsonl"
    _write_jsonl(table, rows)
    description = tmp_path / f"{out.name}.yaml"
    description.write_text(
        f"name: {out.name}\ntemplate: {template}\ntool: pandas\n"
        'tool_version: "1"\n'
    )
    main(
        ["annotate", "import", "--structures-dir", str(step)]
        + ["--table", str(table), "--description", str(description)]
        + ["--output-dir", str(out), "--table-format", "jsonl"]
    )


def _contact_pairs(contacts):
    rows = _jsonl(contacts / "layer.jsonl")
    return sorted({(r["chain_id_1"], r["chain_id_2"]) for r in rows})


def _verdicts(pairs):
    """Label the first pair `bad`, the rest `ok`."""

    return [
        {"entry_id": "1a2c", "chain_id_1": c1, "chain_id_2": c2}
        | {"label": "bad" if i == 0 else "ok"}
        for i, (c1, c2) in enumerate(pairs)
    ]


def _filter_inputs(tmp_path, labels_step=None):
    """A contact layer of 1a2c and a labels layer marking its first
    interface `bad`; returns (step, contacts, labels)."""

    step = _step_folder(tmp_path)
    (step / "curation.json").write_text("{}")
    contacts = tmp_path / "annotations" / "interfaces"
    _compute(step, contacts, "--layer", "chain_interfaces")
    labels = tmp_path / "annotations" / "verdict"
    _import_labels(
        tmp_path,
        labels_step or step,
        _verdicts(_contact_pairs(contacts)),
        labels,
    )
    return step, contacts, labels


def _filter(contacts, labels, out, *extra):
    main(
        ["annotate", "filter-contacts", "--contacts", str(contacts)]
        + ["--labels", str(labels), "--output-dir", str(out)]
        + ["--table-format", "jsonl", *extra]
    )


def test_annotate_filter_contacts_writes_layer(tmp_path):
    _, contacts, labels = _filter_inputs(tmp_path)
    out = tmp_path / "filtered"
    _filter(contacts, labels, out, "--drop", "bad")

    rows = read_records(contacts / "layer.jsonl", ResidueContactRow)
    label_rows = read_records(labels / "layer.jsonl", LayerRow)
    written = read_records(out / "layer.jsonl", ResidueContactRow)
    assert written == drop_labelled_interfaces(rows, label_rows, {"bad"})
    bad = set(_contact_pairs(contacts)[0])
    assert all({r.chain_id_1, r.chain_id_2} != bad for r in written)
    assert len(written) < len(rows)

    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "annotation_filter"
    assert summary["policy"]["name"] == "filtered_contacts"
    assert summary["policy"]["template"] == "residue_contacts"
    assert summary["policy"]["parameters"] == {
        "contacts_layer": "chain_interfaces",
        "labels_layer": "verdict",
        "labels_dir": str(labels),
        "drop": ["bad"],
    }
    assert summary["previous_summary"] == str(contacts / "layer.json")
    n = len(_contact_pairs(contacts))
    assert summary["counts"] == {
        "interfaces_dropped": 1,
        "interfaces_in": n,
        "interfaces_unlabelled": 0,
        "labels_unmatched": 0,
        "rows_in": len(rows),
        "rows_out": len(written),
    }


def test_annotate_filter_contacts_output_feeds_records(tmp_path):
    _, contacts, labels = _filter_inputs(tmp_path)
    out = tmp_path / "filtered"
    _filter(contacts, labels, out, "--drop", "bad")
    rows = read_records(out / "layer.jsonl", ResidueContactRow)
    pairs = {
        frozenset((r.chain_id_1, r.chain_id_2))
        for r in extract_interface_records(rows)
    }
    assert pairs
    assert frozenset(_contact_pairs(contacts)[0]) not in pairs


def test_annotate_filter_contacts_chains(tmp_path):
    _, contacts, labels = _filter_inputs(tmp_path)
    first = tmp_path / "first"
    _filter(contacts, labels, first, "--drop", "bad")
    second = tmp_path / "second"
    _filter(first, labels, second, "--drop", "bad")
    assert (second / "layer.jsonl").read_text() == (
        first / "layer.jsonl"
    ).read_text()
    counts = json.loads((second / "layer.json").read_text())["counts"]
    assert counts["interfaces_dropped"] == 0
    assert counts["labels_unmatched"] == 1


def test_annotate_filter_contacts_rejects_wrong_labels_template(tmp_path):
    step, contacts, _ = _filter_inputs(tmp_path)
    chain_labels = tmp_path / "chain_verdict"
    _import_labels(
        tmp_path,
        step,
        [{"entry_id": "1a2c", "chain_id": "A", "label": "bad"}],
        chain_labels,
        template="chain_labels",
    )
    out = tmp_path / "filtered"
    with pytest.raises(SystemExit, match="needs interface_labels"):
        _filter(contacts, chain_labels, out, "--drop", "bad")
    assert not out.exists()


def test_annotate_filter_contacts_rejects_unknown_drop(tmp_path):
    _, contacts, labels = _filter_inputs(tmp_path)
    out = tmp_path / "filtered"
    with pytest.raises(SystemExit, match="typo"):
        _filter(contacts, labels, out, "--drop", "typo")
    assert not out.exists()


@pytest.mark.parametrize("which", ["contacts", "labels"])
def test_annotate_filter_contacts_refuses_input_as_output(tmp_path, which):
    _, contacts, labels = _filter_inputs(tmp_path)
    target = contacts if which == "contacts" else labels
    before = (target / "layer.jsonl").read_text()
    with pytest.raises(SystemExit, match=f"is the --{which} folder"):
        _filter(contacts, labels, target, "--drop", "bad", "--overwrite")
    assert (target / "layer.jsonl").read_text() == before


def test_annotate_filter_contacts_rejects_different_origins(tmp_path):
    other = tmp_path / "other"
    shutil.copytree(FIXTURES_DIR / "a2", other / "structures" / "a2")
    (other / "curation.json").write_text("{}")
    _, contacts, labels = _filter_inputs(tmp_path, labels_step=other)
    out = tmp_path / "filtered"
    with pytest.raises(SystemExit, match="same structures folder"):
        _filter(contacts, labels, out, "--drop", "bad")
    assert not out.exists()


def test_annotate_filter_contacts_warns_without_origin(tmp_path, capsys):
    other = tmp_path / "other"
    (other / "structures" / "a2").mkdir(parents=True)
    shutil.copy(FIXTURES_DIR / "a2" / "1a2c.cif", other / "structures" / "a2")
    _, contacts, labels = _filter_inputs(tmp_path, labels_step=other)
    out = tmp_path / "filtered"
    _filter(contacts, labels, out, "--drop", "bad")
    assert "warning" in capsys.readouterr().err
    assert (out / "layer.json").exists()


def test_annotate_filter_contacts_empty_contacts(tmp_path):
    step, _, labels = _filter_inputs(tmp_path)
    empty = tmp_path / "empty"
    empty.mkdir()
    (empty / "layer.json").write_text(
        json.dumps(
            {
                "step": "annotation_compute",
                "input": str(step),
                "previous_summary": None,
                "policy": {"name": "x", "template": "residue_contacts"},
            }
        )
    )
    out = tmp_path / "filtered"
    _filter(empty, labels, out, "--drop", "bad")
    assert sorted(p.name for p in out.iterdir()) == ["layer.json"]
    counts = json.loads((out / "layer.json").read_text())["counts"]
    n = len(read_records(labels / "layer.jsonl", LayerRow))
    assert counts == {
        "interfaces_dropped": 0,
        "interfaces_in": 0,
        "interfaces_unlabelled": 0,
        "labels_unmatched": n,
        "rows_in": 0,
        "rows_out": 0,
    }


def test_annotate_filter_contacts_from_another_cwd(
    tmp_path, monkeypatch, capsys
):
    monkeypatch.chdir(tmp_path)
    step = _step_folder(Path("."))
    (step / "curation.json").write_text("{}")
    contacts = Path("annotations") / "interfaces"
    _compute(step, contacts, "--layer", "chain_interfaces")
    labels = Path("annotations") / "verdict"
    _import_labels(Path("."), step, _verdicts(_contact_pairs(contacts)), labels)
    _filter(contacts, labels, Path("first"), "--drop", "bad")

    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    capsys.readouterr()
    out = tmp_path / "second"
    _filter(tmp_path / "first", tmp_path / labels, out, "--drop", "bad")
    assert "warning" in capsys.readouterr().err
    assert (out / "layer.json").exists()


def test_layer_origin_stops_on_a_loop(tmp_path):
    for name, other in (("a", "b"), ("b", "a")):
        (tmp_path / name).mkdir()
        (tmp_path / name / "layer.json").write_text(
            json.dumps(
                {"previous_summary": str(tmp_path / other / "layer.json")}
            )
        )
    assert layer_origin(tmp_path / "a") is None


def _kinds_clusters(tmp_path, groups):
    path = tmp_path / "clusters.json"
    path.write_text(
        json.dumps([{"components": g, "n_components": len(g)} for g in groups])
    )
    return path


def _labels_layer(tmp_path, template, rows):
    folder = tmp_path / "labels"
    folder.mkdir()
    (folder / "layer.json").write_text(
        json.dumps({"policy": {"name": "kinds", "template": template}})
    )
    if rows is not None:
        _write_jsonl(folder / "layer.jsonl", rows)
    return folder


_KIND_GROUPS = [
    [f"1A{i}_A" for i in range(6)],
    [f"1B{i}_A" for i in range(4)],
    ["1H0_A", "1H0_B"],
    ["1H1_A"],
    ["1H2_A"],
    ["1S0_A"],
    ["1S1_A"],
]


def _kind_rows(groups, hetero_entries):
    return [
        {
            "entry_id": item.split("_")[0],
            "chain_id": item.split("_")[1],
            "label": "heteromer"
            if item.split("_")[0] in hetero_entries
            else "homomer",
        }
        for g in groups
        for item in g
    ]


def _partition(tmp_path, clusters, labels, *extra):
    return main(
        [
            "partition",
            "--clusters",
            str(clusters),
            "--labels",
            str(labels),
            "--output",
            str(tmp_path / "out" / "splits.json"),
            *extra,
        ]
    )


def test_partition_labels_writes_counts(tmp_path):
    clusters = _kinds_clusters(tmp_path, _KIND_GROUPS)
    rows = _kind_rows(_KIND_GROUPS, {"1H0", "1H1", "1H2"})
    labels = _labels_layer(tmp_path, "chain_labels", rows)
    _partition(tmp_path, clusters, labels)
    splits = json.loads((tmp_path / "out" / "splits.json").read_text())
    assert set(splits) == {"train", "val", "test"}
    prov = json.loads(
        (tmp_path / "out" / "partition_provenance.json").read_text()
    )
    assert prov["missing_labels"] == {}
    assert prov["labels_layer"] == str(labels / "layer.json")
    for counts in prov["label_counts"].values():
        assert set(counts) == {"homomer", "heteromer"}


def test_partition_labels_warns_missing(tmp_path, capsys):
    clusters = _kinds_clusters(tmp_path, _KIND_GROUPS)
    rows = _kind_rows(_KIND_GROUPS, {"1H0"})
    labels = _labels_layer(tmp_path, "chain_labels", rows)
    _partition(tmp_path, clusters, labels)
    err = capsys.readouterr().err
    assert err.count("has no item labelled heteromer") == 2


def test_partition_labels_wrong_level(tmp_path, capsys):
    clusters = _kinds_clusters(tmp_path, _KIND_GROUPS)
    labels = _labels_layer(
        tmp_path,
        "interface_labels",
        [
            {
                "entry_id": "1H0",
                "chain_id_1": "A",
                "chain_id_2": "B",
                "label": "heteromer",
            }
        ],
    )
    _partition(tmp_path, clusters, labels)
    assert "unlabelled items: 16 of 16" in capsys.readouterr().out


def test_partition_labels_empty_layer(tmp_path, capsys):
    clusters = _kinds_clusters(tmp_path, _KIND_GROUPS)
    labels = _labels_layer(tmp_path, "chain_labels", None)
    _partition(tmp_path, clusters, labels)
    assert "unlabelled items: 16 of 16" in capsys.readouterr().out


def test_partition_labels_rejects_template(tmp_path):
    clusters = _kinds_clusters(tmp_path, _KIND_GROUPS)
    labels = _labels_layer(tmp_path, "residue_values", [])
    with pytest.raises(SystemExit, match="residue_values"):
        _partition(tmp_path, clusters, labels)
    assert not (tmp_path / "out" / "splits.json").exists()


def test_partition_labels_rejects_no_keep(tmp_path):
    clusters = _kinds_clusters(tmp_path, _KIND_GROUPS)
    labels = _labels_layer(tmp_path, "chain_labels", [])
    with pytest.raises(SystemExit, match="no-keep-similar-items"):
        _partition(tmp_path, clusters, labels, "--no-keep-similar-items")
    assert not (tmp_path / "out" / "splits.json").exists()


def _label(contacts, out, *extra):
    main(
        ["annotate", "label-interfaces", "--contacts", str(contacts)]
        + ["--output-dir", str(out), "--table-format", "jsonl", *extra]
    )


def _label_inputs(tmp_path, skip=()):
    """1a2c step folder with a chains table, and its contact layer."""

    step = _step_folder(tmp_path)
    contacts = tmp_path / "interfaces"
    _compute(step, contacts, "--layer", "chain_interfaces")
    chains = sorted({c for pair in _contact_pairs(contacts) for c in pair})
    rows = [
        {"entry_id": "1A2C", "chain_id": c, "entity_id": "1"}
        | {"residue_count": 1, "atom_count": 1, "seqres_length": 100}
        for c in chains
        if c not in skip
    ]
    (step / "chains.jsonl").write_text(
        "".join(json.dumps(r) + "\n" for r in rows)
    )
    return step, contacts


def test_annotate_label_interfaces_writes_layer(tmp_path):
    step, contacts = _label_inputs(tmp_path)
    out = tmp_path / "labels"
    _label(contacts, out, "--min-residues-per-chain", "5")

    rows = read_records(contacts / "layer.jsonl", ResidueContactRow)
    lengths = {
        ("1A2C", c): 100 for pair in _contact_pairs(contacts) for c in pair
    }
    written = read_records(out / "layer.jsonl", LayerRow)
    expected = label_interfaces(rows, lengths, 5, None)
    assert written
    assert [r.model_dump() for r in written] == [
        r.model_dump() for r in expected
    ]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "annotation_label"
    assert summary["policy"]["name"] == "contact_labels"
    assert summary["policy"]["template"] == "interface_labels"
    assert summary["policy"]["parameters"] == {
        "contacts_layer": "chain_interfaces",
        "structures_dir": str(step.resolve()),
        "min_residues_per_chain": 5,
        "max_total_residues": None,
    }
    labels = [r.label for r in written]
    assert summary["counts"] == {
        "rows_in": len(rows),
        "interfaces": len(written),
        "ok": labels.count("ok"),
        "few_contacts": labels.count("few_contacts"),
        "too_large": 0,
        "unknown_length": 0,
    }
    assert summary["previous_summary"] == str(contacts / "layer.json")


def test_annotate_label_interfaces_feeds_filter(tmp_path):
    _, contacts = _label_inputs(tmp_path)
    labels = tmp_path / "labels"
    _label(contacts, labels, "--min-residues-per-chain", "1000")
    out = tmp_path / "real"
    _filter(contacts, labels, out, "--drop", "few_contacts")
    # every interface dropped: filter-contacts writes no rows table
    assert not (out / "layer.jsonl").exists()
    counts = json.loads((out / "layer.json").read_text())["counts"]
    assert counts["interfaces_dropped"] == counts["interfaces_in"] > 0


def test_annotate_label_interfaces_untraceable_origin_exits(tmp_path):
    _, contacts = _label_inputs(tmp_path)
    summary = json.loads((contacts / "layer.json").read_text())
    summary["input"] = str(tmp_path / "gone")
    (contacts / "layer.json").write_text(json.dumps(summary))
    with pytest.raises(SystemExit, match="relative to where"):
        _label(contacts, tmp_path / "labels")
    assert not (tmp_path / "labels").exists()


def test_annotate_label_interfaces_missing_chain_has_unknown_length(
    tmp_path,
):
    _, probe = _label_inputs(tmp_path / "probe")
    first = _contact_pairs(probe)[0][0]
    _, contacts = _label_inputs(tmp_path, skip={first})
    out = tmp_path / "labels"
    _label(contacts, out, "--max-total-residues", "150")
    written = read_records(out / "layer.jsonl", LayerRow)
    unknown = [r for r in written if first in (r.chain_id_1, r.chain_id_2)]
    assert unknown and all(r.total_residues is None for r in unknown)
    assert all(r.label == "ok" for r in unknown)
    assert all(r.label == "too_large" for r in written if r not in unknown)
    counts = json.loads((out / "layer.json").read_text())["counts"]
    assert counts["unknown_length"] == len(unknown)


def test_annotate_label_interfaces_without_chains_table_exits(tmp_path):
    step, contacts = _label_inputs(tmp_path)
    (step / "chains.jsonl").unlink()
    with pytest.raises(SystemExit, match="pandora curate"):
        _label(contacts, tmp_path / "labels")
    assert not (tmp_path / "labels").exists()


def _gaps(contacts, out):
    main(
        ["annotate", "interface-gaps", "--contacts", str(contacts)]
        + ["--output-dir", str(out), "--table-format", "jsonl"]
    )


def test_annotate_interface_gaps_writes_layer(tmp_path):
    step = _step_folder(tmp_path)
    contacts = tmp_path / "interfaces"
    _compute(step, contacts, "--layer", "chain_interfaces")
    out = tmp_path / "gaps"
    _gaps(contacts, out)

    written = read_records(out / "layer.jsonl", LayerRow)
    # 1a2c: heavy chain B (auth H) has a missing loop inside both of
    # its interfaces with the light chain and with hirudin (D/J).
    assert [
        (r.chain_id_1, r.chain_id_2, r.gaps_1, r.gaps_2, r.label)
        for r in written
    ] == [
        ("A", "B", 0, 7, "interface_gap"),
        ("B", "C", 0, 0, "ok"),
        ("B", "D", 7, 0, "interface_gap"),
    ]
    summary = json.loads((out / "layer.json").read_text())
    assert summary["policy"]["name"] == "interface_gaps"
    assert summary["policy"]["template"] == "interface_labels"
    assert summary["policy"]["parameters"] == {
        "contacts_layer": "chain_interfaces",
        "structures_dir": str(step.resolve()),
        "missing_residue_definition": "incomplete_backbone",
    }
    assert summary["counts"] == {
        "rows_in": 120,
        "interfaces": 3,
        "ok": 1,
        "interface_gap": 2,
        "unknown": 0,
        "chains_unmeasured": 0,
    }

    kept = tmp_path / "nogap"
    _filter(contacts, out, kept, "--drop", "interface_gap")
    assert {tuple(sorted(pair)) for pair in _contact_pairs(kept)} == {
        ("B", "C")
    }


def test_annotate_label_interfaces_output_is_contacts_exits(tmp_path):
    _, contacts = _label_inputs(tmp_path)
    before = (contacts / "layer.json").read_text()
    with pytest.raises(SystemExit, match="--contacts folder"):
        _label(contacts, contacts)
    assert (contacts / "layer.json").read_text() == before


def _iface_inputs(tmp_path, clusters=None):
    """1a2c step folder whose chains table has sequences that fit its
    contact layer, plus a clusters file (default: all chains in one)."""

    import gemmi

    step = _step_folder(tmp_path)
    contacts = tmp_path / "interfaces"
    _compute(step, contacts, "--layer", "chain_interfaces")
    letters: dict[str, dict[int, str]] = {}
    for r in _jsonl(contacts / "layer.jsonl"):
        for side in ("1", "2"):
            seq_id = r[f"seq_id_{side}"]
            if seq_id is None:
                continue
            info = gemmi.find_tabulated_residue(r[f"comp_id_{side}"])
            code = info.one_letter_code.strip().upper() if info else ""
            letters.setdefault(r[f"chain_id_{side}"], {})[seq_id] = code or "X"
    seqs = {
        c: "".join(pos.get(i, "A") for i in range(1, max(pos) + 1))
        for c, pos in letters.items()
    }
    rows = [
        {"entry_id": "1A2C", "chain_id": c, "entity_id": "1"}
        | {"residue_count": 1, "atom_count": 1, "sequence": seq}
        for c, seq in seqs.items()
    ]
    _write_jsonl(step / "chains.jsonl", rows)
    if clusters is None:
        clusters = [[f"1A2C_{c}" for c in sorted(seqs)]]
    clusters_file = tmp_path / "seq90.json"
    clusters_file.write_text(
        json.dumps(
            [{"components": c, "n_components": len(c)} for c in clusters]
        )
    )
    return step, contacts, clusters_file, seqs


def _iface_sim(contacts, clusters_file, out, *extra):
    main(
        ["similarity", "interfaces", "--contacts", str(contacts)]
        + ["--chain-clusters", str(clusters_file)]
        + ["--output-dir", str(out), *extra]
    )


def _iface_cluster(layer, contacts, output):
    main(
        ["cluster", "--similarity", str(layer), "--items", "interface"]
        + ["--contacts", str(contacts), "--min-score", "0"]
        + ["--output", str(output)]
    )
    return [c["components"] for c in json.loads(output.read_text())]


def test_similarity_interfaces_writes_layer(tmp_path, capsys):
    step, contacts, clusters_file, seqs = _iface_inputs(tmp_path)
    out = tmp_path / "iface_jaccard"
    _iface_sim(contacts, clusters_file, out, "--min-identity", "0")

    rows = read_records(contacts / "layer.jsonl", ResidueContactRow)
    clusters = [
        SimilarityCluster.model_validate(c)
        for c in json.loads(clusters_file.read_text())
    ]
    sequences = {("1A2C", c): seq for c, seq in seqs.items()}
    expected, counts = interface_similarity(rows, clusters, sequences, 0.0)
    written = _parquet_rows(out / "layer.parquet")
    assert written
    assert written == [r.model_dump() for r in expected]

    summary = json.loads((out / "layer.json").read_text())
    assert summary["step"] == "similarity_interfaces"
    assert summary["policy"]["name"] == "iface_jaccard"
    assert summary["policy"]["template"] == "similarity"
    assert summary["policy"]["items"] == "interface"
    assert summary["policy"]["parameters"] == {
        "contacts_layer": "chain_interfaces",
        "structures_dir": str(step.resolve()),
        "chain_clusters": str(clusters_file),
        "min_identity": 0.0,
    }
    assert summary["counts"] == dict(
        sorted(({"rows_in": len(rows)} | counts).items())
    )
    assert summary["previous_summary"] == str(contacts / "layer.json")
    assert f"scored: {counts['pairs_compared']} interface pairs" in (
        capsys.readouterr().out
    )


def test_similarity_interfaces_feeds_cluster(tmp_path):
    _, contacts, clusters_file, _ = _iface_inputs(tmp_path)
    out = tmp_path / "iface_jaccard"
    _iface_sim(contacts, clusters_file, out, "--min-identity", "0")
    groups = _iface_cluster(out, contacts, tmp_path / "dimers.json")
    pairs = _contact_pairs(contacts)
    expected = sorted({"1A2C_" + "_".join(sorted(p)) for p in pairs})
    assert [sorted(g) for g in groups] == [expected]


def test_similarity_interfaces_no_pairs(tmp_path):
    probe = _iface_inputs(tmp_path / "probe")[3]
    singletons = [[f"1A2C_{c}"] for c in sorted(probe)]
    _, contacts, clusters_file, _ = _iface_inputs(tmp_path, singletons)
    out = tmp_path / "iface_jaccard"
    _iface_sim(contacts, clusters_file, out)
    assert (out / "layer.json").exists()
    assert not (out / "layer.parquet").exists()
    groups = _iface_cluster(out, contacts, tmp_path / "dimers.json")
    assert groups and all(len(g) == 1 for g in groups)


def test_similarity_interfaces_untraceable_origin_exits(tmp_path):
    _, contacts, clusters_file, _ = _iface_inputs(tmp_path)
    summary = json.loads((contacts / "layer.json").read_text())
    summary["input"] = str(tmp_path / "gone")
    (contacts / "layer.json").write_text(json.dumps(summary))
    with pytest.raises(SystemExit, match="relative to where"):
        _iface_sim(contacts, clusters_file, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_similarity_interfaces_no_chains_table_exits(tmp_path):
    step, contacts, clusters_file, _ = _iface_inputs(tmp_path)
    (step / "chains.jsonl").unlink()
    with pytest.raises(SystemExit, match="pandora curate"):
        _iface_sim(contacts, clusters_file, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_similarity_interfaces_output_is_contacts_exits(tmp_path):
    _, contacts, clusters_file, _ = _iface_inputs(tmp_path)
    before = (contacts / "layer.json").read_text()
    with pytest.raises(SystemExit, match="--contacts folder"):
        _iface_sim(contacts, clusters_file, contacts)
    assert (contacts / "layer.json").read_text() == before


def test_similarity_interfaces_wrong_level_clusters_exits(tmp_path):
    _, contacts, clusters_file, _ = _iface_inputs(tmp_path, [["1A2C"]])
    with pytest.raises(SystemExit, match="chain items"):
        _iface_sim(contacts, clusters_file, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def _contained_inputs(tmp_path, clusters=None):
    """A curate-like step folder (1AAA: A, res 2.0; 2BBB: A+B, res 1.0;
    3CCC: A, res 1.5) with dummy structure files, and a clusters.json
    (default: [1AAA_A, 2BBB_A, 2BBB_B], [3CCC_A])."""

    step = tmp_path / "curation"
    entries = []
    chains = []
    for entry, chain_ids, res in (
        ("1AAA", ["A"], 2.0),
        ("2BBB", ["A", "B"], 1.0),
        ("3CCC", ["A"], 1.5),
    ):
        path = f"structures/{entry[1:3].lower()}/{entry.lower()}.cif"
        (step / path).parent.mkdir(parents=True, exist_ok=True)
        (step / path).write_text(f"data_{entry}\n")
        entries.append(
            {
                "entry_id": entry,
                "structure_path": path,
                "resolution": res,
                "chain_count": len(chain_ids),
                "policy_id": "p",
                "policy_name": "p",
                "policy_version": "1",
                "curated_at": "2026-11-10T00:00:00+00:00",
            }
        )
        chains += [
            {
                "entry_id": entry,
                "chain_id": c,
                "entity_id": "1",
                "residue_count": 10,
                "atom_count": 80,
            }
            for c in chain_ids
        ]
    _write_jsonl(step / "entries.jsonl", entries)
    _write_jsonl(step / "chains.jsonl", chains)
    if clusters is None:
        clusters = [["1AAA_A", "2BBB_A", "2BBB_B"], ["3CCC_A"]]
    clusters_file = tmp_path / "clusters.json"
    clusters_file.write_text(
        json.dumps(
            [{"components": c, "n_components": len(c)} for c in clusters]
        )
    )
    return step, clusters_file


def _dedup_contained(step, clusters, out, *extra):
    main(
        ["dedup", "--contained", "--input-dir", str(step)]
        + ["--clusters", str(clusters), "--output-dir", str(out)]
        + ["--table-format", "jsonl", *extra]
    )


def _jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def test_dedup_contained_writes_step_folder(tmp_path):
    step, clusters = _contained_inputs(tmp_path)
    out = tmp_path / "uncontained"
    _dedup_contained(step, clusters, out)
    entries = [row["entry_id"] for row in _jsonl(out / "entries.jsonl")]
    assert entries == ["2BBB", "3CCC"]
    assert len(_jsonl(out / "chains.jsonl")) == 3
    assert _jsonl(out / "contained.jsonl") == [
        {"entry_id": "1AAA", "contained_in": "2BBB", "equal": False}
    ]
    assert (out / "structures/bb/2bbb.cif").exists()
    assert not (out / "structures/aa/1aaa.cif").exists()
    assert (out / "failures.jsonl").read_text() == ""
    summary = json.loads((out / "dedup.json").read_text())
    assert summary["counts"]["entries_contained"] == 1
    assert summary["counts"]["entries_kept"] == 2
    assert summary["counts"]["chains_kept"] == 3
    assert summary["counts"]["chains_unclustered"] == 0


def test_dedup_contained_counts_unclustered_chains(tmp_path, capsys):
    step, clusters = _contained_inputs(
        tmp_path, clusters=[["1AAA_A", "2BBB_A", "2BBB_B"]]
    )
    out = tmp_path / "uncontained"
    _dedup_contained(step, clusters, out)
    summary = json.loads((out / "dedup.json").read_text())
    assert summary["counts"]["chains_unclustered"] == 1
    entries = [row["entry_id"] for row in _jsonl(out / "entries.jsonl")]
    assert "3CCC" in entries
    assert "warning" in capsys.readouterr().err


def test_dedup_contained_disable_keeps_everything(tmp_path):
    step, clusters = _contained_inputs(tmp_path)
    out = tmp_path / "uncontained"
    _dedup_contained(step, clusters, out, "--disable")
    assert len(_jsonl(out / "entries.jsonl")) == 3
    assert [r["entry_id"] for r in _jsonl(out / "contained.jsonl")] == ["1AAA"]
    assert (out / "structures/aa/1aaa.cif").exists()


def test_dedup_contained_needs_clusters(tmp_path):
    with pytest.raises(SystemExit, match="--clusters"):
        main(
            ["dedup", "--contained", "--input-dir", str(tmp_path)]
            + ["--output-dir", str(tmp_path / "o")]
        )


def test_dedup_contained_rejects_unknown_rank_key(tmp_path):
    step, clusters = _contained_inputs(tmp_path)
    policy = tmp_path / "dedup.yaml"
    policy.write_text("contained_rank:\n  - key: n_atoms\n")
    with pytest.raises(SystemExit, match="n_atoms"):
        _dedup_contained(
            step, clusters, tmp_path / "o", "--policy", str(policy)
        )


def test_dedup_contained_uses_policy_rank(tmp_path):
    # Tie 1AAA vs 2BBB (both one chain in the same cluster); rank on
    # resolution max keeps the worse one (1AAA, 2.0).
    step, clusters = _contained_inputs(
        tmp_path, clusters=[["1AAA_A", "2BBB_A"], ["2BBB_B"], ["3CCC_A"]]
    )
    rows = _jsonl(step / "chains.jsonl")
    _write_jsonl(
        step / "chains.jsonl", [r for r in rows if r["chain_id"] != "B"]
    )
    policy = tmp_path / "dedup.yaml"
    policy.write_text("contained_rank:\n  - key: resolution\n    order: max\n")
    out = tmp_path / "o"
    _dedup_contained(step, clusters, out, "--policy", str(policy))
    assert _jsonl(out / "contained.jsonl") == [
        {"entry_id": "2BBB", "contained_in": "1AAA", "equal": True}
    ]


def test_dedup_contained_refuses_existing_output(tmp_path):
    step, clusters = _contained_inputs(tmp_path)
    out = tmp_path / "uncontained"
    _dedup_contained(step, clusters, out)
    with pytest.raises(SystemExit):
        _dedup_contained(step, clusters, out)
    _dedup_contained(step, clusters, out, "--overwrite")
    assert (out / "dedup.json").exists()


def _rep_inputs(tmp_path, missing, hits, clusters=None):
    """Step folder (chains 104M_A, 112M_A, 112M_B copy of A, 118L_A with
    `missing` residues each), clusters.json (default: one of all four)
    and a chain similarity layer of `hits` pairs."""

    step = tmp_path / "step"
    step.mkdir()
    chains = [("104M", "A", None), ("112M", "A", None)]
    chains += [("112M", "B", "A"), ("118L", "A", None)]
    _write_jsonl(
        step / "chains.jsonl",
        [
            {
                "entry_id": e,
                "chain_id": c,
                "copy_of": copy,
                "missing_n_term": 0,
                "missing_c_term": 0,
                "missing_middle": missing[f"{e}_{c}"],
            }
            for e, c, copy in chains
        ],
    )
    if clusters is None:
        clusters = [["104M_A", "112M_A", "112M_B", "118L_A"]]
    clusters_file = tmp_path / "clusters.json"
    clusters_file.write_text(
        json.dumps(
            [{"components": c, "n_components": len(c)} for c in clusters]
        )
    )
    rows = [{"item_id_1": a, "item_id_2": b, "score": 0.95} for a, b in hits]
    layer = _similarity_layer(tmp_path / "layer", rows, "chain")
    return step, clusters_file, layer


def _representatives(step, clusters, out, *extra):
    main(
        ["representatives", "--clusters", str(clusters)]
        + ["--structures-dir", str(step), "--output", str(out), *extra]
    )
    return json.loads(out.read_text())


_MISSING = {"104M_A": 2, "112M_A": 2, "112M_B": 0, "118L_A": 0}
_ALL_HITS = [("104M_A", "112M_A"), ("104M_A", "118L_A"), ("112M_A", "118L_A")]


def test_representatives_default_policy(tmp_path):
    step, clusters, layer = _rep_inputs(tmp_path, _MISSING, _ALL_HITS)
    out = tmp_path / "reps" / "representatives.json"
    reps = _representatives(
        step, clusters, out, "--similarity", str(layer), "--min-score", "0.9"
    )
    assert reps == [
        {
            "representative": "118L_A",
            "members": ["104M_A", "112M_A", "112M_B", "118L_A"],
            "n_members": 4,
        }
    ]
    prov = json.loads(
        (out.parent / "representatives_provenance.json").read_text()
    )
    assert prov["n_clusters"] == 1 and prov["items"] == "chain"
    assert prov["n_members"] == 4 and len(prov["layers"]) == 1


def test_representatives_n_hits_breaks_tie(tmp_path):
    missing = dict.fromkeys(_MISSING, 0)
    hits = [("104M_A", "112M_A"), ("112M_A", "118L_A")]
    step, clusters, layer = _rep_inputs(tmp_path, missing, hits)
    reps = _representatives(
        step,
        clusters,
        tmp_path / "r.json",
        "--similarity",
        str(layer),
        "--min-score",
        "0.9",
    )
    assert reps[0]["representative"] == "112M_A"


def test_representatives_drops_chains_outside_step_folder(tmp_path):
    step, clusters, layer = _rep_inputs(
        tmp_path,
        _MISSING,
        _ALL_HITS,
        clusters=[["104M_A", "112M_A", "112M_B", "118L_A", "999X_A"]]
        + [["999Y_A"]],
    )
    out = tmp_path / "r.json"
    reps = _representatives(
        step, clusters, out, "--similarity", str(layer), "--min-score", "0.9"
    )
    assert len(reps) == 1 and "999X_A" not in reps[0]["members"]
    prov = json.loads(
        (tmp_path / "representatives_provenance.json").read_text()
    )
    assert prov["n_members_dropped"] == 2
    assert prov["n_clusters_dropped"] == 1


def test_representatives_n_hits_needs_similarity(tmp_path):
    step, clusters, _ = _rep_inputs(tmp_path, _MISSING, _ALL_HITS)
    with pytest.raises(SystemExit, match="n_hits"):
        _representatives(step, clusters, tmp_path / "r.json")


def test_representatives_policy_without_n_hits_needs_no_layer(tmp_path):
    step, clusters, _ = _rep_inputs(tmp_path, _MISSING, _ALL_HITS)
    policy = tmp_path / "reps.yaml"
    policy.write_text("rank:\n  - key: missing_residues\n")
    reps = _representatives(
        step, clusters, tmp_path / "r.json", "--policy", str(policy)
    )
    assert reps[0]["representative"] == "118L_A"


def test_representatives_interface_items_need_contacts(tmp_path):
    step, clusters, _ = _rep_inputs(tmp_path, _MISSING, _ALL_HITS)
    with pytest.raises(SystemExit, match="--contacts"):
        _representatives(
            step, clusters, tmp_path / "r.json", "--items", "interface"
        )


def _iface_rep_inputs(tmp_path, copies=()):
    """1a2c contact layer, its chains table (`copies` marked as symmetry
    copies) and one cluster of all its interfaces, largest first."""

    step, contacts, _, seqs = _iface_inputs(tmp_path)
    _write_jsonl(
        step / "chains.jsonl",
        [
            {"entry_id": "1A2C", "chain_id": c, "resolution": 2.0}
            | {"copy_of": "X" if c in copies else None}
            | dict.fromkeys(("missing_n_term", "missing_c_term"), 0)
            | {"missing_middle": 0}
            for c in seqs
        ],
    )
    rows = read_records(contacts / "layer.jsonl", ResidueContactRow)
    labelled = label_interfaces(rows, {})
    labelled.sort(key=lambda r: -(r.residues_1 + r.residues_2))
    ids = [f"1A2C_{r.chain_id_1}_{r.chain_id_2}" for r in labelled]
    clusters = tmp_path / "dimers.json"
    clusters.write_text(
        json.dumps([{"components": sorted(ids), "n_components": len(ids)}])
    )
    return contacts, clusters, labelled


def _iface_reps(contacts, clusters, out, *extra):
    main(
        ["representatives", "--items", "interface"]
        + ["--clusters", str(clusters), "--contacts", str(contacts)]
        + ["--output", str(out), *extra]
    )
    return json.loads(out.read_text())


def test_representatives_interfaces_pick_largest_interface(tmp_path):
    contacts, clusters, labelled = _iface_rep_inputs(tmp_path)
    best = labelled[0]
    assert best.residues_1 + best.residues_2 > (
        labelled[1].residues_1 + labelled[1].residues_2
    )
    out = tmp_path / "reps.json"
    reps = _iface_reps(contacts, clusters, out)
    assert reps[0]["representative"] == (
        f"1A2C_{best.chain_id_1}_{best.chain_id_2}"
    )
    assert reps[0]["n_members"] == len(labelled)
    prov = json.loads(
        (tmp_path / "representatives_provenance.json").read_text()
    )
    assert prov["items"] == "interface"
    assert prov["contacts"] == str(contacts)
    assert prov["policy"]["rank"][0]["key"] == "interface_residues"


def test_representatives_interfaces_skip_symmetry_copies(tmp_path):
    _, _, labelled = _iface_rep_inputs(tmp_path / "probe")
    copy = labelled[0].chain_id_1
    contacts, clusters, _ = _iface_rep_inputs(tmp_path, copies={copy})
    reps = _iface_reps(contacts, clusters, tmp_path / "reps.json")
    # The largest interface without the copy, not the lowest such id.
    rest = [r for r in labelled if copy not in (r.chain_id_1, r.chain_id_2)]
    best = f"1A2C_{rest[0].chain_id_1}_{rest[0].chain_id_2}"
    assert best != min(f"1A2C_{r.chain_id_1}_{r.chain_id_2}" for r in rest)
    assert reps[0]["representative"] == best


def test_representatives_interfaces_n_hits_from_interface_layer(tmp_path):
    contacts, clusters, labelled = _iface_rep_inputs(tmp_path)
    a, b = (f"1A2C_{r.chain_id_1}_{r.chain_id_2}" for r in labelled[-2:])
    layer = _similarity_layer(
        tmp_path / "layer",
        [{"item_id_1": a, "item_id_2": b, "score": 0.5}],
        "interface",
    )
    policy = tmp_path / "rank.yaml"
    policy.write_text("rank:\n  - key: n_hits\n    order: max\n")
    reps = _iface_reps(
        contacts,
        clusters,
        tmp_path / "reps.json",
        "--policy",
        str(policy),
        "--similarity",
        str(layer),
        "--min-score",
        "0.1",
    )
    assert reps[0]["representative"] == min(a, b)


def test_representatives_interfaces_reject_structures_dir(tmp_path):
    contacts, clusters, _ = _iface_rep_inputs(tmp_path)
    with pytest.raises(SystemExit, match="not --structures-dir"):
        _iface_reps(
            contacts,
            clusters,
            tmp_path / "reps.json",
            "--structures-dir",
            str(tmp_path),
        )


def test_representatives_rejects_unknown_rank_key(tmp_path):
    step, clusters, _ = _rep_inputs(tmp_path, _MISSING, _ALL_HITS)
    policy = tmp_path / "reps.yaml"
    policy.write_text("rank:\n  - key: bogus\n")
    with pytest.raises(SystemExit, match="bogus"):
        _representatives(
            step, clusters, tmp_path / "r.json", "--policy", str(policy)
        )


def test_layer_on_dedup_step_folder_links_its_summary(tmp_path):
    args = _layer_inputs(tmp_path)
    curated = tmp_path / "curated"
    (curated / "curation.json").rename(curated / "dedup.json")
    out = tmp_path / "layer"
    main(["annotate", "import", *args, "--output-dir", str(out)])
    summary = json.loads((out / "layer.json").read_text())
    assert summary["previous_summary"].endswith("dedup.json")


def test_plain_dedup_rejects_contained_only_flags(tmp_path):
    for flag in (["--clusters", "c.json"], ["--policy", "p.yaml"]):
        with pytest.raises(SystemExit, match="--contained"):
            main(
                ["dedup", "--input-dir", str(tmp_path)]
                + ["--output-dir", str(tmp_path / "o"), *flag]
            )


def _rebalance_inputs(tmp_path, homomers=False, groups=None):
    """1a2c contact layer, chain clusters (every chain alone, or all in
    one with `homomers`), dimer groups (default: one of every dimer)
    and an empty interface layer."""

    contacts, _, labelled = _iface_rep_inputs(tmp_path)
    dimers = sorted(f"1A2C_{r.chain_id_1}_{r.chain_id_2}" for r in labelled)
    chains = sorted(
        {f"1A2C_{c}" for r in labelled for c in (r.chain_id_1, r.chain_id_2)}
    )
    chain_groups = [chains] if homomers else [[c] for c in chains]
    seq90 = tmp_path / "seq90_clusters.json"
    seq90.write_text(
        json.dumps(
            [{"components": c, "n_components": len(c)} for c in chain_groups]
        )
    )
    groups_file = tmp_path / "groups.json"
    groups_file.write_text(
        json.dumps(
            [
                {"components": g, "n_components": len(g)}
                for g in (groups or [dimers])
            ]
        )
    )
    iface = _similarity_layer(tmp_path / "iface", None, "interface")
    return contacts, seq90, groups_file, iface, dimers


def _rebalance(contacts, seq90, groups, iface, out, *extra):
    main(
        ["rebalance", "--clusters", str(groups), "--contacts", str(contacts)]
        + ["--chain-clusters", str(seq90), "--similarity", str(iface)]
        + ["--min-score", "0.5", "--output-dir", str(out)]
        + ["--table-format", "jsonl", *extra]
    )


def test_rebalance_writes_positives_and_dropped(tmp_path):
    contacts, seq90, groups, iface, dimers = _rebalance_inputs(tmp_path)
    out = tmp_path / "positives"
    _rebalance(contacts, seq90, groups, iface, out)

    positives = _jsonl(out / "positives.jsonl")
    dropped = (
        _jsonl(out / "dropped.jsonl") if len(positives) < len(dimers) else []
    )
    assert (
        sorted([p["dimer"] for p in positives] + [d["dimer"] for d in dropped])
        == dimers
    )
    assert {p["kind"] for p in positives} == {"heteromer"}
    # No two positives share a chain (sharing a chain makes them similar).
    sides = [(p["chain_id_1"], p["chain_id_2"]) for p in positives]
    flat = [c for pair in sides for c in pair]
    assert len(flat) == len(set(flat))
    summary = json.loads((out / "rebalance.json").read_text())
    assert summary["step"] == "rebalance"
    assert summary["previous_summary"] == str(contacts / "layer.json")
    assert summary["counts"]["groups"] == 1
    assert summary["counts"]["dimers_in"] == len(dimers)
    assert summary["counts"]["positives:heteromer"] == len(positives)
    assert summary["policy"]["rank"]["items"] == "interface"


def test_rebalance_homomer_only_group_keeps_one(tmp_path):
    contacts, seq90, groups, iface, dimers = _rebalance_inputs(
        tmp_path, homomers=True
    )
    out = tmp_path / "positives"
    _rebalance(contacts, seq90, groups, iface, out)

    positives = _jsonl(out / "positives.jsonl")
    assert [p["kind"] for p in positives] == ["homomer"]
    dropped = _jsonl(out / "dropped.jsonl")
    assert {d["reason"] for d in dropped} == {"not_most_connected"}
    assert len(dropped) == len(dimers) - 1


def test_rebalance_chain_missing_from_chain_clusters_exits(tmp_path):
    contacts, seq90, groups, iface, _ = _rebalance_inputs(tmp_path)
    clusters = json.loads(seq90.read_text())
    missing = clusters.pop()["components"][0]
    seq90.write_text(json.dumps(clusters))
    with pytest.raises(SystemExit, match=missing):
        _rebalance(contacts, seq90, groups, iface, tmp_path / "out")


def test_rebalance_dimer_missing_from_contacts_exits(tmp_path):
    contacts, seq90, _, iface, dimers = _rebalance_inputs(tmp_path, groups=None)
    groups = tmp_path / "groups2.json"
    groups.write_text(
        json.dumps([{"components": [*dimers, "9ZZZ_A_B"], "n_components": 99}])
    )
    with pytest.raises(SystemExit, match="9ZZZ_A_B"):
        _rebalance(contacts, seq90, groups, iface, tmp_path / "out")


def test_rebalance_refuses_existing_output(tmp_path):
    contacts, seq90, groups, iface, _ = _rebalance_inputs(tmp_path)
    out = tmp_path / "positives"
    _rebalance(contacts, seq90, groups, iface, out)
    with pytest.raises(SystemExit, match="overwrite"):
        _rebalance(contacts, seq90, groups, iface, out)
    _rebalance(contacts, seq90, groups, iface, out, "--overwrite")
    assert (out / "rebalance.json").exists()


def test_rebalance_inputs_mismatch_exits(tmp_path):
    # Each dimer in its own group, but dimers sharing a chain are similar.
    contacts, seq90, _, iface, dimers = _rebalance_inputs(tmp_path)
    groups = tmp_path / "split.json"
    groups.write_text(
        json.dumps([{"components": [d], "n_components": 1} for d in dimers])
    )
    with pytest.raises(SystemExit, match="different groups"):
        _rebalance(contacts, seq90, groups, iface, tmp_path / "out")


_OUTSIDE = [f"{n}XXX_A" for n in range(2, 10)]


def _negatives_inputs(tmp_path):
    """A rebalance step folder over 1a2c, representatives (every chain
    its own, plus eight outside proteins) and an empty TM layer."""

    contacts, seq90, groups, iface, _ = _rebalance_inputs(tmp_path)
    positives = tmp_path / "positives"
    _rebalance(contacts, seq90, groups, iface, positives)
    chains = [c for g in json.loads(seq90.read_text()) for c in g["components"]]
    reps = _chain_reps(
        tmp_path / "reps.json",
        {c: c for c in [*chains, *_OUTSIDE]},
    )
    tm = _similarity_layer(tmp_path / "tm", None, "chain")
    return positives, reps, tm, contacts


def _negatives(positives, reps, tm, contacts, out, *extra):
    main(
        ["negatives", "--positives", str(positives)]
        + ["--representatives", str(reps), "--similarity", str(tm)]
        + ["--min-score", "0.5", "--contacts", str(contacts)]
        + ["--output-dir", str(out), "--table-format", "jsonl", *extra]
    )


def test_negatives_writes_pairs_and_summary(tmp_path):
    positives, reps, tm, contacts = _negatives_inputs(tmp_path)
    out = tmp_path / "negatives"
    _negatives(positives, reps, tm, contacts, out, "--seed", "3")

    n_positives = len(_jsonl(positives / "positives.jsonl"))
    pairs = _jsonl(out / "negatives.jsonl")
    summary = json.loads((out / "negatives.json").read_text())
    assert summary["step"] == "negatives"
    assert summary["policy"]["seed"] == 3
    assert summary["policy"]["hit_filter"]["min_score"] == 0.5
    assert summary["previous_summary"] == str(positives / "rebalance.json")
    rows = _jsonl(positives / "positives.jsonl")
    in_positives = {
        f"1A2C_{r[side]}" for r in rows for side in ("chain_id_1", "chain_id_2")
    }
    counts = summary["counts"]
    assert counts["positives"] == n_positives
    n_reps = len(json.loads(reps.read_text()))
    assert counts["pool"] == n_reps - len(in_positives)
    assert counts["negatives"] == len(pairs)
    assert counts["negatives"] + counts["unpaired"] == n_positives
    assert len(summary["policy"]["unpaired"]) == counts["unpaired"]
    assert all(p["chain_1"] in in_positives for p in pairs)
    assert not {p["chain_2"] for p in pairs} & in_positives
    assert len({p["chain_2"] for p in pairs}) == len(pairs)


def test_negatives_same_seed_reproduces(tmp_path):
    positives, reps, tm, contacts = _negatives_inputs(tmp_path)
    _negatives(positives, reps, tm, contacts, tmp_path / "n1")
    _negatives(positives, reps, tm, contacts, tmp_path / "n2")
    assert _jsonl(tmp_path / "n1" / "negatives.jsonl") == _jsonl(
        tmp_path / "n2" / "negatives.jsonl"
    )


def test_negatives_needs_chain_layer(tmp_path):
    positives, reps, _, contacts = _negatives_inputs(tmp_path)
    iface = _similarity_layer(tmp_path / "iface2", None, "interface")
    with pytest.raises(SystemExit, match="chain layer"):
        _negatives(positives, reps, iface, contacts, tmp_path / "out")


def test_negatives_empty_contacts_layer(tmp_path):
    positives, reps, tm, _ = _negatives_inputs(tmp_path)
    empty = tmp_path / "no_contacts"
    empty.mkdir()
    (empty / "layer.json").write_text(
        json.dumps({"policy": {"template": "residue_contacts"}})
    )
    out = tmp_path / "negatives"
    _negatives(positives, reps, tm, empty, out)
    summary = json.loads((out / "negatives.json").read_text())
    assert summary["counts"]["unpaired"] == 0


def test_negatives_needs_rebalance_folder(tmp_path):
    _, reps, tm, contacts = _negatives_inputs(tmp_path)
    with pytest.raises(SystemExit, match="pandora rebalance"):
        _negatives(tmp_path, reps, tm, contacts, tmp_path / "out")


def test_negatives_refuses_existing_output(tmp_path):
    positives, reps, tm, contacts = _negatives_inputs(tmp_path)
    out = tmp_path / "negatives"
    _negatives(positives, reps, tm, contacts, out)
    with pytest.raises(SystemExit, match="overwrite"):
        _negatives(positives, reps, tm, contacts, out)


def test_rebalance_passes_interface_residues_to_coverage_filter(
    tmp_path, monkeypatch
):
    # A chain layer filtered by min_interface_coverage, as `cluster` takes.
    contacts, seq90, groups, _, _ = _rebalance_inputs(tmp_path)
    row = {"item_id_1": "1A2C_A", "item_id_2": "1A2C_B", "score": 0.9}
    row |= {"start_1": 1, "end_1": 50, "start_2": 1, "end_2": 50}
    tm = _similarity_layer(tmp_path / "tm", [row], "chain")
    captured = {}

    def fake_edges(batches, hit_filter, interface_residues=None):
        captured["residues"] = interface_residues
        raise _Captured

    monkeypatch.setattr("pandora.cli.similarity.similarity_edges", fake_edges)
    hit_filter = tmp_path / "filter.yaml"
    hit_filter.write_text("min_interface_coverage: 0.5\n")
    with pytest.raises(_Captured):
        main(
            ["rebalance", "--clusters", str(groups)]
            + ["--contacts", str(contacts), "--chain-clusters", str(seq90)]
            + ["--similarity", str(tm), "--hit-filter", str(hit_filter)]
            + ["--output-dir", str(tmp_path / "out")]
        )

    assert captured["residues"]
    assert "1A2C_A" in captured["residues"]


def test_rebalance_writes_groups_and_kinds(tmp_path):
    contacts, seq90, groups, iface, _ = _rebalance_inputs(tmp_path)
    out = tmp_path / "positives"
    _rebalance(contacts, seq90, groups, iface, out)

    positives = _jsonl(out / "positives.jsonl")
    kept = json.loads((out / "groups.json").read_text())
    assert sorted(d for g in kept for d in g["components"]) == sorted(
        p["dimer"] for p in positives
    )
    assert all(g["n_components"] == len(g["components"]) for g in kept)
    kinds = _jsonl(out / "kinds" / "layer.jsonl")
    assert {
        f"{r['entry_id']}_{r['chain_id_1']}_{r['chain_id_2']}": r["label"]
        for r in kinds
    } == {p["dimer"]: p["kind"] for p in positives}
    summary = json.loads((out / "kinds" / "layer.json").read_text())
    assert summary["policy"]["template"] == "interface_labels"
    assert summary["previous_summary"] == str(out / "rebalance.json")

    split = tmp_path / "splits.json"
    main(
        ["partition", "--clusters", str(out / "groups.json")]
        + ["--labels", str(out / "kinds"), "--output", str(split)]
    )
    assert split.exists()


def test_rebalance_overwrite_clears_groups_and_kinds(tmp_path):
    contacts, seq90, groups, iface, _ = _rebalance_inputs(tmp_path)
    out = tmp_path / "positives"
    _rebalance(contacts, seq90, groups, iface, out)
    (out / "kinds" / "stale.txt").write_text("x")
    _rebalance(contacts, seq90, groups, iface, out, "--overwrite")
    assert not (out / "kinds" / "stale.txt").exists()
    assert (out / "groups.json").exists()


def test_partition_stratify_needs_labels(tmp_path):
    clusters = _kinds_clusters(tmp_path, _KIND_GROUPS)
    with pytest.raises(SystemExit, match="--labels"):
        main(
            ["partition", "--clusters", str(clusters), "--stratify"]
            + ["--output", str(tmp_path / "s.json")]
        )


def test_partition_stratify_records_flag(tmp_path):
    clusters = _kinds_clusters(tmp_path, _KIND_GROUPS)
    rows = _kind_rows(_KIND_GROUPS, {"1H0", "1H1", "1H2"})
    labels = _labels_layer(tmp_path, "chain_labels", rows)
    _partition(tmp_path, clusters, labels, "--stratify")
    prov = json.loads(
        (tmp_path / "out" / "partition_provenance.json").read_text()
    )
    assert prov["stratify"] is True


def _assemble_inputs(tmp_path):
    """Positives, negatives (seed 0), stratified splits and the inputs
    `negatives` used."""

    positives, reps, tm, contacts = _negatives_inputs(tmp_path)
    negatives = tmp_path / "negatives"
    _negatives(positives, reps, tm, contacts, negatives)
    splits = tmp_path / "splits.json"
    main(
        ["partition", "--clusters", str(positives / "groups.json")]
        + ["--labels", str(positives / "kinds"), "--stratify"]
        + ["--output", str(splits)]
    )
    return positives, negatives, splits, reps, tm


def _assemble(positives, negatives, splits, reps, tm, out, *extra):
    main(
        ["assemble", "--positives", str(positives)]
        + ["--negatives", str(negatives), "--splits", str(splits)]
        + ["--representatives", str(reps), "--similarity", str(tm)]
        + ["--min-score", "0.5", "--output-dir", str(out)]
        + ["--table-format", "jsonl", *extra]
    )


def test_assemble_writes_pairs_and_summary(tmp_path, capsys):
    positives, negatives, splits, reps, tm = _assemble_inputs(tmp_path)
    out = tmp_path / "dataset"
    _assemble(positives, negatives, splits, reps, tm, out)

    rows = _jsonl(out / "pairs.jsonl")
    n_pos = len(_jsonl(positives / "positives.jsonl"))
    n_neg = len(_jsonl(negatives / "negatives.jsonl"))
    assert [r["label"] for r in rows] == [1] * n_pos + [0] * n_neg
    split_of = {
        d: name
        for name, members in json.loads(splits.read_text()).items()
        for d in members
    }
    assert all(r["split"] == split_of[r["positive"]] for r in rows)
    summary = json.loads((out / "assemble.json").read_text())
    assert summary["step"] == "assemble"
    assert summary["previous_summary"] == str(negatives / "negatives.json")
    counts = summary["counts"]
    assert sum(v for k, v in counts.items() if k.startswith("pairs:")) == len(
        rows
    )
    assert counts["negatives_c_leak_positive"] == 0
    assert counts["negatives_c_leak_negative"] == 0
    assert counts["positives_unpaired"] == n_pos - n_neg
    assert "assemble: warning" not in capsys.readouterr().err
    # Every split and label has a count, zero included.
    for split in ("train", "val", "test"):
        for label in (0, 1):
            assert f"pairs:{split}:{label}" in counts


def test_assemble_warns_on_other_representatives(tmp_path, capsys):
    positives, negatives, splits, reps, tm = _assemble_inputs(tmp_path)
    other = tmp_path / "other_reps.json"
    other.write_text(reps.read_text())
    _assemble(positives, negatives, splits, other, tm, tmp_path / "d")
    assert "assemble: warning" in capsys.readouterr().err


def test_assemble_refuses_existing_output(tmp_path):
    positives, negatives, splits, reps, tm = _assemble_inputs(tmp_path)
    out = tmp_path / "dataset"
    _assemble(positives, negatives, splits, reps, tm, out)
    with pytest.raises(SystemExit, match="overwrite"):
        _assemble(positives, negatives, splits, reps, tm, out)


def test_assemble_same_file_other_spelling_no_warning(tmp_path, capsys):
    positives, negatives, splits, reps, tm = _assemble_inputs(tmp_path)
    (tmp_path / "sub").mkdir()
    spelled = tmp_path / "sub" / ".." / reps.name
    _assemble(positives, negatives, splits, spelled, tm, tmp_path / "d")
    assert "assemble: warning" not in capsys.readouterr().err
