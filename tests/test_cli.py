import gzip
import json
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

from pandora.cli._io import read_structures_dir, write_structures_dir
from pandora.cli.app import _build_parser, main
from pandora.datasets import extract_interface_records
from pandora.export import read_records
from pandora.ingestion import read_entry_paths
from pandora.parsing import mmcif_to_structure
from pandora.schemas.annotation import LayerRow, ResidueContactRow
from pandora.schemas.provenance import CanonicalisationEntryRecord

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
    assert provenance["hit_filter"]["min_score"] == 0.99
    assert provenance["hit_filter"]["min_coverage"] == 0.5
    assert provenance["similarity_layer"] == str(layer / "layer.json")
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
    assert rows["A_2"]["copy_of"] == "A"
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
    with pytest.raises(_Captured):
        main(
            ["cluster", "--similarity", str(layer), "--min-score", "0.5"]
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
