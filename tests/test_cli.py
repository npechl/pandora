import json
from pathlib import Path

from pandora.cli.app import _build_parser, main

FIXTURES_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"
CANON_POLICY = (
    Path(__file__).parent.parent / "datasets" / "canonicalisation.yaml"
)

_SUBCOMMANDS = [
    "fetch",
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
