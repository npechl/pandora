from pathlib import Path

import pytest

from pandora.export import structure_to_mmcif, write_json, write_records
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import ChainRecord

MMCIF_PATH = (
    Path(__file__).parent.parent / "datasets" / "dev" / "mmcif" / "1ayi.cif"
)


def test_structure_to_mmcif_round_trips_core_fields(tmp_path):
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))

    out_path = structure_to_mmcif(structure, tmp_path / "1ayi.out.cif")
    reparsed, _, status = mmcif_to_structure(str(out_path))

    assert status in ("success", "warning")
    assert reparsed.entry_id == structure.entry_id
    assert len(reparsed.atoms) == len(structure.atoms)
    assert {e.id for e in reparsed.entities} == {
        e.id for e in structure.entities
    }
    assert len(reparsed.assemblies) == len(structure.assemblies)


def test_write_json_round_trips(tmp_path):
    
    MMCIF_PATH = (
        Path(__file__).parent.parent / "datasets" / "dev" / "mmcif" / "1ayi.cif"
    )

    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))

    out_path = write_json(structure.entry, tmp_path / "entry.json")

    assert out_path.exists()
    loaded = type(structure.entry).model_validate_json(out_path.read_text())
    assert loaded == structure.entry


def test_write_records_jsonl_and_json(tmp_path):
    records = [
        ChainRecord(
            entry_id="1ayi",
            chain_id="A",
            entity_id="1",
            residue_count=10,
            atom_count=80,
        )
    ]

    jsonl_path = write_records(records, tmp_path / "chains.jsonl")
    json_path = write_records(records, tmp_path / "chains.json")

    assert jsonl_path.read_text().strip().count("\n") == 0
    reloaded = ChainRecord.model_validate_json(jsonl_path.read_text().strip())
    assert reloaded == records[0]

    import json

    assert json.loads(json_path.read_text()) == [records[0].model_dump()]


def test_write_records_rejects_empty_and_bad_suffix(tmp_path):
    with pytest.raises(ValueError):
        write_records([], tmp_path / "empty.json")

    record = ChainRecord(
        entry_id="1ayi",
        chain_id="A",
        entity_id="1",
        residue_count=1,
        atom_count=1,
    )
    with pytest.raises(ValueError):
        write_records([record], tmp_path / "chains.csv")


def test_write_records_parquet(tmp_path):
    pytest.importorskip("pandas")
    pytest.importorskip("pyarrow")
    pd = pytest.importorskip("pandas")

    record = ChainRecord(
        entry_id="1ayi",
        chain_id="A",
        entity_id="1",
        residue_count=1,
        atom_count=1,
    )
    out_path = write_records([record], tmp_path / "chains.parquet")

    df = pd.read_parquet(out_path)
    assert df.iloc[0]["chain_id"] == "A"
