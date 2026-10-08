from pathlib import Path

import pytest

from pandora.canonicalisation import canonicalise_structure
from pandora.export import (
    export_chain_mmcif,
    structure_to_mmcif,
    write_json,
    write_records,
)
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    AssemblyRules,
    IncompleteChainRules,
    MissingDataRules,
    canonicalisationPolicy,
)
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


def test_export_chain_mmcif_keeps_only_that_chain(tmp_path):
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    all_chain_ids = {a.label_asym_id for a in structure.atoms}
    assert "B" in all_chain_ids and len(all_chain_ids) > 1

    out_path = export_chain_mmcif(structure, "B", tmp_path / "1ayi_B.cif")
    reparsed, _, status = mmcif_to_structure(str(out_path))

    assert status in ("success", "warning")
    assert {a.label_asym_id for a in reparsed.atoms} == {"B"}
    assert len(reparsed.atoms) == sum(
        1 for a in structure.atoms if a.label_asym_id == "B"
    )


def test_write_json_round_trips(tmp_path):

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


@pytest.mark.xfail(
    strict=True,
    reason="raw categories such as _pdbx_sifts_xref_db keep references to "
    "residues canonicalisation removed, so gemmi can't read the file back",
)
def test_truncated_structure_round_trips_through_mmcif(tmp_path):
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH.parent / "10mv.cif"))
    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        missing_data_rules=MissingDataRules(
            incomplete_chains=IncompleteChainRules(
                strategy="truncate_to_complete_regions"
            )
        ),
    )
    canonical, _, _ = canonicalise_structure(structure, policy)
    path = tmp_path / "10mv.cif"
    structure_to_mmcif(canonical, str(path))
    reparsed, _, status = mmcif_to_structure(str(path))
    assert reparsed is not None, status


@pytest.mark.xfail(
    strict=True,
    reason="raw categories such as _pdbx_sifts_xref_db keep references to "
    "chains outside the expanded assembly, so gemmi can't read the file back",
)
def test_expanded_assembly_round_trips_through_mmcif(tmp_path):
    # 1a08: assembly 1 is chains A, B (+ ligands); C, D, G, H are dropped.
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH.parent / "1a08.cif"))
    policy = canonicalisationPolicy(
        policy_id="p",
        policy_name="p",
        policy_version="1.0.0",
        assembly_rules=AssemblyRules(
            strategy="standardize_biological_assembly"
        ),
    )
    canonical, _, _ = canonicalise_structure(structure, policy)
    path = tmp_path / "1a08.cif"
    structure_to_mmcif(canonical, str(path))
    reparsed, _, status = mmcif_to_structure(str(path))
    assert reparsed is not None, status
