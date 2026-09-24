from pathlib import Path

import gemmi
import pytest

from pandora.datasets import extract_chain_records
from pandora.parsing import mmcif_to_structure

MMCIF_PATH = (
    Path(__file__).parent.parent / "datasets" / "dev" / "mmcif" / "104m.cif"
)


def test_parses_core_fields_against_gemmi():
    structure, diag, status = mmcif_to_structure(str(MMCIF_PATH))

    assert status == "success"
    assert not diag.errors and not diag.warnings
    assert structure.entry_id == "104M"
    expected = gemmi.read_structure(str(MMCIF_PATH))[0].count_atom_sites()
    assert len(structure.atoms) == expected
    assert [(e.id, e.type) for e in structure.entities] == [
        ("1", "polymer"),
        ("2", "non-polymer"),
        ("3", "non-polymer"),
        ("4", "non-polymer"),
        ("5", "water"),
    ]


def test_raw_keeps_untyped_categories_but_not_atom_site():
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))

    assert "_atom_site" not in structure.raw
    assert structure.raw["_exptl"][0]["method"] is not None
    assert len(structure.raw["_entity"]) == 5


def test_missing_model_falls_back_with_warning():
    structure, diag, status = mmcif_to_structure(str(MMCIF_PATH), model_num=99)

    assert status == "warning"
    assert [w.code for w in diag.warnings] == ["MODEL_NOT_FOUND"]
    assert structure is not None


@pytest.mark.parametrize(
    ("path", "code"), [("", "EMPTY_CONTENT"), ("missing.cif", "PARSE_ERROR")]
)
def test_unreadable_input_fails(path, code):
    structure, diag, status = mmcif_to_structure(path)

    assert structure is None
    assert status == "failed"
    assert diag.errors[0].code == code


@pytest.mark.xfail(
    strict=True, reason="_cs() keeps CIF quote delimiters in string values"
)
def test_quoted_values_are_unquoted():
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))

    assert structure.entry.title == (
        "SPERM WHALE MYOGLOBIN N-BUTYL ISOCYANIDE AT PH 7.0"
    )


@pytest.mark.xfail(
    strict=True, reason="_cs() keeps ';' text-field delimiters in sequences"
)
def test_chain_sequence_has_no_text_field_delimiters():
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))

    (chain,) = extract_chain_records(structure)
    assert chain.sequence.startswith("VLSEGEWQL")
    assert chain.sequence.endswith("GYQG")
