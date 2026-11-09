from collections import Counter
from itertools import combinations
from pathlib import Path

import pytest
from pydantic import ValidationError

from pandora.annotations.entry import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    polymer_asym_ids,
)
from pandora.annotations.layers import check_layer_columns, resolve_layer_rows
from pandora.datasets import extract_interface_records
from pandora.parsing import mmcif_to_structure
from pandora.schemas.annotation import (
    LayerDescription,
    LayerRow,
    ResidueContactRow,
)
from pandora.schemas.dataset import InterfaceRecord
from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
    EntityRecord,
    EntryRecord,
    Structure,
)

MMCIF_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"
HYDROGENS = {"H", "D"}


def _atom(
    *,
    group_PDB="ATOM",
    id=1,
    label_asym_id="A",
    label_seq_id=1,
    label_comp_id="ALA",
    auth_seq_id="1",
    entity_id="1",
    atom_name="CA",
    element="C",
    alt_id=None,
    x=0.0,
    y=0.0,
    z=0.0,
) -> AtomSiteRecord:
    return AtomSiteRecord(
        group_PDB=group_PDB,
        id=id,
        type_symbol=element,
        label_atom_id=atom_name,
        label_alt_id=alt_id,
        label_comp_id=label_comp_id,
        label_asym_id=label_asym_id,
        label_entity_id=entity_id,
        label_seq_id=label_seq_id,
        Cartn_x=x,
        Cartn_y=y,
        Cartn_z=z,
        occupancy=1.0,
        B_iso_or_equiv=20.0,
        auth_seq_id=auth_seq_id,
        auth_comp_id=label_comp_id,
        auth_asym_id=label_asym_id,
        auth_atom_id=atom_name,
    )


def _structure(atoms, entities, asym_units) -> Structure:
    return Structure(
        entry_id="test",
        entry=EntryRecord(id="test"),
        entities=entities,
        asym_units=asym_units,
        atoms=atoms,
    )


def _load(entry_id: str) -> Structure:
    structure, _, _ = mmcif_to_structure(
        str(MMCIF_DIR / entry_id[1:3] / f"{entry_id}.cif")
    )
    return structure


def _reference_pairs(structure, cutoff, keep=lambda atom: True):
    """Check-every-atom-pair reference: {(chain_1, chain_2): {pairs}}."""

    chains = sorted(polymer_asym_ids(structure))
    by_chain = {
        chain: [
            a for a in structure.atoms if a.label_asym_id == chain and keep(a)
        ]
        for chain in chains
    }
    result = {}
    for chain_1, chain_2 in combinations(chains, 2):
        pairs = {
            (
                f"{a.label_asym_id}:{a.label_seq_id}",
                f"{b.label_asym_id}:{b.label_seq_id}",
            )
            for a in by_chain[chain_1]
            for b in by_chain[chain_2]
            if (a.Cartn_x - b.Cartn_x) ** 2
            + (a.Cartn_y - b.Cartn_y) ** 2
            + (a.Cartn_z - b.Cartn_z) ** 2
            <= cutoff * cutoff
        }
        if pairs:
            result[(chain_1, chain_2)] = pairs
    return result


def _layer_pairs(rows):
    """{(chain_1, chain_2): {residue pairs}} from chain-interface rows."""

    pairs = {}
    for r in rows:
        pairs.setdefault((r.chain_id_1, r.chain_id_2), set()).add(
            (f"{r.chain_id_1}:{r.seq_id_1}", f"{r.chain_id_2}:{r.seq_id_2}")
        )
    return pairs


def _two_chain_structure(atoms):
    """Two polymer chains A (entity 1) and B (entity 2) around atoms."""

    return _structure(
        atoms,
        [
            EntityRecord(id="1", type="polymer"),
            EntityRecord(id="2", type="polymer"),
        ],
        [AsymRecord(id="A", entity_id="1"), AsymRecord(id="B", entity_id="2")],
    )


def test_chain_interfaces_finds_contact_across_a_grid_cell_boundary():
    # Contact detection buckets atoms into cutoff-sized grid cells and
    # only checks the 27 cells around each atom. Place the contacting
    # pair straddling a cell boundary (cutoff=4.0 -> cell size 4.0) so a
    # same-cell-only search would miss it, plus an atom far enough away
    # that it must not show up as a contact.
    atoms = [
        _atom(id=1, label_asym_id="A", label_seq_id=1, entity_id="1", x=3.9),
        _atom(id=2, label_asym_id="B", label_seq_id=1, entity_id="2", x=4.1),
        _atom(id=3, label_asym_id="B", label_seq_id=2, entity_id="2", x=100.0),
    ]
    entities = [
        EntityRecord(id="1", type="polymer"),
        EntityRecord(id="2", type="polymer"),
    ]
    asym_units = [
        AsymRecord(id="A", entity_id="1"),
        AsymRecord(id="B", entity_id="2"),
    ]
    structure = _structure(atoms, entities, asym_units)

    rows = annotate_chain_interfaces(structure, distance_cutoff=4.0)

    assert _layer_pairs(rows) == {("A", "B"): {("A:1", "B:1")}}


def test_ligand_contacts_reports_nearest_polymer_residue():
    atoms = [
        _atom(id=1, label_asym_id="A", label_seq_id=1, entity_id="1", x=0.0),
        _atom(id=2, label_asym_id="A", label_seq_id=2, entity_id="1", x=10.0),
        _atom(
            id=3,
            group_PDB="HETATM",
            label_asym_id="B",
            label_seq_id=None,
            label_comp_id="ZN",
            auth_seq_id="101",
            entity_id="2",
            x=1.0,
        ),
    ]
    entities = [
        EntityRecord(id="1", type="polymer"),
        EntityRecord(id="2", type="non-polymer"),
    ]
    asym_units = [
        AsymRecord(id="A", entity_id="1"),
        AsymRecord(id="B", entity_id="2"),
    ]
    structure = _structure(atoms, entities, asym_units)

    (row,) = annotate_ligand_contacts(structure, distance_cutoff=4.0)

    assert (row.chain_id_1, row.seq_id_1) == ("A", 1)
    assert (row.chain_id_2, row.seq_id_2) == ("B", None)
    assert (row.auth_seq_id_2, row.comp_id_2) == ("101", "ZN")
    assert row.distance == 1.0


def test_ligand_contacts_fixture():
    rows = annotate_ligand_contacts(_load("1a2c"))

    assert {
        (r.chain_id_2, r.auth_chain_id_2, r.auth_seq_id_2, r.comp_id_2)
        for r in rows
    } == {("E", "H", "626", "NA")}
    assert len(rows) == 2


def test_ligand_contacts_hetatm_in_polymer_chain_is_not_a_ligand():
    # A: polymer with a HETATM modified residue (MSE) at seq 2.
    # B: a second polymer chain. C: a real ligand near A:2.
    atoms = [
        _atom(id=1, label_asym_id="A", label_seq_id=1, entity_id="1", x=0.0),
        _atom(
            id=2,
            group_PDB="HETATM",
            label_asym_id="A",
            label_seq_id=2,
            label_comp_id="MSE",
            entity_id="1",
            x=10.0,
        ),
        _atom(id=3, label_asym_id="B", label_seq_id=1, entity_id="2", x=11.0),
        _atom(
            id=4,
            group_PDB="HETATM",
            label_asym_id="C",
            label_seq_id=None,
            label_comp_id="ZN",
            auth_seq_id="101",
            entity_id="3",
            x=9.0,
        ),
    ]
    entities = [
        EntityRecord(id="1", type="polymer"),
        EntityRecord(id="2", type="polymer"),
        EntityRecord(id="3", type="non-polymer"),
    ]
    asym_units = [
        AsymRecord(id=c, entity_id=e)
        for c, e in (("A", "1"), ("B", "2"), ("C", "3"))
    ]

    rows = annotate_ligand_contacts(_structure(atoms, entities, asym_units))

    assert {r.chain_id_2 for r in rows} == {"C"}  # MSE is not a ligand
    assert ("A", 2) in {(r.chain_id_1, r.seq_id_1) for r in rows}
    assert extract_interface_records(rows) == []


def test_chain_interfaces_match_reference_on_fixture():
    structure = _load("1a7f")

    layer = annotate_chain_interfaces(
        structure, distance_cutoff=4.0, atom_set="all"
    )

    assert _layer_pairs(layer) == _reference_pairs(structure, 4.0)


def test_chain_interfaces_skip_far_chains_and_near_miss_boxes():
    # A-B touch. C is far from both (its box is skipped). D's box
    # overlaps A's grown box but no atom pair is within the cutoff.
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1", x=0.0),
        _atom(id=2, label_asym_id="A", label_seq_id=2, entity_id="1", y=10.0),
        _atom(id=3, label_asym_id="B", entity_id="2", x=3.0),
        _atom(id=4, label_asym_id="C", entity_id="3", x=100.0),
        _atom(id=5, label_asym_id="D", entity_id="4", x=3.0, y=6.0),
    ]
    entities = [EntityRecord(id=str(n), type="polymer") for n in range(1, 5)]
    asym_units = [
        AsymRecord(id=chain, entity_id=str(n))
        for n, chain in enumerate("ABCD", start=1)
    ]
    structure = _structure(atoms, entities, asym_units)

    layer = annotate_chain_interfaces(structure, distance_cutoff=4.0)

    assert _layer_pairs(layer) == _reference_pairs(structure, 4.0)
    assert _layer_pairs(layer) == {("A", "B"): {("A:1", "B:1")}}


def test_interface_records_residue_pairs_agree_with_sides():
    rows = annotate_chain_interfaces(_load("1a02"), distance_cutoff=4.0)
    records = extract_interface_records(rows)

    assert records
    for record in records:
        pairs = record.residue_pairs
        assert pairs == sorted(pairs)
        assert sorted({a for a, _ in pairs}) == (
            record.interface_residues_chain_1
        )
        assert sorted({b for _, b in pairs}) == (
            record.interface_residues_chain_2
        )


def test_chain_interfaces_altloc_atoms_give_one_residue_pair():
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1", alt_id="A", x=0.0),
        _atom(id=2, label_asym_id="A", entity_id="1", alt_id="B", x=0.5),
        _atom(id=3, label_asym_id="B", entity_id="2", x=3.0),
    ]

    rows = annotate_chain_interfaces(_two_chain_structure(atoms))

    assert len(rows) == 1
    assert _layer_pairs(rows) == {("A", "B"): {("A:1", "B:1")}}


@pytest.mark.parametrize(
    ("atom_set", "keep"),
    [
        ("all", lambda a: True),
        ("heavy", lambda a: a.type_symbol.upper() not in HYDROGENS),
        ("backbone", lambda a: a.label_atom_id in {"N", "CA", "C", "O"}),
        ("ca", lambda a: a.label_atom_id == "CA"),
    ],
)
def test_chain_interfaces_atom_sets_match_reference(atom_set, keep):
    structure = _load("1a7f")  # two protein chains, explicit hydrogens

    layer = annotate_chain_interfaces(structure, atom_set=atom_set)

    assert _layer_pairs(layer) == _reference_pairs(structure, 4.0, keep)


def test_chain_interfaces_heavy_is_subset_of_all():
    structure = _load("1a7f")

    heavy = _layer_pairs(annotate_chain_interfaces(structure))
    every = _layer_pairs(annotate_chain_interfaces(structure, atom_set="all"))

    assert heavy
    for chain_pair, pairs in heavy.items():
        assert pairs <= every[chain_pair]


def test_chain_interfaces_nucleic_backbone_matches_reference():
    structure = _load("1a02")  # DNA chains A, B; protein chains C, D, E
    nucleic = {"O5'", "C5'", "C4'", "C3'", "O3'"}
    protein = {"N", "CA", "C", "O"}

    layer = annotate_chain_interfaces(structure, atom_set="backbone")

    def keep(atom):
        names = nucleic if atom.label_asym_id in {"A", "B"} else protein
        return atom.label_atom_id in names

    assert _layer_pairs(layer) == _reference_pairs(structure, 4.0, keep)


def test_chain_interfaces_backbone_and_ca_ignore_side_chains():
    # Only a side-chain CB of A is within 4 A of B's CA.
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1", x=0.0),
        _atom(id=2, label_asym_id="A", entity_id="1", atom_name="CB", x=10.0),
        _atom(id=3, label_asym_id="B", entity_id="2", x=13.0),
    ]
    structure = _two_chain_structure(atoms)

    assert annotate_chain_interfaces(structure)
    for atom_set in ("backbone", "ca"):
        assert annotate_chain_interfaces(structure, atom_set=atom_set) == []


def test_chain_interfaces_ca_ignores_other_backbone_atoms():
    # Only the N atoms are within 4 A; the CAs are 20 A apart.
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1", x=0.0),
        _atom(
            id=2,
            label_asym_id="A",
            entity_id="1",
            atom_name="N",
            element="N",
            x=10.0,
        ),
        _atom(
            id=3,
            label_asym_id="B",
            entity_id="2",
            atom_name="N",
            element="N",
            x=13.0,
        ),
        _atom(id=4, label_asym_id="B", entity_id="2", x=20.0),
    ]
    structure = _two_chain_structure(atoms)

    backbone = annotate_chain_interfaces(structure, atom_set="backbone")
    ca = annotate_chain_interfaces(structure, atom_set="ca")

    assert _layer_pairs(backbone) == {("A", "B"): {("A:1", "B:1")}}
    assert ca == []


def test_chain_interfaces_polymer_types_filters_chains():
    structure = _load("1a02")

    layer = annotate_chain_interfaces(
        structure, polymer_types=["polypeptide(L)"]
    )

    chains = {chain for r in layer for chain in (r.chain_id_1, r.chain_id_2)}
    assert chains == {"C", "D", "E"}


def test_chain_interfaces_empty_polymer_types_means_no_chains():
    layer = annotate_chain_interfaces(_load("1a02"), polymer_types=[])

    assert layer == []


def test_chain_interfaces_ca_skips_chains_without_selected_atoms():
    # DNA chains have no CA atoms; they must be skipped, not crash.
    layer = annotate_chain_interfaces(_load("1a02"), atom_set="ca")

    for row in layer:
        assert row.chain_id_1 not in {"A", "B"}
        assert row.chain_id_2 not in {"A", "B"}


def test_chain_interfaces_rejects_unknown_atom_set():
    with pytest.raises(ValueError, match="atom_set"):
        annotate_chain_interfaces(_load("1a7f"), atom_set="CA")


def test_chain_interfaces_rows_on_fixture():
    rows = annotate_chain_interfaces(_load("1a2c"))

    counts = Counter((r.chain_id_1, r.chain_id_2) for r in rows)
    assert counts == {("A", "B"): 74, ("B", "C"): 21, ("B", "D"): 25}
    first = rows[0]
    assert (first.auth_chain_id_1, first.auth_chain_id_2) == ("L", "H")
    assert first.comp_id_1 and first.auth_seq_id_1
    assert all(0 < r.distance <= 4.0 for r in rows)


def test_chain_interfaces_distance_is_closest_atom_pair():
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1", x=0.0),
        _atom(id=2, label_asym_id="A", entity_id="1", atom_name="CB", x=1.0),
        _atom(id=3, label_asym_id="B", entity_id="2", x=3.5),
    ]

    rows = annotate_chain_interfaces(_two_chain_structure(atoms))

    assert [r.distance for r in rows] == [2.5]


def test_interface_records_from_rows_match_fixture():
    rows = annotate_chain_interfaces(_load("1a2c"))

    records = extract_interface_records(
        rows, distance_cutoff=4.0, atom_set="heavy"
    )

    assert [
        (r.chain_id_1, r.chain_id_2, len(r.residue_pairs)) for r in records
    ] == [("A", "B", 74), ("B", "C", 21), ("B", "D", 25)]
    assert (records[0].distance_cutoff, records[0].atom_set) == (4.0, "heavy")
    assert records[0].contact_count == len(
        records[0].interface_residues_chain_1
    ) + len(records[0].interface_residues_chain_2)


def test_interface_records_unordered_repeats_ligands_same_chain():
    base = {"entry_id": "X"}
    rows = [
        ResidueContactRow(
            **base,
            chain_id_1="A",
            seq_id_1=5,
            chain_id_2="B",
            seq_id_2=7,
            contact_type="hbond",
        ),
        ResidueContactRow(
            **base,
            chain_id_1="B",
            seq_id_1=7,
            chain_id_2="A",
            seq_id_2=5,
            contact_type="hydrophobic",
        ),
        ResidueContactRow(  # ligand
            **base, chain_id_1="A", seq_id_1=5, chain_id_2="E", seq_id_2=None
        ),
        ResidueContactRow(  # same chain
            **base, chain_id_1="A", seq_id_1=5, chain_id_2="A", seq_id_2=9
        ),
    ]

    (record,) = extract_interface_records(rows)

    assert (record.chain_id_1, record.chain_id_2) == ("A", "B")
    assert record.residue_pairs == [("A:5", "B:7")]
    assert record.contact_count == 2
    assert record.distance_cutoff is None
    assert record.atom_set is None


def test_interface_record_without_atom_set_loads_as_all():
    # Records exported before atom_set existed were all-atom contacts.
    old = {
        "entry_id": "1A3N",
        "chain_id_1": "A",
        "chain_id_2": "B",
        "distance_cutoff": 4.0,
        "interface_residues_chain_1": ["A:1"],
        "interface_residues_chain_2": ["B:1"],
        "contact_count": 2,
    }

    record = InterfaceRecord.model_validate(old)

    assert record.atom_set == "all"
    assert record.residue_pairs == []
    with pytest.raises(ValidationError):
        InterfaceRecord.model_validate({**old, "atom_set": "bogus"})


def test_chain_interfaces_rejects_string_polymer_types():
    with pytest.raises(TypeError, match="polymer_types"):
        annotate_chain_interfaces(_load("1a02"), polymer_types="polypeptide(L)")


# Layer templates -------------------------------------------------------


def test_layer_description_rejects_unknown_template_and_fields():
    with pytest.raises(ValidationError):
        LayerDescription(name="x", template="nope", tool="t")
    with pytest.raises(ValidationError):
        LayerDescription(name="x", template="free_form", tool="t", extra=1)


def test_layer_row_keeps_any_columns():
    assert LayerRow(a=1, b="x").model_dump() == {"a": 1, "b": "x"}


def test_check_layer_columns_accepts_either_naming():
    check_layer_columns(
        {"entry_id", "chain_id", "seq_id", "value"}, "residue_values"
    )
    check_layer_columns(
        {"entry_id", "auth_chain_id", "auth_seq_id", "value"},
        "residue_values",
    )
    check_layer_columns(
        {
            "entry_id",
            "chain_id_1",
            "seq_id_1",
            "auth_chain_id_2",
            "auth_seq_id_2",
        },
        "residue_contacts",
    )
    check_layer_columns({"entry_id"}, "free_form")


@pytest.mark.parametrize(
    "columns, template, missing",
    [
        ({"chain_id", "seq_id", "value"}, "residue_values", "entry_id"),
        ({"entry_id", "chain_id", "seq_id"}, "residue_values", "value"),
        ({"entry_id", "chain_id", "label"}, "residue_labels", "seq_id"),
        ({"entry_id", "chain_id_1", "seq_id_1"}, "residue_contacts", "_2"),
        (
            {"entry_id", "chain_id_1", "chain_id_2"},
            "interface_labels",
            "label",
        ),
    ],
)
def test_check_layer_columns_rejects_missing(columns, template, missing):
    with pytest.raises(ValueError, match=missing):
        check_layer_columns(columns, template)


_1A2C = mmcif_to_structure(str(MMCIF_DIR / "a2" / "1a2c.cif"))[0]
S = {"1A2C": _1A2C}


def test_resolve_author_naming_fills_pandora_keys():
    m, u = resolve_layer_rows(
        [
            {
                "entry_id": "1a2c",
                "auth_chain_id": "L",
                "auth_seq_id": "1",
                "ins_code": "G",
                "value": 0.5,
            }
        ],
        "residue_values",
        S,
    )
    assert u == []
    assert m == [
        {
            "entry_id": "1A2C",
            "chain_id": "A",
            "seq_id": 2,
            "auth_chain_id": "L",
            "auth_seq_id": "1",
            "ins_code": "G",
            "comp_id": "PHE",
            "value": 0.5,
        }
    ]


def test_resolve_pandora_naming_fills_author_keys():
    m, _ = resolve_layer_rows(
        [{"entry_id": "1A2C", "chain_id": "A", "seq_id": 1, "value": 1.0}],
        "residue_values",
        S,
    )
    row = m[0]
    assert (
        row["auth_chain_id"],
        row["auth_seq_id"],
        row["ins_code"],
        row["comp_id"],
    ) == ("L", "1", "H", "THR")


def test_resolve_parquet_floats_and_ins_code_spellings():
    rows = [
        {"entry_id": "1A2C", "chain_id": "A", "seq_id": 3.0, "value": 1},
        {
            "entry_id": "1A2C",
            "auth_chain_id": "L",
            "auth_seq_id": 471.0,
            "ins_code": "?",
            "value": 1,
        },
    ]
    m, u = resolve_layer_rows(rows, "residue_values", S)
    assert u == [] and m[0]["seq_id"] == 3
    assert (m[1]["chain_id"], m[1]["seq_id"]) == ("F", None)


@pytest.mark.parametrize(
    "row, reason",
    [
        ({"entry_id": "9XYZ", "chain_id": "A", "seq_id": 1}, "entry_not_found"),
        ({"entry_id": "1A2C", "chain_id": "A"}, "missing_key"),
        ({"entry_id": "1A2C", "chain_id": "Z", "seq_id": 1}, "chain_not_found"),
        (
            {"entry_id": "1A2C", "chain_id": "A", "seq_id": 9999},
            "residue_not_found",
        ),
        (
            {
                "entry_id": "1A2C",
                "auth_chain_id": "H",
                "auth_seq_id": "147",
                "ins_code": "A",
            },
            "no_coordinates",
        ),
        (
            {"entry_id": "1A2C", "chain_id": "B", "seq_id": 148},
            "no_coordinates",
        ),
        (
            {
                "entry_id": "1A2C",
                "chain_id": "A",
                "seq_id": 1,
                "auth_chain_id": "L",
                "auth_seq_id": "1",
                "ins_code": "G",
            },
            "naming_conflict",
        ),
        ({"entry_id": "1A2C", "chain_id": "A", "seq_id": "x"}, "bad_value"),
    ],
)
def test_resolve_unmatched_reasons(row, reason):
    m, u = resolve_layer_rows([{**row, "label": "x"}], "residue_labels", S)
    assert m == [] and u == [{**row, "label": "x", "reason": reason}]


def test_resolve_ambiguous_author_id():
    # strip insertion codes, as preserve_insertion_codes=False does (#125)
    atoms = [
        a.model_copy(update={"pdbx_PDB_ins_code": None}) for a in _1A2C.atoms
    ]
    s = _1A2C.model_copy(update={"atoms": atoms})
    _, u = resolve_layer_rows(
        [
            {
                "entry_id": "1A2C",
                "auth_chain_id": "L",
                "auth_seq_id": "1",
                "label": "x",
            }
        ],
        "residue_labels",
        {"1A2C": s},
    )
    assert u[0]["reason"] == "ambiguous_author_id"


def test_resolve_does_not_mutate_rows():
    row = {"entry_id": "1a2c", "chain_id": "A", "seq_id": 1, "label": "x"}
    resolve_layer_rows([row], "residue_labels", S)
    assert row == {
        "entry_id": "1a2c",
        "chain_id": "A",
        "seq_id": 1,
        "label": "x",
    }


def test_resolve_contacts_both_sides_and_repeats_kept():
    row = {
        "entry_id": "1A2C",
        "chain_id_1": "A",
        "seq_id_1": 1,
        "auth_chain_id_2": "L",
        "auth_seq_id_2": "1",
        "ins_code_2": "G",
        "contact_type": "hbond",
    }
    m, u = resolve_layer_rows(
        [row, {**row, "contact_type": "vdw"}], "residue_contacts", S
    )
    assert u == [] and len(m) == 2
    assert (m[0]["auth_chain_id_1"], m[0]["chain_id_2"], m[0]["seq_id_2"]) == (
        "L",
        "A",
        2,
    )


def test_resolve_contacts_reports_first_failing_side():
    _, u = resolve_layer_rows(
        [
            {
                "entry_id": "1A2C",
                "chain_id_1": "A",
                "seq_id_1": 1,
                "chain_id_2": "Z",
                "seq_id_2": 1,
            }
        ],
        "residue_contacts",
        S,
    )
    assert u[0]["reason"] == "chain_not_found"


def test_resolve_chain_author_naming_picks_polymer():
    m, u = resolve_layer_rows(
        [{"entry_id": "1A2C", "auth_chain_id": "L", "label": "light"}],
        "chain_labels",
        S,
    )
    assert u == [] and m[0]["chain_id"] == "A"


def test_resolve_interface_labels():
    m, _ = resolve_layer_rows(
        [
            {
                "entry_id": "1A2C",
                "chain_id_1": "A",
                "chain_id_2": "B",
                "label": "physiological",
            }
        ],
        "interface_labels",
        S,
    )
    assert (m[0]["auth_chain_id_1"], m[0]["auth_chain_id_2"]) == ("L", "H")


@pytest.mark.parametrize("value", [None, float("nan"), True, "high"])
def test_resolve_bad_value(value):
    _, u = resolve_layer_rows(
        [{"entry_id": "1A2C", "chain_id": "A", "value": value}],
        "chain_values",
        S,
    )
    assert u[0]["reason"] == "bad_value"


def test_resolve_duplicate_key_keeps_first():
    rows = [
        {"entry_id": "1A2C", "chain_id": "A", "seq_id": 1, "label": "a"},
        {
            "entry_id": "1A2C",
            "auth_chain_id": "L",
            "auth_seq_id": "1",
            "ins_code": "H",
            "label": "b",
        },
    ]
    m, u = resolve_layer_rows(rows, "residue_labels", S)
    assert [r["label"] for r in m] == ["a"]
    assert u[0]["reason"] == "duplicate_key"


def test_resolve_free_form():
    rows = [
        {"entry_id": "1a2c", "note": "x"},
        {"entry_id": "1a2c", "note": "y"},
        {"entry_id": "1A2C", "chain_id": "Z", "note": "z"},
    ]
    m, u = resolve_layer_rows(rows, "free_form", S)
    assert [r["note"] for r in m] == ["x", "y"]
    assert u[0]["reason"] == "chain_not_found"


def test_resolve_no_coordinates_without_poly_seq_scheme():
    # canonicalisation drops _pdbx_poly_seq_scheme when residues change;
    # _entity_poly_seq still says B 148 is in the sequence
    raw = {k: v for k, v in _1A2C.raw.items() if k != "_pdbx_poly_seq_scheme"}
    s = _1A2C.model_copy(update={"raw": raw})
    _, u = resolve_layer_rows(
        [{"entry_id": "1A2C", "chain_id": "B", "seq_id": 148, "label": "x"}],
        "residue_labels",
        {"1A2C": s},
    )
    assert u[0]["reason"] == "no_coordinates"


def test_residue_contact_row_reads_parquet_shapes():
    row = ResidueContactRow.model_validate(
        {
            "entry_id": "1A2C",
            "chain_id_1": "A",
            "seq_id_1": 42.0,
            "ins_code_1": "",
            "chain_id_2": "E",
            "seq_id_2": None,
            "ins_code_2": "?",
            "distance": "3.5",
            "tool_score": 7,
        }
    )
    assert row.seq_id_1 == 42
    assert row.ins_code_1 is None and row.ins_code_2 is None
    assert row.distance == 3.5
    assert row.model_dump()["tool_score"] == 7


def test_resolve_fills_comp_id():
    matched, _ = resolve_layer_rows(
        [
            {
                "entry_id": "1A2C",
                "chain_id_1": "A",
                "seq_id_1": 1,
                "chain_id_2": "B",
                "seq_id_2": 1,
            }
        ],
        "residue_contacts",
        S,
    )
    assert (matched[0]["comp_id_1"], matched[0]["comp_id_2"]) == (
        "THR",
        "ILE",
    )


@pytest.mark.parametrize("distance", ["near", True])
def test_resolve_contacts_bad_distance(distance):
    _, unmatched = resolve_layer_rows(
        [
            {
                "entry_id": "1A2C",
                "chain_id_1": "A",
                "seq_id_1": 1,
                "chain_id_2": "B",
                "seq_id_2": 1,
                "distance": distance,
            }
        ],
        "residue_contacts",
        S,
    )
    assert unmatched[0]["reason"] == "bad_value"
