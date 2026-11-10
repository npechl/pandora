from collections import Counter
from itertools import combinations
from pathlib import Path

import pytest
from pydantic import ValidationError

from pandora.annotations.counts import (
    count_chain_contacts,
    count_residue_contacts,
    label_interfaces,
)
from pandora.annotations.entry import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    polymer_asym_ids,
)
from pandora.annotations.filters import (
    count_interface_labels,
    drop_labelled_interfaces,
)
from pandora.annotations.layers import check_layer_columns, resolve_layer_rows
from pandora.datasets import extract_interface_records
from pandora.parsing import mmcif_to_structure
from pandora.schemas.annotation import (
    LayerDescription,
    LayerRow,
    ResidueContactRow,
    SimilarityRow,
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


def test_resolve_microheterogeneity_is_one_residue():
    # A:5 has altloc A = SER and altloc B = THR at one position.
    atoms = [
        _atom(
            id=1,
            label_asym_id="A",
            label_seq_id=5,
            auth_seq_id="5",
            label_comp_id="SER",
            alt_id="A",
            entity_id="1",
        ),
        _atom(
            id=2,
            label_asym_id="A",
            label_seq_id=5,
            auth_seq_id="5",
            label_comp_id="THR",
            alt_id="B",
            entity_id="1",
        ),
        _atom(
            id=3,
            label_asym_id="B",
            label_seq_id=7,
            auth_seq_id="7",
            entity_id="2",
            x=3.0,
        ),
    ]
    structure = _structure(
        atoms,
        [
            EntityRecord(id="1", type="polymer"),
            EntityRecord(id="2", type="polymer"),
        ],
        [
            AsymRecord(id="A", entity_id="1", auth_id="A"),
            AsymRecord(id="B", entity_id="2", auth_id="B"),
        ],
    )
    row = {
        "entry_id": "test",
        "auth_chain_id_1": "A",
        "auth_seq_id_1": "5",
        "auth_chain_id_2": "B",
        "auth_seq_id_2": "7",
    }

    matched, unmatched = resolve_layer_rows(
        [row], "residue_contacts", {"test": structure}
    )

    assert unmatched == []
    assert matched[0]["comp_id_1"] == "SER"  # first atom, as the built-ins
    (built_in,) = annotate_chain_interfaces(structure)
    assert built_in.comp_id_1 == "SER"


def test_ligand_contacts_keep_insertion_codes_apart():
    # Waters 301 and 301A in one chain, each near a different residue.
    atoms = [
        _atom(id=1, label_asym_id="A", label_seq_id=1, entity_id="1", x=0.0),
        _atom(id=2, label_asym_id="A", label_seq_id=2, entity_id="1", x=50.0),
        _atom(
            id=3,
            group_PDB="HETATM",
            label_asym_id="W",
            label_seq_id=None,
            label_comp_id="HOH",
            auth_seq_id="301",
            entity_id="2",
            x=1.0,
        ),
        _atom(
            id=4,
            group_PDB="HETATM",
            label_asym_id="W",
            label_seq_id=None,
            label_comp_id="HOH",
            auth_seq_id="301",
            entity_id="2",
            x=51.0,
        ),
    ]
    atoms[3] = atoms[3].model_copy(update={"pdbx_PDB_ins_code": "A"})
    structure = _structure(
        atoms,
        [
            EntityRecord(id="1", type="polymer"),
            EntityRecord(id="2", type="water"),
        ],
        [AsymRecord(id="A", entity_id="1"), AsymRecord(id="W", entity_id="2")],
    )

    rows = annotate_ligand_contacts(structure, include_waters=True)

    assert {(r.seq_id_1, r.auth_seq_id_2, r.ins_code_2) for r in rows} == {
        (1, "301", None),
        (2, "301", "A"),
    }


def test_similarity_description_needs_items():
    with pytest.raises(ValidationError, match="items"):
        LayerDescription(name="g", template="similarity", tool="R")
    description = LayerDescription(
        name="g", template="similarity", tool="R", items="chain"
    )
    assert description.items == "chain"


def test_items_rejected_for_other_templates():
    with pytest.raises(ValidationError, match="items"):
        LayerDescription(
            name="c", template="residue_contacts", tool="x", items="chain"
        )


def test_similarity_row_coerces_int_score():
    row = SimilarityRow(
        item_id_1="1ABC_A", item_id_2="2XYZ_B", score=3, extra_col="k"
    )
    assert row.score == 3.0
    assert isinstance(row.score, float)
    assert row.model_dump()["extra_col"] == "k"


def test_annotate_import_refuses_similarity_template():
    with pytest.raises(ValueError, match="pandora similarity import"):
        check_layer_columns({"item_id_1", "item_id_2", "score"}, "similarity")


def _contact(c1, s1, c2, s2, **kw) -> ResidueContactRow:
    fields = {"entry_id": "1abc"}
    for side, chain, seq in (("1", c1, s1), ("2", c2, s2)):
        fields |= {
            f"chain_id_{side}": chain,
            f"seq_id_{side}": seq,
            f"auth_chain_id_{side}": chain,
            f"auth_seq_id_{side}": str(seq),
            f"comp_id_{side}": "ALA",
        }
    return ResidueContactRow(**(fields | kw))


def test_count_residue_contacts_counts_distinct_partners():
    rows = [
        _contact("A", 1, "B", 5, contact_type="hbond"),
        _contact("A", 1, "B", 5, contact_type="salt_bridge"),
        _contact("B", 5, "A", 1),
        _contact("A", 1, "B", 6),
    ]
    out = {
        (r.chain_id, r.seq_id): r.value for r in count_residue_contacts(rows)
    }
    assert out == {("A", 1): 2, ("B", 5): 1, ("B", 6): 1}


def test_count_residue_contacts_keeps_both_namings():
    row = count_residue_contacts([_contact("A", 1, "B", 5)])[0]
    assert row.model_dump() == {
        "entry_id": "1abc",
        "chain_id": "A",
        "seq_id": 1,
        "auth_chain_id": "A",
        "auth_seq_id": "1",
        "ins_code": None,
        "comp_id": "ALA",
        "value": 1,
    }


def test_count_residue_contacts_ligands_keyed_by_author_naming():
    rows = [
        _contact("A", 1, "C", None, auth_seq_id_2="301"),
        _contact("A", 2, "C", None, auth_seq_id_2="302"),
    ]
    ligands = [r for r in count_residue_contacts(rows) if r.chain_id == "C"]
    assert [(r.seq_id, r.auth_seq_id, r.value) for r in ligands] == [
        (None, "301", 1),
        (None, "302", 1),
    ]


def test_count_residue_contacts_self_pair_adds_nothing():
    assert count_residue_contacts([_contact("A", 1, "A", 1)]) == []
    assert count_chain_contacts([_contact("A", 1, "A", 1)]) == []


def test_count_chain_contacts_counts_distinct_residues():
    rows = [
        _contact("A", 1, "B", 5),
        _contact("A", 1, "B", 6),
        _contact("A", 2, "B", 5),
    ]
    out = {r.chain_id: r.value for r in count_chain_contacts(rows)}
    assert out == {"A": 2, "B": 2}
    assert count_chain_contacts(rows)[0].model_dump() == {
        "entry_id": "1abc",
        "chain_id": "A",
        "auth_chain_id": "A",
        "value": 2,
    }


def test_count_layers_pass_their_template_check():
    rows = [_contact("A", 1, "B", 5)]
    check_layer_columns(
        count_residue_contacts(rows)[0].model_dump(), "residue_values"
    )
    check_layer_columns(
        count_chain_contacts(rows)[0].model_dump(), "chain_values"
    )


def test_count_functions_do_not_mutate_input():
    rows = [_contact("A", 1, "B", 5)]
    before = [r.model_copy(deep=True) for r in rows]
    count_residue_contacts(rows)
    count_chain_contacts(rows)
    assert rows == before


def _label(c1, c2, label, entry_id="1abc") -> LayerRow:
    return LayerRow(
        entry_id=entry_id, chain_id_1=c1, chain_id_2=c2, label=label
    )


def _pairs(rows):
    return [(r.chain_id_1, r.seq_id_1, r.chain_id_2, r.seq_id_2) for r in rows]


def test_drop_labelled_interfaces_drops_both_orders():
    rows = [_contact("A", 1, "B", 5), _contact("B", 6, "A", 2)]
    rows.append(_contact("A", 1, "C", 3))
    kept = drop_labelled_interfaces(
        rows, [_label("B", "A", "low_hbond")], {"low_hbond"}
    )
    assert _pairs(kept) == [("A", 1, "C", 3)]


def test_drop_labelled_interfaces_keeps_order_and_unlabelled():
    rows = [
        _contact("A", 1, "C", 1),
        _contact("A", 2, "B", 2),
        _contact("B", 3, "C", 3),
    ]
    kept = drop_labelled_interfaces(rows, [_label("A", "B", "bad")], {"bad"})
    assert _pairs(kept) == [("A", 1, "C", 1), ("B", 3, "C", 3)]


def test_drop_labelled_interfaces_any_label_drops():
    rows = [_contact("A", 1, "B", 2), _contact("A", 1, "C", 2)]
    labels = [_label("A", "B", "ok"), _label("B", "A", "bad")]
    kept = drop_labelled_interfaces(rows, labels, {"bad"})
    assert _pairs(kept) == [("A", 1, "C", 2)]


def test_drop_labelled_interfaces_self_pair_and_ligand():
    rows = [_contact("A", 1, "A", 9), _contact("A", 1, "J", None)]
    labels = [_label("A", "A", "x"), _label("A", "J", "y")]
    assert _pairs(drop_labelled_interfaces(rows, labels, {"x"})) == [
        ("A", 1, "J", None)
    ]
    assert _pairs(drop_labelled_interfaces(rows, labels, {"y"})) == [
        ("A", 1, "A", 9)
    ]


def test_drop_labelled_interfaces_entry_case():
    rows = [_contact("A", 1, "B", 2), _contact("A", 1, "C", 2)]
    labels = [_label("A", "B", "bad", entry_id="1ABC")]
    kept = drop_labelled_interfaces(rows, labels, {"bad"})
    assert _pairs(kept) == [("A", 1, "C", 2)]


def test_drop_labelled_interfaces_unknown_label_raises():
    rows = [_contact("A", 1, "B", 2)]
    with pytest.raises(ValueError, match="typo"):
        drop_labelled_interfaces(
            rows, [_label("A", "B", "bad")], {"bad", "typo"}
        )


def test_count_interface_labels():
    rows = [
        _contact("A", 1, "B", 1),
        _contact("B", 2, "A", 2),
        _contact("A", 1, "C", 1),
        _contact("B", 1, "C", 1),
    ]
    labels = [
        _label("A", "B", "bad"),
        _label("A", "C", "ok"),
        _label("C", "D", "bad"),
    ]
    assert count_interface_labels(rows, labels, {"bad"}) == {
        "interfaces_in": 3,
        "interfaces_dropped": 1,
        "interfaces_unlabelled": 1,
        "labels_unmatched": 1,
    }


def test_filter_functions_do_not_mutate_input():
    rows = [_contact("A", 1, "B", 2), _contact("A", 1, "C", 2)]
    labels = [_label("A", "B", "bad")]
    before = [r.model_dump() for r in rows + labels]
    drop_labelled_interfaces(rows, labels, {"bad"})
    count_interface_labels(rows, labels, {"bad"})
    assert [r.model_dump() for r in rows + labels] == before


def test_drop_labelled_interfaces_accepts_resolved_dicts():
    rows = [_contact("A", 1, "B", 2), _contact("A", 1, "C", 2)]
    labels = [_label("A", "B", "bad").model_dump()]
    kept = drop_labelled_interfaces(rows, labels, {"bad"})
    assert _pairs(kept) == [("A", 1, "C", 2)]
    assert (
        count_interface_labels(rows, labels, {"bad"})["interfaces_dropped"] == 1
    )


def _lengths(**chains):
    return {("1ABC", c): n for c, n in chains.items()}


def test_label_interfaces_counts_each_side():
    rows = [
        _contact("B", 1, "A", 1),
        _contact("A", 2, "B", 1),
        _contact("A", 2, "B", 2),
        _contact("A", 1, "B", 1),  # repeat of the first, other order
    ]
    (row,) = label_interfaces(rows, _lengths(A=10, B=20))
    assert row.model_dump() == {
        "entry_id": "1abc",
        "chain_id_1": "A",
        "chain_id_2": "B",
        "auth_chain_id_1": "A",
        "auth_chain_id_2": "B",
        "residues_1": 2,
        "residues_2": 2,
        "length_1": 10,
        "length_2": 20,
        "total_residues": 30,
        "label": "ok",
    }


@pytest.mark.parametrize(
    "min_res, max_total, label",
    [
        (None, None, "ok"),
        (3, None, "few_contacts"),
        (2, None, "ok"),
        (None, 29, "too_large"),
        (None, 30, "ok"),
        (3, 29, "too_large"),  # both apply: too_large wins
    ],
)
def test_label_interfaces_rules(min_res, max_total, label):
    rows = [_contact("A", 1, "B", 1), _contact("A", 2, "B", 2)]
    (row,) = label_interfaces(rows, _lengths(A=10, B=20), min_res, max_total)
    assert row.label == label


def test_label_interfaces_unknown_length_is_never_too_large():
    rows = [_contact("A", 1, "B", 1)]
    (row,) = label_interfaces(rows, _lengths(A=10), max_total_residues=1)
    assert (row.length_2, row.total_residues, row.label) == (None, None, "ok")


def test_label_interfaces_entry_case_insensitive_and_sorted():
    rows = [_contact("C", 1, "B", 1), _contact("A", 1, "B", 1)]
    out = label_interfaces(rows, _lengths(A=1, B=1, C=1))
    assert [(r.chain_id_1, r.chain_id_2) for r in out] == [
        ("A", "B"),
        ("B", "C"),
    ]
    assert all(r.total_residues == 2 for r in out)


def test_label_interfaces_passes_template_check_and_keeps_input():
    rows = [_contact("A", 1, "B", 1)]
    before = [r.model_copy() for r in rows]
    out = label_interfaces(rows, _lengths(A=1, B=1))
    check_layer_columns(out[0].model_dump().keys(), "interface_labels")
    assert rows == before


def test_label_interfaces_orders_reversed_pair_and_its_counts():
    rows = [_contact("B", s, "A", 1) for s in (1, 2, 3)]
    (row,) = label_interfaces(rows, _lengths(A=10, B=20))
    assert (row.chain_id_1, row.residues_1, row.length_1) == ("A", 1, 10)
    assert (row.chain_id_2, row.residues_2, row.length_2) == ("B", 3, 20)
