from itertools import combinations
from pathlib import Path

import pytest
from pydantic import ValidationError

from pandora.annotations.entry import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    polymer_asym_ids,
)
from pandora.datasets import extract_interface_records
from pandora.parsing import mmcif_to_structure
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
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / f"{entry_id}.cif"))
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


def _layer_pairs(layer):
    """{(chain_1, chain_2): {residue pairs}} from a chain_interfaces layer."""

    return {
        (i["chain_id_1"], i["chain_id_2"]): {
            tuple(pair) for pair in i["residue_pairs"]
        }
        for i in layer.data["interfaces"]
    }


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

    layer = annotate_chain_interfaces(structure, distance_cutoff=4.0)

    interfaces = layer.data["interfaces"]
    assert len(interfaces) == 1
    assert interfaces[0]["interface_residues_chain_1"] == ["A:1"]
    assert interfaces[0]["interface_residues_chain_2"] == ["B:1"]
    assert layer.parameters == {
        "distance_cutoff": 4.0,
        "atom_set": "heavy",
        "polymer_types": None,
    }


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

    layer = annotate_ligand_contacts(structure, distance_cutoff=4.0)

    ligands = layer.data["ligands"]
    assert len(ligands) == 1
    contacts = ligands[0]["contacts"]
    assert len(contacts) == 1
    assert contacts[0]["label_seq_id"] == 1
    assert contacts[0]["distance"] == 1.0
    assert layer.parameters == {
        "distance_cutoff": 4.0,
        "include_waters": False,
    }


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


def test_chain_interfaces_residue_pairs_agree_with_sides():
    layer = annotate_chain_interfaces(_load("1a02"), distance_cutoff=4.0)

    assert layer.data["interfaces"]
    for interface in layer.data["interfaces"]:
        pairs = interface["residue_pairs"]
        assert pairs == sorted(pairs)
        assert (
            sorted({a for a, _ in pairs})
            == (interface["interface_residues_chain_1"])
        )
        assert (
            sorted({b for _, b in pairs})
            == (interface["interface_residues_chain_2"])
        )


def test_chain_interfaces_altloc_atoms_give_one_residue_pair():
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1", alt_id="A", x=0.0),
        _atom(id=2, label_asym_id="A", entity_id="1", alt_id="B", x=0.5),
        _atom(id=3, label_asym_id="B", entity_id="2", x=3.0),
    ]

    layer = annotate_chain_interfaces(_two_chain_structure(atoms))

    assert layer.data["interfaces"][0]["residue_pairs"] == [["A:1", "B:1"]]


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

    assert annotate_chain_interfaces(structure).data["interfaces"]
    for atom_set in ("backbone", "ca"):
        layer = annotate_chain_interfaces(structure, atom_set=atom_set)
        assert layer.data["interfaces"] == []


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

    assert backbone.data["interfaces"][0]["residue_pairs"] == [["A:1", "B:1"]]
    assert ca.data["interfaces"] == []


def test_chain_interfaces_polymer_types_filters_chains():
    structure = _load("1a02")

    layer = annotate_chain_interfaces(
        structure, polymer_types=["polypeptide(L)"]
    )

    chains = {
        chain
        for i in layer.data["interfaces"]
        for chain in (i["chain_id_1"], i["chain_id_2"])
    }
    assert chains == {"C", "D", "E"}


def test_chain_interfaces_empty_polymer_types_means_no_chains():
    layer = annotate_chain_interfaces(_load("1a02"), polymer_types=[])

    assert layer.data["interfaces"] == []


def test_chain_interfaces_ca_skips_chains_without_selected_atoms():
    # DNA chains have no CA atoms; they must be skipped, not crash.
    layer = annotate_chain_interfaces(_load("1a02"), atom_set="ca")

    for interface in layer.data["interfaces"]:
        assert interface["chain_id_1"] not in {"A", "B"}
        assert interface["chain_id_2"] not in {"A", "B"}


def test_chain_interfaces_rejects_unknown_atom_set():
    with pytest.raises(ValueError, match="atom_set"):
        annotate_chain_interfaces(_load("1a7f"), atom_set="CA")


def test_chain_interfaces_rebuild_from_recorded_parameters():
    structure = _load("1a02")
    polymer_types = ["polypeptide(L)"]

    layer = annotate_chain_interfaces(
        structure,
        distance_cutoff=8.0,
        atom_set="ca",
        polymer_types=polymer_types,
    )
    polymer_types.append("polyribonucleotide")  # caller mutates its list
    rebuilt = annotate_chain_interfaces(structure, **layer.parameters)

    assert layer.parameters == {
        "distance_cutoff": 8.0,
        "atom_set": "ca",
        "polymer_types": ["polypeptide(L)"],
    }
    assert layer.method == "pandora.basic.distance_cutoff_contacts.v2"
    assert layer.data["interfaces"]  # non-empty, so equality means something
    assert rebuilt.data == layer.data


def test_interface_records_carry_atom_set_and_residue_pairs():
    structure = _load("1a02")

    settings = {
        "distance_cutoff": 8.0,
        "atom_set": "ca",
        "polymer_types": ["polypeptide(L)"],
    }

    records = extract_interface_records(structure, **settings)
    layer = annotate_chain_interfaces(structure, **settings)

    assert len(records) == len(layer.data["interfaces"]) > 0
    for record, interface in zip(records, layer.data["interfaces"]):
        assert record.atom_set == "ca"
        assert record.residue_pairs == [
            tuple(pair) for pair in interface["residue_pairs"]
        ]
        assert record.chain_id_1 in {"C", "D", "E"}


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
