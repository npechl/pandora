from itertools import combinations
from pathlib import Path

from pandora.annotations.entry import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    polymer_asym_ids,
)
from pandora.parsing import mmcif_to_structure
from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
    EntityRecord,
    EntryRecord,
    Structure,
)

MMCIF_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"


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
    assert layer.parameters == {"distance_cutoff": 4.0}


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

    layer = annotate_chain_interfaces(structure, distance_cutoff=4.0)

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
