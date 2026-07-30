from pandora.annotations.entry import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
)
from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
    EntityRecord,
    EntryRecord,
    Structure,
)


def _atom(
    *,
    group_PDB="ATOM",
    id=1,
    label_asym_id="A",
    label_seq_id=1,
    label_comp_id="ALA",
    auth_seq_id="1",
    entity_id="1",
    x=0.0,
    y=0.0,
    z=0.0,
) -> AtomSiteRecord:
    return AtomSiteRecord(
        group_PDB=group_PDB,
        id=id,
        type_symbol="C",
        label_atom_id="CA",
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
        auth_atom_id="CA",
    )


def _structure(atoms, entities, asym_units) -> Structure:
    return Structure(
        entry_id="test",
        entry=EntryRecord(id="test"),
        entities=entities,
        asym_units=asym_units,
        atoms=atoms,
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
