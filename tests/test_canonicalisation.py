from pathlib import Path

import pytest

from pandora.canonicalisation.altlocs import _resolve_altlocs
from pandora.canonicalisation.assemblies import _normalize_assemblies
from pandora.canonicalisation.canonicalise import canonicalise_structure
from pandora.canonicalisation.chain_ids import (
    _apply_chain_map,
    _normalize_chain_ids,
)
from pandora.canonicalisation.entities import _normalize_entities
from pandora.canonicalisation.ligands import filter_ligands
from pandora.canonicalisation.missing_data import (
    _handle_incomplete_chains,
    _handle_missing_atoms,
    _handle_missing_residues,
)
from pandora.canonicalisation.residues import _normalize_residue_numbering
from pandora.canonicalisation.validation import _validate
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    AltlocRules,
    AssemblyChainCopy,
    AssemblyRules,
    EntityRules,
    IdentifierRules,
    LigandRules,
    MissingAtomsRules,
    MissingDataRules,
    MissingResiduesRules,
    ResidueNumberingRules,
    ValidationRules,
    canonicalisationPolicy,
    ChainIdRules,
    IncompleteChainRules,
)
from pandora.schemas.common import Diagnostic, DiagnosticBundle
from pandora.schemas.structure import (
    AsymRecord,
    AtomSiteRecord,
    ConfRecord,
    ConnPartner,
    ConnRecord,
    EntityPolyRecord,
    EntityRecord,
    SheetStrandRecord,
    SSRecord,
)

MMCIF_PATH = (
    Path(__file__).parent.parent / "datasets" / "dev" / "mmcif" / "1ayi.cif"
)


def _atom(
    *,
    group_PDB="ATOM",
    id=1,
    label_asym_id="A",
    label_seq_id=1,
    label_comp_id="ALA",
    label_alt_id=None,
    auth_seq_id="1",
    ins_code=None,
    entity_id="1",
    atom_id="CA",
    occupancy=1.0,
    b_iso=20.0,
    x=0.0,
    y=0.0,
    z=0.0,
) -> AtomSiteRecord:
    return AtomSiteRecord(
        group_PDB=group_PDB,
        id=id,
        type_symbol="C",
        label_atom_id=atom_id,
        label_alt_id=label_alt_id,
        label_comp_id=label_comp_id,
        label_asym_id=label_asym_id,
        label_entity_id=entity_id,
        label_seq_id=label_seq_id,
        pdbx_PDB_ins_code=ins_code,
        Cartn_x=x,
        Cartn_y=y,
        Cartn_z=z,
        occupancy=occupancy,
        B_iso_or_equiv=b_iso,
        auth_seq_id=auth_seq_id,
        auth_comp_id=label_comp_id,
        auth_asym_id=label_asym_id,
        auth_atom_id=atom_id,
    )


# chain_ids -------------------------------------------------------------


def test_normalize_chain_ids_remap_assigns_sequential_letters():
    asym_units = [
        AsymRecord(id="C", entity_id="1", auth_id="C"),
        AsymRecord(id="A", entity_id="2", auth_id="A"),
    ]
    chain_map, mapping = _normalize_chain_ids(asym_units, "remap", record=True)

    assert chain_map == {"C": "A", "A": "B"}
    assert [item.canonical_chain_id for item in mapping.items] == ["A", "B"]


def test_apply_chain_map_updates_atoms_and_asyms():
    atoms = [_atom(label_asym_id="C")]
    asym_units = [AsymRecord(id="C", entity_id="1", auth_id="C")]
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    structure = structure.model_copy(
        update={"atoms": atoms, "asym_units": asym_units}
    )

    new_atoms, new_asyms, _ = _apply_chain_map(
        atoms, asym_units, structure, {"C": "A"}
    )

    assert new_atoms[0].label_asym_id == "A"
    assert new_asyms[0].id == "A"


# residues ----------------------------------------------------------------


def test_renumber_gives_each_distinct_residue_sequential_ids():
    atoms = [
        _atom(id=1, label_seq_id=5, auth_seq_id="5", atom_id="N"),
        _atom(id=2, label_seq_id=5, auth_seq_id="5", atom_id="CA"),
        _atom(id=3, label_seq_id=9, auth_seq_id="9", atom_id="N"),
    ]
    result, mapping = _normalize_residue_numbering(
        atoms, {}, "renumber", preserve_insertion_codes=True, record=True
    )

    assert [a.label_seq_id for a in result] == [1, 1, 2]
    assert len(mapping.items) == 2


def test_use_auth_seq_parses_digits_from_auth_seq_id():
    atoms = [_atom(label_seq_id=1, auth_seq_id="42A")]
    result, _ = _normalize_residue_numbering(
        atoms, {}, "use_auth_seq", preserve_insertion_codes=False, record=False
    )

    assert result[0].label_seq_id == 42
    assert result[0].pdbx_PDB_ins_code is None


def test_preserve_drops_insertion_codes_when_not_requested():
    atoms = [_atom(ins_code="A")]
    result, _ = _normalize_residue_numbering(
        atoms, {}, "preserve", preserve_insertion_codes=False, record=False
    )

    assert result[0].pdbx_PDB_ins_code is None


# altlocs -------------------------------------------------------------------


def test_select_best_occupancy_keeps_higher_occupancy_atom():
    atoms = [
        _atom(id=1, label_alt_id="A", occupancy=0.4, atom_id="CA"),
        _atom(id=2, label_alt_id="B", occupancy=0.6, atom_id="CA"),
    ]
    rules = AltlocRules(strategy="select_best_occupancy")

    result, mapping = _resolve_altlocs(atoms, rules)

    assert len(result) == 1
    assert result[0].id == 2
    assert result[0].label_alt_id is None
    assert mapping.items[0].selected_altloc == "B"


def test_altloc_resolution_does_not_merge_distinct_hetatm_residues():
    # Two separate waters: both HETATM, both label_seq_id=None (the field
    # is only populated for polymer atoms), disambiguated only by
    # auth_seq_id. Regression for a bug where these collapsed into one
    # altloc group.
    atoms = [
        _atom(
            id=1,
            group_PDB="HETATM",
            label_asym_id="B",
            label_seq_id=None,
            label_comp_id="HOH",
            auth_seq_id="101",
            label_alt_id="A",
            occupancy=0.5,
        ),
        _atom(
            id=2,
            group_PDB="HETATM",
            label_asym_id="B",
            label_seq_id=None,
            label_comp_id="HOH",
            auth_seq_id="101",
            label_alt_id="B",
            occupancy=0.5,
        ),
        _atom(
            id=3,
            group_PDB="HETATM",
            label_asym_id="B",
            label_seq_id=None,
            label_comp_id="HOH",
            auth_seq_id="102",
            label_alt_id="A",
            occupancy=1.0,
        ),
    ]
    rules = AltlocRules(strategy="select_best_occupancy")

    result, _ = _resolve_altlocs(atoms, rules)

    assert {a.auth_seq_id for a in result} == {"101", "102"}
    assert len(result) == 2


def test_altloc_resolution_picks_a_winner_for_compositional_disorder():
    # mmCIF also uses label_alt_id for *compositional* disorder: two
    # different residue identities modeled at the same physical site
    # (real example: PDB 1AL4, altloc A = FVA, altloc B = QIL at chain A
    # seq_id 1). Grouping by (chain, auth_seq_id, label_seq_id, comp_id)
    # would treat these as two unrelated single-altloc residues, leave
    # both in the output, and produce an unresolved residue-number
    # collision downstream even though the file already resolved it via
    # label_alt_id.
    atoms = [
        _atom(
            id=1,
            label_seq_id=1,
            label_comp_id="FVA",
            label_alt_id="A",
            occupancy=0.7,
            atom_id="CA",
        ),
        _atom(
            id=2,
            label_seq_id=1,
            label_comp_id="QIL",
            label_alt_id="B",
            occupancy=0.3,
            atom_id="CA",
        ),
    ]
    rules = AltlocRules(strategy="select_best_occupancy")

    result, mapping = _resolve_altlocs(atoms, rules)

    assert len(result) == 1
    assert result[0].label_comp_id == "FVA"
    assert result[0].label_alt_id is None
    assert mapping.items[0].residue_id == "1:FVA"
    assert mapping.items[0].available_altlocs == ["A", "B"]


# missing_data ----------------------------------------------------------


def test_drop_partial_residue_removes_residue_missing_backbone_atoms():
    atoms = [
        _atom(id=1, label_seq_id=1, atom_id="N"),
        _atom(id=2, label_seq_id=1, atom_id="CA"),
        # residue 1 is missing C/O -> dropped
        _atom(id=3, label_seq_id=2, atom_id="N"),
        _atom(id=4, label_seq_id=2, atom_id="CA"),
        _atom(id=5, label_seq_id=2, atom_id="C"),
        _atom(id=6, label_seq_id=2, atom_id="O"),
    ]
    diagnostics = DiagnosticBundle()
    rules = MissingAtomsRules(strategy="drop_partial_residue")

    result = _handle_missing_atoms(atoms, rules, diagnostics, "test", [])

    assert {a.label_seq_id for a in result} == {2}
    assert diagnostics.warnings


def test_missing_atoms_does_not_collide_distinct_hetatm_residues():
    # Two waters sharing label_seq_id=None must not be treated as the
    # same residue when checking for missing backbone atoms.
    atoms = [
        _atom(
            id=1,
            group_PDB="HETATM",
            label_seq_id=None,
            label_comp_id="HOH",
            auth_seq_id="101",
            atom_id="O",
        ),
        _atom(
            id=2,
            group_PDB="HETATM",
            label_seq_id=None,
            label_comp_id="HOH",
            auth_seq_id="102",
            atom_id="O",
        ),
    ]
    diagnostics = DiagnosticBundle()
    rules = MissingAtomsRules(strategy="drop_partial_residue")

    result = _handle_missing_atoms(atoms, rules, diagnostics, "test", [])

    # Neither is a polymer ATOM record, so nothing should be dropped or
    # flagged regardless of grouping.
    assert len(result) == 2
    assert not diagnostics.warnings


def test_drop_chain_segment_removes_chain_with_a_gap():
    atoms = [
        _atom(id=1, label_asym_id="A", label_seq_id=1),
        _atom(id=2, label_asym_id="A", label_seq_id=3),  # gap at 2
        _atom(id=3, label_asym_id="B", label_seq_id=1),
        _atom(id=4, label_asym_id="B", label_seq_id=2),
    ]
    diagnostics = DiagnosticBundle()
    rules = MissingResiduesRules(strategy="drop_chain_segment")

    result = _handle_missing_residues(atoms, rules, diagnostics, "test")

    assert {a.label_asym_id for a in result} == {"B"}


def test_incomplete_chains_truncate_keeps_longest_contiguous_run():
    atoms = [
        _atom(id=1, label_asym_id="A", label_seq_id=1),
        _atom(id=2, label_asym_id="A", label_seq_id=2),
        _atom(id=3, label_asym_id="A", label_seq_id=5),
        _atom(id=4, label_asym_id="A", label_seq_id=6),
        _atom(id=5, label_asym_id="A", label_seq_id=7),
    ]
    asym_units = [AsymRecord(id="A", entity_id="1")]
    rules = IncompleteChainRules(strategy="truncate_to_complete_regions")

    result, _ = _handle_incomplete_chains(atoms, asym_units, rules)

    assert {a.label_seq_id for a in result} == {5, 6, 7}


# validation --------------------------------------------------------------


def test_validate_does_not_flag_a_normal_multi_atom_residue():
    atoms = [
        _atom(id=1, label_seq_id=1, label_comp_id="ALA", atom_id="N"),
        _atom(id=2, label_seq_id=1, label_comp_id="ALA", atom_id="CA"),
    ]
    diagnostics = DiagnosticBundle()

    status = _validate(atoms, [], ValidationRules(), diagnostics, "test")

    assert status == "success"
    assert not diagnostics.errors


def test_validate_flags_a_genuine_residue_number_collision_once():
    atoms = [
        _atom(id=1, label_seq_id=1, label_comp_id="ALA", atom_id="N"),
        _atom(id=2, label_seq_id=1, label_comp_id="ALA", atom_id="CA"),
        # Same seq_id, different residue identity -> collision.
        _atom(id=3, label_seq_id=1, label_comp_id="GLY", atom_id="N"),
        _atom(id=4, label_seq_id=1, label_comp_id="GLY", atom_id="CA"),
    ]
    diagnostics = DiagnosticBundle()
    rules = ValidationRules(fail_on_unresolved_issues=False)

    status = _validate(atoms, [], rules, diagnostics, "test")

    assert status == "warning"
    assert len(diagnostics.errors) == 1
    assert diagnostics.errors[0].code == "RESIDUE_NUMBER_COLLISION"


def test_validate_fails_when_fail_on_unresolved_issues():
    atoms = [
        _atom(id=1, label_seq_id=1, label_comp_id="ALA"),
        _atom(id=2, label_seq_id=1, label_comp_id="GLY"),
    ]
    diagnostics = DiagnosticBundle()
    rules = ValidationRules(fail_on_unresolved_issues=True)

    status = _validate(atoms, [], rules, diagnostics, "test")

    assert status == "failed"


def test_validate_strict_promotes_warnings_to_errors():
    diagnostics = DiagnosticBundle()
    diagnostics.warnings.append(
        Diagnostic(code="SOME_WARNING", severity="warning", message="x")
    )
    rules = ValidationRules(strictness="strict", fail_on_unresolved_issues=True)

    status = _validate([], [], rules, diagnostics, "test")

    assert status == "failed"
    assert any(e.code == "SOME_WARNING" for e in diagnostics.errors)


# entities ------------------------------------------------------------------


def test_merge_equivalent_entities_by_canonical_sequence():
    entities = [
        EntityRecord(
            id="1",
            type="polymer",
            poly=EntityPolyRecord(pdbx_seq_one_letter_code_can="ACDE"),
        ),
        EntityRecord(
            id="2",
            type="polymer",
            poly=EntityPolyRecord(pdbx_seq_one_letter_code_can="ACDE"),
        ),
    ]
    asym_units = [
        AsymRecord(id="A", entity_id="1"),
        AsymRecord(id="B", entity_id="2"),
    ]
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1"),
        _atom(id=2, label_asym_id="B", entity_id="2"),
    ]
    rules = EntityRules(strategy="merge_equivalent_entities")

    new_entities, new_asyms, new_atoms, mapping = _normalize_entities(
        entities, asym_units, atoms, rules, record=True
    )

    assert len(new_entities) == 1
    assert {a.entity_id for a in new_asyms} == {new_entities[0].id}
    assert {a.label_entity_id for a in new_atoms} == {new_entities[0].id}
    assert sorted(mapping.items[0].original_entity_ids) == ["1", "2"]


# ligands ---------------------------------------------------------------


def testfilter_ligands_drops_waters_when_keep_waters_false():
    entities = [
        EntityRecord(id="1", type="polymer"),
        EntityRecord(id="2", type="water"),
    ]
    asym_units = [
        AsymRecord(id="A", entity_id="1"),
        AsymRecord(id="B", entity_id="2"),
    ]
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1"),
        _atom(
            id=2,
            group_PDB="HETATM",
            label_asym_id="B",
            entity_id="2",
            label_seq_id=None,
            label_comp_id="HOH",
        ),
    ]
    rules = LigandRules(strategy="filter", keep_waters=False)
    diagnostics = DiagnosticBundle()

    new_atoms, new_asyms = filter_ligands(
        atoms, asym_units, entities, rules, diagnostics, "test"
    )

    assert {a.label_asym_id for a in new_atoms} == {"A"}
    assert {a.id for a in new_asyms} == {"A"}


def testfilter_ligands_annotate_only_removes_but_records_diagnostic():
    entities = [
        EntityRecord(id="1", type="polymer"),
        EntityRecord(id="2", type="non-polymer", pdbx_description="ZINC ION"),
    ]
    asym_units = [
        AsymRecord(id="A", entity_id="1"),
        AsymRecord(id="B", entity_id="2"),
    ]
    atoms = [
        _atom(id=1, label_asym_id="A", entity_id="1"),
        _atom(
            id=2,
            group_PDB="HETATM",
            label_asym_id="B",
            entity_id="2",
            label_seq_id=None,
            label_comp_id="ZN",
        ),
    ]
    rules = LigandRules(strategy="annotate_only")
    diagnostics = DiagnosticBundle()

    new_atoms, new_asyms = filter_ligands(
        atoms, asym_units, entities, rules, diagnostics, "test"
    )

    assert {a.label_asym_id for a in new_atoms} == {"A"}
    assert diagnostics.warnings[0].code == "LIGAND_ANNOTATED_ONLY"


# assemblies --------------------------------------------------------------


def test_select_first_assembly_keeps_only_the_first():
    from pandora.schemas.structure import AssemblyRecord

    assemblies = [AssemblyRecord(id="1"), AssemblyRecord(id="2")]
    rules = AssemblyRules(strategy="select_first_assembly")
    diagnostics = DiagnosticBundle()

    result, _, _, _, _, _ = _normalize_assemblies(
        assemblies,
        [],
        [],
        [],
        SSRecord(),
        rules,
        "preserve",
        False,
        diagnostics,
        "test",
    )

    assert [a.id for a in result] == ["1"]


def test_standardize_biological_assembly_noop_without_assemblies():
    rules = AssemblyRules(strategy="standardize_biological_assembly")
    diagnostics = DiagnosticBundle()
    atoms = [_atom(label_asym_id="A")]
    asym_units = [AsymRecord(id="A", entity_id="1")]

    result, new_atoms, new_asyms, new_conns, new_ss, _ = _normalize_assemblies(
        [],
        atoms,
        asym_units,
        [],
        SSRecord(),
        rules,
        "preserve",
        False,
        diagnostics,
        "test",
    )

    assert result == []
    assert new_atoms == atoms
    assert new_asyms == asym_units
    assert new_conns == []
    assert new_ss == SSRecord()


def test_standardize_biological_assembly_identity_only_is_unchanged():
    from pandora.schemas.structure import (
        AssemblyGenRecord,
        AssemblyOperRecord,
        AssemblyRecord,
    )

    identity_op = AssemblyOperRecord(
        id="1",
        matrix=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        vector=[0.0, 0.0, 0.0],
    )
    assemblies = [
        AssemblyRecord(
            id="1",
            author_determined=True,
            generators=[
                AssemblyGenRecord(
                    assembly_id="1", oper_expression="1", asym_id_list=["A"]
                )
            ],
            operators=[identity_op],
        )
    ]
    atoms = [_atom(label_asym_id="A")]
    asym_units = [AsymRecord(id="A", entity_id="1")]
    rules = AssemblyRules(strategy="standardize_biological_assembly")
    diagnostics = DiagnosticBundle()

    result, new_atoms, new_asyms, _, _, mapping = _normalize_assemblies(
        assemblies,
        atoms,
        asym_units,
        [],
        SSRecord(),
        rules,
        "preserve",
        True,
        diagnostics,
        "test",
    )

    assert [a.label_asym_id for a in new_atoms] == ["A"]
    assert [a.id for a in new_asyms] == ["A"]
    assert result[0].generators == []
    assert mapping.items[0].chain_copies == []


def test_standardize_biological_assembly_expands_symmetry_copy():
    from pandora.schemas.structure import (
        AssemblyGenRecord,
        AssemblyOperRecord,
        AssemblyRecord,
    )

    translate_op = AssemblyOperRecord(
        id="2",
        matrix=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        vector=[10.0, 0.0, 0.0],
    )
    assemblies = [
        AssemblyRecord(
            id="1",
            author_determined=True,
            oligomeric_count=1,
            generators=[
                AssemblyGenRecord(
                    assembly_id="1", oper_expression="2", asym_id_list=["A"]
                )
            ],
            operators=[translate_op],
        )
    ]
    atoms = [_atom(label_asym_id="A", x=1.0, y=2.0, z=3.0)]
    asym_units = [AsymRecord(id="A", entity_id="1")]
    rules = AssemblyRules(strategy="standardize_biological_assembly")
    diagnostics = DiagnosticBundle()

    result, new_atoms, new_asyms, _, _, mapping = _normalize_assemblies(
        assemblies,
        atoms,
        asym_units,
        [],
        SSRecord(),
        rules,
        "preserve",
        True,
        diagnostics,
        "test",
    )

    assert {a.label_asym_id for a in new_atoms} == {"A", "A_2"}
    copy_atom = next(a for a in new_atoms if a.label_asym_id == "A_2")
    assert (copy_atom.Cartn_x, copy_atom.Cartn_y, copy_atom.Cartn_z) == (
        11.0,
        2.0,
        3.0,
    )
    assert {a.id for a in new_asyms} == {"A", "A_2"}
    assert result[0].oligomeric_count == 2
    assert mapping.items[0].chain_copies == [
        AssemblyChainCopy(
            canonical_chain_id="A_2", source_chain_id="A", operator_id="2"
        )
    ]


def _assembly_with_symmetry_copy(*, asym_id_list):
    from pandora.schemas.structure import (
        AssemblyGenRecord,
        AssemblyOperRecord,
        AssemblyRecord,
    )

    translate_op = AssemblyOperRecord(
        id="2",
        matrix=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        vector=[10.0, 0.0, 0.0],
    )
    return [
        AssemblyRecord(
            id="1",
            author_determined=True,
            generators=[
                AssemblyGenRecord(
                    assembly_id="1",
                    oper_expression="2",
                    asym_id_list=asym_id_list,
                )
            ],
            operators=[translate_op],
        )
    ]


def test_standardize_biological_assembly_expands_intra_chain_connection():
    assemblies = _assembly_with_symmetry_copy(asym_id_list=["A"])
    atoms = [_atom(label_asym_id="A", x=1.0, y=2.0, z=3.0)]
    asym_units = [AsymRecord(id="A", entity_id="1")]
    # a disulfide-style bond within chain A -- rigid-body copies of A
    # (i.e. A_2) preserve this bond exactly, so it should get duplicated.
    conn = ConnRecord(
        id="disulf1",
        conn_type_id="disulf",
        ptnr1=ConnPartner(label_asym_id="A", label_comp_id="CYS"),
        ptnr2=ConnPartner(label_asym_id="A", label_comp_id="CYS"),
        pdbx_dist_value=2.05,
    )
    rules = AssemblyRules(strategy="standardize_biological_assembly")
    diagnostics = DiagnosticBundle()

    _, _, _, new_conns, _, _ = _normalize_assemblies(
        assemblies,
        atoms,
        asym_units,
        [conn],
        SSRecord(),
        rules,
        "preserve",
        False,
        diagnostics,
        "test",
    )

    assert conn in new_conns  # original, untouched
    copy = next(c for c in new_conns if c is not conn)
    assert copy.ptnr1.label_asym_id == "A_2"
    assert copy.ptnr2.label_asym_id == "A_2"
    assert copy.pdbx_dist_value == 2.05  # rigid transform preserves distance


def test_standardize_biological_assembly_drops_inter_chain_connection_when_only_one_side_copied():
    # only chain A is in the generator's asym_id_list, so B never gets a
    # copy -- an A-B connection has no valid single-operator remap and
    # must not be fabricated.
    assemblies = _assembly_with_symmetry_copy(asym_id_list=["A"])
    atoms = [
        _atom(label_asym_id="A", x=1.0, y=2.0, z=3.0),
        _atom(label_asym_id="B", id=2, x=4.0, y=5.0, z=6.0),
    ]
    asym_units = [
        AsymRecord(id="A", entity_id="1"),
        AsymRecord(id="B", entity_id="2"),
    ]
    conn = ConnRecord(
        id="metalc1",
        conn_type_id="metalc",
        ptnr1=ConnPartner(label_asym_id="A", label_comp_id="CYS"),
        ptnr2=ConnPartner(label_asym_id="B", label_comp_id="ZN"),
    )
    rules = AssemblyRules(strategy="standardize_biological_assembly")
    diagnostics = DiagnosticBundle()

    _, _, _, new_conns, _, _ = _normalize_assemblies(
        assemblies,
        atoms,
        asym_units,
        [conn],
        SSRecord(),
        rules,
        "preserve",
        False,
        diagnostics,
        "test",
    )

    assert new_conns == [conn]  # no fabricated A_2-B connection


def test_standardize_biological_assembly_expands_inter_chain_connection_when_both_sides_copied_together():
    # both A and B are copied under the same operator (listed together in
    # the same generator), so A-B's connection has a valid remap: A_2-B_2.
    assemblies = _assembly_with_symmetry_copy(asym_id_list=["A", "B"])
    atoms = [
        _atom(label_asym_id="A", x=1.0, y=2.0, z=3.0),
        _atom(label_asym_id="B", id=2, x=4.0, y=5.0, z=6.0),
    ]
    asym_units = [
        AsymRecord(id="A", entity_id="1"),
        AsymRecord(id="B", entity_id="2"),
    ]
    conn = ConnRecord(
        id="metalc1",
        conn_type_id="metalc",
        ptnr1=ConnPartner(label_asym_id="A", label_comp_id="CYS"),
        ptnr2=ConnPartner(label_asym_id="B", label_comp_id="ZN"),
    )
    rules = AssemblyRules(strategy="standardize_biological_assembly")
    diagnostics = DiagnosticBundle()

    _, _, _, new_conns, _, _ = _normalize_assemblies(
        assemblies,
        atoms,
        asym_units,
        [conn],
        SSRecord(),
        rules,
        "preserve",
        False,
        diagnostics,
        "test",
    )

    copy = next(c for c in new_conns if c is not conn)
    assert copy.ptnr1.label_asym_id == "A_2"
    assert copy.ptnr2.label_asym_id == "B_2"


def test_standardize_biological_assembly_expands_secondary_structure():
    assemblies = _assembly_with_symmetry_copy(asym_id_list=["A"])
    atoms = [_atom(label_asym_id="A", x=1.0, y=2.0, z=3.0)]
    asym_units = [AsymRecord(id="A", entity_id="1")]
    helix = ConfRecord(
        id="H1",
        conf_type_id="HELX_P",
        beg_label_asym_id="A",
        beg_label_seq_id=1,
        end_label_asym_id="A",
        end_label_seq_id=10,
    )
    strand = SheetStrandRecord(
        sheet_id="S1",
        id="1",
        beg_label_asym_id="A",
        beg_label_seq_id=20,
        end_label_asym_id="A",
        end_label_seq_id=25,
    )
    ss = SSRecord(conf_records=[helix], sheet_strands=[strand])
    rules = AssemblyRules(strategy="standardize_biological_assembly")
    diagnostics = DiagnosticBundle()

    _, _, _, _, new_ss, _ = _normalize_assemblies(
        assemblies,
        atoms,
        asym_units,
        [],
        ss,
        rules,
        "preserve",
        False,
        diagnostics,
        "test",
    )

    assert helix in new_ss.conf_records
    helix_copy = next(c for c in new_ss.conf_records if c is not helix)
    assert helix_copy.beg_label_asym_id == "A_2"
    assert helix_copy.end_label_asym_id == "A_2"
    assert helix_copy.beg_label_seq_id == 1  # residue range unchanged

    assert strand in new_ss.sheet_strands
    strand_copy = next(s for s in new_ss.sheet_strands if s is not strand)
    assert strand_copy.beg_label_asym_id == "A_2"
    assert strand_copy.end_label_asym_id == "A_2"


def test_standardize_biological_assembly_pdbe_source_falls_back_with_warning():
    from pandora.schemas.structure import AssemblyRecord

    assemblies = [AssemblyRecord(id="1", software_determined=True)]
    rules = AssemblyRules(
        strategy="standardize_biological_assembly",
        preferred_assembly_source="pdbe",
    )
    diagnostics = DiagnosticBundle()

    result, _, _, _, _, _ = _normalize_assemblies(
        assemblies,
        [],
        [],
        [],
        SSRecord(),
        rules,
        "preserve",
        False,
        diagnostics,
        "test",
    )

    assert result[0].id == "1"
    assert diagnostics.warnings[0].code == "ASSEMBLY_SOURCE_UNAVAILABLE"


# canonicalise_structure (end to end) ------------------------------------


def test_canonicalise_structure_records_transforms_for_non_preserve_rules():
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    policy = canonicalisationPolicy(
        policy_id="test",
        policy_name="test",
        policy_version="1.0.0",
        identifier_rules=IdentifierRules(
            chain_id=ChainIdRules(strategy="remap"),
            residue_numbering=ResidueNumberingRules(strategy="renumber"),
        ),
        altloc_rules=AltlocRules(strategy="select_best_occupancy"),
    )

    canonical, mappings, provenance = canonicalise_structure(structure, policy)

    assert "chain_id:remap" in provenance.transforms
    assert "residue_numbering:renumber" in provenance.transforms
    assert canonical.entry_id == structure.entry_id
    assert mappings.chain_id_mapping.items


def test_canonicalise_structure_raises_on_unresolved_validation_failure():
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    # Force a residue-number collision: renumber every atom to seq_id 1
    # while keeping distinct comp_ids, then validate strictly.
    colliding_atoms = [
        a.model_copy(update={"label_seq_id": 1, "label_comp_id": "ALA"})
        if i % 2 == 0
        else a.model_copy(update={"label_seq_id": 1, "label_comp_id": "GLY"})
        for i, a in enumerate(structure.atoms[:4])
    ]
    structure = structure.model_copy(update={"atoms": colliding_atoms})
    policy = canonicalisationPolicy(
        policy_id="strict",
        policy_name="strict",
        policy_version="1.0.0",
        validation_rules=ValidationRules(fail_on_unresolved_issues=True),
    )

    with pytest.raises(ValueError):
        canonicalise_structure(structure, policy)


MMCIF_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"


def _load(entry_id: str):
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / f"{entry_id}.cif"))
    return structure


def _lenient_policy(**rules) -> canonicalisationPolicy:
    return canonicalisationPolicy(
        policy_id="test",
        policy_name="test",
        policy_version="1.0.0",
        validation_rules=ValidationRules(fail_on_unresolved_issues=False),
        **rules,
    )


def test_use_auth_chain_id_keeps_ligand_and_water_asym_units_distinct():
    # 104m: one protein chain plus ligands and waters, all with auth "A".
    structure = _load("104m")
    policy = _lenient_policy(
        identifier_rules=IdentifierRules(
            chain_id=ChainIdRules(strategy="use_auth_chain_id")
        )
    )

    canonical, _, _ = canonicalise_structure(structure, policy)

    asym_ids = [asym.id for asym in canonical.asym_units]
    assert len(asym_ids) == len(structure.asym_units)
    assert len(set(asym_ids)) == len(asym_ids)
    assert asym_ids[0] == "A"
    assert len(canonical.atoms) == len(structure.atoms)


def test_filter_ligands_does_not_treat_phosphate_ligand_as_ion():
    # 1bvi: guanosine-2'-monophosphate (2GP) and a calcium ion.
    structure = _load("1bvi")
    rules = LigandRules(
        strategy="filter",
        keep_waters=False,
        keep_ions=False,
        keep_nonpolymer_ligands=True,
    )

    atoms, _ = filter_ligands(
        structure.atoms,
        structure.asym_units,
        structure.entities,
        rules,
        DiagnosticBundle(),
        structure.entry_id,
    )

    het_comp_ids = {a.label_comp_id for a in atoms if a.group_PDB == "HETATM"}
    assert het_comp_ids == {"2GP"}


def test_drop_partial_residue_keeps_complete_nucleotides():
    # 11kb: protein plus a DNA chain; nucleotides have no N/CA/C/O.
    structure = _load("11kb")
    dna_entities = {
        e.id
        for e in structure.entities
        if e.poly and "nucleotide" in e.poly.type
    }
    policy = _lenient_policy(
        missing_data_rules=MissingDataRules(
            missing_atoms=MissingAtomsRules(strategy="drop_partial_residue")
        ),
        altloc_rules=AltlocRules(strategy="preserve"),
    )

    canonical, _, _ = canonicalise_structure(structure, policy)

    def n_dna_atoms(s):
        return sum(a.label_entity_id in dna_entities for a in s.atoms)

    assert n_dna_atoms(structure) > 0
    assert n_dna_atoms(canonical) == n_dna_atoms(structure)


def test_renumber_does_not_hide_sequence_gaps():
    # 1a3e: polymer chain B has a gap in label_seq_id; A and C don't.
    structure = _load("1a3e")
    policy = _lenient_policy(
        identifier_rules=IdentifierRules(
            residue_numbering=ResidueNumberingRules(strategy="renumber")
        ),
        missing_data_rules=MissingDataRules(
            missing_residues=MissingResiduesRules(strategy="drop_chain_segment")
        ),
    )

    canonical, _, _ = canonicalise_structure(structure, policy)

    polymer_chains = {
        a.label_asym_id for a in canonical.atoms if a.group_PDB == "ATOM"
    }
    assert polymer_chains == {"A", "C"}
