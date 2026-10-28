from itertools import permutations, product
from pathlib import Path

import pytest

from pandora.canonicalisation.altlocs import _resolve_altlocs
from pandora.canonicalisation.assemblies import (
    _normalize_assemblies,
    _select_biological_assembly,
)
from pandora.canonicalisation.canonicalise import canonicalise_structure
from pandora.canonicalisation.chain_ids import (
    _apply_chain_map,
    _normalize_chain_ids,
)
from pandora.canonicalisation.entities import _normalize_entities
from pandora.canonicalisation.ligands import filter_ligands
from pandora.canonicalisation.missing_data import (
    _drop_emptied_records,
    _handle_incomplete_chains,
    _handle_missing_atoms,
    _handle_missing_residues,
)
from pandora.canonicalisation.modified_residues import _map_modified_residues
from pandora.canonicalisation.residues import _normalize_residue_numbering
from pandora.canonicalisation.validation import _validate
from pandora.export import structure_to_mmcif
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    AltlocRules,
    AssemblyChainCopy,
    AssemblyRules,
    ChainIdRules,
    EntityRules,
    IdentifierRules,
    IncompleteChainRules,
    LigandRules,
    MissingAtomsRules,
    MissingDataRules,
    MissingResiduesRules,
    ModifiedResidueRules,
    ResidueNumberingRules,
    ValidationRules,
    canonicalisationPolicy,
    canonicalisationProvenanceRules,
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
    Structure,
)

MMCIF_PATH = (
    Path(__file__).parent.parent
    / "datasets"
    / "dev"
    / "mmcif"
    / "ay"
    / "1ayi.cif"
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


def test_normalize_chain_ids_remap_continues_past_two_letters():
    asym_units = [
        AsymRecord(id=f"X{i}", entity_id="1") for i in range(26 + 26**2 + 2)
    ]
    chain_map, _ = _normalize_chain_ids(asym_units, "remap", record=False)

    ids = list(chain_map.values())
    assert len(set(ids)) == len(ids)
    assert ids[26 + 26**2 - 1 :] == ["ZZ", "AAA", "AAB"]


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


def test_renumbering_leaves_ligands_without_seq_id():
    structure = _load("10tm")  # CLR ligands have no label_seq_id
    ligand_ids = {a.id for a in structure.atoms if a.label_seq_id is None}
    assert ligand_ids
    for strategy in ("renumber", "use_auth_seq"):
        result, mapping = _normalize_residue_numbering(
            structure.atoms, {}, strategy, True, record=True
        )
        assert {a.id for a in result if a.label_seq_id is None} == ligand_ids
        assert all(i.original_seq_id is not None for i in mapping.items)


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


def test_drop_emptied_records_reports_emptied_chains():
    atoms = [_atom(label_asym_id="A", entity_id="1")]
    asym_units = [
        AsymRecord(id="A", entity_id="1"),
        AsymRecord(id="B", entity_id="2"),
    ]
    entities = [
        EntityRecord(id="1", type="polymer"),
        EntityRecord(id="2", type="polymer"),
    ]
    diagnostics = DiagnosticBundle()

    kept_asyms, kept_entities = _drop_emptied_records(
        atoms, asym_units, entities, diagnostics, "TEST"
    )

    assert [u.id for u in kept_asyms] == ["A"]
    assert [e.id for e in kept_entities] == ["1"]
    [warning] = diagnostics.warnings
    assert warning.code == "CHAIN_EMPTIED"
    assert warning.context["asym_id"] == "B"


@pytest.mark.parametrize(
    ("entry_id", "policy_update"),
    [
        (
            "1p58",
            {
                "missing_data_rules": MissingDataRules(
                    missing_atoms=MissingAtomsRules(
                        strategy="drop_partial_residue"
                    )
                )
            },
        ),
        (
            "1ayi",
            {
                "ligand_rules": LigandRules(
                    strategy="filter",
                    keep_waters=False,
                    keep_ions=False,
                    keep_nonpolymer_ligands=False,
                )
            },
        ),
    ],
)
def test_canonicalise_structure_drops_asyms_and_entities_without_atoms(
    entry_id, policy_update
):
    structure, _, _ = mmcif_to_structure(
        str(MMCIF_DIR / entry_id[1:3] / f"{entry_id}.cif")
    )
    policy = canonicalisationPolicy(
        policy_id="test",
        policy_name="test",
        policy_version="1.0.0",
        **policy_update,
    )

    canonical, _, _ = canonicalise_structure(structure, policy)

    assert len(canonical.entities) < len(structure.entities)
    assert {u.id for u in canonical.asym_units} == {
        a.label_asym_id for a in canonical.atoms
    }
    assert {e.id for e in canonical.entities} == {
        a.label_entity_id for a in canonical.atoms
    }


# validation --------------------------------------------------------------


@pytest.fixture
def compositional_disorder(tmp_path: Path) -> Structure:
    """Parse a synthetic position with two explicitly labelled residues."""

    path = tmp_path / "alternates.cif"
    path.write_text(
        """data_synthetic_alternates
_entry.id synthetic_alternates
loop_
_entity.id
_entity.type
_entity.pdbx_description
_entity.formula_weight
_entity.src_method
1 polymer 'Synthetic format fixture' ? syn
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code
_entity_poly.pdbx_seq_one_letter_code_can
_entity_poly.pdbx_strand_id
1 'polypeptide(L)' A A A
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 ALA y
1 1 GLY y
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N  A ALA A 1 1 ? 0 0 0 0.25 20 1 ALA A N  1
ATOM 2 C CA A ALA A 1 1 ? 1 0 0 0.25 20 1 ALA A CA 1
ATOM 3 C C  A ALA A 1 1 ? 2 0 0 0.25 20 1 ALA A C  1
ATOM 4 O O  A ALA A 1 1 ? 3 0 0 0.25 20 1 ALA A O  1
ATOM 5 N N  B GLY A 1 1 ? 0 1 0 0.75 20 1 GLY A N  1
ATOM 6 C CA B GLY A 1 1 ? 1 1 0 0.75 20 1 GLY A CA 1
ATOM 7 C C  B GLY A 1 1 ? 2 1 0 0.75 20 1 GLY A C  1
ATOM 8 O O  B GLY A 1 1 ? 3 1 0 0.75 20 1 GLY A O  1
"""
    )
    structure, diagnostics, status = mmcif_to_structure(str(path))
    assert status == "success"
    assert not diagnostics.errors
    assert structure is not None
    assert len(structure.atoms) == 8
    assert {a.group_PDB for a in structure.atoms} == {"ATOM"}
    assert {
        (a.label_seq_id, a.label_comp_id, a.label_alt_id)
        for a in structure.atoms
    } == {(1, "ALA", "A"), (1, "GLY", "B")}
    return structure


def test_preserve_compositional_disorder_roundtrip(
    compositional_disorder: Structure, tmp_path: Path
) -> None:
    original = compositional_disorder.model_dump()
    policy = canonicalisationPolicy(
        policy_id="preserve-alternates",
        policy_name="Preserve alternates",
        policy_version="1",
        altloc_rules=AltlocRules(strategy="preserve"),
        provenance_rules=canonicalisationProvenanceRules(
            emit_canonicalisation_report=True
        ),
    )
    canonical, mappings, provenance = canonicalise_structure(
        compositional_disorder, policy
    )
    assert canonical is not compositional_disorder
    assert canonical.atoms == compositional_disorder.atoms
    assert mappings.altloc_selection_mapping.items == []
    assert provenance.report == {"warnings": 0, "errors": 0}

    path = structure_to_mmcif(canonical, tmp_path / "canonical.cif")
    reparsed, diagnostics, status = mmcif_to_structure(str(path))
    assert status == "success"
    assert not diagnostics.errors
    assert reparsed is not None
    assert reparsed.atoms == canonical.atoms
    recanonical, _, _ = canonicalise_structure(reparsed, policy)
    assert recanonical.atoms == canonical.atoms
    assert compositional_disorder.model_dump() == original


def test_default_selection_of_compositional_disorder_is_unchanged(
    compositional_disorder: Structure,
) -> None:
    original = compositional_disorder.model_dump()
    policy = canonicalisationPolicy(
        policy_id="default", policy_name="Default", policy_version="1"
    )
    assert policy.altloc_rules.strategy == "select_best_occupancy"
    canonical, mappings, _ = canonicalise_structure(
        compositional_disorder, policy
    )
    assert [a.id for a in canonical.atoms] == [5, 6, 7, 8]
    assert {(a.label_comp_id, a.label_alt_id) for a in canonical.atoms} == {
        ("GLY", None)
    }
    [selection] = mappings.altloc_selection_mapping.items
    assert selection.selected_altloc == "B"
    assert selection.available_altlocs == ["A", "B"]
    assert selection.selection_reason == "best_occupancy"
    assert compositional_disorder.model_dump() == original


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize(
    "altlocs", [("A", "A"), *product((None, "", ".", "?"), (None, "A"))]
)
def test_preserve_still_rejects_conflicting_residues(
    compositional_disorder: Structure,
    altlocs: tuple[str | None, str | None],
    reverse: bool,
) -> None:
    atoms = [
        atom.model_copy(update={"label_alt_id": altlocs[i // 4]})
        for i, atom in enumerate(compositional_disorder.atoms)
    ]
    structure = compositional_disorder.model_copy(
        update={"atoms": atoms[::-1] if reverse else atoms}
    )
    original = structure.model_dump()
    policy = canonicalisationPolicy(
        policy_id="preserve",
        policy_name="Preserve",
        policy_version="1",
        altloc_rules=AltlocRules(strategy="preserve"),
    )
    with pytest.raises(ValueError, match="Validation of canonicalised"):
        canonicalise_structure(structure, policy)
    assert structure.model_dump() == original


@pytest.mark.parametrize(
    ("residues", "collision"),
    [
        ([("ALA", "A"), ("GLY", "B"), ("SER", "C")], False),
        ([("ALA", "A"), ("GLY", "B"), ("SER", "B")], True),
        ([("ALA", "A"), ("ALA", "B"), ("GLY", "B")], True),
        ([("ALA", None), ("ALA", "A"), ("ALA", "B")], False),
        ([("ALA", None), ("ALA", "A"), ("GLY", "B")], True),
        *[
            ([("ALA", left), ("GLY", right)], True)
            for left, right in product((None, "", ".", "?"), repeat=2)
        ],
    ],
)
def test_validate_alternate_residue_collisions_in_every_order(
    residues: list[tuple[str, str | None]], collision: bool
) -> None:
    for order in permutations(residues):
        # Repeated atoms of the same residue must not add more errors.
        atoms = [
            _atom(
                id=i * 2 + j + 1,
                label_comp_id=comp,
                label_alt_id=alt,
                atom_id=name,
            )
            for i, (comp, alt) in enumerate(order)
            for j, name in enumerate(("N", "CA"))
        ]
        original = [a.model_dump() for a in atoms]
        diagnostics = DiagnosticBundle()
        status = _validate(atoms, [], ValidationRules(), diagnostics, "test")
        assert status == ("failed" if collision else "success"), order
        assert [e.code for e in diagnostics.errors] == (
            ["RESIDUE_NUMBER_COLLISION"] if collision else []
        ), order
        assert [a.model_dump() for a in atoms] == original


@pytest.mark.parametrize(
    "other_slot",
    [{"label_asym_id": "B"}, {"label_seq_id": 2}, {"pdbx_PDB_ins_code": "A"}],
)
def test_validate_keeps_chain_sequence_and_insertion_slots_separate(
    other_slot: dict[str, str | int],
) -> None:
    atoms = [
        _atom(label_comp_id="ALA", label_alt_id="A"),
        _atom(id=2, label_comp_id="GLY", label_alt_id="A").model_copy(
            update=other_slot
        ),
    ]
    for order in (atoms, atoms[::-1]):
        diagnostics = DiagnosticBundle()
        assert (
            _validate(order, [], ValidationRules(), diagnostics, "test")
            == "success"
        )
        assert not diagnostics.errors


@pytest.mark.parametrize("count", [128, 256])
def test_validate_many_named_alternates_uses_linear_label_comparisons(
    count: int,
) -> None:
    class CountedAltloc(str):
        comparisons = 0

        def __eq__(self, other: object) -> bool:
            type(self).comparisons += 1
            return super().__eq__(other)

        __hash__ = str.__hash__

    # model_copy retains the counting str subclass used by this control.
    atoms = [
        _atom(id=i + 1, label_comp_id=f"RES{i}").model_copy(
            update={"label_alt_id": CountedAltloc(f"alt{i}")}
        )
        for i in range(count)
    ]
    diagnostics = DiagnosticBundle()

    status = _validate(atoms, [], ValidationRules(), diagnostics, "test")

    assert status == "success"
    assert not diagnostics.errors
    # This deterministic bound catches scans of all earlier labels without
    # relying on machine speed or a fragile wall-clock threshold.
    assert CountedAltloc.comparisons <= 8 * count


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

    new_entities, new_asyms, new_atoms, mapping, _ = _normalize_entities(
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

    # The generator lists A under operator 2 only, so the untransformed A
    # isn't part of the assembly.
    assert {a.label_asym_id for a in new_atoms} == {"A_2"}
    copy_atom = next(a for a in new_atoms if a.label_asym_id == "A_2")
    assert (copy_atom.Cartn_x, copy_atom.Cartn_y, copy_atom.Cartn_z) == (
        11.0,
        2.0,
        3.0,
    )
    assert {a.id for a in new_asyms} == {"A_2"}
    # The deposited count already describes the whole assembly.
    assert result[0].oligomeric_count == 1
    assert mapping.items[0].chain_copies == [
        AssemblyChainCopy(
            canonical_chain_id="A_2", source_chain_id="A", operator_id="2"
        )
    ]


def _assembly_with_identity_and_copy(*, identity_chains, copied_chains):
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
                    oper_expression="1",
                    asym_id_list=identity_chains,
                ),
                AssemblyGenRecord(
                    assembly_id="1",
                    oper_expression="2",
                    asym_id_list=copied_chains,
                ),
            ],
            operators=[identity_op, translate_op],
        )
    ]


def test_standardize_biological_assembly_drops_chains_outside_assembly():
    # B isn't in any generator: its atoms, asym unit, its connection to A
    # and its secondary structure all go.
    assemblies = _assembly_with_identity_and_copy(
        identity_chains=["A"], copied_chains=[]
    )
    atoms = [
        _atom(label_asym_id="A", x=1.0, y=2.0, z=3.0),
        _atom(label_asym_id="B", id=2, x=4.0, y=5.0, z=6.0),
    ]
    asym_units = [
        AsymRecord(id="A", entity_id="1"),
        AsymRecord(id="B", entity_id="1"),
    ]
    conn = ConnRecord(
        id="disulf1",
        conn_type_id="disulf",
        ptnr1=ConnPartner(label_asym_id="A", label_comp_id="CYS"),
        ptnr2=ConnPartner(label_asym_id="B", label_comp_id="CYS"),
    )
    helix_b = ConfRecord(
        id="HELX_P1",
        conf_type_id="HELX_P",
        beg_label_asym_id="B",
        beg_label_seq_id=1,
        end_label_asym_id="B",
        end_label_seq_id=5,
    )

    _, new_atoms, new_asyms, new_conns, new_ss, _ = _normalize_assemblies(
        assemblies,
        atoms,
        asym_units,
        [conn],
        SSRecord(conf_records=[helix_b]),
        AssemblyRules(strategy="standardize_biological_assembly"),
        "preserve",
        True,
        DiagnosticBundle(),
        "test",
    )

    assert [a.label_asym_id for a in new_atoms] == ["A"]
    assert [a.id for a in new_asyms] == ["A"]
    assert new_conns == []
    assert new_ss.conf_records == []
    assert atoms[1].label_asym_id == "B"  # input untouched


def test_standardize_biological_assembly_expands_intra_chain_connection():
    assemblies = _assembly_with_identity_and_copy(
        identity_chains=["A"], copied_chains=["A"]
    )
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
    # A and B are in the assembly at identity, but only A also gets a copy
    # under operator 2 -- an A_2-B connection has no valid single-operator
    # remap and must not be fabricated.
    assemblies = _assembly_with_identity_and_copy(
        identity_chains=["A", "B"], copied_chains=["A"]
    )
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
    assemblies = _assembly_with_identity_and_copy(
        identity_chains=["A", "B"], copied_chains=["A", "B"]
    )
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
    assemblies = _assembly_with_identity_and_copy(
        identity_chains=["A"], copied_chains=["A"]
    )
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


def test_first_assembly_source_ignores_author_flag():
    from pandora.schemas.structure import AssemblyRecord

    assemblies = [
        AssemblyRecord(id="1", software_determined=True),
        AssemblyRecord(id="2", author_determined=True),
    ]
    picked = _select_biological_assembly(
        assemblies, "first", DiagnosticBundle(), "test"
    )
    assert picked.id == "1"


def test_author_assembly_source_still_prefers_author_flag():
    from pandora.schemas.structure import AssemblyRecord

    assemblies = [
        AssemblyRecord(id="1", software_determined=True),
        AssemblyRecord(id="2", author_determined=True),
    ]
    picked = _select_biological_assembly(
        assemblies, "author", DiagnosticBundle(), "test"
    )
    assert picked.id == "2"


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


@pytest.mark.parametrize("record", [True, False])
def test_canonicalise_structure_honours_record_original_assembly_mapping(
    record,
):
    structure, _, _ = mmcif_to_structure(str(MMCIF_PATH))
    policy = canonicalisationPolicy(
        policy_id="test",
        policy_name="test",
        policy_version="1.0.0",
        assembly_rules=AssemblyRules(record_original_assembly_mapping=record),
    )

    _, mappings, _ = canonicalise_structure(structure, policy)

    assert bool(mappings.assembly_mapping.items) is record


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
    structure, _, _ = mmcif_to_structure(
        str(MMCIF_DIR / entry_id[1:3] / f"{entry_id}.cif")
    )
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


def test_modified_residues_preserved_by_default():
    structure = _load("1b6w")
    canonical, mappings, _ = canonicalise_structure(
        structure, _lenient_policy()
    )
    assert any(a.label_comp_id == "MSE" for a in canonical.atoms)
    assert mappings.modified_residue_mapping.items == []


def test_mse_mapped_to_met():
    structure = _load("1b6w")
    policy = _lenient_policy(
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent")
    )
    canonical, mappings, provenance = canonicalise_structure(structure, policy)

    assert not any(a.label_comp_id == "MSE" for a in canonical.atoms)
    residue_35 = [
        a
        for a in canonical.atoms
        if a.label_asym_id == "A" and a.label_seq_id == 35
    ]
    assert residue_35
    assert all(
        a.label_comp_id == "MET" and a.auth_comp_id == "MET" for a in residue_35
    )
    assert {
        (i.chain_id, i.seq_id, i.original_comp_id, i.parent_comp_id)
        for i in mappings.modified_residue_mapping.items
    } == {("A", 1, "MSE", "MET"), ("A", 35, "MSE", "MET")}
    assert "modified_residues:map_to_parent" in provenance.transforms
    # input untouched
    assert any(a.label_comp_id == "MSE" for a in structure.atoms)


def test_mse_selenium_becomes_sulfur_delta():
    se = _atom(label_comp_id="MSE", atom_id="SE").model_copy(
        update={"type_symbol": "SE", "auth_atom_id": "SE"}
    )
    atoms, mapping = _map_modified_residues(
        [se],
        [],
        ModifiedResidueRules(strategy="map_to_parent"),
        DiagnosticBundle(),
        "test",
    )
    assert atoms[0].label_comp_id == "MET"
    assert atoms[0].label_atom_id == "SD"
    assert atoms[0].auth_atom_id == "SD"
    assert atoms[0].type_symbol == "S"
    assert len(mapping.items) == 1
    assert se.label_atom_id == "SE"


def test_listed_modified_residue_mapped_from_mod_residue_table():
    structure = _load("1a08")
    diagnostics = DiagnosticBundle()
    atoms, mapping = _map_modified_residues(
        list(structure.atoms),
        structure.raw["_pdbx_struct_mod_residue"],
        ModifiedResidueRules(strategy="map_to_parent", comp_ids=["FTY"]),
        diagnostics,
        structure.entry_id,
    )
    assert not any(a.label_comp_id == "FTY" for a in atoms)
    assert {
        (i.chain_id, i.seq_id, i.parent_comp_id) for i in mapping.items
    } == {("B", 2, "TYR"), ("D", 2, "TYR")}
    # ACE and DIP were not asked for: untouched, no warning
    assert any(a.label_comp_id == "ACE" for a in atoms)
    assert diagnostics.warnings == []


def test_unmapped_modified_residue_gets_warning_not_guess():
    structure = _load("1a08")
    diagnostics = DiagnosticBundle()
    atoms, _ = _map_modified_residues(
        list(structure.atoms),
        structure.raw["_pdbx_struct_mod_residue"],
        ModifiedResidueRules(strategy="map_to_parent", comp_ids=[]),
        diagnostics,
        structure.entry_id,
    )
    assert not any(a.label_comp_id == "FTY" for a in atoms)
    assert any(a.label_comp_id == "ACE" for a in atoms)
    unmapped = {
        d.context["comp_id"]
        for d in diagnostics.warnings
        if d.code == "MODIFIED_RESIDUE_UNMAPPED"
    }
    assert unmapped == {"ACE", "DIP"}


def test_assembly_expansion_keeps_only_assembly_chains():
    # 13dg: assembly 1 is chains A-D under operators 1, 2, 3; E-H belong
    # to assembly 2.
    policy = _lenient_policy(
        assembly_rules=AssemblyRules(strategy="standardize_biological_assembly")
    )
    canonical, _, _ = canonicalise_structure(_load("13dg"), policy)
    source_chains = {u.id.split("_")[0] for u in canonical.asym_units}
    assert source_chains == {"A", "B", "C", "D"}
    assert len(canonical.asym_units) == 12
    assert {a.label_asym_id for a in canonical.atoms} == {
        u.id for u in canonical.asym_units
    }
    assert canonical.assemblies[0].oligomeric_count == 12


def test_policy_that_keeps_residues_keeps_raw_categories():
    # Altloc selection and MSE mapping change no chain or residue ids.
    # (MSE mapping renames SEQRES mon_ids; see
    # test_mse_mapping_updates_bonds_and_seqres.)
    structure = _load("1b6w")
    policy = _lenient_policy(
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent")
    )
    canonical, _, provenance = canonicalise_structure(structure, policy)
    assert canonical.raw.keys() == structure.raw.keys()
    for category, rows in structure.raw.items():
        if category != "_entity_poly_seq":
            assert canonical.raw[category] == rows
    assert provenance.dropped_raw_categories == []


def test_renumbering_drops_raw_categories_indexed_by_residue():
    structure = _load("10mv")
    policy = _lenient_policy(
        identifier_rules=IdentifierRules(
            residue_numbering=ResidueNumberingRules(strategy="renumber")
        )
    )
    canonical, _, provenance = canonicalise_structure(structure, policy)

    dropped = set(provenance.dropped_raw_categories)
    assert {"_pdbx_sifts_xref_db", "_pdbx_poly_seq_scheme"} <= dropped
    assert dropped.isdisjoint(canonical.raw)
    # SEQRES and entry-level data don't index chains or residues.
    assert "_entity_poly_seq" in canonical.raw
    assert "_exptl" in canonical.raw
    assert "_pdbx_sifts_xref_db" in structure.raw  # input untouched


def test_entity_merge_rewrites_raw_entity_ids():
    # 1a3n: polymer entities 1 (alpha) and 2 (beta). Give entity 2 the
    # same canonical sequence so merge_equivalent_entities folds it into 1.
    structure = _load("1a3n")
    alpha = next(e for e in structure.entities if e.id == "1")
    seq = alpha.poly.pdbx_seq_one_letter_code_can
    entities = [
        e.model_copy(
            update={
                "poly": e.poly.model_copy(
                    update={"pdbx_seq_one_letter_code_can": seq}
                )
            }
        )
        if e.id == "2"
        else e
        for e in structure.entities
    ]
    structure = structure.model_copy(update={"entities": entities})
    policy = _lenient_policy(
        entity_rules=EntityRules(strategy="merge_equivalent_entities")
    )

    canonical, _, _ = canonicalise_structure(structure, policy)

    entity_ids = {e.id for e in canonical.entities}
    rows = canonical.raw["_entity_poly_seq"]
    assert {r["entity_id"] for r in rows} <= entity_ids
    assert rows == [
        r for r in structure.raw["_entity_poly_seq"] if r["entity_id"] == "1"
    ]


def test_mse_mapping_updates_bonds_and_seqres():
    structure = _load("1b6w")  # MSE at A1 and A35, entity 1
    bond = ConnRecord(
        id="test",
        conn_type_id="covale",
        ptnr1=ConnPartner(
            label_asym_id="A",
            label_comp_id="MSE",
            label_seq_id=35,
            label_atom_id="SE",
        ),
        ptnr2=ConnPartner(
            label_asym_id="A",
            label_comp_id="MSE",
            label_seq_id=1,
            label_atom_id="SE",
        ),
    )
    structure = structure.model_copy(update={"connections": [bond]})
    policy = _lenient_policy(
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent")
    )

    canonical, _, _ = canonicalise_structure(structure, policy)

    [conn] = canonical.connections
    for side in (conn.ptnr1, conn.ptnr2):
        assert (side.label_comp_id, side.label_atom_id) == ("MET", "SD")
    seqres = {
        int(r["num"]): r["mon_id"]
        for r in canonical.raw["_entity_poly_seq"]
        if r["entity_id"] == "1"
    }
    assert seqres[1] == seqres[35] == "MET"
    # input untouched
    assert structure.connections[0].ptnr1.label_comp_id == "MSE"
    assert any(r["mon_id"] == "MSE" for r in structure.raw["_entity_poly_seq"])


def test_residue_mapping_records_chain_ids_before_and_after_remap():
    # 10tm's asym order is A, C, B, D, so remap swaps B and C.
    structure = _load("10tm")
    policy = _lenient_policy(
        identifier_rules=IdentifierRules(
            chain_id=ChainIdRules(strategy="remap"),
            residue_numbering=ResidueNumberingRules(strategy="renumber"),
        )
    )
    _, mappings, _ = canonicalise_structure(structure, policy)

    chain_pairs = {
        (i.original_chain_id, i.canonical_chain_id)
        for i in mappings.chain_id_mapping.items
    }
    assert ("C", "B") in chain_pairs
    # Each entry must name a residue that exists under its original chain
    # (C has 282 residues, B only 5, so a swapped chain id can't pass).
    original_residues = {
        (a.label_asym_id, a.label_seq_id, a.auth_seq_id)
        for a in structure.atoms
    }
    for item in mappings.residue_number_mapping.items:
        assert (item.original_chain_id, item.canonical_chain_id) in chain_pairs
        assert (
            item.original_chain_id,
            item.original_seq_id,
            item.original_auth_seq_id,
        ) in original_residues
