from pathlib import Path

from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import chain_completeness
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    IdentifierRules,
    ModifiedResidueRules,
    ResidueNumberingRules,
    ValidationRules,
    canonicalisationPolicy,
)

MMCIF_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"


def _load(entry_id: str):
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / f"{entry_id}.cif"))
    return structure


def _canonical(structure, **rules):
    policy = canonicalisationPolicy(
        policy_id="t",
        policy_name="t",
        policy_version="1",
        validation_rules=ValidationRules(fail_on_unresolved_issues=False),
        **rules,
    )
    canonical, _, _ = canonicalise_structure(structure, policy)
    return canonical


def _counts(records):
    return {
        r.chain_id: (
            r.seqres_length,
            r.missing_n_term,
            r.missing_c_term,
            r.missing_middle,
        )
        for r in records
    }


def test_completeness_counts_tails_and_middle():
    records, diagnostics = chain_completeness(_load("10mv"))
    assert _counts(records) == {"A": (275, 21, 18, 5)}
    assert records[0].entry_id == "10MV"
    assert diagnostics.warnings == []


def test_completeness_reports_every_polymer_chain():
    records, _ = chain_completeness(_load("1aui"))
    assert _counts(records) == {
        "A": (521, 13, 35, 95),
        "B": (169, 4, 0, 0),
    }


def test_incomplete_backbone_counts_ca_only_residues_as_missing():
    structure = _load("1p58")
    unobserved, _ = chain_completeness(structure, "unobserved")
    backbone, _ = chain_completeness(structure, "incomplete_backbone")
    assert _counts(unobserved)["A"] == (495, 0, 0, 5)
    assert _counts(backbone)["A"] == (495, 495, 0, 0)


def test_chain_without_seqres_gets_diagnostic_not_record():
    structure = _load("10mv")
    raw = {k: v for k, v in structure.raw.items() if k != "_entity_poly_seq"}
    records, diagnostics = chain_completeness(
        structure.model_copy(update={"raw": raw})
    )
    assert records == []
    assert [d.code for d in diagnostics.warnings] == ["NO_SEQRES"]


def test_renumbered_chain_gets_seqres_mismatch():
    canonical = _canonical(
        _load("10mv"),
        identifier_rules=IdentifierRules(
            residue_numbering=ResidueNumberingRules(strategy="renumber")
        ),
    )
    records, diagnostics = chain_completeness(canonical)
    assert records == []
    assert [d.code for d in diagnostics.warnings] == ["SEQRES_MISMATCH"]


def test_completeness_accepts_mapped_modified_residue():
    canonical = _canonical(
        _load("1b6w"),
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
    )
    records, diagnostics = chain_completeness(canonical)
    assert _counts(records) == {"A": (69, 0, 1, 0)}
    assert diagnostics.warnings == []
