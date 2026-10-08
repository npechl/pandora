from pathlib import Path

from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import curate_structure, deduplicate_structures
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import (
    ModifiedResidueRules,
    ValidationRules,
    canonicalisationPolicy,
)
from pandora.schemas.dataset import (
    ContentRules,
    DatasetCurationPolicy,
    DeduplicationRules,
    OrganismRules,
    QualityRules,
)
from pandora.schemas.metadata import MetadataRecord, QualityRecord

MMCIF_DIR = Path(__file__).parent.parent / "datasets" / "dev" / "mmcif"


def _load(entry_id: str):
    structure, _, _ = mmcif_to_structure(str(MMCIF_DIR / f"{entry_id}.cif"))
    return structure


def _policy(**overrides) -> DatasetCurationPolicy:
    return DatasetCurationPolicy(
        policy_id="test",
        policy_name="test",
        policy_version="1.0.0",
        **overrides,
    )


def test_curate_structure_pass_and_fail():
    structure = _load("1ayi")
    metadata = collect_metadata(structure)
    policy = _policy()

    curated, exclusions, provenance = curate_structure(
        structure, metadata, policy
    )
    assert curated is not None
    assert exclusions == []
    assert provenance.policy_id == "test"
    assert provenance.policy_version == "1.0.0"

    tight_policy = _policy(
        quality_rules=QualityRules(
            max_resolution=metadata.quality.resolution - 0.01
        )
    )
    curated, exclusions, provenance = curate_structure(
        structure, metadata, tight_policy
    )
    assert curated is None
    assert [e.reason_code for e in exclusions] == ["RESOLUTION_THRESHOLD"]
    # provenance is populated on exclusion too, so callers always know
    # which policy made the call.
    assert provenance.policy_id == "test"


def test_missing_metadata_treated_as_unknown_not_skipped():
    structure = _load("1ayi")

    curated, exclusions, _ = curate_structure(structure, None, _policy())
    assert curated is not None
    assert exclusions == []

    # No metadata means resolution is unknown, same as a null resolution —
    # excluded by default, retained only if null_resolution_behavior="include".
    curated, exclusions, _ = curate_structure(
        structure, None, _policy(quality_rules=QualityRules(max_resolution=0.1))
    )
    assert curated is None
    assert [e.reason_code for e in exclusions] == ["NULL_RESOLUTION"]

    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(
            quality_rules=QualityRules(
                max_resolution=0.1, null_resolution_behavior="include"
            )
        ),
    )
    assert curated is not None
    assert exclusions == []


def test_min_polymer_chains_excludes_too_few_chains():
    structure = _load("1ayi")  # one polymer chain (A) + water (B)

    policy = _policy(quality_rules=QualityRules(min_polymer_chains=2))
    curated, exclusions, _ = curate_structure(structure, None, policy)
    assert curated is None
    assert [e.reason_code for e in exclusions] == ["TOO_FEW_CHAINS"]

    policy = _policy(quality_rules=QualityRules(min_polymer_chains=1))
    curated, exclusions, _ = curate_structure(structure, None, policy)
    assert curated is not None
    assert exclusions == []


def test_organism_filter_requires_taxonomy():
    structure = _load("1ayi")
    metadata = collect_metadata(structure)
    policy = _policy(organism_rules=OrganismRules(include_taxa=["9999999"]))

    curated, exclusions, _ = curate_structure(structure, metadata, policy)

    assert curated is None
    assert [e.reason_code for e in exclusions] in (
        ["MISSING_TAXONOMY"],
        ["ORGANISM_EXCLUDED"],
    )


def test_content_rules_strip_ligands():
    structure = _load("1ayi")
    has_hetatm = any(a.group_PDB == "HETATM" for a in structure.atoms)
    assert has_hetatm, "fixture must contain ligand atoms for this test"

    policy = _policy(
        content_rules=ContentRules(
            keep_ligands=False, keep_waters=False, keep_ions=False
        )
    )

    curated, exclusions, _ = curate_structure(structure, None, policy)

    assert exclusions == []
    assert not any(a.group_PDB == "HETATM" for a in curated.atoms)


def test_deduplicate_structures():
    structure = _load("1ayi")

    retained, removed, provenance = deduplicate_structures(
        [structure, structure], DeduplicationRules(enabled=True)
    )
    assert len(retained) == 1
    assert len(removed) == 1
    assert removed[0].reason_code == "DUPLICATE"
    assert provenance.enabled is True
    assert provenance.duplicates_found == 1

    retained, removed, provenance = deduplicate_structures(
        [structure, structure], DeduplicationRules(enabled=False)
    )
    assert len(retained) == 2
    assert removed == []
    assert provenance.enabled is False
    assert provenance.duplicates_found == 0


if __name__ == "__main__":
    test_curate_structure_pass_and_fail()
    test_missing_metadata_treated_as_unknown_not_skipped()
    test_min_polymer_chains_excludes_too_few_chains()
    test_organism_filter_requires_taxonomy()
    test_content_rules_strip_ligands()
    test_deduplicate_structures()
    print("ok")


def test_experimental_method_filter_matches_mmcif_method_values():
    structure = _load("104m")  # X-ray diffraction
    metadata = collect_metadata(structure)
    policy = _policy(
        quality_rules=QualityRules(
            include_experimental_methods=["X-RAY DIFFRACTION"]
        )
    )

    curated, exclusions, _ = curate_structure(structure, metadata, policy)

    assert exclusions == []
    assert curated is not None


def test_entry_exclusion_has_no_chain_id():
    structure = _load("1ayi")
    _, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(max_resolution=0.1)),
    )
    assert [(e.reason_code, e.chain_id) for e in exclusions] == [
        ("NULL_RESOLUTION", None)
    ]


def _meta(entry_id: str, **quality) -> MetadataRecord:
    return MetadataRecord(entry_id=entry_id, quality=QualityRecord(**quality))


def _codes(exclusions):
    return [e.reason_code for e in exclusions]


def test_cryo_em_entry_passes_resolution_rule_with_real_metadata():
    # Regression: cryo-EM resolution used to be None -> NULL_RESOLUTION.
    structure = _load("22jy")
    metadata = collect_metadata(structure)
    curated, exclusions, _ = curate_structure(
        structure,
        metadata,
        _policy(quality_rules=QualityRules(max_resolution=3.5)),
    )
    assert curated is not None and exclusions == []

    strict = QualityRules(
        max_resolution=3.5,
        max_resolution_by_method={"ELECTRON MICROSCOPY": 2.0},
    )
    curated, exclusions, _ = curate_structure(
        structure, metadata, _policy(quality_rules=strict)
    )
    assert curated is None
    assert _codes(exclusions) == ["RESOLUTION_THRESHOLD"]


def test_resolution_by_method_falls_back_to_max_resolution():
    structure = _load("1ayi")
    rules = QualityRules(
        max_resolution=3.0,
        max_resolution_by_method={"Electron Microscopy": 2.0},
    )
    meta = _meta(
        "1ayi", experimental_method="X-ray diffraction", resolution=2.4
    )
    curated, exclusions, _ = curate_structure(
        structure, meta, _policy(quality_rules=rules)
    )
    assert curated is not None and exclusions == []


def test_resolution_by_method_strictest_wins():
    structure = _load("1ayi")
    rules = QualityRules(
        max_resolution_by_method={
            " x-ray diffraction ": 2.5,
            "NEUTRON DIFFRACTION": 2.0,
        }
    )
    meta = _meta(
        "1ayi",
        experimental_method="X-ray diffraction; Neutron diffraction",
        resolution=2.4,
    )
    _, exclusions, _ = curate_structure(
        structure, meta, _policy(quality_rules=rules)
    )
    assert _codes(exclusions) == ["RESOLUTION_THRESHOLD"]


def test_rfactor_rules():
    structure = _load("1ayi")
    meta = _meta(
        "1ayi",
        experimental_method="X-ray diffraction",
        r_free=0.30,
        r_work=0.20,
        r_merge=0.12,
    )
    cases = [
        (QualityRules(max_r_free=0.25), ["RFREE_THRESHOLD"]),
        (QualityRules(max_r_free_gap=0.07), ["RFREE_GAP_THRESHOLD"]),
        (QualityRules(max_r_sym=0.10), ["RSYM_THRESHOLD"]),  # via r_merge
        (
            QualityRules(max_r_free=0.35, max_r_free_gap=0.11, max_r_sym=0.2),
            [],
        ),
    ]
    for rules, expected in cases:
        _, exclusions, _ = curate_structure(
            structure, meta, _policy(quality_rules=rules)
        )
        assert _codes(exclusions) == expected, rules


def test_rfactor_rules_skip_cryo_em():
    structure = _load("1ayi")
    meta = _meta(
        "1ayi", experimental_method="Electron Microscopy", resolution=2.0
    )
    rules = QualityRules(max_r_free=0.25, null_rfactor_behavior="exclude")
    curated, exclusions, _ = curate_structure(
        structure, meta, _policy(quality_rules=rules)
    )
    assert curated is not None and exclusions == []


def test_null_rfactor_behavior():
    structure = _load("1ayi")
    meta = _meta("1ayi", experimental_method="X-ray diffraction", r_free=None)
    curated, _, _ = curate_structure(
        structure, meta, _policy(quality_rules=QualityRules(max_r_free=0.25))
    )
    assert curated is not None  # default: include
    _, exclusions, _ = curate_structure(
        structure,
        meta,
        _policy(
            quality_rules=QualityRules(
                max_r_free=0.25, null_rfactor_behavior="exclude"
            )
        ),
    )
    assert _codes(exclusions) == ["NULL_RFACTOR"]


def test_rfactor_rules_skip_entries_without_metadata():
    structure = _load("1ayi")
    rules = QualityRules(max_r_free=0.25, null_rfactor_behavior="exclude")
    curated, exclusions, _ = curate_structure(
        structure, None, _policy(quality_rules=rules)
    )
    assert curated is not None and exclusions == []


def test_nonstandard_residues_exclude_entry_unless_allowed():
    structure = _load("1a08")
    _, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(quality_rules=QualityRules(exclude_nonstandard_residues=True)),
    )
    assert _codes(exclusions) == ["NONSTANDARD_RESIDUE"]
    assert "FTY" in exclusions[0].message

    curated, exclusions, _ = curate_structure(
        structure,
        None,
        _policy(
            quality_rules=QualityRules(
                exclude_nonstandard_residues=True,
                allowed_nonstandard_residues=["ACE", "DIP", "FTY"],
            )
        ),
    )
    assert curated is not None and exclusions == []


def test_mapped_mse_is_standard_for_curation():
    canon_policy = canonicalisationPolicy(
        policy_id="t",
        policy_name="t",
        policy_version="1",
        validation_rules=ValidationRules(fail_on_unresolved_issues=False),
        modified_residue_rules=ModifiedResidueRules(strategy="map_to_parent"),
    )
    canonical, _, _ = canonicalise_structure(_load("1b6w"), canon_policy)
    curated, exclusions, _ = curate_structure(
        canonical,
        None,
        _policy(quality_rules=QualityRules(exclude_nonstandard_residues=True)),
    )
    assert curated is not None and exclusions == []


def test_max_atoms():
    structure = _load("1ayi")  # 704 atoms
    _, exclusions, _ = curate_structure(
        structure, None, _policy(quality_rules=QualityRules(max_atoms=700))
    )
    assert _codes(exclusions) == ["TOO_MANY_ATOMS"]
    curated, _, _ = curate_structure(
        structure, None, _policy(quality_rules=QualityRules(max_atoms=704))
    )
    assert curated is not None
