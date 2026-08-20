from pathlib import Path

from pandora.datasets import curate_structure, deduplicate_structures
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import (
    ContentRules,
    DatasetCurationPolicy,
    DeduplicationRules,
    OrganismRules,
    QualityRules,
)

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

    curated, exclusion, provenance = curate_structure(
        structure, metadata, policy
    )
    assert curated is not None
    assert exclusion is None
    assert provenance.policy_id == "test"
    assert provenance.policy_version == "1.0.0"

    tight_policy = _policy(
        quality_rules=QualityRules(
            max_resolution=metadata.quality.resolution - 0.01
        )
    )
    curated, exclusion, provenance = curate_structure(
        structure, metadata, tight_policy
    )
    assert curated is None
    assert exclusion.reason_code == "RESOLUTION_THRESHOLD"
    # provenance is populated on exclusion too, so callers always know
    # which policy made the call.
    assert provenance.policy_id == "test"


def test_missing_metadata_treated_as_unknown_not_skipped():
    structure = _load("1ayi")

    curated, exclusion, _ = curate_structure(structure, None, _policy())
    assert curated is not None
    assert exclusion is None

    # No metadata means resolution is unknown, same as a null resolution —
    # excluded by default, retained only if null_resolution_behavior="include".
    curated, exclusion, _ = curate_structure(
        structure, None, _policy(quality_rules=QualityRules(max_resolution=0.1))
    )
    assert curated is None
    assert exclusion.reason_code == "NULL_RESOLUTION"

    curated, exclusion, _ = curate_structure(
        structure,
        None,
        _policy(
            quality_rules=QualityRules(
                max_resolution=0.1, null_resolution_behavior="include"
            )
        ),
    )
    assert curated is not None
    assert exclusion is None


def test_organism_filter_requires_taxonomy():
    structure = _load("1ayi")
    metadata = collect_metadata(structure)
    policy = _policy(organism_rules=OrganismRules(include_taxa=["9999999"]))

    curated, exclusion, _ = curate_structure(structure, metadata, policy)

    assert curated is None
    assert exclusion.reason_code in ("MISSING_TAXONOMY", "ORGANISM_EXCLUDED")


def test_content_rules_strip_ligands():
    structure = _load("1ayi")
    has_hetatm = any(a.group_PDB == "HETATM" for a in structure.atoms)
    assert has_hetatm, "fixture must contain ligand atoms for this test"

    policy = _policy(
        content_rules=ContentRules(
            keep_ligands=False, keep_waters=False, keep_ions=False
        )
    )

    curated, exclusion, _ = curate_structure(structure, None, policy)

    assert exclusion is None
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
    test_organism_filter_requires_taxonomy()
    test_content_rules_strip_ligands()
    test_deduplicate_structures()
    print("ok")
