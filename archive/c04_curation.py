"""Component 04 — Dataset Curation: public functions."""
from __future__ import annotations

from datetime import datetime, timezone

from archive.schemas.c03_metadata import AnnotatedStructureWithPlugins
from archive.schemas.c04_curation import (
    ChainDataset,
    ChainDatasetCounts,
    ChainDatasetProvenance,
    Dataset,
    DatasetCounts,
    DatasetCurationPolicy,
    DatasetProvenance,
    DeduplicationReport,
    InterfaceDataset,
    InterfaceDatasetCounts,
    InterfaceDatasetProvenance,
    ResidueDataset,
    ResidueDatasetCounts,
    ResidueDatasetProvenance,
)
from pandora.schemas.common import AppliedPolicyRef


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _policy_ref(policy: DatasetCurationPolicy) -> AppliedPolicyRef:
    return AppliedPolicyRef(
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
    )


# ── Public API ────────────────────────────────────────────────────────────────

def build_dataset(
    annotated_structures: list[AnnotatedStructureWithPlugins],
    policy: DatasetCurationPolicy,
    dataset_id: str = "pandora-dataset",
    dataset_name: str = "Pandora Dataset",
    dataset_version: str = "1.0.0",
) -> Dataset:
    """Apply selection, filtering, and deduplication to produce a structure-level Dataset."""
    # TODO: implement apply_selection_rules → apply_quality_filters →
    #   apply_content_filters → apply_organism_filters → deduplicate_dataset
    #   per policy sub-rules
    return Dataset(
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        dataset_version=dataset_version,
        structures=annotated_structures,
        counts=DatasetCounts(
            total_input=len(annotated_structures),
            total_selected=len(annotated_structures),
            total_excluded=0,
            total_duplicates_removed=0,
        ),
        deduplication_report=DeduplicationReport(
            enabled=policy.deduplication_rules.enabled,
            strategy=policy.deduplication_rules.strategy,
        ),
        applied_policy=_policy_ref(policy),
        provenance=DatasetProvenance(
            created_at=_now_iso(),
            source_count=len(annotated_structures),
        ),
    )


def extract_chains(dataset: Dataset, policy: DatasetCurationPolicy) -> ChainDataset:
    """Extract chain-level records from a structure-level Dataset."""
    # TODO: implement — iterate AnnotatedStructureWithPlugins, build ChainRecord
    #   per policy.extraction_rules.chain_extraction_rules
    now = _now_iso()
    return ChainDataset(
        dataset_id=f"{dataset.dataset_id}:chains",
        dataset_name=f"{dataset.dataset_name} — Chains",
        dataset_version=dataset.dataset_version,
        chains=[],
        source_dataset_id=dataset.dataset_id,
        counts=ChainDatasetCounts(
            total_structures_input=len(dataset.structures),
            total_chains_extracted=0,
            total_chains_excluded=0,
        ),
        applied_policy=_policy_ref(policy),
        provenance=ChainDatasetProvenance(
            created_at=now,
            source_dataset_id=dataset.dataset_id,
        ),
    )


def extract_interfaces(dataset: Dataset, policy: DatasetCurationPolicy) -> InterfaceDataset:
    """Extract interface-level records from a structure-level Dataset."""
    # TODO: implement — detect interfaces using BSA/contact criteria per
    #   policy.extraction_rules.interface_extraction_rules (requires FreeSASA)
    now = _now_iso()
    return InterfaceDataset(
        dataset_id=f"{dataset.dataset_id}:interfaces",
        dataset_name=f"{dataset.dataset_name} — Interfaces",
        dataset_version=dataset.dataset_version,
        interfaces=[],
        source_dataset_id=dataset.dataset_id,
        counts=InterfaceDatasetCounts(
            total_structures_input=len(dataset.structures),
            total_interfaces_extracted=0,
            total_interfaces_excluded=0,
        ),
        applied_policy=_policy_ref(policy),
        provenance=InterfaceDatasetProvenance(
            created_at=now,
            source_dataset_id=dataset.dataset_id,
        ),
    )


def extract_residues(
    source: Dataset | ChainDataset,
    policy: DatasetCurationPolicy,
) -> ResidueDataset:
    """Extract residue-level records from a structure- or chain-level dataset."""
    # TODO: implement — build ResidueRecord per
    #   policy.extraction_rules.residue_extraction_rules
    now = _now_iso()
    if isinstance(source, Dataset):
        source_granularity = "structure"
        source_count = len(source.structures)
    else:
        source_granularity = "chain"
        source_count = len(source.chains)
    return ResidueDataset(
        dataset_id=f"{source.dataset_id}:residues",
        dataset_name=f"{source.dataset_name} — Residues",
        dataset_version=source.dataset_version,
        residues=[],
        source_dataset_id=source.dataset_id,
        source_granularity=source_granularity,
        counts=ResidueDatasetCounts(
            total_source_units_input=source_count,
            total_residues_extracted=0,
            total_residues_excluded=0,
        ),
        applied_policy=_policy_ref(policy),
        provenance=ResidueDatasetProvenance(
            created_at=now,
            source_dataset_id=source.dataset_id,
        ),
    )
