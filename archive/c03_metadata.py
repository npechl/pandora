"""Component 03 — Metadata & Annotation: public functions."""

from __future__ import annotations

from datetime import datetime, timezone

from archive.schemas.c02_canonicalization import CanonicalStructureResult
from archive.schemas.c03_metadata import (
    AnnotatedStructureProvenance,
    AnnotatedStructureWithPlugins,
    AnnotationPluginPolicy,
    MetadataAndAnnotationBatchResult,
    MetadataAndAnnotationBatchResultItem,
    MetadataAndAnnotationBatchSummary,
    MetadataAnnotatedStructure,
    MetadataAnnotatedStructureProvenance,
    MetadataAnnotations,
    MetadataIntegrationPolicy,
    MetadataRecord,
    MetadataRetrievalStatus,
    Plugin,
)
from pandora.schemas.common import AppliedPolicyRef


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _policy_ref(policy: MetadataIntegrationPolicy) -> AppliedPolicyRef:
    return AppliedPolicyRef(
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
    )


# ── Public API ────────────────────────────────────────────────────────────────


def retrieve_metadata(
    entry_id: str, policy: MetadataIntegrationPolicy
) -> MetadataRecord:
    # TODO: implement — fetch from PDBe, SIFTS, UniProt, taxonomy APIs per
    #   policy.include_sources and policy.include_categories
    return MetadataRecord(
        entry_id=entry_id,
        retrieval_status=MetadataRetrievalStatus(),
    )


def attach_metadata(
    canonical_structure_result: CanonicalStructureResult,
    policy: MetadataIntegrationPolicy,
) -> MetadataAnnotatedStructure:
    """Retrieve metadata and attach it to a canonical structure result."""
    # TODO: implement full metadata retrieval and merging
    record = retrieve_metadata(canonical_structure_result.entry_id, policy)
    annotations = MetadataAnnotations(
        archive_metadata=record.archive_metadata,
        biological_mappings=record.biological_mappings,
        structural_annotations=record.structural_annotations,
        provenance_metadata=record.provenance_metadata,
    )
    return MetadataAnnotatedStructure(
        canonical_structure_result=canonical_structure_result,
        metadata_annotations=annotations,
        applied_metadata_policy=_policy_ref(policy),
        provenance=MetadataAnnotatedStructureProvenance(
            retrieved_at=_now_iso()
        ),
    )


def attach_plugins(
    metadata_annotated_structure: MetadataAnnotatedStructure,
    plugin_policy: AnnotationPluginPolicy,
    plugins: list[Plugin],
) -> AnnotatedStructureWithPlugins:
    """Execute annotation plugins and collect derived annotation layers."""
    # TODO: implement — execute each plugin, collect AnnotationLayer results,
    #   handle plugin failures per plugin_policy.execution_rules
    return AnnotatedStructureWithPlugins(
        canonical_structure_result=metadata_annotated_structure.canonical_structure_result,
        metadata_annotations=metadata_annotated_structure.metadata_annotations,
        derived_annotations=[],
        applied_metadata_policy=metadata_annotated_structure.applied_metadata_policy,
        applied_plugins=[],
        provenance=AnnotatedStructureProvenance(retrieved_at=_now_iso()),
    )


def attach_metadata_many(
    canonical_results: list[CanonicalStructureResult],
    metadata_policy: MetadataIntegrationPolicy,
    plugin_policy: AnnotationPluginPolicy | None = None,
    plugins: list[Plugin] | None = None,
    mode: str = "sequential",
) -> MetadataAndAnnotationBatchResult:
    """Attach metadata (and optionally plugins) to a batch of canonical results."""
    # TODO: implement parallel mode via concurrent.futures
    items: list[MetadataAndAnnotationBatchResultItem] = []
    for result in canonical_results:
        if result.status == "failed":
            items.append(
                MetadataAndAnnotationBatchResultItem(
                    entry_id=result.entry_id,
                    status="failed",
                    annotated_structure=None,
                )
            )
            continue
        meta = attach_metadata(result, metadata_policy)
        if plugin_policy is not None and plugins:
            annotated = attach_plugins(meta, plugin_policy, plugins)
        else:
            annotated = AnnotatedStructureWithPlugins(
                canonical_structure_result=meta.canonical_structure_result,
                metadata_annotations=meta.metadata_annotations,
                applied_metadata_policy=meta.applied_metadata_policy,
                provenance=AnnotatedStructureProvenance(
                    retrieved_at=_now_iso()
                ),
            )
        items.append(
            MetadataAndAnnotationBatchResultItem(
                entry_id=result.entry_id,
                status=result.status,
                annotated_structure=annotated,
            )
        )
    summary = MetadataAndAnnotationBatchSummary(
        total=len(items),
        success=sum(1 for i in items if i.status == "success"),
        warning=sum(1 for i in items if i.status == "warning"),
        failed=sum(1 for i in items if i.status == "failed"),
    )
    return MetadataAndAnnotationBatchResult(
        mode=mode, summary=summary, results=items
    )
