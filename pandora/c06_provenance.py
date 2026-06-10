"""Component 06 — Provenance & Reproducibility: public functions."""
from __future__ import annotations

import uuid
from datetime import datetime, timezone

import pandora
from pandora.schemas.c05_splitting import LeakageSafeDataset
from pandora.schemas.c06_provenance import (
    ArtifactChecksums,
    ArtifactProvenance,
    DatasetSummary,
    ExportPolicy,
    ManifestChecksums,
    PandoraArtifact,
    PandoraManifest,
    PolicyProvenanceRecord,
    ProvenanceBundle,
    ProvenancePolicy,
    ReproducibilityReport,
    ReproducibilityReportSummary,
    SourceReleaseProvenanceRecord,
)
from pandora.schemas.common import AppliedPolicyRef


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _policy_ref(policy: ProvenancePolicy) -> AppliedPolicyRef:
    return AppliedPolicyRef(
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
    )


# ── Public API ────────────────────────────────────────────────────────────────

def assemble_provenance(
    artifact_id: str,
    leakage_safe_dataset: LeakageSafeDataset,
    policy: ProvenancePolicy,
) -> ProvenanceBundle:
    # TODO: implement — traverse the full LeakageSafeDataset object hierarchy,
    #   collect per-stage provenance into PipelineProvenance sub-objects
    return ProvenanceBundle(
        source_release_provenance=SourceReleaseProvenanceRecord(),
        policy_provenance=PolicyProvenanceRecord(),
    )


def generate_manifest(
    artifact_id: str,
    leakage_safe_dataset: LeakageSafeDataset,
    provenance: ProvenanceBundle,
    policy: ProvenancePolicy,
) -> PandoraManifest:
    # TODO: implement — compute real SHA-256 checksums for split artifacts
    ps = leakage_safe_dataset.partition_summary
    summary = DatasetSummary(
        dataset_id=leakage_safe_dataset.dataset_id,
        dataset_version=leakage_safe_dataset.dataset_version,
        granularity=leakage_safe_dataset.granularity,
        total_items=ps.train_count + ps.validation_count + ps.test_count,
        train_count=ps.train_count,
        validation_count=ps.validation_count,
        test_count=ps.test_count,
        train_fraction_achieved=ps.train_fraction_achieved,
        validation_fraction_achieved=ps.validation_fraction_achieved,
        test_fraction_achieved=ps.test_fraction_achieved,
    )
    return PandoraManifest(
        manifest_id=str(uuid.uuid4()),
        manifest_format="json",
        pandora_version=pandora.__version__,
        generated_at=_now_iso(),
        artifact_id=artifact_id,
        dataset_summary=summary,
        checksums=ManifestChecksums(),
    )


def finalize_artifact(
    artifact_id: str,
    leakage_safe_dataset: LeakageSafeDataset,
    provenance_policy: ProvenancePolicy,
    export_policy: ExportPolicy | None = None,
    artifact_name: str | None = None,
) -> PandoraArtifact:
    """Assemble provenance, generate manifest, and return a sealed PandoraArtifact."""
    # TODO: implement export per export_policy (YAML/JSON manifest, lineage graph,
    #   checksum bundle) per spec Section 5
    provenance = assemble_provenance(artifact_id, leakage_safe_dataset, provenance_policy)
    manifest = generate_manifest(
        artifact_id, leakage_safe_dataset, provenance, provenance_policy
    )
    report = ReproducibilityReport(
        report_id=str(uuid.uuid4()),
        artifact_id=artifact_id,
        summary=ReproducibilityReportSummary(
            pipeline_steps=6,
            source_count=0,
            policy_count=0,
            plugin_count=0,
        ),
    )
    return PandoraArtifact(
        artifact_id=artifact_id,
        artifact_name=artifact_name,
        leakage_safe_dataset=leakage_safe_dataset,
        provenance_bundle=provenance,
        manifest=manifest,
        checksums=ArtifactChecksums(),
        reproducibility_report=report,
        applied_policy=_policy_ref(provenance_policy),
        provenance=ArtifactProvenance(
            generated_at=_now_iso(),
            pandora_version=pandora.__version__,
        ),
    )
