"""Component 02 — Canonical Structure: public functions."""

from __future__ import annotations

from datetime import datetime, timezone

from archive.schemas.c01_ingestion import MmCIFIngestionResult, ParallelOptions
from archive.schemas.c02_canonicalization import (
    AltlocSelectionMapping,
    AssemblyMapping,
    CanonicalizationBatchResult,
    CanonicalizationBatchSummary,
    CanonicalizationPolicy,
    CanonicalMappings,
    CanonicalStructure,
    CanonicalStructureProvenance,
    CanonicalStructureResult,
    ChainIdMapping,
    EntityMapping,
    ResidueNumberMapping,
)
from pandora.schemas.common import AppliedPolicyRef, DiagnosticBundle


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _policy_ref(policy: CanonicalizationPolicy) -> AppliedPolicyRef:
    return AppliedPolicyRef(
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
    )


# ── Public API ────────────────────────────────────────────────────────────────


def canonicalize_structure(
    ingestion_result: MmCIFIngestionResult,
    policy: CanonicalizationPolicy,
) -> CanonicalStructureResult:
    """Convert a parsed structure into a canonical structure per policy.

    Full workflow (spec Section 5.1):
    normalize_chain_ids → normalize_residue_numbering → normalize_assemblies
    → handle_missing_data → resolve_altlocs → normalize_entities → filter_ligands
    → validate_canonical_structure
    """
    # TODO: implement each sub-step above
    prov = CanonicalStructureProvenance(
        provider=ingestion_result.provenance.provider,
        source_uri=ingestion_result.provenance.source_uri,
        retrieved_at=ingestion_result.provenance.retrieved_at,
    )
    if ingestion_result.status == "failed":
        return CanonicalStructureResult(
            entry_id=ingestion_result.entry_id,
            status="failed",
            canonical_structure=None,
            canonical_mappings=None,
            applied_policy=_policy_ref(policy),
            provenance=prov,
        )

    ps = ingestion_result.parsed_structure  # type: ignore[union-attr]
    canonical = CanonicalStructure(
        atoms=ps.atoms,
        residues=ps.residues,
        chains=ps.chains,
        entities=ps.entities,
        assemblies=ps.assemblies,
        ligands=ps.ligands,
    )
    mappings = CanonicalMappings(
        chain_id_mapping=ChainIdMapping(),
        residue_number_mapping=ResidueNumberMapping(),
        assembly_mapping=AssemblyMapping(),
        entity_mapping=EntityMapping(),
        altloc_selection_mapping=AltlocSelectionMapping(),
    )
    prov.canonicalized_at = _now_iso()
    return CanonicalStructureResult(
        entry_id=ingestion_result.entry_id,
        status=ingestion_result.status,
        canonical_structure=canonical,
        canonical_mappings=mappings,
        applied_policy=_policy_ref(policy),
        provenance=prov,
    )


def validate_canonical_structure(
    canonical_structure: CanonicalStructure,
    canonical_mappings: CanonicalMappings,
    policy: CanonicalizationPolicy,
) -> tuple[str, DiagnosticBundle]:
    # TODO: implement V1 validation rules per spec Section 5.2
    #   (CANONICAL_CHAIN_ID_COLLISION, RESIDUE_NUMBER_COLLISION, ...)
    return "valid", DiagnosticBundle()


def canonicalize_many_structures(
    ingestion_results: list[MmCIFIngestionResult],
    policy: CanonicalizationPolicy,
    mode: str = "sequential",
    parallel_options: ParallelOptions | None = None,
) -> CanonicalizationBatchResult:
    """Canonicalize a batch of ingestion results (sequential or parallel)."""
    # TODO: implement parallel mode via concurrent.futures
    results = [canonicalize_structure(r, policy) for r in ingestion_results]
    summary = CanonicalizationBatchSummary(
        total=len(results),
        success=sum(1 for r in results if r.status == "success"),
        warning=sum(1 for r in results if r.status == "warning"),
        failed=sum(1 for r in results if r.status == "failed"),
    )
    return CanonicalizationBatchResult(
        mode=mode, summary=summary, results=results
    )
