from __future__ import annotations

from pandora._util import now_iso
from pandora.canonicalisation import filter_ligands
from pandora.datasets.records import extract_chain_records
from pandora.schemas.canonicalisation import LigandRules
from pandora.schemas.common import DiagnosticBundle
from pandora.schemas.dataset import (
    ContentRules,
    CurationProvenance,
    DatasetCurationPolicy,
    DeduplicationProvenance,
    DeduplicationRules,
    ExclusionRecord,
    OrganismRules,
    QualityRules,
)
from pandora.schemas.metadata import MetadataRecord
from pandora.schemas.structure import Structure


def _check_quality(
    structure: Structure, metadata: MetadataRecord | None, rules: QualityRules
) -> ExclusionRecord | None:
    """ExclusionRecord if structure fails the resolution/method/chain-length
    quality rules, else None."""

    quality = metadata.quality if metadata else None
    resolution = quality.resolution if quality else None

    if resolution is None:
        if (
            rules.max_resolution is not None
            and rules.null_resolution_behavior == "exclude"
        ):
            return ExclusionRecord(
                entry_id=structure.entry_id,
                reason_code="NULL_RESOLUTION",
                message="resolution is null and "
                "null_resolution_behavior='exclude'",
            )
    elif rules.max_resolution is not None and resolution > rules.max_resolution:
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="RESOLUTION_THRESHOLD",
            message=f"resolution {resolution} exceeds "
            f"max_resolution {rules.max_resolution}",
        )

    method = quality.experimental_method if quality else None
    # Case/whitespace-insensitive: raw mmCIF tokens and policy values may
    # differ in casing (e.g. "X-RAY DIFFRACTION" vs "x-ray diffraction").
    normalized_method = method.strip().upper() if method else None
    include_methods = {
        m.strip().upper() for m in rules.include_experimental_methods
    }
    exclude_methods = {
        m.strip().upper() for m in rules.exclude_experimental_methods
    }
    if (include_methods and normalized_method not in include_methods) or (
        normalized_method in exclude_methods
    ):
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="METHOD_EXCLUDED",
            message=f"experimental_method={method!r} excluded by policy",
        )

    if rules.min_chain_length is not None:
        chain_lengths = [
            record.residue_count for record in extract_chain_records(structure)
        ]
        if max(chain_lengths, default=0) < rules.min_chain_length:
            return ExclusionRecord(
                entry_id=structure.entry_id,
                reason_code="CHAIN_TOO_SHORT",
                message="no polymer chain reaches min_chain_length="
                f"{rules.min_chain_length}",
            )

    return None


def _check_organism(
    structure: Structure, metadata: MetadataRecord | None, rules: OrganismRules
) -> ExclusionRecord | None:
    """ExclusionRecord if structure's taxonomy is excluded by the organism
    rules, else None."""

    if not rules.include_taxa and not rules.exclude_taxa:
        return None

    taxa = {
        str(taxon.ncbi_taxon_id)
        for taxon in (metadata.taxonomies if metadata else [])
        if taxon.ncbi_taxon_id is not None
    }
    if not taxa:
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="MISSING_TAXONOMY",
            message="organism filter is active but no taxonomy "
            "metadata is available",
        )

    if taxa & set(rules.exclude_taxa) or (
        rules.include_taxa and not taxa & set(rules.include_taxa)
    ):
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="ORGANISM_EXCLUDED",
            message=f"taxa={sorted(taxa)} excluded by organism_rules",
        )

    return None


def _apply_content_rules(
    structure: Structure, rules: ContentRules
) -> Structure:
    """Filter ligands/waters/ions out of structure per the content rules."""

    if rules.keep_ligands and rules.keep_waters and rules.keep_ions:
        return structure

    ligand_rules = LigandRules(
        strategy="filter",
        keep_waters=rules.keep_waters,
        keep_ions=rules.keep_ions,
        keep_nonpolymer_ligands=rules.keep_ligands,
    )
    atoms, asym_units = filter_ligands(
        list(structure.atoms),
        list(structure.asym_units),
        structure.entities,
        ligand_rules,
        DiagnosticBundle(),
        structure.entry_id,
    )
    return structure.model_copy(
        update={"atoms": atoms, "asym_units": asym_units}
    )


def deduplicate_structures(
    structures: list[Structure], rules: DeduplicationRules
) -> tuple[list[Structure], list[ExclusionRecord], DeduplicationProvenance]:
    """Remove structures sharing the same `Structure.entry_id`.

    Args:
        structures: Structures to deduplicate, in order. The first
            occurrence of each duplicate entry_id is kept.
        rules: Ignored (returns `structures` unchanged) if
            `rules.enabled` is False.

    Returns:
        `(retained, removed, provenance)` — the deduplicated
        structures, an `ExclusionRecord` (reason_code="DUPLICATE") per
        structure removed, and a record of how many duplicates were found.
    """

    if not rules.enabled:
        return (
            structures,
            [],
            DeduplicationProvenance(deduplicated_at=now_iso(), enabled=False),
        )

    seen: set[str] = set()
    retained: list[Structure] = []
    removed: list[ExclusionRecord] = []
    for structure in structures:
        if structure.entry_id in seen:
            removed.append(
                ExclusionRecord(
                    entry_id=structure.entry_id,
                    reason_code="DUPLICATE",
                    message="duplicate by entry_id",
                )
            )
            continue
        seen.add(structure.entry_id)
        retained.append(structure)

    provenance = DeduplicationProvenance(
        deduplicated_at=now_iso(),
        enabled=True,
        duplicates_found=len(removed),
    )
    return retained, removed, provenance


def curate_structure(
    structure: Structure,
    metadata: MetadataRecord | None,
    policy: DatasetCurationPolicy,
) -> tuple[Structure | None, ExclusionRecord | None, CurationProvenance]:
    """Apply quality, organism, and content rules to one structure.

    Args:
        structure: The canonical structure to curate.
        metadata: The structure's `MetadataRecord`, if available.
            Without it, quality/organism checks that depend on missing
            data (e.g. null resolution) apply their configured default
            rather than being skipped.
        policy: The curation policy governing every rule applied.

    Returns:
        `(curated_structure, None, provenance)` if the structure
        passes selection (with content rules applied), or
        `(None, exclusion, provenance)` if it was excluded. Provenance
        (which policy id/version ran, and when) is always populated,
        regardless of outcome.
    """

    provenance = CurationProvenance(
        curated_at=now_iso(),
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
    )

    exclusion = _check_quality(
        structure, metadata, policy.quality_rules
    ) or _check_organism(structure, metadata, policy.organism_rules)
    if exclusion is not None:
        return None, exclusion, provenance
    return (
        _apply_content_rules(structure, policy.content_rules),
        None,
        provenance,
    )
