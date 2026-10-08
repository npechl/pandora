from __future__ import annotations

from pandora._util import now_iso
from pandora.annotations.entry import polymer_asym_ids
from pandora.canonicalisation import filter_ligands
from pandora.canonicalisation.modified_residues import STANDARD_RESIDUES
from pandora.datasets.completeness import chain_completeness
from pandora.datasets.records import extract_chain_records
from pandora.schemas.canonicalisation import LigandRules
from pandora.schemas.common import DiagnosticBundle
from pandora.schemas.dataset import (
    ContentRules,
    CurationProvenance,
    DatasetCurationPolicy,
    DeduplicationProvenance,
    DeduplicationRules,
    ExclusionReason,
    ExclusionRecord,
    OrganismRules,
    QualityRules,
)
from pandora.schemas.metadata import MetadataRecord, QualityRecord
from pandora.schemas.structure import Structure


def _norm_method(value: str) -> str:
    """Method name upper-cased with runs of whitespace collapsed."""

    return " ".join(value.split()).upper()


def _upper_set(values: list[str]) -> set[str]:
    """Normalised (`_norm_method`) copies of values."""

    return {_norm_method(v) for v in values}


def _entry_methods(quality: QualityRecord | None) -> set[str]:
    """The entry's experimental methods, split on ';' and normalised."""

    method = quality.experimental_method if quality else None
    return _upper_set([m for m in (method or "").split(";") if m.strip()])


def _resolution_limit(methods: set[str], rules: QualityRules) -> float | None:
    """The strictest per-method limit that applies, else max_resolution."""

    by_method = {
        _norm_method(m): limit
        for m, limit in rules.max_resolution_by_method.items()
    }
    limits = [by_method[m] for m in methods if m in by_method]
    return min(limits) if limits else rules.max_resolution


def _check_rfactors(
    entry_id: str,
    quality: QualityRecord | None,
    methods: set[str],
    rules: QualityRules,
) -> ExclusionRecord | None:
    """ExclusionRecord if an R-factor rule fails, for rfactor_methods only."""

    if not methods & _upper_set(rules.rfactor_methods):
        return None
    r_free = quality.r_free if quality else None
    r_work = quality.r_work if quality else None
    r_sym = None
    if quality:
        r_sym = quality.r_sym if quality.r_sym is not None else quality.r_merge
    gap = (
        abs(r_free - r_work)
        if r_free is not None and r_work is not None
        else None
    )
    checks = [
        (rules.max_r_free, r_free, "RFREE_THRESHOLD", "r_free"),
        (rules.max_r_free_gap, gap, "RFREE_GAP_THRESHOLD", "|r_free - r_work|"),
        (rules.max_r_sym, r_sym, "RSYM_THRESHOLD", "r_sym"),
    ]
    for limit, value, code, label in checks:
        if limit is None:
            continue
        if value is None:
            if rules.null_rfactor_behavior == "exclude":
                return ExclusionRecord(
                    entry_id=entry_id,
                    reason_code="NULL_RFACTOR",
                    message=f"{label} is null and "
                    "null_rfactor_behavior='exclude'",
                )
            continue
        if value > limit:
            return ExclusionRecord(
                entry_id=entry_id,
                reason_code=code,
                message=f"{label} {value:.3f} exceeds {limit}",
            )
    return None


def _check_composition(
    structure: Structure, rules: QualityRules
) -> ExclusionRecord | None:
    """ExclusionRecord if the entry has non-standard residues or too many
    atoms, else None."""

    if rules.exclude_nonstandard_residues:
        allowed = STANDARD_RESIDUES | set(rules.allowed_nonstandard_residues)
        found = sorted(
            {
                a.label_comp_id
                for a in structure.atoms
                if a.label_seq_id is not None and a.label_comp_id not in allowed
            }
        )
        if found:
            return ExclusionRecord(
                entry_id=structure.entry_id,
                reason_code="NONSTANDARD_RESIDUE",
                message=f"non-standard polymer residues {found}",
            )
    if rules.max_atoms is not None and len(structure.atoms) > rules.max_atoms:
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="TOO_MANY_ATOMS",
            message=f"{len(structure.atoms)} atoms exceed "
            f"max_atoms={rules.max_atoms}",
        )
    return None


def _check_quality(
    structure: Structure, metadata: MetadataRecord | None, rules: QualityRules
) -> ExclusionRecord | None:
    """ExclusionRecord if structure fails the resolution/method/R-factor
    quality rules, else None."""

    quality = metadata.quality if metadata else None
    resolution = quality.resolution if quality else None
    methods = _entry_methods(quality)
    limit = _resolution_limit(methods, rules)

    if resolution is None:
        if limit is not None and rules.null_resolution_behavior == "exclude":
            return ExclusionRecord(
                entry_id=structure.entry_id,
                reason_code="NULL_RESOLUTION",
                message="resolution is null and "
                "null_resolution_behavior='exclude'",
            )
    elif limit is not None and resolution > limit:
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="RESOLUTION_THRESHOLD",
            message=f"resolution {resolution} exceeds limit {limit}",
        )

    method = quality.experimental_method if quality else None
    # Case/whitespace-insensitive, per method of a multi-method entry:
    # raw mmCIF tokens and policy values may differ in casing (e.g.
    # "X-RAY DIFFRACTION" vs "x-ray diffraction").
    include_methods = _upper_set(rules.include_experimental_methods)
    exclude_methods = _upper_set(rules.exclude_experimental_methods)
    if (include_methods and not methods & include_methods) or (
        methods & exclude_methods
    ):
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="METHOD_EXCLUDED",
            message=f"experimental_method={method!r} excluded by policy",
        )

    rfactor = _check_rfactors(structure.entry_id, quality, methods, rules)
    if rfactor is not None:
        return rfactor

    return None


def _chain_exclusions(
    structure: Structure, rules: QualityRules
) -> list[ExclusionRecord]:
    """One ExclusionRecord per polymer chain failing a chain rule."""

    failed: dict[str, ExclusionRecord] = {}

    def fail(chain_id: str, code: ExclusionReason, message: str) -> None:
        failed.setdefault(
            chain_id,
            ExclusionRecord(
                entry_id=structure.entry_id,
                chain_id=chain_id,
                reason_code=code,
                message=message,
            ),
        )

    if rules.max_chain_length is not None:
        for chain in extract_chain_records(structure):
            if chain.residue_count > rules.max_chain_length:
                fail(
                    chain.chain_id,
                    "CHAIN_TOO_LONG",
                    f"{chain.residue_count} residues exceed "
                    f"max_chain_length={rules.max_chain_length}",
                )

    max_tail = rules.max_missing_tail_fraction
    max_middle = rules.max_missing_middle_fraction
    if max_tail is not None or max_middle is not None:
        records, _ = chain_completeness(
            structure, rules.missing_residue_definition
        )
        by_chain = {r.chain_id: r for r in records}
        for chain_id in sorted(polymer_asym_ids(structure)):
            c = by_chain.get(chain_id)
            if c is None:
                fail(
                    chain_id,
                    "NO_SEQRES",
                    "completeness could not be measured against SEQRES",
                )
                continue
            tails = c.missing_n_term + c.missing_c_term
            tail_fraction = tails / c.seqres_length
            core = c.seqres_length - tails
            middle_fraction = c.missing_middle / core if core else 0.0
            if max_tail is not None and tail_fraction > max_tail:
                fail(
                    chain_id,
                    "MISSING_TAILS",
                    f"missing tails {tail_fraction:.3f} exceed "
                    f"max_missing_tail_fraction={max_tail}",
                )
            elif max_middle is not None and middle_fraction > max_middle:
                fail(
                    chain_id,
                    "MISSING_MIDDLE",
                    f"missing middle {middle_fraction:.3f} exceeds "
                    f"max_missing_middle_fraction={max_middle}",
                )

    return [failed[c] for c in sorted(failed)]


def _drop_chains(structure: Structure, chain_ids: set[str]) -> Structure:
    """Copy of structure without the atoms and asym units of chain_ids."""

    return structure.model_copy(
        update={
            "atoms": [
                a for a in structure.atoms if a.label_asym_id not in chain_ids
            ],
            "asym_units": [
                u for u in structure.asym_units if u.id not in chain_ids
            ],
        }
    )


def _check_chain_counts(
    structure: Structure, rules: QualityRules, chains_dropped: bool
) -> ExclusionRecord | None:
    """ExclusionRecord if too few or too short chains are left, else None."""

    if chains_dropped and not polymer_asym_ids(structure):
        return ExclusionRecord(
            entry_id=structure.entry_id,
            reason_code="NO_CHAINS_LEFT",
            message="every polymer chain was removed by the chain rules",
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

    if rules.min_polymer_chains is not None:
        chain_count = len(polymer_asym_ids(structure))
        if chain_count < rules.min_polymer_chains:
            return ExclusionRecord(
                entry_id=structure.entry_id,
                reason_code="TOO_FEW_CHAINS",
                message=f"{chain_count} polymer chain(s) found, fewer than "
                f"min_polymer_chains={rules.min_polymer_chains}",
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
) -> tuple[Structure | None, list[ExclusionRecord], CurationProvenance]:
    """Apply quality, organism, and content rules to one structure.

    Args:
        structure: The canonical structure to curate.
        metadata: The structure's `MetadataRecord`, if available.
            Without it, quality/organism checks that depend on missing
            data (e.g. null resolution) apply their configured default
            rather than being skipped.
        policy: The curation policy governing every rule applied.

    Returns:
        `(curated_structure, exclusions, provenance)`. When the entry
        is kept, `curated_structure` has content rules applied and
        failing chains removed, and `exclusions` holds one record per
        removed chain. When the entry is excluded, `curated_structure`
        is None and the last record in `exclusions` has
        `chain_id=None`. Provenance (which policy id/version ran, and
        when) is always populated.
    """

    provenance = CurationProvenance(
        curated_at=now_iso(),
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
    )

    rules = policy.quality_rules
    exclusion = (
        _check_quality(structure, metadata, rules)
        or _check_organism(structure, metadata, policy.organism_rules)
        or _check_composition(structure, rules)
    )
    if exclusion is not None:
        return None, [exclusion], provenance

    records = _chain_exclusions(structure, rules)
    if records:
        structure = _drop_chains(structure, {r.chain_id for r in records})
    exclusion = _check_chain_counts(structure, rules, bool(records))
    if exclusion is not None:
        return None, [*records, exclusion], provenance

    return (
        _apply_content_rules(structure, policy.content_rules),
        records,
        provenance,
    )
