from __future__ import annotations

from pathlib import Path

from pandora.canonicalisation import canonicalise_structure
from pandora.datasets.curation import curate_structure, deduplicate_structures
from pandora.ingestion.cache import find_cached
from pandora.ingestion.mmcif import fetch_mmcif
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.provenance.manifest import (
    build_dataset_manifest,
    build_provenance_bundle,
)
from pandora.schemas.canonicalisation import canonicalisationProvenance
from pandora.schemas.dataset import (
    DatasetCurationPolicy,
    DeduplicationRules,
    ExclusionRecord,
)
from pandora.schemas.ingestion import FetchOptions, IngestionProvenance
from pandora.schemas.provenance import DatasetManifest
from pandora.schemas.structure import Structure


def reproduce_dataset(
    manifest: DatasetManifest,
    output_dir: str | Path,
    *,
    fetch_options: FetchOptions | None = None,
) -> tuple[dict[str, Structure], DatasetManifest]:
    """Re-run the pipeline described by a `DatasetManifest` from scratch.

    Re-fetches every structure via its recorded `ingestion` provenance,
    then replays canonicalisation, curation and deduplication as the
    manifest describes them. Similarity, clustering, partitioning and
    annotation layers are not replayed: the new manifest has no
    `clustering`, `partition` or `splits`, and its bundles carry no
    annotations.

    This is a best-effort re-run, not a guaranteed byte-identical
    rebuild: source data can change upstream. The returned manifest is the
    record of what *this* run actually produced — diff it against the
    input manifest to see what changed.

    Args:
        manifest: A `DatasetManifest`, as returned by
            `build_dataset_manifest`.
        output_dir: Directory to fetch raw mmCIF files into.
        fetch_options: Passed through to every `fetch_mmcif` call.

    Returns:
        `(structures, new_manifest)` — the reproduced, retained
        structures keyed by entry_id, and a fresh `DatasetManifest`
        describing this run.

    Raises:
        ValueError: A structure's `ProvenanceBundle` has no `ingestion`
            provenance to fetch from.
    """

    output_dir = Path(output_dir)
    structures, ingestion_provenance, canonicalisation_provenance = (
        _reproduce_structures(manifest, output_dir, fetch_options)
    )

    excluded: list[ExclusionRecord] = []
    if manifest.curation_policy is not None:
        structures, excluded = _reproduce_curation(
            structures, manifest.curation_policy
        )

    dedup_prov = None
    if manifest.deduplication is not None and manifest.deduplication.enabled:
        retained, removed, dedup_prov = deduplicate_structures(
            list(structures.values()), DeduplicationRules(enabled=True)
        )
        excluded.extend(removed)
        structures = {s.entry_id: s for s in retained}

    new_structures = {
        entry_id: build_provenance_bundle(
            structures[entry_id],
            ingestion=ingestion_provenance.get(entry_id),
            canonicalisation=canonicalisation_provenance.get(entry_id),
        )
        for entry_id in structures
    }

    new_manifest = build_dataset_manifest(
        dataset_id=manifest.dataset_id,
        dataset_name=manifest.dataset_name,
        dataset_version=manifest.dataset_version,
        curation_policy=manifest.curation_policy,
        canonicalisation_policy=manifest.canonicalisation_policy,
        excluded=excluded,
        deduplication=dedup_prov,
        structures=list(new_structures.values()),
    )
    return structures, new_manifest


def _reproduce_structures(
    manifest: DatasetManifest,
    output_dir: Path,
    fetch_options: FetchOptions | None,
) -> tuple[
    dict[str, Structure],
    dict[str, IngestionProvenance],
    dict[str, canonicalisationProvenance],
]:
    """Re-fetch, parse and (if the manifest has a policy) canonicalise
    every structure; each result is keyed by entry_id."""

    structures: dict[str, Structure] = {}
    ingestion_provenance: dict[str, IngestionProvenance] = {}
    canonicalisation_provenance: dict[str, canonicalisationProvenance] = {}

    for bundle in manifest.structures:
        if bundle.ingestion is None:
            raise ValueError(
                f"cannot reproduce {bundle.entry_id!r}: its ProvenanceBundle "
                "has no ingestion provenance to fetch from"
            )
        ingestion_prov = fetch_mmcif(
            bundle.entry_id,
            bundle.ingestion.provider,
            bundle.ingestion.source_uri,
            output_dir,
            fetch_options,
        )
        cached_path = find_cached(bundle.entry_id, output_dir)
        structure, _, _ = mmcif_to_structure(str(cached_path))
        if manifest.canonicalisation_policy is not None:
            structure, _, canon_prov = canonicalise_structure(
                structure, manifest.canonicalisation_policy
            )
            canonicalisation_provenance[structure.entry_id] = canon_prov
        # Keyed by structure.entry_id (not bundle.entry_id) from here on,
        # matching curate_structure/deduplicate_structures/ExclusionRecord,
        # which all key by Structure.entry_id too — the two only differ if
        # the source data's own entry_id changed since the manifest was
        # built, which is exactly the kind of drift reproduction tolerates.
        ingestion_provenance[structure.entry_id] = ingestion_prov
        structures[structure.entry_id] = structure

    return structures, ingestion_provenance, canonicalisation_provenance


def _reproduce_curation(
    structures: dict[str, Structure], policy: DatasetCurationPolicy
) -> tuple[dict[str, Structure], list[ExclusionRecord]]:
    """Structures that pass the curation policy (curated), and the
    exclusion records of the rest."""

    kept: dict[str, Structure] = {}
    excluded: list[ExclusionRecord] = []
    for entry_id, structure in structures.items():
        metadata = collect_metadata(structure)
        curated, records, _ = curate_structure(structure, metadata, policy)
        excluded.extend(records)
        if curated is not None:
            kept[entry_id] = curated
    return kept, excluded
