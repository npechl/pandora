from __future__ import annotations

from pathlib import Path

from pandora.annotations.entry import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    annotate_structure_counts,
)
from pandora.annotations.pairwise import annotate_pairwise_sequence_identity
from pandora.canonicalisation import canonicalise_structure
from pandora.datasets.curation import curate_structure, deduplicate_structures
from pandora.datasets.records import entry_sequences
from pandora.export.mmcif import structure_to_mmcif
from pandora.ingestion.cache import find_cached
from pandora.ingestion.mmcif import fetch_mmcif
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.provenance.manifest import (
    build_dataset_manifest,
    build_provenance_bundle,
)
from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.canonicalisation import canonicalisationProvenance
from pandora.schemas.dataset import DeduplicationRules, ExclusionRecord
from pandora.schemas.ingestion import FetchOptions, IngestionProvenance
from pandora.schemas.provenance import DatasetManifest
from pandora.schemas.similarity import SimilarityRelationship
from pandora.schemas.structure import Structure
from pandora.similarity.clustering import cluster_similar_items
from pandora.similarity.partition import partition_dataset
from pandora.similarity.sequence import compute_sequence_similarity
from pandora.similarity.structure import compute_structure_similarity

ENTRY_ANNOTATION_DISPATCH = {
    "structure_counts": annotate_structure_counts,
    "ligand_contacts": annotate_ligand_contacts,
    "chain_interfaces": annotate_chain_interfaces,
}
PAIRWISE_ANNOTATION_DISPATCH = {
    "pairwise_sequence_identity": annotate_pairwise_sequence_identity,
}


def _reproduce_similarity(
    method_engine: str,
    method_parameters: dict,
    structures: dict[str, Structure],
    output_dir: Path,
) -> list[SimilarityRelationship]:
    if method_engine == "MMseqs2":
        return compute_sequence_similarity(
            entry_sequences(structures), **method_parameters
        )
    if method_engine == "Foldseek":
        paths: dict[str, Path] = {}
        for entry_id, structure in structures.items():
            paths[entry_id] = structure_to_mmcif(
                structure, output_dir / f"{entry_id}.reproduced.cif"
            )
        return compute_structure_similarity(paths, **method_parameters)
    raise ValueError(
        f"cannot reproduce similarity network: unknown engine {method_engine!r} "
        "(only 'MMseqs2' and 'Foldseek' can be auto-reproduced)"
    )


def _reproduce_annotations(
    manifest: DatasetManifest, structures: dict[str, Structure]
) -> dict[str, list[AnnotationLayer]]:
    layers_by_entry: dict[str, list[AnnotationLayer]] = {
        entry_id: [] for entry_id in structures
    }
    seen_pairs: set[frozenset[str]] = set()

    for bundle in manifest.structures:
        for record in bundle.annotations:
            if record.layer_type in ENTRY_ANNOTATION_DISPATCH:
                if bundle.entry_id not in structures:
                    continue
                try:
                    layer = ENTRY_ANNOTATION_DISPATCH[record.layer_type](
                        structures[bundle.entry_id], **record.parameters
                    )
                except TypeError as exc:
                    raise TypeError(
                        f"cannot reproduce annotation {record.layer_type!r} "
                        f"for entry {bundle.entry_id!r}: {exc}"
                    ) from exc
                layers_by_entry[bundle.entry_id].append(layer)
            elif record.layer_type in PAIRWISE_ANNOTATION_DISPATCH:
                pair = frozenset(record.target_ids)
                if len(pair) != 2 or pair in seen_pairs:
                    continue
                left_id, right_id = record.target_ids
                if left_id not in structures or right_id not in structures:
                    continue
                seen_pairs.add(pair)
                layer = PAIRWISE_ANNOTATION_DISPATCH[record.layer_type](
                    structures[left_id], structures[right_id]
                )
                layers_by_entry[left_id].append(layer)
                layers_by_entry[right_id].append(layer)

    return layers_by_entry


def reproduce_dataset(
    manifest: DatasetManifest,
    output_dir: str | Path,
    *,
    fetch_options: FetchOptions | None = None,
) -> tuple[dict[str, Structure], DatasetManifest]:
    """Re-run the pipeline described by a `DatasetManifest` from scratch.

    Re-fetches every structure via its recorded `ingestion` provenance,
    then replays canonicalisation, curation, deduplication, similarity/
    clustering, partitioning, and annotation exactly as the manifest
    describes them. Every step is optional in the same way it is in
    `build_dataset_manifest` (e.g. no `clustering` means no similarity
    network is (re)computed).

    This is a best-effort re-run, not a guaranteed byte-identical
    rebuild: source data can change upstream, and external similarity
    tools can drift between versions. The returned manifest is the
    record of what *this* run actually produced — diff it against the
    input manifest to see what changed.

    Args:
        manifest: A `DatasetManifest`, as returned by
            `build_dataset_manifest`.
        output_dir: Directory to fetch raw mmCIF files (and, for a
            Foldseek re-run, write canonicalised structures) into.
        fetch_options: Passed through to every `fetch_mmcif` call.

    Returns:
        `(structures, new_manifest)` — the reproduced, retained
        structures keyed by entry_id, and a fresh `DatasetManifest`
        describing this run.

    Raises:
        ValueError: A structure's `ProvenanceBundle` has no `ingestion`
            provenance to fetch from, or `clustering.similarity_method`
            names an engine that can't be auto-reproduced.
    """

    output_dir = Path(output_dir)
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

    excluded: list[ExclusionRecord] = []
    if manifest.curation_policy is not None:
        for entry_id in list(structures):
            metadata = collect_metadata(structures[entry_id])
            curated, exclusion, _ = curate_structure(
                structures[entry_id], metadata, manifest.curation_policy
            )
            if curated is None:
                excluded.append(exclusion)
                del structures[entry_id]
            else:
                structures[entry_id] = curated

    dedup_prov = None
    if manifest.deduplication is not None and manifest.deduplication.enabled:
        retained, removed, dedup_prov = deduplicate_structures(
            list(structures.values()), DeduplicationRules(enabled=True)
        )
        excluded.extend(removed)
        structures = {s.entry_id: s for s in retained}

    cluster_prov = None
    partition_prov = None
    splits: dict[str, list[str]] = {}
    if manifest.clustering is not None:
        method = manifest.clustering.similarity_method
        if method is None:
            raise ValueError(
                "cannot reproduce clustering: the original manifest did not "
                "record which similarity method produced the network"
            )
        relationships = _reproduce_similarity(
            method.engine, method.parameters or {}, structures, output_dir
        )
        clusters, cluster_prov = cluster_similar_items(
            list(structures), relationships, manifest.clustering.threshold
        )
        if manifest.partition is not None:
            splits, partition_prov = partition_dataset(
                clusters,
                manifest.partition.pct_train,
                manifest.partition.pct_val,
                manifest.partition.pct_test,
                keep_similar_items=manifest.partition.keep_similar_items,
            )

    annotation_layers = _reproduce_annotations(manifest, structures)
    new_structures = {
        entry_id: build_provenance_bundle(
            structures[entry_id],
            ingestion=ingestion_provenance.get(entry_id),
            canonicalisation=canonicalisation_provenance.get(entry_id),
            annotations=annotation_layers.get(entry_id, []),
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
        clustering=cluster_prov,
        partition=partition_prov,
        splits=splits,
        structures=list(new_structures.values()),
    )
    return structures, new_manifest
