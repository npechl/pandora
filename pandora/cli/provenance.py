"""CLI `manifest` and `reproduce`: dataset manifests and rebuilds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pandora.cli._io import (
    load_annotations,
    load_json_dict,
    load_json_list,
    load_json_model,
    load_yaml_model,
    read_structures_dir,
    write_structures_dir,
)
from pandora.export import write_json
from pandora.ingestion.policy import load_policy
from pandora.provenance.manifest import (
    build_dataset_manifest,
    build_provenance_bundle,
)
from pandora.schemas.canonicalisation import (
    canonicalisationProvenance,
)
from pandora.schemas.dataset import (
    DatasetCurationPolicy,
    DeduplicationProvenance,
    ExclusionRecord,
)
from pandora.schemas.ingestion import FetchOptions, IngestionProvenance
from pandora.schemas.provenance import DatasetManifest
from pandora.schemas.similarity import (
    ClusteringProvenance,
    PartitionProvenance,
)


def _cmd_manifest(args: argparse.Namespace) -> None:
    """Handle the `manifest` subcommand: assemble a DatasetManifest from
    prior-stage outputs."""

    structures = read_structures_dir(Path(args.structures_dir))

    canonicalisation_policy = (
        load_policy(args.canonicalisation_policy)
        if args.canonicalisation_policy
        else None
    )
    curation_policy = (
        load_yaml_model(DatasetCurationPolicy, Path(args.curation_policy))
        if args.curation_policy
        else None
    )
    ingestion_prov = (
        load_json_dict(IngestionProvenance, Path(args.ingestion_provenance))
        if args.ingestion_provenance
        else {}
    )
    canon_prov = (
        load_json_dict(
            canonicalisationProvenance, Path(args.canonicalisation_provenance)
        )
        if args.canonicalisation_provenance
        else {}
    )
    annotations = (
        load_annotations(Path(args.annotations)) if args.annotations else {}
    )
    exclusions: list[ExclusionRecord] = []
    for path in args.exclusions or []:
        exclusions.extend(load_json_list(ExclusionRecord, Path(path)))

    dedup_prov = (
        load_json_model(DeduplicationProvenance, Path(args.dedup_provenance))
        if args.dedup_provenance
        else None
    )
    cluster_prov = (
        load_json_model(ClusteringProvenance, Path(args.cluster_provenance))
        if args.cluster_provenance
        else None
    )
    partition_prov = (
        load_json_model(PartitionProvenance, Path(args.partition_provenance))
        if args.partition_provenance
        else None
    )
    splits = json.loads(Path(args.splits).read_text()) if args.splits else {}

    bundles = [
        build_provenance_bundle(
            structure,
            ingestion=ingestion_prov.get(entry_id),
            canonicalisation=canon_prov.get(entry_id),
            annotations=annotations.get(entry_id, []),
        )
        for entry_id, structure in structures.items()
    ]

    manifest = build_dataset_manifest(
        dataset_id=args.dataset_id,
        dataset_name=args.dataset_name,
        dataset_version=args.dataset_version,
        curation_policy=curation_policy,
        canonicalisation_policy=canonicalisation_policy,
        excluded=exclusions,
        deduplication=dedup_prov,
        clustering=cluster_prov,
        partition=partition_prov,
        splits=splits,
        structures=bundles,
    )
    write_json(manifest, Path(args.output))
    print(f"manifest for {len(bundles)} structures -> {args.output}")


def _cmd_reproduce(args: argparse.Namespace) -> None:
    """Handle the `reproduce` subcommand: replay a DatasetManifest from
    scratch."""

    from pandora.provenance.reproduce import reproduce_dataset

    manifest = DatasetManifest.model_validate_json(
        Path(args.manifest).read_text()
    )
    output_dir = Path(args.output_dir)
    fetch_options = FetchOptions(use_cache=not args.no_cache)

    structures, new_manifest = reproduce_dataset(
        manifest, output_dir, fetch_options=fetch_options
    )

    write_structures_dir(structures, output_dir / "structures")
    write_json(new_manifest, output_dir / "reproduced_manifest.json")
    print(f"reproduced {len(structures)} structures -> {output_dir}")


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    p = subparsers.add_parser(
        "manifest", help="Assemble a DatasetManifest from prior-stage outputs."
    )
    p.add_argument("--dataset-id", required=True)
    p.add_argument("--dataset-name", required=True)
    p.add_argument("--dataset-version", required=True)
    p.add_argument("--structures-dir", required=True)
    p.add_argument("--canonicalisation-policy")
    p.add_argument("--curation-policy")
    p.add_argument("--ingestion-provenance")
    p.add_argument("--canonicalisation-provenance")
    p.add_argument("--exclusions", nargs="+")
    p.add_argument("--dedup-provenance")
    p.add_argument("--cluster-provenance")
    p.add_argument("--partition-provenance")
    p.add_argument("--splits")
    p.add_argument("--annotations")
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_manifest)

    p = subparsers.add_parser(
        "reproduce", help="Replay a DatasetManifest from scratch."
    )
    p.add_argument("--manifest", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--no-cache", action="store_true")
    p.set_defaults(func=_cmd_reproduce)
