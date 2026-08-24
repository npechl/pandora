from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

import yaml
from pydantic import BaseModel

from pandora.annotations import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    annotate_pairwise_sequence_identity,
    annotate_structure_counts,
)
from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import (
    curate_structure,
    deduplicate_structures,
    entry_sequences,
)
from pandora.export import structure_to_mmcif, write_json, write_records
from pandora.ingestion.policy import load_policy
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.provenance.manifest import (
    build_dataset_manifest,
    build_provenance_bundle,
)
from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.canonicalisation import canonicalisationProvenance
from pandora.schemas.dataset import (
    DatasetCurationPolicy,
    DeduplicationProvenance,
    DeduplicationRules,
    ExclusionRecord,
)
from pandora.schemas.ingestion import FetchOptions, IngestionProvenance
from pandora.schemas.provenance import DatasetManifest
from pandora.schemas.similarity import (
    ClusteringProvenance,
    PartitionProvenance,
    SimilarityCluster,
    SimilarityRelationship,
)
from pandora.schemas.structure import Structure
from pandora.similarity import (
    cluster_similar_items,
    compute_sequence_similarity,
    compute_structure_similarity,
    partition_dataset,
)

_ENTRY_ANNOTATORS = {
    "structure_counts": annotate_structure_counts,
    "ligand_contacts": annotate_ligand_contacts,
    "chain_interfaces": annotate_chain_interfaces,
}


# I/O helpers -----------------------------------------------------------


def _read_structures_dir(input_dir: Path) -> dict[str, Structure]:
    """Parse every `*.cif` in `input_dir`, keyed by `Structure.entry_id`.

    Files that fail to parse are skipped with a warning rather than
    aborting the whole batch.
    """

    structures: dict[str, Structure] = {}
    for path in sorted(input_dir.glob("*.cif")):
        structure, _, status = mmcif_to_structure(str(path))
        if structure is None:
            print(f"warning: failed to parse {path}: {status}", file=sys.stderr)
            continue
        structures[structure.entry_id] = structure
    return structures


def _write_structures_dir(
    structures: dict[str, Structure], output_dir: Path
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for entry_id, structure in structures.items():
        structure_to_mmcif(structure, output_dir / f"{entry_id.lower()}.cif")


def _write_json_dict(data: dict[str, BaseModel], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({k: v.model_dump() for k, v in data.items()}, indent=2),
        encoding="utf-8",
    )


def _load_json_model(model_cls: type[BaseModel], path: Path) -> BaseModel:
    return model_cls.model_validate(json.loads(path.read_text()))


def _load_json_list(model_cls: type[BaseModel], path: Path) -> list:
    data = json.loads(path.read_text())
    return [model_cls.model_validate(item) for item in data]


def _load_json_dict(model_cls: type[BaseModel], path: Path) -> dict:
    data = json.loads(path.read_text())
    return {k: model_cls.model_validate(v) for k, v in data.items()}


def _load_annotations(path: Path) -> dict[str, list[AnnotationLayer]]:
    data = json.loads(path.read_text())
    return {
        entry_id: [AnnotationLayer.model_validate(item) for item in layers]
        for entry_id, layers in data.items()
    }


def _load_yaml_model(model_cls: type[BaseModel], path: Path) -> BaseModel:
    with path.open() as stream:
        data = yaml.safe_load(stream)
    return model_cls.model_validate(data)


# Subcommands -------------------------------------------------------------


def _cmd_fetch(args: argparse.Namespace) -> None:
    from pandora.ingestion.mmcif import fetch_mmcif

    output_dir = Path(args.output_dir)
    fetch_options = FetchOptions(use_cache=not args.no_cache)

    provenance: dict[str, IngestionProvenance] = {}
    for entry_id in args.entry_ids:
        try:
            provenance[entry_id.upper()] = fetch_mmcif(
                entry_id, args.provider, None, output_dir, fetch_options
            )
        except (RuntimeError, ValueError) as exc:
            if not args.allow_partial:
                raise
            print(
                f"warning: failed to fetch {entry_id}: {exc}", file=sys.stderr
            )

    _write_json_dict(provenance, output_dir / "ingestion_provenance.json")
    print(
        f"fetched {len(provenance)}/{len(args.entry_ids)} entries -> {output_dir}"
    )


def _cmd_canonicalise(args: argparse.Namespace) -> None:
    output_dir = Path(args.output_dir)
    policy = load_policy(args.policy)
    structures = _read_structures_dir(Path(args.input_dir))

    canonical: dict[str, Structure] = {}
    provenance: dict[str, canonicalisationProvenance] = {}
    for entry_id, structure in structures.items():
        canonical_structure, _, prov = canonicalise_structure(structure, policy)
        canonical[entry_id] = canonical_structure
        provenance[entry_id] = prov

    _write_structures_dir(canonical, output_dir)
    _write_json_dict(
        provenance, output_dir / "canonicalisation_provenance.json"
    )
    print(f"canonicalised {len(canonical)} structures -> {output_dir}")


def _cmd_curate(args: argparse.Namespace) -> None:
    output_dir = Path(args.output_dir)
    policy = _load_yaml_model(DatasetCurationPolicy, Path(args.policy))
    structures = _read_structures_dir(Path(args.input_dir))

    retained: dict[str, Structure] = {}
    provenance = {}
    exclusions: list[ExclusionRecord] = []
    for entry_id, structure in structures.items():
        metadata = collect_metadata(structure)
        curated, exclusion, prov = curate_structure(structure, metadata, policy)
        provenance[entry_id] = prov
        if curated is None:
            exclusions.append(exclusion)
        else:
            retained[entry_id] = curated

    _write_structures_dir(retained, output_dir)
    _write_json_dict(provenance, output_dir / "curation_provenance.json")
    if exclusions:
        write_records(exclusions, output_dir / "curation_exclusions.json")
    print(
        f"curated: {len(retained)} retained, {len(exclusions)} excluded "
        f"-> {output_dir}"
    )


def _cmd_dedup(args: argparse.Namespace) -> None:
    output_dir = Path(args.output_dir)
    structures = _read_structures_dir(Path(args.input_dir))
    rules = DeduplicationRules(enabled=not args.disable)

    retained, removed, prov = deduplicate_structures(
        list(structures.values()), rules
    )

    _write_structures_dir({s.entry_id: s for s in retained}, output_dir)
    write_json(prov, output_dir / "dedup_provenance.json")
    if removed:
        write_records(removed, output_dir / "dedup_exclusions.json")
    print(
        f"dedup: {len(retained)} retained, {len(removed)} removed "
        f"-> {output_dir}"
    )


def _cmd_similarity(args: argparse.Namespace) -> None:
    output = Path(args.output)

    if args.engine == "mmseqs2":
        structures = _read_structures_dir(Path(args.input_dir))
        sequences = entry_sequences(structures)
        kwargs = {"mmseqs_bin": args.mmseqs_bin}
        if args.sensitivity is not None:
            kwargs["sensitivity"] = args.sensitivity
        relationships = compute_sequence_similarity(sequences, **kwargs)
    else:
        kwargs = {
            "foldseek_bin": args.foldseek_bin,
            "alignment_type": args.alignment_type,
        }
        if args.sensitivity is not None:
            kwargs["sensitivity"] = args.sensitivity
        relationships = compute_structure_similarity(
            Path(args.input_dir), **kwargs
        )

    # Zero hits above the tool's reporting threshold is a legitimate
    # result (e.g. an all-vs-all search over unrelated sequences), not
    # an error — write_records() refuses empty lists, so write directly.
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps([r.model_dump() for r in relationships], indent=2),
        encoding="utf-8",
    )
    print(f"computed {len(relationships)} relationships -> {output}")


def _cmd_cluster(args: argparse.Namespace) -> None:
    input_dir = Path(args.input_dir)
    relationships = _load_json_list(
        SimilarityRelationship, Path(args.relationships)
    )
    item_ids = sorted(path.stem.upper() for path in input_dir.glob("*.cif"))

    clusters, prov = cluster_similar_items(
        item_ids, relationships, args.threshold
    )

    output = Path(args.output)
    write_records(clusters, output)
    write_json(prov, output.parent / "cluster_provenance.json")
    print(f"{len(clusters)} clusters at threshold={args.threshold} -> {output}")


def _cmd_partition(args: argparse.Namespace) -> None:
    clusters = _load_json_list(SimilarityCluster, Path(args.clusters))
    splits, prov = partition_dataset(
        clusters,
        args.pct_train,
        args.pct_val,
        args.pct_test,
        keep_similar_items=not args.no_keep_similar_items,
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(splits, indent=2), encoding="utf-8")
    write_json(prov, output.parent / "partition_provenance.json")
    print(f"split sizes: {prov.split_sizes} -> {output}")


def _cmd_annotate(args: argparse.Namespace) -> None:
    output_dir = Path(args.output_dir)
    structures = _read_structures_dir(Path(args.input_dir))

    layers_by_entry: dict[str, list[AnnotationLayer]] = {
        entry_id: [] for entry_id in structures
    }
    for layer_name in args.layers:
        annotator = _ENTRY_ANNOTATORS[layer_name]
        for entry_id, structure in structures.items():
            layers_by_entry[entry_id].append(annotator(structure))

    if args.pairwise:
        for left_id, right_id in combinations(sorted(structures), 2):
            layer = annotate_pairwise_sequence_identity(
                structures[left_id], structures[right_id]
            )
            layers_by_entry[left_id].append(layer)
            layers_by_entry[right_id].append(layer)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "annotations.json").write_text(
        json.dumps(
            {
                entry_id: [layer.model_dump() for layer in layers]
                for entry_id, layers in layers_by_entry.items()
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"annotated {len(structures)} entries -> {output_dir}")


def _cmd_manifest(args: argparse.Namespace) -> None:
    structures = _read_structures_dir(Path(args.structures_dir))

    canonicalisation_policy = (
        load_policy(args.canonicalisation_policy)
        if args.canonicalisation_policy
        else None
    )
    curation_policy = (
        _load_yaml_model(DatasetCurationPolicy, Path(args.curation_policy))
        if args.curation_policy
        else None
    )
    ingestion_prov = (
        _load_json_dict(IngestionProvenance, Path(args.ingestion_provenance))
        if args.ingestion_provenance
        else {}
    )
    canon_prov = (
        _load_json_dict(
            canonicalisationProvenance, Path(args.canonicalisation_provenance)
        )
        if args.canonicalisation_provenance
        else {}
    )
    annotations = (
        _load_annotations(Path(args.annotations)) if args.annotations else {}
    )
    exclusions: list[ExclusionRecord] = []
    for path in args.exclusions or []:
        exclusions.extend(_load_json_list(ExclusionRecord, Path(path)))

    dedup_prov = (
        _load_json_model(DeduplicationProvenance, Path(args.dedup_provenance))
        if args.dedup_provenance
        else None
    )
    cluster_prov = (
        _load_json_model(ClusteringProvenance, Path(args.cluster_provenance))
        if args.cluster_provenance
        else None
    )
    partition_prov = (
        _load_json_model(PartitionProvenance, Path(args.partition_provenance))
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
    from pandora.provenance.reproduce import reproduce_dataset

    manifest = DatasetManifest.model_validate_json(
        Path(args.manifest).read_text()
    )
    output_dir = Path(args.output_dir)
    fetch_options = FetchOptions(use_cache=not args.no_cache)

    structures, new_manifest = reproduce_dataset(
        manifest, output_dir, fetch_options=fetch_options
    )

    _write_structures_dir(structures, output_dir / "structures")
    write_json(new_manifest, output_dir / "reproduced_manifest.json")
    print(f"reproduced {len(structures)} structures -> {output_dir}")


def _cmd_export(args: argparse.Namespace) -> None:
    structure, _, status = mmcif_to_structure(args.input)
    if structure is None:
        raise SystemExit(f"failed to parse {args.input}: {status}")

    output = Path(args.output)
    if output.suffix == ".json":
        write_json(structure, output)
    else:
        structure_to_mmcif(structure, output)
    print(f"exported -> {output}")


# Argument parser -----------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pandora", description="Pandora dataset-curation pipeline CLI."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p = subparsers.add_parser("fetch", help="Fetch raw mmCIF files.")
    p.add_argument("entry_ids", nargs="+")
    p.add_argument("--provider", choices=["pdbe", "pdb"], default="pdbe")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--no-cache", action="store_true")
    p.add_argument("--allow-partial", action="store_true")
    p.set_defaults(func=_cmd_fetch)

    p = subparsers.add_parser(
        "canonicalise", help="Parse + canonicalise a directory of mmCIF files."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument(
        "--policy", required=True, help="canonicalisationPolicy YAML"
    )
    p.add_argument("--output-dir", required=True)
    p.set_defaults(func=_cmd_canonicalise)

    p = subparsers.add_parser(
        "curate", help="Apply quality/organism/content rules to a directory."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument("--policy", required=True, help="DatasetCurationPolicy YAML")
    p.add_argument("--output-dir", required=True)
    p.set_defaults(func=_cmd_curate)

    p = subparsers.add_parser(
        "dedup", help="Remove structures sharing the same entry_id."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument(
        "--disable", action="store_true", help="record but keep dupes"
    )
    p.set_defaults(func=_cmd_dedup)

    p = subparsers.add_parser(
        "similarity", help="All-vs-all sequence or structure similarity."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument("--engine", choices=["mmseqs2", "foldseek"], required=True)
    p.add_argument("--mmseqs-bin", default="mmseqs")
    p.add_argument("--foldseek-bin", default="foldseek")
    p.add_argument("--sensitivity", type=float, default=None)
    p.add_argument("--alignment-type", type=int, default=2)
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_similarity)

    p = subparsers.add_parser(
        "cluster",
        help="Connected-component clustering of a relationship network.",
    )
    p.add_argument(
        "--input-dir", required=True, help="dir of the clustered *.cif"
    )
    p.add_argument("--relationships", required=True)
    p.add_argument("--threshold", type=float, default=0.9)
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_cluster)

    p = subparsers.add_parser(
        "partition", help="Leakage-safe train/val/test split of clusters."
    )
    p.add_argument("--clusters", required=True)
    p.add_argument("--pct-train", type=float, default=0.6)
    p.add_argument("--pct-val", type=float, default=0.2)
    p.add_argument("--pct-test", type=float, default=0.2)
    p.add_argument("--no-keep-similar-items", action="store_true")
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_partition)

    p = subparsers.add_parser(
        "annotate", help="Entry-level and/or pairwise annotation layers."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument(
        "--layers",
        nargs="+",
        choices=list(_ENTRY_ANNOTATORS),
        default=list(_ENTRY_ANNOTATORS),
    )
    p.add_argument("--pairwise", action="store_true")
    p.add_argument("--output-dir", required=True)
    p.set_defaults(func=_cmd_annotate)

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

    p = subparsers.add_parser(
        "export", help="Convert one mmCIF file to mmCIF or JSON."
    )
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=_cmd_export)

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    args.func(args)


def app() -> None:
    """Pandora CLI entry point."""

    main()
