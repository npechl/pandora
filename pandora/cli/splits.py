"""CLI `splits assign`: leakage-safe train/val/test split of clusters."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pandora._util import now_iso
from pandora.cli._io import (
    SUMMARY,
    add_out,
    check_out,
    find_table,
    input_summary,
    load_json_list,
    prepare_output_dir,
    read_summary,
    resolve_file,
    resolve_table,
    write_summary,
)
from pandora.export import read_records
from pandora.provenance import build_step_summary
from pandora.schemas.annotation import LayerRow
from pandora.schemas.similarity import SimilarityCluster
from pandora.similarity import partition_dataset
from pandora.similarity.layer import item_labels


def _partition_labels(
    command: str, source: str | None
) -> tuple[dict[str, str] | None, Path | None]:
    """The item labels and table of a labels layer, or (None, None);
    exits on a bad layer."""

    if source is None:
        return None, None
    path = Path(source)
    if path.is_dir():
        table = find_table(path, "layer")  # an empty layer has no table
        folder = path
    else:
        table = resolve_table(command, path, "layer")
        folder = path.parent
    summary = read_summary(folder)
    if summary is None:
        raise SystemExit(
            f"{command}: {source} needs a labels layer folder "
            "(with summary.json), to tell chain from interface labels"
        )
    template = (summary.get("policy") or {}).get("template")
    try:
        rows = read_records(table, LayerRow) if table else []
        return item_labels(rows, template), table
    except ValueError as exc:
        raise SystemExit(f"{command}: {source}: {exc}") from exc


def _cmd_assign(args: argparse.Namespace) -> None:
    """Handle `splits assign`: leakage-safe train/val/test split of
    clusters."""

    command = "splits assign"
    started_at = now_iso()
    source = args.labels or args.stratify
    labels, table = _partition_labels(command, source)
    clusters_file = resolve_file(command, Path(args.clusters), "clusters.json")
    out = Path(args.out)
    check_out(command, clusters_file, out)
    clusters = load_json_list(SimilarityCluster, clusters_file)
    try:
        splits, prov = partition_dataset(
            clusters,
            args.pct_train,
            args.pct_val,
            args.pct_test,
            labels=labels,
            stratify=args.stratify is not None,
        )
    except ValueError as exc:
        raise SystemExit(f"{command}: {exc}") from exc
    if table is not None:
        prov = prov.model_copy(update={"labels_layer": str(table.resolve())})

    prepare_output_dir(command, out, ["splits.json", SUMMARY], args.overwrite)
    (out / "splits.json").write_text(
        json.dumps(splits, indent=2), encoding="utf-8"
    )
    write_summary(
        out,
        build_step_summary(
            "partition",
            prov,
            input=clusters_file.resolve(),
            output=out,
            started_at=started_at,
            counts={f"items:{k}": n for k, n in prov.split_sizes.items()},
            previous_summary=input_summary(clusters_file.parent),
        ),
    )
    print(f"split sizes: {prov.split_sizes} -> {out}")
    if labels is not None:
        items = [i for c in clusters for i in c.components]
        unlabelled = sum(i not in labels for i in items)
        print(f"unlabelled items: {unlabelled} of {len(items)}")
        for split, missing in prov.missing_labels.items():
            print(
                f"{command}: warning: {split} has no item labelled "
                f"{', '.join(missing)}",
                file=sys.stderr,
            )


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add the `splits` group to the `pandora` parser."""

    group = subparsers.add_parser(
        "splits", help="Assign clusters to train/val/test splits."
    )
    commands = group.add_subparsers(dest="splits_command", required=True)
    p = commands.add_parser(
        "assign",
        help="Leakage-safe train/val/test split: each cluster goes whole "
        "to one split.",
    )
    p.add_argument("clusters", help="clusters folder or clusters.json")
    add_out(p, tables=False)
    p.add_argument(
        "--pct-train",
        type=float,
        default=0.6,
        help="fraction of items in train (default 0.6)",
    )
    p.add_argument(
        "--pct-val",
        type=float,
        default=0.2,
        help="fraction of items in val (default 0.2)",
    )
    p.add_argument(
        "--pct-test",
        type=float,
        default=0.2,
        help="fraction of items in test (default 0.2)",
    )
    labels = p.add_mutually_exclusive_group()
    labels.add_argument(
        "--labels",
        help="chain_labels or interface_labels layer folder: place "
        "clusters so every label reaches every split",
    )
    labels.add_argument(
        "--stratify",
        help="chain_labels or interface_labels layer folder: split each "
        "label's clusters separately, keeping the label mix in every split",
    )
    p.set_defaults(func=_cmd_assign)
