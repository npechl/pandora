"""CLI `annotate`: entry-level and pairwise annotation layers."""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

from pandora.annotations import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    annotate_pairwise_sequence_identity,
    annotate_structure_counts,
)
from pandora.cli._io import (
    read_structures_dir,
)
from pandora.schemas.annotation import AnnotationLayer

_ENTRY_ANNOTATORS = {
    "structure_counts": annotate_structure_counts,
    "ligand_contacts": annotate_ligand_contacts,
    "chain_interfaces": annotate_chain_interfaces,
}


def _cmd_annotate(args: argparse.Namespace) -> None:
    """Handle the `annotate` subcommand: compute entry-level and/or
    pairwise annotation layers."""

    output_dir = Path(args.output_dir)
    structures = read_structures_dir(Path(args.input_dir))

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


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

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
