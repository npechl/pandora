"""File helpers shared by the CLI subcommands."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml
from pydantic import BaseModel

from pandora.export import structure_to_mmcif
from pandora.parsing import mmcif_to_structure
from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.structure import Structure


def read_structures_dir(input_dir: Path) -> dict[str, Structure]:
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


def write_structures_dir(
    structures: dict[str, Structure], output_dir: Path
) -> None:
    """Write each structure to `<entry_id>.cif` (lowercase) in output_dir."""

    output_dir.mkdir(parents=True, exist_ok=True)
    for entry_id, structure in structures.items():
        structure_to_mmcif(structure, output_dir / f"{entry_id.lower()}.cif")


def write_json_dict(data: dict[str, BaseModel], path: Path) -> None:
    """Write a dict of pydantic models as one pretty-printed JSON object,
    keyed by dict key."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({k: v.model_dump() for k, v in data.items()}, indent=2),
        encoding="utf-8",
    )


def load_json_model(model_cls: type[BaseModel], path: Path) -> BaseModel:
    """Load and validate a single model instance from a JSON file."""

    return model_cls.model_validate(json.loads(path.read_text()))


def load_json_list(model_cls: type[BaseModel], path: Path) -> list:
    """Load and validate a JSON array as a list of model instances."""

    data = json.loads(path.read_text())
    return [model_cls.model_validate(item) for item in data]


def load_json_dict(model_cls: type[BaseModel], path: Path) -> dict:
    """Load and validate a JSON object as a dict of model instances."""

    data = json.loads(path.read_text())
    return {k: model_cls.model_validate(v) for k, v in data.items()}


def load_annotations(path: Path) -> dict[str, list[AnnotationLayer]]:
    """Load a JSON object of entry_id -> list of AnnotationLayer."""

    data = json.loads(path.read_text())
    return {
        entry_id: [AnnotationLayer.model_validate(item) for item in layers]
        for entry_id, layers in data.items()
    }


def load_yaml_model(model_cls: type[BaseModel], path: Path) -> BaseModel:
    """Load and validate a single model instance from a YAML file."""

    with path.open() as stream:
        data = yaml.safe_load(stream)
    return model_cls.model_validate(data)


def load_interface_residues(path: Path) -> dict[str, set[int]]:
    """Load a JSON object of item_id -> [residue position, ...] as
    `compute_structure_similarity(interface_residues=...)` expects."""

    data = json.loads(path.read_text())
    return {item_id: set(positions) for item_id, positions in data.items()}


def load_pairs(path: Path) -> list[tuple[str, str]]:
    """Load a JSON array of [item_id_1, item_id_2] pairs."""

    data = json.loads(path.read_text())
    return [(item_1, item_2) for item_1, item_2 in data]
