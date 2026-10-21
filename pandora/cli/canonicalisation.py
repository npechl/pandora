"""CLI `canonicalise`: parse and canonicalise a directory of mmCIF files."""

from __future__ import annotations

import argparse
from pathlib import Path

from pandora.canonicalisation import canonicalise_structure
from pandora.cli._io import (
    read_structures_dir,
    write_json_dict,
    write_structures_dir,
)
from pandora.ingestion.policy import load_policy
from pandora.schemas.canonicalisation import (
    CanonicalMappings,
    canonicalisationProvenance,
)
from pandora.schemas.structure import Structure


def _cmd_canonicalise(args: argparse.Namespace) -> None:
    """Handle the `canonicalise` subcommand: parse and canonicalise a
    directory of mmCIF files."""

    output_dir = Path(args.output_dir)
    policy = load_policy(args.policy)
    structures = read_structures_dir(Path(args.input_dir))

    canonical: dict[str, Structure] = {}
    provenance: dict[str, canonicalisationProvenance] = {}
    mappings: dict[str, CanonicalMappings] = {}
    for entry_id, structure in structures.items():
        canonical_structure, entry_mappings, prov = canonicalise_structure(
            structure, policy
        )
        canonical[entry_id] = canonical_structure
        provenance[entry_id] = prov
        mappings[entry_id] = entry_mappings

    write_structures_dir(canonical, output_dir)
    write_json_dict(provenance, output_dir / "canonicalisation_provenance.json")
    write_json_dict(mappings, output_dir / "canonicalisation_mappings.json")
    print(f"canonicalised {len(canonical)} structures -> {output_dir}")


def register(subparsers: argparse._SubParsersAction) -> None:
    """Add this module's subcommands to the `pandora` parser."""

    p = subparsers.add_parser(
        "canonicalise", help="Parse + canonicalise a directory of mmCIF files."
    )
    p.add_argument("--input-dir", required=True)
    p.add_argument(
        "--policy", required=True, help="canonicalisationPolicy YAML"
    )
    p.add_argument("--output-dir", required=True)
    p.set_defaults(func=_cmd_canonicalise)
