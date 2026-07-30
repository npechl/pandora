"""Multi-chain, ligand-bound structure walkthrough: both annotation
layers, dataset record extraction, every export format, and a full
provenance bundle — run against 1a3n from datasets/dev/mmcif/.

examples/overview.py stays deliberately small (two monomers, no
export). This exercises the stages that only do something interesting
on a richer structure: 1a3n (deoxyhemoglobin — 2 alpha + 2 beta
chains, 4 heme groups) for chain-chain interfaces and ligand pockets,
plus every pandora.export format and a pandora.provenance bundle tying
the run together.

Ingestion is skipped, as in overview.py: the file is already on disk.
"""

from pathlib import Path

from pandora.annotations import (
    annotate_chain_interfaces,
    annotate_ligand_contacts,
    annotate_structure_counts,
)
from pandora.canonicalisation import canonicalise_structure
from pandora.datasets import extract_chain_records, extract_interface_records
from pandora.export import structure_to_mmcif, write_json, write_records
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.provenance import build_provenance_bundle
from pandora.schemas.canonicalisation import (
    AltlocRules,
    canonicalisationPolicy,
    EntityRules,
    LigandRules,
)

MMCIF_PATH = Path("./datasets/dev/mmcif/1a3n.cif")
OUTPUT_DIR = Path("./datasets/output/dev/")
DISTANCE_CUTOFF = 4.0

# Unlike overview.py's policy, this keeps chain_id/residue_numbering at
# their "preserve" default rather than remap/renumber. 1a3n carries SIFTS
# and validation categories in Structure.raw (_pdbx_sifts_xref_db,
# _pdbx_validate_*, ...) that cross-reference label_asym_id/label_seq_id;
# structure_to_mmcif() re-emits those verbatim (see its docstring — the
# round-trip is lossy but not renumbering-aware), so a policy that
# remaps/renumbers desyncs them from the rewritten _atom_site and the
# reparse below fails gemmi's cross-checks. keep_waters is the only
# ligand filter: heme is a non-polymer ligand, not water or an ion, so
# keep_nonpolymer_ligands (default True) leaves it in the canonical atoms
# for annotate_ligand_contacts to find.
policy = canonicalisationPolicy(
    policy_id="complex-preserve",
    policy_name="Preserve Identifiers",
    policy_version="1.0.0",
    altloc_rules=AltlocRules(strategy="select_best_occupancy"),
    entity_rules=EntityRules(strategy="merge_equivalent_entities"),
    ligand_rules=LigandRules(strategy="filter", keep_waters=False),
)

# 1. Parsing + canonicalisation -------------------------------------------
structure, diagnostics, status = mmcif_to_structure(str(MMCIF_PATH))
print(f"parsed 1a3n: status={status} atoms={len(structure.atoms)} chains={len(structure.asym_units)}")

canonical, mappings, canon_prov = canonicalise_structure(structure, policy)
print(f"canonicalised: transforms={canon_prov.transforms}")

# 2. Metadata ---------------------------------------------------------------
metadata = collect_metadata(canonical)
ligand_names = [ligand.comp_id for ligand in metadata.ligands]
print(f"metadata: title={metadata.entry.title!r} ligands={ligand_names}")

# 3. Annotations — structure counts, chain-chain interfaces, ligand pockets
counts = annotate_structure_counts(canonical)
print(
    f"structure_counts: {counts.data['asym_unit_count']} chains, "
    f"{counts.data['residue_count']} residues"
)

interfaces = annotate_chain_interfaces(
    canonical, distance_cutoff=DISTANCE_CUTOFF
)
for interface in interfaces.data["interfaces"]:
    print(
        f"  interface {interface['chain_id_1']}-{interface['chain_id_2']}: "
        f"{interface['contact_count']} contacting residues"
    )

contacts = annotate_ligand_contacts(canonical, distance_cutoff=DISTANCE_CUTOFF)
for ligand in contacts.data["ligands"]:
    print(
        f"  ligand {ligand['ligand_comp_id']} in chain "
        f"{ligand['ligand_asym_id']}: {ligand['contact_count']} contacting "
        "residues"
    )

# 4. Dataset curation — reshape the canonical structure into flat records --
chain_records = extract_chain_records(canonical)
interface_records = extract_interface_records(
    canonical, distance_cutoff=DISTANCE_CUTOFF
)
print(
    f"records: {len(chain_records)} chains, {len(interface_records)} interfaces"
)

# 5. Export — mmCIF round-trip, flat records, provenance JSON --------------
mmcif_out = structure_to_mmcif(canonical, OUTPUT_DIR / "1a3n.canonical.cif")
reparsed, _, reparse_status = mmcif_to_structure(str(mmcif_out))
print(
    f"round-tripped mmCIF: status={reparse_status} atoms={len(reparsed.atoms)}"
)

write_records(chain_records, OUTPUT_DIR / "chains.jsonl")
write_records(interface_records, OUTPUT_DIR / "interfaces.jsonl")
try:
    write_records(chain_records, OUTPUT_DIR / "chains.parquet")
except RuntimeError as exc:
    print(f"skipped parquet export: {exc}")

# 6. Provenance — every stage's provenance in one bundle --------------------
bundle = build_provenance_bundle(
    canonical,
    canonicalisation=canon_prov,
    metadata=metadata,
    annotations=[counts, interfaces, contacts],
)
write_json(bundle, OUTPUT_DIR / "provenance.json")
print(f"provenance checksum: {bundle.checksums.structure_checksum[:12]}...")

print(f"\nall output written to {OUTPUT_DIR}")
