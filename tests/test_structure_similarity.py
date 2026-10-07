from __future__ import annotations

import stat

from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.structure import AtomSiteRecord, EntryRecord, Structure
from pandora.similarity.hits import HIT_COLUMNS, _interface_coverage
from pandora.similarity.structure import (
    compute_structure_similarity,
    chain_item_id,
    interface_residues_from_annotation,
    residue_positions,
)


def test_interface_coverage_full_and_partial_overlap() -> None:
    assert _interface_coverage({10, 11, 12}, start=1, end=20) == 1.0
    assert _interface_coverage({10, 11, 12, 13}, start=1, end=11) == 0.5


def test_interface_coverage_no_overlap() -> None:
    assert _interface_coverage({1, 2, 3}, start=10, end=20) == 0.0


def test_interface_coverage_empty_residues() -> None:
    assert _interface_coverage(set(), start=1, end=20) == 0.0


def _atom(chain_id: str, seq_id: int) -> AtomSiteRecord:
    return AtomSiteRecord(
        group_PDB="ATOM",
        id=1,
        type_symbol="C",
        label_atom_id="CA",
        label_comp_id="ALA",
        label_asym_id=chain_id,
        label_entity_id="1",
        label_seq_id=seq_id,
        Cartn_x=0.0,
        Cartn_y=0.0,
        Cartn_z=0.0,
        occupancy=1.0,
        B_iso_or_equiv=0.0,
    )


def _structure_with_gap() -> Structure:
    # chain A: residues 1,2,3 resolved, 4,5 missing (no density), 6,7
    # resolved -> label_seq_id 6/7 land at file positions 4/5, not 6/7.
    atoms = [_atom("A", i) for i in (1, 2, 3, 6, 7)]
    atoms += [_atom("B", i) for i in (1, 2, 3)]
    return Structure(
        entry_id="ENTRY1", entry=EntryRecord(id="ENTRY1"), atoms=atoms
    )


def test_residue_positions_reindexes_around_gaps() -> None:
    structure = _structure_with_gap()
    assert residue_positions(structure, "A") == {1: 1, 2: 2, 3: 3, 6: 4, 7: 5}
    assert residue_positions(structure, "B") == {1: 1, 2: 2, 3: 3}


def test_chain_item_id() -> None:
    assert chain_item_id("ENTRY1", "A") == "ENTRY1_A"


def test_interface_residues_from_annotation_maps_to_file_positions() -> None:
    structure = _structure_with_gap()
    layer = AnnotationLayer(
        layer_name="Chain-chain interfaces",
        layer_type="chain_interfaces",
        scope="interface",
        method="pandora.basic.distance_cutoff_contacts.v1",
        target_ids=["ENTRY1"],
        data={
            "distance_cutoff": 4.0,
            "interfaces": [
                {
                    "chain_id_1": "A",
                    "chain_id_2": "B",
                    "interface_residues_chain_1": ["A:6", "A:7"],
                    "interface_residues_chain_2": ["B:2"],
                    "contact_count": 3,
                }
            ],
        },
    )

    result = interface_residues_from_annotation(
        {"ENTRY1": structure}, {"ENTRY1": layer}
    )

    # A:6/A:7 are label_seq_id's, but land at file positions 4/5 (the gap
    # shift) -- this is exactly the mismatch compute_structure_similarity's
    # interface_residues needs positions for, not label_seq_id's.
    assert result == {
        "ENTRY1_A": {4, 5},
        "ENTRY1_B": {2},
    }


FAKE_FOLDSEEK = """#!/bin/sh
if [ "$1" = version ]; then echo 10.941cd33; exit 0; fi
echo "$@" > "$(dirname "$0")/args.txt"
printf 'A\\tB\\t0.5\\t100\\t0.9\\t0.9\\t0.6\\t0.6\\t0.6\\t1\\t100\\t1\\t100\\n' > "$4"
"""


def _fake_foldseek(directory):
    foldseek = directory / "foldseek"
    foldseek.write_text(FAKE_FOLDSEEK)
    foldseek.chmod(foldseek.stat().st_mode | stat.S_IEXEC)
    return foldseek


def test_structure_search_passes_caps_and_keeps_hits(tmp_path):
    foldseek = _fake_foldseek(tmp_path)
    struct_dir = tmp_path / "structures"
    struct_dir.mkdir()
    hits_path = tmp_path / "hits.tsv"

    search = compute_structure_similarity(
        struct_dir,
        hits_path,
        foldseek_bin=str(foldseek),
        alignment_type=1,
        max_seqs=2000,
        exhaustive_search=True,
    )

    args = (tmp_path / "args.txt").read_text().split()
    assert args[args.index("--max-seqs") + 1] == "2000"
    assert args[args.index("--exhaustive-search") + 1] == "1"
    assert args[args.index("--alignment-type") + 1] == "1"
    assert args[args.index("--format-output") + 1] == ",".join(
        HIT_COLUMNS["Foldseek"]
    )
    assert hits_path.read_text().count("\n") == 1
    assert search.engine == "Foldseek"
    assert search.version == "10.941cd33"
    assert search.columns == HIT_COLUMNS["Foldseek"]
    assert search.parameters == {
        "foldseek_bin": str(foldseek),
        "sensitivity": 9.5,
        "alignment_type": 1,
        "max_seqs": 2000,
        "exhaustive_search": True,
        "foldseek_options": [],
    }


def test_structure_search_omits_exhaustive_flag_by_default(tmp_path):
    foldseek = _fake_foldseek(tmp_path)
    struct_dir = tmp_path / "structures"
    struct_dir.mkdir()

    compute_structure_similarity(
        struct_dir, tmp_path / "hits.tsv", foldseek_bin=str(foldseek)
    )

    args = (tmp_path / "args.txt").read_text().split()
    assert "--exhaustive-search" not in args
    assert args[args.index("--max-seqs") + 1] == "1000"


def test_structure_search_records_absolute_hits_path(tmp_path, monkeypatch):
    foldseek = _fake_foldseek(tmp_path)
    (tmp_path / "structures").mkdir()
    monkeypatch.chdir(tmp_path)

    search = compute_structure_similarity(
        "structures", "hits.tsv", foldseek_bin=str(foldseek)
    )

    assert search.hits_path == str((tmp_path / "hits.tsv").resolve())


if __name__ == "__main__":
    test_interface_coverage_full_and_partial_overlap()
    test_interface_coverage_no_overlap()
    test_interface_coverage_empty_residues()
    test_residue_positions_reindexes_around_gaps()
    test_chain_item_id()
    test_interface_residues_from_annotation_maps_to_file_positions()
    print("ok")
