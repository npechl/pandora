from __future__ import annotations

import tempfile
from pathlib import Path

from pandora.schemas.dataset import ChainRecord
from pandora.similarity.sequence import run_sequence_similarity_engine

SEQUENCES = {
    "1abc_A": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR",
    "1abc_B": "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR",
    "2xyz_A": "GSSGSSGPSSGGSSGSSGPSSGGSSGSSGPSSGGSSGSSGPSSGGSSGSSGPSSGGSSGSSGPSSGGSSGSSGPSSGG",
}


def test_identical_sequences_score_one() -> None:
    relationships = run_sequence_similarity_engine(SEQUENCES)
    by_pair = {(r.source_id, r.target_id): r for r in relationships}

    identical = by_pair[("1abc_A", "1abc_B")]
    assert identical.score == 1.0
    assert identical.identity == 1.0
    assert identical.method.engine == "MMseqs2"


def test_unrelated_sequence_pairs_are_ordered_and_absent_or_low() -> None:
    relationships = run_sequence_similarity_engine(SEQUENCES)
    for r in relationships:
        assert r.source_id < r.target_id
        if {r.source_id, r.target_id} != {"1abc_A", "1abc_B"}:
            assert r.score < 1.0


def test_directory_of_fasta_files_input() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        for seq_id, sequence in SEQUENCES.items():
            (directory / f"{seq_id}.fasta").write_text(
                f">{seq_id}\n{sequence}\n"
            )

        relationships = run_sequence_similarity_engine(directory)
        by_pair = {(r.source_id, r.target_id): r for r in relationships}
        assert by_pair[("1abc_A", "1abc_B")].score == 1.0


def test_unknown_binary_raises() -> None:
    try:
        run_sequence_similarity_engine(
            SEQUENCES, mmseqs_bin="not-a-real-binary"
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("expected RuntimeError for missing binary")


def test_chain_record_list_input() -> None:
    records = [
        ChainRecord(
            entry_id=seq_id.rsplit("_", 1)[0],
            chain_id=seq_id.rsplit("_", 1)[1],
            entity_id="1",
            sequence=sequence,
            residue_count=len(sequence),
            atom_count=len(sequence),
        )
        for seq_id, sequence in SEQUENCES.items()
    ] + [
        ChainRecord(
            entry_id="3noseq",
            chain_id="A",
            entity_id="1",
            sequence=None,
            residue_count=0,
            atom_count=0,
        )
    ]

    relationships = run_sequence_similarity_engine(records)
    by_pair = {(r.source_id, r.target_id): r for r in relationships}
    assert by_pair[("1abc_A", "1abc_B")].score == 1.0
    assert all("3noseq" not in pair for pair in by_pair)


if __name__ == "__main__":
    test_identical_sequences_score_one()
    test_unrelated_sequence_pairs_are_ordered_and_absent_or_low()
    test_directory_of_fasta_files_input()
    test_unknown_binary_raises()
    test_chain_record_list_input()
    print("ok")
