import stat

from pandora.schemas.dataset import ChainRecord
from pandora.similarity.sequence import compute_sequence_similarity

# Stands in for `mmseqs`: records the query FASTA and writes a canned
# easy-search result (self-hit, plus a<->b hit in both directions).
FAKE_MMSEQS = """#!/bin/sh
if [ "$1" = version ]; then echo 15.6f452; exit 0; fi
cp "$2" "$(dirname "$0")/query.fasta"
printf 'a_A\\ta_A\\t1.0\\t10\\t1.0\\t1.0\\n' > "$4"
printf 'a_A\\tb_A\\t0.5\\t10\\t0.9\\t0.8\\n' >> "$4"
printf 'b_A\\ta_A\\t0.7\\t10\\t0.6\\t0.9\\n' >> "$4"
"""


def test_best_hit_per_pair_from_chain_records(tmp_path):
    mmseqs = tmp_path / "mmseqs"
    mmseqs.write_text(FAKE_MMSEQS)
    mmseqs.chmod(mmseqs.stat().st_mode | stat.S_IEXEC)
    records = [
        ChainRecord(
            entry_id=entry_id,
            chain_id="A",
            entity_id="1",
            sequence=sequence,
            residue_count=3,
            atom_count=3,
        )
        for entry_id, sequence in [("a", "MKV"), ("b", "MKL"), ("c", None)]
    ]

    (hit,) = compute_sequence_similarity(records, mmseqs_bin=str(mmseqs))

    assert (tmp_path / "query.fasta").read_text() == ">a_A\nMKV\n>b_A\nMKL\n"
    assert (hit.source_id, hit.target_id) == ("a_A", "b_A")
    assert hit.score == 0.7
    assert hit.coverage == 0.6
    assert hit.method.engine == "MMseqs2"
    assert hit.method.version == "15.6f452"
