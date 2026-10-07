import stat

from pandora.schemas.dataset import ChainRecord
from pandora.similarity.hits import HIT_COLUMNS
from pandora.similarity.sequence import compute_sequence_similarity

# Stands in for `mmseqs`: records its arguments and the query FASTA, and
# writes a canned easy-search result (a self-hit plus an a->b hit).
FAKE_MMSEQS = """#!/bin/sh
if [ "$1" = version ]; then echo 15.6f452; exit 0; fi
echo "$@" > "$(dirname "$0")/args.txt"
cp "$2" "$(dirname "$0")/query.fasta"
printf 'a_A\\ta_A\\t1.0\\t10\\t1.0\\t1.0\\n' > "$4"
printf 'a_A\\tb_A\\t0.5\\t10\\t0.9\\t0.8\\n' >> "$4"
"""


def _fake_mmseqs(directory):
    mmseqs = directory / "mmseqs"
    mmseqs.write_text(FAKE_MMSEQS)
    mmseqs.chmod(mmseqs.stat().st_mode | stat.S_IEXEC)
    return mmseqs


def test_search_keeps_hit_file_and_records_parameters(tmp_path):
    mmseqs = _fake_mmseqs(tmp_path)
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
    hits_path = tmp_path / "out" / "hits.tsv"  # parent doesn't exist yet

    search = compute_sequence_similarity(
        records, hits_path, mmseqs_bin=str(mmseqs), max_seqs=500
    )

    assert (tmp_path / "query.fasta").read_text() == ">a_A\nMKV\n>b_A\nMKL\n"
    assert hits_path.read_text().count("\n") == 2
    assert search.engine == "MMseqs2"
    assert search.version == "15.6f452"
    assert search.hits_path == str(hits_path)
    assert search.columns == HIT_COLUMNS["MMseqs2"]
    assert search.origin == "computed"
    assert search.searched_at is not None
    assert search.parameters == {
        "mmseqs_bin": str(mmseqs),
        "sensitivity": 5.7,
        "max_seqs": 500,
        "mmseqs_options": [],
    }
    args = (tmp_path / "args.txt").read_text().split()
    assert args[args.index("--max-seqs") + 1] == "500"
    assert args[args.index("--format-output") + 1] == ",".join(
        HIT_COLUMNS["MMseqs2"]
    )
