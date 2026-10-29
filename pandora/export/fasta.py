from __future__ import annotations

import re
from pathlib import Path

_WHITESPACE = re.compile(r"\s+")


def write_fasta(sequences: dict[str, str], path: str | Path) -> Path:
    """Write sequences to a FASTA file, one record per sequence.

    Whitespace inside a sequence is removed, so sequences copied from
    mmCIF (which wraps long one-letter codes) are written on one line.

    Args:
        sequences: Record id to sequence, written in dict order.
        path: Destination file path; parent directories are created
            if missing.

    Returns:
        The `Path` the FASTA file was written to.
    """

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        for seq_id, sequence in sequences.items():
            handle.write(f">{seq_id}\n{_WHITESPACE.sub('', sequence)}\n")
    return out_path
