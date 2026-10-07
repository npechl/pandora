from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from pandora._util import now_iso
from pandora.schemas.dataset import ChainRecord
from pandora.schemas.similarity import SimilaritySearch
from pandora.similarity.hits import HIT_COLUMNS

_FASTA_GLOBS = ("*.fasta", "*.fa", "*.fna", "*.faa")


def _mmseqs_version(mmseqs_bin: str) -> str | None:
    """Installed MMseqs2 version string, or None if it can't be determined."""

    result = subprocess.run(
        [mmseqs_bin, "version"], capture_output=True, text=True, check=False
    )
    return result.stdout.strip() or None


def _write_fasta(sequences: dict[str, str], path: Path) -> None:
    """Write sequences as a FASTA file at path."""

    with path.open("w") as handle:
        for seq_id, sequence in sequences.items():
            handle.write(f">{seq_id}\n{sequence}\n")


def _concat_fasta_dir(directory: Path, path: Path) -> None:
    """Concatenate every FASTA file in directory into a single file at path."""

    fasta_files = sorted(
        f for pattern in _FASTA_GLOBS for f in directory.glob(pattern)
    )
    if not fasta_files:
        raise ValueError(
            f"no FASTA files ({'/'.join(_FASTA_GLOBS)}) found in {directory}"
        )

    with path.open("w") as out:
        for fasta_file in fasta_files:
            text = fasta_file.read_text()
            out.write(text if text.endswith("\n") else text + "\n")


def compute_sequence_similarity(
    sequences: dict[str, str] | list[ChainRecord] | str | Path,
    hits_path: str | Path,
    *,
    mmseqs_bin: str = "mmseqs",
    sensitivity: float = 5.7,
    max_seqs: int = 300,
    tmp_dir: str | Path | None = None,
    mmseqs_options: list[str] | None = None,
) -> SimilaritySearch:
    """All-vs-all sequence similarity via MMseqs2 `easy-search`.

    Writes MMseqs2's hits to `hits_path` as a TSV with the columns in
    `HIT_COLUMNS["MMseqs2"]` and keeps it there; nothing is loaded into
    memory. Read it with `iter_edges()` / `cluster_similar_items()`.

    Args:
        sequences: Mapping of item id -> sequence, a list of `ChainRecord`
            (keyed as "{entry_id}_{chain_id}", records with no sequence are
            skipped), or a path to a directory of existing FASTA files to
            run similarity over directly.
        hits_path: Where to write the hit TSV. Parent directories are
            created.
        mmseqs_bin: Path or name of the MMseqs2 binary.
        sensitivity: MMseqs2 `-s` sensitivity value.
        max_seqs: MMseqs2 `--max-seqs`: hits kept per query after the
            prefilter. The default (MMseqs2's own) can miss similar
            pairs in large families; raise it for leakage control.
        tmp_dir: Working directory for FASTA/temporary files. None =
            system temp.
        mmseqs_options: Additional options passed to MMseqs2.

    Returns:
        A `SimilaritySearch` pointing at `hits_path`. Its `parameters`
        are this function's keyword arguments, so the search can be
        re-run with `compute_sequence_similarity(seqs, path,
        **search.parameters)`.

    Raises:
        RuntimeError: `mmseqs_bin` is not on PATH.
        ValueError: `sequences` is a directory with no FASTA files.
        subprocess.CalledProcessError: MMseqs2 failed.
    """

    if shutil.which(mmseqs_bin) is None:
        raise RuntimeError(
            f"mmseqs binary {mmseqs_bin!r} not found (required for "
            "sequence similarity)"
        )

    if isinstance(sequences, list):
        sequences = {
            f"{record.entry_id}_{record.chain_id}": record.sequence
            for record in sequences
            if record.sequence is not None
        }

    options = list(mmseqs_options or [])
    hits = Path(hits_path)
    hits.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=tmp_dir) as work_dir:
        work = Path(work_dir)
        fasta_path = work / "sequences.fasta"
        if isinstance(sequences, dict):
            _write_fasta(sequences, fasta_path)
        else:
            _concat_fasta_dir(Path(sequences), fasta_path)

        subprocess.run(
            [
                mmseqs_bin,
                "easy-search",
                str(fasta_path),
                str(fasta_path),
                str(hits),
                str(work / "tmp"),
                "-s",
                str(sensitivity),
                "--max-seqs",
                str(max_seqs),
                "--format-output",
                ",".join(HIT_COLUMNS["MMseqs2"]),
                "-v",
                "1",
            ]
            + options,
            check=True,
            capture_output=True,
            text=True,
        )

    return SimilaritySearch(
        engine="MMseqs2",
        version=_mmseqs_version(mmseqs_bin),
        hits_path=str(hits),
        columns=list(HIT_COLUMNS["MMseqs2"]),
        parameters={
            "mmseqs_bin": mmseqs_bin,
            "sensitivity": sensitivity,
            "max_seqs": max_seqs,
            "mmseqs_options": options,
        },
        searched_at=now_iso(),
    )
