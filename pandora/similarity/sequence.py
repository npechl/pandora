from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from pandora.schemas.dataset import ChainRecord
from pandora.schemas.similarity import (
    SimilarityMethod,
    SimilarityRelationship,
)

_OUTPUT_COLUMNS = "query,target,fident,alnlen,qcov,tcov"
_FASTA_GLOBS = ("*.fasta", "*.fa", "*.fna", "*.faa")


def _mmseqs_version(mmseqs_bin: str) -> str | None:
    result = subprocess.run(
        [mmseqs_bin, "version"], capture_output=True, text=True, check=False
    )
    return result.stdout.strip() or None


def _write_fasta(sequences: dict[str, str], path: Path) -> None:
    with path.open("w") as handle:
        for seq_id, sequence in sequences.items():
            handle.write(f">{seq_id}\n{sequence}\n")


def _concat_fasta_dir(directory: Path, path: Path) -> None:
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
    *,
    mmseqs_bin: str = "mmseqs",
    sensitivity: float = 5.7,
    tmp_dir: str | Path | None = None,
    mmseqs_options: list[str] = [],
) -> list[SimilarityRelationship]:
    """All-vs-all sequence similarity via MMseqs2 `easy-search`.

    Args:
        sequences: Mapping of item id -> sequence, a list of `ChainRecord`
            (keyed as "{entry_id}_{chain_id}", records with no sequence are
            skipped), or a path to a directory of existing FASTA files to
            run similarity over directly.
        mmseqs_bin: Path or name of the MMseqs2 binary.
        sensitivity: MMseqs2 `-s` sensitivity value.
        tmp_dir: Working directory for FASTA/result files. None = system temp.
        mmseqs_options: addition options passed to MMseqs2 implementation.

    Returns:
        One `SimilarityRelationship` per unordered pair of items with a hit,
        `source_id < target_id`. Unthresholded — callers filter by score
        when building a similarity network.
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

    with tempfile.TemporaryDirectory(dir=tmp_dir) as work_dir:
        work = Path(work_dir)
        fasta_path = work / "sequences.fasta"
        result_path = work / "result.m8"
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
                str(result_path),
                str(work / "tmp"),
                "-s",
                str(sensitivity),
                "--format-output",
                _OUTPUT_COLUMNS,
                "-v 1",
            ]
            + mmseqs_options,
            check=True,
            capture_output=True,
            text=True,
        )

        version = _mmseqs_version(mmseqs_bin)
        best_hits: dict[tuple[str, str], tuple[float, float]] = {}
        for line in result_path.read_text().splitlines():
            query, target, fident, _alnlen, qcov, tcov = line.split("\t")
            if query == target:
                continue

            source_id, target_id = sorted((query, target))
            score = float(fident)
            coverage = min(float(qcov), float(tcov))

            pair = (source_id, target_id)
            current = best_hits.get(pair)
            if current is None or score > current[0]:
                best_hits[pair] = (score, coverage)

    return [
        SimilarityRelationship(
            source_id=source_id,
            target_id=target_id,
            similarity_type="sequence_similarity",
            score=score,
            coverage=coverage,
            identity=score,
            method=SimilarityMethod(
                engine="MMseqs2",
                version=version,
                parameters={"sensitivity": sensitivity, "binary": mmseqs_bin},
            ),
        )
        for (source_id, target_id), (score, coverage) in sorted(
            best_hits.items()
        )
    ]
