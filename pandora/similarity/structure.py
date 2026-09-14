from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from pandora.schemas.similarity import SimilarityMethod, SimilarityRelationship

_OUTPUT_COLUMNS = (
    "query,target,fident,alnlen,qcov,tcov,alntmscore,qstart,qend,tstart,tend"
)


def _foldseek_version(foldseek_bin: str) -> str | None:
    """Installed Foldseek version string, or None if it can't be determined."""

    result = subprocess.run(
        [foldseek_bin, "version"], capture_output=True, text=True, check=False
    )
    return result.stdout.strip() or None


def _structure_suffix(path: Path) -> str:
    """path's file suffix, treating a trailing .gz as part of it
    (e.g. ".cif.gz")."""

    suffixes = path.suffixes
    if len(suffixes) > 1 and suffixes[-1] == ".gz":
        return "".join(suffixes[-2:])
    return path.suffix


def _link_structures(
    structures: dict[str, str | Path], directory: Path
) -> None:
    """Symlink (or copy, if symlinking fails) each structure file into
    directory, named by its item id."""

    for item_id, path in structures.items():
        src = Path(path)
        dest = directory / f"{item_id}{_structure_suffix(src)}"
        try:
            dest.symlink_to(src.resolve())
        except OSError:
            shutil.copy(src, dest)


def _interface_coverage(residues: set[int], start: int, end: int) -> float:
    """Fraction of residues (1-indexed positions) that fall within the
    inclusive [start, end] alignment range. 0.0 if residues is empty."""

    if not residues:
        return 0.0
    return sum(1 for r in residues if start <= r <= end) / len(residues)


def compute_structure_similarity(
    structures: dict[str, str | Path] | str | Path,
    *,
    interface_residues: dict[str, set[int]] | None = None,
    foldseek_bin: str = "foldseek",
    sensitivity: float = 9.5,
    alignment_type: int = 2,
    tmp_dir: str | Path | None = None,
    foldseek_options: list[str] | None = None,
) -> list[SimilarityRelationship]:
    """All-vs-all structural similarity via Foldseek `easy-search`.

    Args:
        structures: Mapping of item id -> structure file path (PDB/mmCIF,
            optionally gzipped), or a path to a directory of existing
            structure files to run similarity over directly (ids are then
            Foldseek's own file-derived names).
        interface_residues: Optional mapping of item id -> 1-indexed
            interface residue positions, as Foldseek numbers residues for
            that item's structure file (identity mapping from
            `label_seq_id` only holds for a single-chain, gap-free file —
            see `docs/usage/similarity.md`). When given, each resulting
            relationship's `interface_coverage` is the fraction of the
            (source or target) item's interface residues falling inside
            the alignment range, rather than Foldseek's whole-chain
            coverage. A pair gets no `interface_coverage` unless both its
            items are present in this mapping.
        foldseek_bin: Path or name of the Foldseek binary.
        sensitivity: Foldseek `-s` sensitivity value.
        alignment_type: Foldseek `--alignment-type` (0: 3Di alignment,
            1: TM-align, 2: 3Di+AA — Foldseek's own default).
        tmp_dir: Working directory for structure/result files. None = system
            temp.
        foldseek_options: additional options passed to the Foldseek binary.

    Returns:
        One `SimilarityRelationship` per unordered pair of items with a hit,
        `source_id < target_id`. Unthresholded — callers filter by score
        when building a similarity network. `score` is the best hit's
        TM-score (`alntmscore`), `identity` its fraction of identical
        aligned residues (`fident`), `coverage` the min of query/target
        whole-chain coverage, `interface_coverage` the min of query/target
        interface-restricted coverage (if `interface_residues` given).
    """

    interface_residues = interface_residues or {}

    if shutil.which(foldseek_bin) is None:
        raise RuntimeError(
            f"foldseek binary {foldseek_bin!r} not found (required for "
            "structure similarity)"
        )

    foldseek_options = foldseek_options or []

    with tempfile.TemporaryDirectory(dir=tmp_dir) as work_dir:
        work = Path(work_dir)
        result_path = work / "result.m8"

        if isinstance(structures, dict):
            struct_dir = work / "structures"
            struct_dir.mkdir()
            _link_structures(structures, struct_dir)
        else:
            struct_dir = Path(structures)

        subprocess.run(
            [
                foldseek_bin,
                "easy-search",
                str(struct_dir),
                str(struct_dir),
                str(result_path),
                str(work / "tmp"),
                "-s",
                str(sensitivity),
                "--alignment-type",
                str(alignment_type),
                "--format-output",
                _OUTPUT_COLUMNS,
                "-v",
                "1",
            ]
            + foldseek_options,
            check=True,
            capture_output=True,
            text=True,
        )

        version = _foldseek_version(foldseek_bin)
        best_hits: dict[
            tuple[str, str], tuple[float, float, float, float | None]
        ] = {}
        for line in result_path.read_text().splitlines():
            (
                query,
                target,
                fident,
                _alnlen,
                qcov,
                tcov,
                alntmscore,
                qstart,
                qend,
                tstart,
                tend,
            ) = line.split("\t")
            if query == target:
                continue

            source_id, target_id = sorted((query, target))
            score = float(alntmscore)
            identity = float(fident)
            coverage = min(float(qcov), float(tcov))

            interface_coverage = None
            if query in interface_residues and target in interface_residues:
                q_iface_cov = _interface_coverage(
                    interface_residues[query], int(qstart), int(qend)
                )
                t_iface_cov = _interface_coverage(
                    interface_residues[target], int(tstart), int(tend)
                )
                interface_coverage = min(q_iface_cov, t_iface_cov)

            pair = (source_id, target_id)
            current = best_hits.get(pair)
            if current is None or score > current[0]:
                best_hits[pair] = (
                    score,
                    identity,
                    coverage,
                    interface_coverage,
                )

    return [
        SimilarityRelationship(
            source_id=source_id,
            target_id=target_id,
            similarity_type="structure_similarity",
            score=score,
            coverage=coverage,
            interface_coverage=interface_coverage,
            identity=identity,
            method=SimilarityMethod(
                engine="Foldseek",
                version=version,
                parameters={
                    "sensitivity": sensitivity,
                    "alignment_type": alignment_type,
                    "foldseek_bin": foldseek_bin,
                },
            ),
        )
        for (source_id, target_id), (
            score,
            identity,
            coverage,
            interface_coverage,
        ) in sorted(best_hits.items())
    ]
