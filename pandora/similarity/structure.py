from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from pandora.schemas.annotation import AnnotationLayer
from pandora.schemas.similarity import SimilarityMethod, SimilarityRelationship
from pandora.schemas.structure import Structure

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


def residue_positions(structure: Structure, chain_id: str) -> dict[int, int]:
    """`label_seq_id -> 1-indexed position` for one chain, in the order
    `export_chain_mmcif()`/`structure_to_mmcif()` write that chain's
    atoms — i.e. the same numbering Foldseek assigns that chain's
    residues when it parses the resulting file.

    Args:
        structure: The structure to read chain_id's atoms from.
        chain_id: The `label_asym_id` of the chain.

    Returns:
        One entry per residue that has atoms for this chain in
        `structure.atoms` (first-appearance order). A residue missing
        from `structure.atoms` (no resolved density) is absent here too
        — it's also absent from the exported file, so there's nothing
        for Foldseek to number.
    """

    positions: dict[int, int] = {}
    for atom in structure.atoms:
        if atom.label_asym_id != chain_id or atom.label_seq_id is None:
            continue
        if atom.label_seq_id not in positions:
            positions[atom.label_seq_id] = len(positions) + 1
    return positions


def chain_item_id(entry_id: str, chain_id: str) -> str:
    """Item id for one chain's exported file — pair with
    `export_chain_mmcif(structure, chain_id, ...)` so the id used for
    `compute_structure_similarity()` matches the file it names."""

    return f"{entry_id}_{chain_id}"


def interface_residues_from_annotation(
    structures: dict[str, Structure],
    interface_layers: dict[str, AnnotationLayer],
) -> dict[str, set[int]]:
    """Bridge `annotate_chain_interfaces()` output into the
    `interface_residues` mapping `compute_structure_similarity()` expects,
    for one exported file per chain (`export_chain_mmcif()`).

    Converts each interface's `label_asym_id:label_seq_id` residue ids to
    Foldseek-aligned positions via `residue_positions()`, so callers never
    have to reason about the numbering mismatch between mmCIF's
    `label_seq_id` and Foldseek's own per-file residue order themselves.

    Args:
        structures: `entry_id -> Structure`, the same structures the
            layers in `interface_layers` were computed from.
        interface_layers: `entry_id -> its "chain_interfaces"
            AnnotationLayer` (from `annotate_chain_interfaces()`).

    Returns:
        `{chain_item_id(entry_id, chain_id): {positions...}}` — one
        entry per chain that appears in at least one interface. Pass
        this straight to `compute_structure_similarity(structures=...,
        interface_residues=...)` where `structures` was built with
        `export_chain_mmcif()` using the same ids.
    """

    result: dict[str, set[int]] = {}
    for entry_id, layer in interface_layers.items():
        structure = structures[entry_id]
        positions_by_chain: dict[str, dict[int, int]] = {}
        for interface in layer.data.get("interfaces", []):
            sides = (
                (
                    interface["chain_id_1"],
                    interface["interface_residues_chain_1"],
                ),
                (
                    interface["chain_id_2"],
                    interface["interface_residues_chain_2"],
                ),
            )
            for chain_id, residue_ids in sides:
                if chain_id not in positions_by_chain:
                    positions_by_chain[chain_id] = residue_positions(
                        structure, chain_id
                    )
                positions = positions_by_chain[chain_id]
                item_id = chain_item_id(entry_id, chain_id)
                for residue_id in residue_ids:
                    label_seq_id = int(residue_id.rsplit(":", 1)[1])
                    if label_seq_id in positions:
                        result.setdefault(item_id, set()).add(
                            positions[label_seq_id]
                        )
    return result


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
