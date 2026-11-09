from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections.abc import Iterable, Mapping
from pathlib import Path

from pandora._util import now_iso
from pandora.datasets.records import interface_residue_pairs
from pandora.schemas.annotation import ResidueContactRow
from pandora.schemas.similarity import SimilaritySearch
from pandora.schemas.structure import Structure
from pandora.similarity.hits import HIT_COLUMNS


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
    structures: Mapping[str, Structure],
    rows: Iterable[ResidueContactRow],
) -> dict[str, set[int]]:
    """Bridge a `residue_contacts` layer into the `interface_residues`
    mapping `iter_edges()` expects, for one exported file per chain
    (`export_chain_mmcif()`).

    Works on any contact layer, Pandora's own
    (`annotate_chain_interfaces()`) or imported. Only contacts between
    polymer residues on different chains count. Each residue's
    `label_seq_id` is converted to its Foldseek-aligned position via
    `residue_positions()`, so callers never have to reason about the
    numbering mismatch between mmCIF's `label_seq_id` and Foldseek's
    own per-file residue order themselves.

    Args:
        structures: `entry_id -> Structure`, the structures the rows
            refer to.
        rows: `residue_contacts` rows.

    Returns:
        `{chain_item_id(entry_id, chain_id): {positions...}}` — one
        entry per chain that appears in at least one interface. Pass
        this to `iter_edges()` / `cluster_similar_items()` as
        `interface_residues`, for a search run on files written by
        `export_chain_mmcif()` with the same ids.

    Raises:
        KeyError: A row's entry is not in `structures`.
    """

    result: dict[str, set[int]] = {}
    positions_by_chain: dict[tuple[str, str], dict[int, int]] = {}
    pairs = interface_residue_pairs(rows)
    for (entry_id, _, _), residue_pairs in sorted(pairs.items()):
        if entry_id not in structures:
            raise KeyError(f"no structure for entry {entry_id!r}")
        for residue_id in {r for pair in residue_pairs for r in pair}:
            chain_id, seq = residue_id.rsplit(":", 1)
            key = (entry_id, chain_id)
            if key not in positions_by_chain:
                positions_by_chain[key] = residue_positions(
                    structures[entry_id], chain_id
                )
            position = positions_by_chain[key].get(int(seq))
            if position is not None:
                result.setdefault(chain_item_id(entry_id, chain_id), set()).add(
                    position
                )
    return result


def compute_structure_similarity(
    structures: dict[str, str | Path] | str | Path,
    hits_path: str | Path,
    *,
    foldseek_bin: str = "foldseek",
    sensitivity: float = 9.5,
    alignment_type: int = 2,
    max_seqs: int = 1000,
    exhaustive_search: bool = False,
    tmp_dir: str | Path | None = None,
    foldseek_options: list[str] | None = None,
) -> SimilaritySearch:
    """All-vs-all structural similarity via Foldseek `easy-search`.

    Writes Foldseek's hits to `hits_path` as a TSV with the columns in
    `HIT_COLUMNS["Foldseek"]` (including `qtmscore`/`ttmscore` and the
    alignment ranges) and keeps it there; nothing is loaded into
    memory. Interface-restricted coverage is computed when reading,
    with `iter_edges(..., interface_residues=...)`.

    Args:
        structures: Mapping of item id -> structure file path (PDB/mmCIF,
            optionally gzipped), or a path to a directory of existing
            structure files (ids are then Foldseek's own file-derived
            names).
        hits_path: Where to write the hit TSV. Parent directories are
            created.
        foldseek_bin: Path or name of the Foldseek binary.
        sensitivity: Foldseek `-s` sensitivity value.
        alignment_type: Foldseek `--alignment-type` (0: 3Di alignment,
            1: TM-align, 2: 3Di+AA — Foldseek's own default).
        max_seqs: Foldseek `--max-seqs`: hits kept per query after the
            prefilter. The default (Foldseek's own) can miss similar
            pairs in large families; raise it for leakage control.
        exhaustive_search: Pass `--exhaustive-search 1` (skip the
            prefilter; slower, finds every pair).
        tmp_dir: Working directory for structure/temporary files.
            None = system temp.
        foldseek_options: Additional options passed to Foldseek.

    Returns:
        A `SimilaritySearch` pointing at `hits_path`. Its `parameters`
        are this function's keyword arguments, so the search can be
        re-run with `compute_structure_similarity(structures, path,
        **search.parameters)`.

    Raises:
        RuntimeError: `foldseek_bin` is not on PATH.
        subprocess.CalledProcessError: Foldseek failed.
    """

    if shutil.which(foldseek_bin) is None:
        raise RuntimeError(
            f"foldseek binary {foldseek_bin!r} not found (required for "
            "structure similarity)"
        )

    options = list(foldseek_options or [])
    hits = Path(hits_path).resolve()
    hits.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=tmp_dir) as work_dir:
        work = Path(work_dir)
        if isinstance(structures, dict):
            struct_dir = work / "structures"
            struct_dir.mkdir()
            _link_structures(structures, struct_dir)
        else:
            struct_dir = Path(structures)

        command = [
            foldseek_bin,
            "easy-search",
            str(struct_dir),
            str(struct_dir),
            str(hits),
            str(work / "tmp"),
            "-s",
            str(sensitivity),
            "--alignment-type",
            str(alignment_type),
            "--max-seqs",
            str(max_seqs),
            "--format-output",
            ",".join(HIT_COLUMNS["Foldseek"]),
            "-v",
            "1",
        ]
        if exhaustive_search:
            command += ["--exhaustive-search", "1"]
        subprocess.run(
            command + options, check=True, capture_output=True, text=True
        )

    return SimilaritySearch(
        engine="Foldseek",
        version=_foldseek_version(foldseek_bin),
        hits_path=str(hits),
        columns=list(HIT_COLUMNS["Foldseek"]),
        parameters={
            "foldseek_bin": foldseek_bin,
            "sensitivity": sensitivity,
            "alignment_type": alignment_type,
            "max_seqs": max_seqs,
            "exhaustive_search": exhaustive_search,
            "foldseek_options": options,
        },
        searched_at=now_iso(),
    )
