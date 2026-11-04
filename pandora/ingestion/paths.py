from __future__ import annotations

from pathlib import Path


def read_entry_paths(path: str | Path) -> list[Path]:
    """Read mmCIF file paths from a text file, one per line.

    Blank lines and lines starting with `#` are skipped. Relative paths
    are resolved against the list file's folder, so the list can move
    together with its files. Paths are not checked for existence.

    Args:
        path: Path to the text file listing the mmCIF files.

    Returns:
        The listed paths, in file order.
    """

    list_path = Path(path)
    lines = (line.strip() for line in list_path.read_text().splitlines())
    return [
        list_path.parent / line
        for line in lines
        if line and not line.startswith("#")
    ]


def mirror_path(root: str | Path, stem: str, suffix: str = ".cif") -> Path:
    """Path of a file in the PDB mirror layout: `root/<mid>/<stem><suffix>`.

    `<mid>` is the 2nd and 3rd characters of `stem`, as in the PDB
    archive (`1abc.cif` lives in `ab/`). The stem's case is kept
    as given. Nothing is created on disk.

    Args:
        root: The mirror's root folder.
        stem: The file name without suffix, e.g. an entry id (`10kz`) or
            an entry-and-chain id (`10kz_A`).
        suffix: The file suffix, including the dot.

    Returns:
        The file's path under `root`.
    """

    return Path(root) / stem[1:3] / f"{stem}{suffix}"


def list_mirror_files(root: str | Path, suffix: str = ".cif") -> list[Path]:
    """List the files of a PDB mirror layout folder, sorted by path.

    Args:
        root: The mirror's root folder, as passed to `mirror_path`.
        suffix: The file suffix to match, including the dot.

    Returns:
        Every `root/<mid>/*<suffix>` file, sorted; empty if `root` does
        not exist.
    """

    return sorted(Path(root).glob(f"*/*{suffix}"))
