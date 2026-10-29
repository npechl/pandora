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
