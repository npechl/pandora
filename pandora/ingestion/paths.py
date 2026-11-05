from __future__ import annotations

import os
import re
from collections.abc import Iterable
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


def write_entry_paths(paths: Iterable[Path], path: str | Path) -> None:
    """Write mmCIF file paths to a text file, one per line.

    The counterpart of `read_entry_paths`. Parent folders are created.

    Args:
        paths: The paths to write, in the order given.
        path: Path of the text file to write.
    """

    list_path = Path(path)
    list_path.parent.mkdir(parents=True, exist_ok=True)
    list_path.write_text("".join(f"{p}\n" for p in paths))


def template_entry_paths(
    root: str | Path, template: str, ids: list[str] | None = None
) -> list[Path]:
    """List an archive's files from a path template, reading as few
    folders as possible.

    `template` is a path relative to `root` with the placeholders `{id}`
    (the entry ID) and `{mid}` (`id[1:3]`, the PDB mirror folder), e.g.
    `{mid}/{id}/clean_mmcif/{id}_updated.cif.gz` for the PDBe release
    tree. The first path segment containing `{id}` is the ID segment.

    With `ids`, the paths are built from them (lowercased) without
    reading the filesystem. Without, the folders above the ID segment
    are listed (`{mid}` matches any folder name), and the names there
    that match the ID segment give the IDs, as spelled on disk; each path
    continues from the folder where its ID was found. Names starting
    with `.` are skipped, and so are files where the template continues
    below the ID. Nothing below the ID segment is read, and no path is
    checked for existence.

    Args:
        root: The archive's root folder.
        template: The path template, relative to `root`.
        ids: Entry IDs to build paths for, or None to list them.

    Returns:
        The paths, sorted.

    Raises:
        ValueError: If no segment of `template` contains `{id}`, or the
            ID segment also contains `{mid}`.
    """

    root = Path(root)
    parts = template.split("/")
    k = next((i for i, part in enumerate(parts) if "{id}" in part), None)
    if k is None:
        raise ValueError(f"template {template!r} has no {{id}}")
    if "{mid}" in parts[k]:
        raise ValueError(f"template segment {parts[k]!r} mixes {{mid}}, {{id}}")
    if ids is not None:
        return sorted(
            root / template.format(id=i.lower(), mid=i.lower()[1:3])
            for i in ids
        )
    dirs = [root]
    for part in parts[:k]:
        pattern = part.replace("{mid}", "*")
        dirs = [
            c
            for d in dirs
            for c in d.glob(pattern)
            if c.is_dir() and not c.name.startswith(".")
        ]
    id_re = re.compile(re.escape(parts[k]).replace(r"\{id\}", "(.+)"))
    rest = parts[k + 1 :]
    found = []
    for d in dirs:
        with os.scandir(d) as names:
            for entry in names:
                if entry.name.startswith(".") or (rest and not entry.is_dir()):
                    continue
                if m := id_re.fullmatch(entry.name):
                    i = m.group(1)
                    below = (r.format(id=i, mid=i[1:3]) for r in rest)
                    found.append(Path(d, entry.name, *below))
    return sorted(found)
