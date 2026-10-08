"""Run the examples on the usage pages and check they print what the
page shows.

An example is a ```python or ```bash block followed directly by a
```text block holding its output. Python runs from the repository root.
Shell blocks run in a temporary directory with `datasets/` linked in,
so `pandora` commands can read the fixtures and write their output
directories freely; they need `jq` on PATH. A ```yaml block whose first
line is `# <path>` must show that file exactly.
"""

import contextlib
import io
import os
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
PAGES = [
    "docs/usage/canonicalisation.md",
    "docs/usage/datasets.md",
    "docs/recipes/benchmark-quality-filters.md",
]
FENCE = re.compile(
    r"^(?P<indent>[ ]*)```(?P<lang>\w*)[^\n]*\n(?P<body>.*?)^(?P=indent)```[ ]*$",
    re.MULTILINE | re.DOTALL,
)


def _examples() -> list[tuple[str, str, str, str]]:
    """(id, language, code, expected output) for every example."""

    found = []
    for page in PAGES:
        text = (ROOT / page).read_text()
        blocks = [
            (m.group("lang"), textwrap.dedent(m.group("body")), m.start())
            for m in FENCE.finditer(text)
        ]
        for (lang, code, start), (next_lang, output, _) in zip(
            blocks, blocks[1:]
        ):
            if lang in ("python", "bash") and next_lang == "text":
                line = text.count("\n", 0, start) + 1
                found.append((f"{page}:{line}", lang, code, output))
    return found


EXAMPLES = _examples()


def test_pages_have_examples():
    assert len(EXAMPLES) > 20


@pytest.mark.parametrize(
    "lang,code,expected",
    [e[1:] for e in EXAMPLES],
    ids=[e[0] for e in EXAMPLES],
)
def test_example_prints_what_the_page_shows(
    lang, code, expected, tmp_path, monkeypatch
):
    if lang == "python":
        monkeypatch.chdir(ROOT)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            exec(compile(code, "<example>", "exec"), {"__name__": "__main__"})
        assert out.getvalue() == expected
        return

    if shutil.which("jq") is None:
        pytest.skip("jq is not on PATH")
    (tmp_path / "datasets").symlink_to(ROOT / "datasets")
    result = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", code],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        env=os.environ,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == expected


def _file_blocks() -> list[tuple[str, str, str]]:
    """(id, path, body) for every ```yaml block whose first line is
    `# <path>`."""

    found = []
    for page in PAGES:
        text = (ROOT / page).read_text()
        for m in FENCE.finditer(text):
            body = textwrap.dedent(m.group("body"))
            first, _, rest = body.partition("\n")
            if m.group("lang") == "yaml" and first.startswith("# "):
                line = text.count("\n", 0, m.start()) + 1
                found.append((f"{page}:{line}", first[2:].strip(), rest))
    return found


FILE_BLOCKS = _file_blocks()


@pytest.mark.parametrize(
    "path,body", [b[1:] for b in FILE_BLOCKS], ids=[b[0] for b in FILE_BLOCKS]
)
def test_shown_file_matches_repository(path, body):
    assert body == (ROOT / path).read_text()
