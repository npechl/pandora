"""Run every docs example script and compare its output with the saved
`.out` file next to it.

The usage pages include both files through `pymdownx.snippets`, so a
passing test means the docs show code that runs and the output it
really prints. After an intended behaviour change, regenerate one with
`uv run python docs/examples/<area>/<name>.py > docs/examples/<area>/<name>.out`.
"""

import runpy
from pathlib import Path

import pytest

ROOT = Path(__file__).parent.parent
EXAMPLES = sorted((ROOT / "docs" / "examples").rglob("*.py"))


def test_examples_exist():
    assert EXAMPLES


@pytest.mark.parametrize(
    "script", EXAMPLES, ids=[str(p.relative_to(ROOT)) for p in EXAMPLES]
)
def test_example_output_matches(script, monkeypatch, capsys):
    monkeypatch.chdir(ROOT)
    runpy.run_path(str(script), run_name="__main__")
    expected = script.with_suffix(".out").read_text()
    assert capsys.readouterr().out == expected
