from __future__ import annotations

from pandora.similarity.structure import _interface_coverage


def test_interface_coverage_full_and_partial_overlap() -> None:
    assert _interface_coverage({10, 11, 12}, start=1, end=20) == 1.0
    assert _interface_coverage({10, 11, 12, 13}, start=1, end=11) == 0.5


def test_interface_coverage_no_overlap() -> None:
    assert _interface_coverage({1, 2, 3}, start=10, end=20) == 0.0


def test_interface_coverage_empty_residues() -> None:
    assert _interface_coverage(set(), start=1, end=20) == 0.0


if __name__ == "__main__":
    test_interface_coverage_full_and_partial_overlap()
    test_interface_coverage_no_overlap()
    test_interface_coverage_empty_residues()
    print("ok")
