from __future__ import annotations

from pandora.schemas.structure import AtomSiteRecord

# A raw category with a column whose name contains one of these refers to
# chains or residues by id.
_INDEX_COLUMN_MARKERS = ("asym_id", "seq_id", "strand_id")


def _residue_keys(atoms: list[AtomSiteRecord]) -> set[tuple[str, int | None]]:
    """(chain, residue) ids present in atoms."""

    return {(a.label_asym_id, a.label_seq_id) for a in atoms}


def _drop_stale_raw_categories(
    raw: dict[str, list[dict[str, str | None]]],
    original_atoms: list[AtomSiteRecord],
    canonical_atoms: list[AtomSiteRecord],
) -> tuple[dict[str, list[dict[str, str | None]]], list[str]]:
    """Raw without chain/residue-indexed categories, if canonicalisation
    changed which chains or residues exist; plus the dropped names."""

    if _residue_keys(original_atoms) == _residue_keys(canonical_atoms):
        return raw, []
    stale = sorted(
        category
        for category, rows in raw.items()
        if any(
            marker in column
            for row in rows[:1]
            for column in row
            for marker in _INDEX_COLUMN_MARKERS
        )
    )
    kept = {c: rows for c, rows in raw.items() if c not in stale}
    return kept, stale
