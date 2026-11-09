"""Contact counts: turn a `residue_contacts` layer into value layers."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from typing import Any

from pandora.schemas.annotation import LayerRow, ResidueContactRow

# (entry_id, chain_id, seq_id, auth_seq_id, ins_code)
_Key = tuple[str, str, int | None, str | None, str | None]


def _side(row: ResidueContactRow, side: str) -> tuple[_Key, dict[str, Any]]:
    """One side's residue key and its both-naming columns."""

    def get(name: str) -> Any:
        return getattr(row, f"{name}_{side}")

    key = (
        row.entry_id,
        get("chain_id"),
        get("seq_id"),
        get("auth_seq_id"),
        get("ins_code"),
    )
    columns = {
        "entry_id": row.entry_id,
        "chain_id": get("chain_id"),
        "seq_id": get("seq_id"),
        "auth_chain_id": get("auth_chain_id"),
        "auth_seq_id": get("auth_seq_id"),
        "ins_code": get("ins_code"),
        "comp_id": get("comp_id"),
    }
    return key, columns


def _sort_key(key: tuple[Any, ...]) -> tuple[Any, ...]:
    """Sort key putting None first."""

    return tuple((value is not None, value) for value in key)


def _partners(
    rows: Iterable[ResidueContactRow],
) -> tuple[dict[_Key, set[_Key]], dict[_Key, dict[str, Any]]]:
    """Each residue's distinct partner residues, and its columns."""

    partners: dict[_Key, set[_Key]] = defaultdict(set)
    columns: dict[_Key, dict[str, Any]] = {}
    for row in rows:
        key_1, columns_1 = _side(row, "1")
        key_2, columns_2 = _side(row, "2")
        if key_1 == key_2:
            continue
        partners[key_1].add(key_2)
        partners[key_2].add(key_1)
        columns.setdefault(key_1, columns_1)
        columns.setdefault(key_2, columns_2)
    return partners, columns


def count_residue_contacts(
    rows: Iterable[ResidueContactRow],
) -> list[LayerRow]:
    """Count each residue's distinct contact partners.

    A residue pair given several times, in either order, counts once.
    A row pairing a residue with itself adds nothing. Residues with no
    partner get no row.

    Args:
        rows: The rows of one `residue_contacts` layer.

    Returns:
        A `residue_values` row per residue on either side, with both
        namings (`chain_id`, `seq_id`, `auth_chain_id`, `auth_seq_id`,
        `ins_code`), `comp_id`, and `value` = its number of distinct
        partner residues, sorted by residue.
    """

    partners, columns = _partners(rows)
    return [
        LayerRow(**columns[key], value=len(partners[key]))
        for key in sorted(partners, key=_sort_key)
    ]


def count_chain_contacts(
    rows: Iterable[ResidueContactRow],
) -> list[LayerRow]:
    """Count each chain's distinct residues in contact.

    Args:
        rows: The rows of one `residue_contacts` layer.

    Returns:
        A `chain_values` row per chain with at least one residue in
        contact: `entry_id`, `chain_id`, `auth_chain_id`, and `value` =
        its number of distinct residues in contact, sorted by chain.
    """

    partners, columns = _partners(rows)
    counts: dict[tuple[str, str], int] = defaultdict(int)
    auth: dict[tuple[str, str], str | None] = {}
    for key in partners:
        chain = (key[0], key[1])
        counts[chain] += 1
        auth.setdefault(chain, columns[key]["auth_chain_id"])
    return [
        LayerRow(
            entry_id=entry_id,
            chain_id=chain_id,
            auth_chain_id=auth[(entry_id, chain_id)],
            value=counts[(entry_id, chain_id)],
        )
        for entry_id, chain_id in sorted(counts)
    ]
