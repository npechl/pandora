"""Contact counts: turn a `residue_contacts` layer into value layers."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
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


def label_interfaces(
    contacts: Iterable[ResidueContactRow],
    lengths: Mapping[tuple[str, str], int],
    min_residues_per_chain: int | None = None,
    max_total_residues: int | None = None,
) -> list[LayerRow]:
    """Label each interface by its contacts per side and its total size.

    An interface is an entry's unordered chain pair. Each side counts
    its distinct residues in contact with the other chain; a residue
    pair given several times, in either order, counts once. The first
    matching label wins:

    1. `too_large`: `max_total_residues` is set and the two chains'
       lengths add up to more (an unknown length never matches);
    2. `few_contacts`: `min_residues_per_chain` is set and the side
       with fewer residues in contact has fewer than that;
    3. `ok` otherwise.

    Args:
        contacts: The rows of one `residue_contacts` layer.
        lengths: Chain length by `(entry_id.upper(), chain_id)`.
            Chains missing here have an unknown length.
        min_residues_per_chain: The fewest residues in contact each
            side needs, or None to skip the rule.
        max_total_residues: The most residues the two chains may have
            together, or None to skip the rule.

    Returns:
        An `interface_labels` row per interface, sorted by entry and
        chain pair, with the smaller chain id in `chain_id_1`:
        `entry_id`, `chain_id_1`, `chain_id_2`, `auth_chain_id_1`,
        `auth_chain_id_2`, `residues_1`, `residues_2`, `length_1`,
        `length_2`, `total_residues` (None if a length is unknown) and
        `label`.
    """

    residues: dict[tuple[str, str, str], set[_Key]] = defaultdict(set)
    entry: dict[str, str] = {}
    auth: dict[tuple[str, str], str | None] = {}
    for row in contacts:
        key_1, columns_1 = _side(row, "1")
        key_2, columns_2 = _side(row, "2")
        upper = row.entry_id.upper()
        entry.setdefault(upper, row.entry_id)
        chain_1, chain_2 = key_1[1], key_2[1]
        residues[(upper, chain_1, chain_2)].add(key_1)
        residues[(upper, chain_2, chain_1)].add(key_2)
        auth.setdefault((upper, chain_1), columns_1["auth_chain_id"])
        auth.setdefault((upper, chain_2), columns_2["auth_chain_id"])

    pairs = sorted({(e, *sorted((a, b))) for e, a, b in residues})
    out = []
    for upper, chain_1, chain_2 in pairs:
        count_1 = len(residues[(upper, chain_1, chain_2)])
        count_2 = len(residues[(upper, chain_2, chain_1)])
        length_1 = lengths.get((upper, chain_1))
        length_2 = lengths.get((upper, chain_2))
        total = (
            None
            if length_1 is None or length_2 is None
            else length_1 + length_2
        )
        if (
            max_total_residues is not None
            and total is not None
            and (total > max_total_residues)
        ):
            label = "too_large"
        elif min_residues_per_chain is not None and (
            min(count_1, count_2) < min_residues_per_chain
        ):
            label = "few_contacts"
        else:
            label = "ok"
        out.append(
            LayerRow(
                entry_id=entry[upper],
                chain_id_1=chain_1,
                chain_id_2=chain_2,
                auth_chain_id_1=auth[(upper, chain_1)],
                auth_chain_id_2=auth[(upper, chain_2)],
                residues_1=count_1,
                residues_2=count_2,
                length_1=length_1,
                length_2=length_2,
                total_residues=total,
                label=label,
            )
        )
    return out
