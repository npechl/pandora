"""Filters: drop labelled interfaces from a `residue_contacts` layer."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Collection, Iterable

from pandora.schemas.annotation import LayerRow, ResidueContactRow

_Interface = tuple[str, frozenset[str]]


def _interface(entry_id: str, chain_1: str, chain_2: str) -> _Interface:
    """An interface key: upper-cased entry id and unordered chain pair."""

    return entry_id.upper(), frozenset((chain_1, chain_2))


def _contact_interface(row: ResidueContactRow) -> _Interface:
    """The interface a contact row belongs to."""

    return _interface(row.entry_id, row.chain_id_1, row.chain_id_2)


def _labels_by_interface(
    labels: Iterable[LayerRow],
) -> dict[_Interface, set[str]]:
    """Every label given to each interface."""

    found: dict[_Interface, set[str]] = defaultdict(set)
    for row in labels:
        key = _interface(row.entry_id, row.chain_id_1, row.chain_id_2)
        found[key].add(row.label)
    return found


def drop_labelled_interfaces(
    contacts: Iterable[ResidueContactRow],
    labels: Iterable[LayerRow],
    drop: Collection[str],
) -> list[ResidueContactRow]:
    """Remove the contact rows of interfaces labelled for dropping.

    An interface is an entry's unordered chain pair, so a label on B–A
    also drops rows written A–B. An interface is dropped if any of its
    label rows carries a `drop` label; interfaces without a label are
    kept.

    Args:
        contacts: The rows of one `residue_contacts` layer.
        labels: The rows of one `interface_labels` layer.
        drop: The label values that mean drop.

    Returns:
        The kept contact rows, in input order.

    Raises:
        ValueError: If a `drop` value is on no label row.
    """

    by_interface = _labels_by_interface(labels)
    unknown = set(drop) - set().union(*by_interface.values())
    if unknown:
        raise ValueError(
            f"no interface has the label(s) {sorted(unknown)}; "
            "check --drop for typos"
        )
    dropped = {key for key, found in by_interface.items() if found & set(drop)}
    return [row for row in contacts if _contact_interface(row) not in dropped]


def count_interface_labels(
    contacts: Iterable[ResidueContactRow],
    labels: Iterable[LayerRow],
    drop: Collection[str],
) -> dict[str, int]:
    """Count how a label layer matches a contact layer's interfaces.

    Args:
        contacts: The rows of one `residue_contacts` layer.
        labels: The rows of one `interface_labels` layer.
        drop: The label values that mean drop.

    Returns:
        `interfaces_in` (distinct interfaces in `contacts`),
        `interfaces_dropped` (of those, with a `drop` label),
        `interfaces_unlabelled` (of those, with no label at all) and
        `labels_unmatched` (labelled interfaces not in `contacts`).
    """

    by_interface = _labels_by_interface(labels)
    present = {_contact_interface(row) for row in contacts}
    return {
        "interfaces_in": len(present),
        "interfaces_dropped": sum(
            bool(by_interface.get(key, set()) & set(drop)) for key in present
        ),
        "interfaces_unlabelled": len(present - by_interface.keys()),
        "labels_unmatched": len(by_interface.keys() - present),
    }
