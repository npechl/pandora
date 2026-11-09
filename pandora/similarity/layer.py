"""Similarity layers: item IDs, the items a layer may name, and matching
a table's rows against them.

A similarity layer is a table of hits between two items (entries,
chains or interfaces, see `SimilarityItems`). Tables can hold billions
of rows, so everything here works on pyarrow record batches.
"""

from __future__ import annotations

from collections.abc import Collection, Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from pandora.schemas.annotation import SimilarityItems

if TYPE_CHECKING:
    import pyarrow as pa

SIMILARITY_REQUIRED = ("item_id_1", "item_id_2", "score")
SIMILARITY_FLOATS = ("score", "identity", "coverage_1", "coverage_2")
SIMILARITY_INTS = ("start_1", "end_1", "start_2", "end_2")

_CHAINS = {"entry": 0, "chain": 1, "interface": 2}
_ITEM_COLUMNS = ("item_id_1", "item_id_2")


@dataclass(frozen=True)
class ItemSet:
    """The items a similarity layer may name, at one level.

    Attributes:
        level: What the items are.
        ids: Lookup key (entry part upper-cased, interface chains
            sorted) -> the item ID as stored.
        entries: Upper-cased entry IDs that have items.
        chains: (upper-cased entry ID, chain ID) pairs, for chain items;
            empty otherwise.
    """

    level: SimilarityItems
    ids: Mapping[str, str]
    entries: frozenset[str]
    chains: frozenset[tuple[str, str]] = frozenset()


def item_id(level: SimilarityItems, entry_id: str, *chain_ids: str) -> str:
    """Build the item ID of an entry, chain or interface.

    Args:
        level: What the item is.
        entry_id: The item's entry ID.
        *chain_ids: No chain for an entry, one for a chain, two for an
            interface (stored sorted).

    Returns:
        `E`, `E_C` or `E_C1_C2`.

    Raises:
        ValueError: The number of chains doesn't fit `level`.
    """

    if len(chain_ids) != _CHAINS[level]:
        raise ValueError(
            f"a {level} item needs {_CHAINS[level]} chain ID(s), got "
            f"{len(chain_ids)}"
        )
    if level == "interface":
        chain_ids = tuple(sorted(chain_ids))
    return "_".join((entry_id, *chain_ids))


def parse_item_id(level: SimilarityItems, item: str) -> tuple[str, ...]:
    """Split an item ID into its entry and chain IDs.

    Splits from the right: chain IDs never contain `_`, extended PDB
    entry IDs (`pdb_00001abc`) do.

    Args:
        level: What the item is.
        item: The item ID.

    Returns:
        `(entry,)`, `(entry, chain)` or `(entry, chain_1, chain_2)`.

    Raises:
        ValueError: The ID has too few parts, or an empty one.
    """

    n_chains = _CHAINS[level]
    parts = tuple(item.rsplit("_", n_chains)) if n_chains else (item,)
    if len(parts) != n_chains + 1 or not all(parts):
        raise ValueError(f"{item!r} is not a {level} item ID")
    return parts


def _key(level: SimilarityItems, parts: tuple[str, ...]) -> str:
    """Lookup key of a parsed item ID."""

    chains = sorted(parts[1:]) if level == "interface" else parts[1:]
    return "_".join((parts[0].upper(), *chains))


def entry_items(entry_ids: Iterable[str]) -> ItemSet:
    """The entry items of the given entries.

    Args:
        entry_ids: Entry IDs, spelled as they should be stored.

    Returns:
        An `ItemSet` of level "entry".
    """

    ids = {entry.upper(): entry for entry in entry_ids}
    return ItemSet("entry", ids, frozenset(ids))


def chain_items(chains: Iterable[tuple[str, str]]) -> ItemSet:
    """The chain items of the given chains.

    Args:
        chains: (entry ID, chain ID) pairs (`label_asym_id`).

    Returns:
        An `ItemSet` of level "chain".
    """

    pairs = list(chains)
    ids = {
        _key("chain", (entry, chain)): item_id("chain", entry, chain)
        for entry, chain in pairs
    }
    return ItemSet(
        "chain",
        ids,
        frozenset(entry.upper() for entry, _ in pairs),
        frozenset((entry.upper(), chain) for entry, chain in pairs),
    )


def interface_items(rows: Iterable[Mapping[str, Any]]) -> ItemSet:
    """The interface items of a `residue_contacts` layer.

    An interface is a pair of chains with at least one row where both
    sides are polymer residues (`seq_id` present) on different chains,
    as in `extract_interface_records`.

    Args:
        rows: Contact rows as mappings with `entry_id`, `chain_id_1`,
            `chain_id_2`, `seq_id_1` and `seq_id_2`.

    Returns:
        An `ItemSet` of level "interface".
    """

    ids: dict[str, str] = {}
    entries: set[str] = set()
    for row in rows:
        entries.add(row["entry_id"].upper())
        if row["seq_id_1"] is None or row["seq_id_2"] is None:
            continue
        chains = (row["chain_id_1"], row["chain_id_2"])
        if chains[0] == chains[1]:
            continue
        stored = item_id("interface", row["entry_id"], *chains)
        ids.setdefault(_key("interface", (row["entry_id"], *chains)), stored)
    return ItemSet("interface", ids, frozenset(entries))


def check_similarity_columns(columns: Collection[str]) -> None:
    """Check that a table's columns can hold `similarity` rows.

    Args:
        columns: The table's column names.

    Raises:
        ValueError: A required column is absent.
    """

    missing = [c for c in SIMILARITY_REQUIRED if c not in columns]
    if missing:
        raise ValueError(
            f"table has no {', '.join(missing)} column(s), which the "
            "similarity template requires"
        )


def match_similarity_batches(
    batches: Iterable[pa.RecordBatch], items: ItemSet
) -> Iterator[tuple[pa.RecordBatch, pa.RecordBatch]]:
    """Link each row of a similarity table to the items it names.

    Self-hits and repeated pairs are kept as given.

    Args:
        batches: The table, as record batches with the similarity
            template's columns (see `check_similarity_columns`).
        items: The items rows may name.

    Returns:
        An iterator of one `(matched, unmatched)` pair per batch. In
        `matched` the item columns hold the stored item IDs, recognised
        number columns are float64 / int64, and other columns are kept
        as given. `unmatched` holds every input column as a string, plus
        `reason`: `bad_item_id`, `entry_not_found`, `chain_not_found`,
        `interface_not_found` (side 1 checked first, then side 2), or
        `bad_value` (a recognised column that isn't a number, or no
        `score`).
    """

    import pyarrow as pa

    stored = pa.array(sorted(set(items.ids.values())), pa.string())
    for batch in batches:
        yield _match_batch(batch, items, stored)


def _match_batch(
    batch: pa.RecordBatch, items: ItemSet, stored: pa.Array
) -> tuple[pa.RecordBatch, pa.RecordBatch]:
    """`match_similarity_batches` for one batch.

    Rows whose item IDs are stored IDs and whose values are numbers pass
    column-wise; only the rest are looked at one by one.
    """

    import pyarrow as pa
    import pyarrow.compute as pc

    n_rows = batch.num_rows
    ids = {
        name: pc.cast(batch.column(name), pa.string()) for name in _ITEM_COLUMNS
    }
    good = pc.fill_null(
        pc.and_(*(pc.is_in(ids[name], value_set=stored) for name in ids)),
        False,
    )
    values_ok = pa.repeat(True, n_rows)
    numbers: dict[str, pa.Array] = {}
    for names, kind in ((SIMILARITY_FLOATS, float), (SIMILARITY_INTS, int)):
        for name in names:
            if name in batch.schema.names:
                array, ok = _numbers(batch.column(name), kind, name == "score")
                numbers[name] = array
                values_ok = pc.and_(values_ok, ok)
    good = pc.and_(good, values_ok)

    slow_matched: list[int] = []
    replaced: dict[str, list[str]] = {name: [] for name in ids}
    reasons: list[str] = []
    for row in pc.indices_nonzero(pc.invert(good)).to_pylist():
        resolved: dict[str, str] = {}
        reason = None
        for name, column in ids.items():
            resolved[name], reason = _resolve(column[row].as_py(), items)
            if reason is not None:
                break
        if reason is None and not values_ok[row].as_py():
            reason = "bad_value"
        if reason is None:
            slow_matched.append(row)
            for name in ids:
                replaced[name].append(resolved[name])
        else:
            reasons.append(reason)

    matched = good
    if slow_matched:
        slow = [False] * n_rows
        for row in slow_matched:
            slow[row] = True
        slow_mask = pa.array(slow)
        matched = pc.or_(good, slow_mask)
        ids = {
            name: pc.replace_with_mask(
                column, slow_mask, pa.array(replaced[name], pa.string())
            )
            for name, column in ids.items()
        }

    columns = [
        ids.get(name, numbers.get(name, batch.column(name)))
        for name in batch.schema.names
    ]
    matched_batch = pa.RecordBatch.from_arrays(
        columns, names=batch.schema.names
    ).filter(matched)
    rest = batch.filter(pc.invert(matched))
    unmatched_batch = pa.RecordBatch.from_arrays(
        [pc.cast(c, pa.string()) for c in rest.columns]
        + [pa.array(reasons, pa.string())],
        names=[*batch.schema.names, "reason"],
    )
    return matched_batch, unmatched_batch


def _numbers(
    column: pa.Array, kind: type, required: bool
) -> tuple[pa.Array, pa.BooleanArray]:
    """`column` as float64/int64, and per row whether its value is usable
    (a number, or empty unless `required`)."""

    import pyarrow as pa
    import pyarrow.compute as pc
    from pyarrow.lib import ArrowInvalid, ArrowNotImplementedError

    target = pa.float64() if kind is float else pa.int64()
    try:
        array = pc.cast(column, target)
        ok = pa.repeat(True, len(column))
    except (ArrowInvalid, ArrowNotImplementedError):
        raw = column.to_pylist()
        values = [_number(v, kind) for v in raw]
        array = pa.array(values, target)
        ok = pa.array(
            [
                v is not None or r is None
                for v, r in zip(values, raw, strict=True)
            ]
        )
    if required:
        ok = pc.and_(ok, pc.is_valid(array))
    return array, ok


def _number(value: Any, kind: type) -> float | int | None:
    """`value` as a float or int, or None if it isn't one."""

    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if kind is int:
        return int(number) if number.is_integer() else None
    return number


def _resolve(raw: str | None, items: ItemSet) -> tuple[str | None, str | None]:
    """(stored item ID, None), or (None, reason) if `raw` names no item."""

    if raw is None:
        return None, "bad_item_id"
    try:
        parts = parse_item_id(items.level, raw)
    except ValueError:
        return None, "bad_item_id"
    stored = items.ids.get(_key(items.level, parts))
    if stored is not None:
        return stored, None
    if parts[0].upper() not in items.entries:
        return None, "entry_not_found"
    return None, f"{items.level}_not_found"
