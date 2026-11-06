"""Import annotation layers produced by external tools.

A layer is a table in one of the `LayerTemplate` shapes. Each row names
an entry and, per side, a chain or residue, in Pandora's naming
(`chain_id`, `seq_id`) or the author's (`auth_chain_id`, `auth_seq_id`,
`ins_code`).
"""

from __future__ import annotations

import math
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from pandora.schemas.annotation import LayerTemplate
from pandora.schemas.structure import Structure

# template -> (sides, level, required value column); free_form infers
# its sides and levels from the columns present.
_TEMPLATES: dict[str, tuple[tuple[str, ...], str, str | None]] = {
    "residue_contacts": (("_1", "_2"), "residue", None),
    "residue_values": (("",), "residue", "value"),
    "residue_labels": (("",), "residue", "label"),
    "chain_values": (("",), "chain", "value"),
    "chain_labels": (("",), "chain", "label"),
    "interface_labels": (("_1", "_2"), "chain", "label"),
}

_KEY_STEMS = ("chain_id", "seq_id", "auth_chain_id", "auth_seq_id", "ins_code")


def _namings(side: str, level: str) -> tuple[tuple[str, ...], ...]:
    """The column sets of Pandora's and the author's naming for a side."""

    if level == "chain":
        return (f"chain_id{side}",), (f"auth_chain_id{side}",)
    return (
        (f"chain_id{side}", f"seq_id{side}"),
        (f"auth_chain_id{side}", f"auth_seq_id{side}"),
    )


def _sides(
    template: LayerTemplate, columns: Collection[str]
) -> list[tuple[str, str]]:
    """(suffix, level) for each side of `template`, given `columns`."""

    if template != "free_form":
        sides, level, _ = _TEMPLATES[template]
        return [(side, level) for side in sides]
    found = []
    for side in ("", "_1", "_2"):
        if not any(f"{stem}{side}" in columns for stem in _KEY_STEMS):
            continue
        residue = {f"seq_id{side}", f"auth_seq_id{side}"} & set(columns)
        found.append((side, "residue" if residue else "chain"))
    return found


def check_layer_columns(
    columns: Collection[str], template: LayerTemplate
) -> None:
    """Check that a table's columns can hold rows of `template`.

    Args:
        columns: The table's column names.
        template: The layer template the table claims to follow.

    Raises:
        ValueError: `entry_id` is absent, a side has neither Pandora's
            nor the author's naming columns complete, or the template's
            required value column is absent.
    """

    if "entry_id" not in columns:
        raise ValueError(
            f"table has no 'entry_id' column, which template "
            f"{template!r} requires"
        )
    for side, level in _sides(template, columns):
        namings = _namings(side, level)
        if not any(set(naming) <= set(columns) for naming in namings):
            options = " or ".join(", ".join(n) for n in namings)
            raise ValueError(
                f"table has no complete key columns for side {side!r} of "
                f"template {template!r}: needs {options}"
            )
    required = _TEMPLATES.get(template, ((), "", None))[2]
    if required is not None and required not in columns:
        raise ValueError(
            f"table has no {required!r} column, which template "
            f"{template!r} requires"
        )


# A residue: (chain_id, seq_id, auth_chain_id, auth_seq_id, ins_code).
_Residue = tuple[str, int | None, str | None, str | None, str | None]
_NO_INS = ("", "?", ".")


def _present(value: Any) -> bool:
    """Whether a cell holds a value (not None, not NaN)."""

    return value is not None and not (
        isinstance(value, float) and math.isnan(value)
    )


def _as_int(value: Any) -> int | None:
    """An integer from an int, an integral float or a digit string."""

    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value) if value.is_integer() else None
    if isinstance(value, str):
        try:
            return int(value.strip())
        except ValueError:
            return None
    return None


def _as_auth(value: Any) -> str:
    """An author residue number as a string (45.0 -> "45")."""

    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _as_ins(value: Any) -> str:
    """An insertion code, "" for none."""

    if not _present(value) or str(value).strip() in _NO_INS:
        return ""
    return str(value).strip()


@dataclass
class _Index:
    """Lookups from either naming to a structure's residues."""

    structure: Structure
    asym_ids: set[str] = field(default_factory=set)
    auth_chains: set[str] = field(default_factory=set)
    by_label: dict[tuple[str, int], _Residue] = field(default_factory=dict)
    by_auth: dict[tuple[str, str, str], list[_Residue]] = field(
        default_factory=dict
    )
    unobserved_label: set[tuple[str, int]] = field(default_factory=set)
    unobserved_auth: set[tuple[str, str, str]] = field(default_factory=set)


def _build_index(structure: Structure) -> _Index:
    """Index a structure's residues by both namings."""

    index = _Index(structure)
    index.asym_ids = {u.id for u in structure.asym_units}
    index.auth_chains = {u.auth_id for u in structure.asym_units if u.auth_id}
    seen: set[_Residue] = set()
    for atom in structure.atoms:
        ins = _as_ins(atom.pdbx_PDB_ins_code)
        residue = (
            atom.label_asym_id,
            atom.label_seq_id,
            atom.auth_asym_id,
            atom.auth_seq_id,
            ins or None,
        )
        if residue in seen:
            continue
        seen.add(residue)
        if atom.label_seq_id is not None:
            index.by_label[(atom.label_asym_id, atom.label_seq_id)] = residue
        if atom.auth_asym_id is not None and atom.auth_seq_id is not None:
            key = (atom.auth_asym_id, atom.auth_seq_id, ins)
            index.by_auth.setdefault(key, []).append(residue)
    for row in structure.raw.get("_pdbx_poly_seq_scheme", []):
        seq_id = _as_int(row.get("seq_id"))
        if row.get("asym_id") and seq_id is not None:
            index.unobserved_label.add((row["asym_id"], seq_id))
        if row.get("pdb_strand_id") and _present(row.get("pdb_seq_num")):
            index.unobserved_auth.add(
                (
                    row["pdb_strand_id"],
                    _as_auth(row["pdb_seq_num"]),
                    _as_ins(row.get("pdb_ins_code")),
                )
            )
    return index


def _by_label(index: _Index, chain: Any, seq: Any) -> _Residue | str:
    """The residue at a Pandora-named position, or a reason."""

    seq_id = _as_int(seq)
    if seq_id is None:
        return "bad_value"
    if chain not in index.asym_ids:
        return "chain_not_found"
    residue = index.by_label.get((chain, seq_id))
    if residue is not None:
        return residue
    if (chain, seq_id) in index.unobserved_label:
        return "no_coordinates"
    return "residue_not_found"


def _by_auth(index: _Index, chain: Any, seq: Any, ins: Any) -> _Residue | str:
    """The residue at an author-named position, or a reason."""

    if chain not in index.auth_chains:
        return "chain_not_found"
    key = (chain, _as_auth(seq), _as_ins(ins))
    hits = index.by_auth.get(key, [])
    if len(hits) > 1:
        return "ambiguous_author_id"
    if hits:
        return hits[0]
    if key in index.unobserved_auth:
        return "no_coordinates"
    return "residue_not_found"


def _resolve_residue(
    index: _Index, row: Mapping[str, Any], side: str
) -> dict[str, Any] | str:
    """The filled key columns of one residue side, or a reason."""

    label = (row.get(f"chain_id{side}"), row.get(f"seq_id{side}"))
    auth = (row.get(f"auth_chain_id{side}"), row.get(f"auth_seq_id{side}"))
    found = []
    if all(_present(v) for v in label):
        found.append(_by_label(index, *label))
    if all(_present(v) for v in auth):
        found.append(_by_auth(index, *auth, row.get(f"ins_code{side}")))
    if not found:
        return "missing_key"
    for result in found:
        if isinstance(result, str):
            return result
    if len(set(found)) > 1:
        return "naming_conflict"
    chain_id, seq_id, auth_chain_id, auth_seq_id, ins_code = found[0]
    return {
        f"chain_id{side}": chain_id,
        f"seq_id{side}": seq_id,
        f"auth_chain_id{side}": auth_chain_id,
        f"auth_seq_id{side}": auth_seq_id,
        f"ins_code{side}": ins_code,
    }


def resolve_layer_rows(
    rows: Sequence[Mapping[str, Any]],
    template: LayerTemplate,
    structures: Mapping[str, Structure],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Link a layer table's rows to the residues and chains they name.

    Each row may name chains and residues in Pandora's naming
    (`chain_id`, `seq_id`) or the author's (`auth_chain_id`,
    `auth_seq_id`, `ins_code`), per side (`_1`/`_2` suffixes for pair
    templates). Entry ids match case-insensitively.

    Args:
        rows: The table's rows, one mapping per row.
        template: The layer template the rows follow.
        structures: The structures the tool ran on, keyed by entry id.

    Returns:
        `(matched, unmatched)`. A matched row is the input row with
        `entry_id` set to the structure's and both namings filled in for
        every side. An unmatched row is the input row plus a `reason`
        (an `UnmatchedReason`). Input rows are not modified.
    """

    by_entry = {s.entry_id.upper(): s for s in structures.values()}
    indexes: dict[str, _Index] = {}
    matched: list[dict[str, Any]] = []
    unmatched: list[dict[str, Any]] = []
    for row in rows:
        result = _resolve_row(row, template, by_entry, indexes)
        if isinstance(result, str):
            unmatched.append({**row, "reason": result})
        else:
            matched.append(result)
    return matched, unmatched


def _resolve_row(
    row: Mapping[str, Any],
    template: LayerTemplate,
    by_entry: Mapping[str, Structure],
    indexes: dict[str, _Index],
) -> dict[str, Any] | str:
    """One row with its keys filled in, or the reason it didn't match."""

    entry = by_entry.get(str(row.get("entry_id", "")).upper())
    if entry is None:
        return "entry_not_found"
    index = indexes.get(entry.entry_id)
    if index is None:
        index = indexes[entry.entry_id] = _build_index(entry)
    filled: dict[str, Any] = {**row, "entry_id": entry.entry_id}
    for side, _level in _sides(template, row):
        keys = _resolve_residue(index, row, side)
        if isinstance(keys, str):
            return keys
        filled.update(keys)
    return filled
