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

from pandora.annotations.entry import polymer_asym_ids
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

# Templates whose rows must be one per key; the rest keep repeats.
_UNIQUE = {
    "residue_values",
    "residue_labels",
    "chain_values",
    "chain_labels",
    "interface_labels",
}


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
    auth_of: dict[str, str | None] = field(default_factory=dict)
    polymers_by_auth: dict[str, list[str]] = field(default_factory=dict)
    by_label: dict[tuple[str, int], _Residue] = field(default_factory=dict)
    by_auth: dict[tuple[str, str, str], list[_Residue]] = field(
        default_factory=dict
    )
    unobserved_label: set[tuple[str, int]] = field(default_factory=set)
    unobserved_auth: set[tuple[str, str, str]] = field(default_factory=set)
    # A residue's name, from its first atom (microheterogeneity: one
    # residue, several altloc residue names).
    comp_of: dict[_Residue, str] = field(default_factory=dict)


def _build_index(structure: Structure) -> _Index:
    """Index a structure's residues by both namings."""

    index = _Index(structure)
    index.asym_ids = {u.id for u in structure.asym_units}
    index.auth_chains = {u.auth_id for u in structure.asym_units if u.auth_id}
    index.auth_of = {u.id: u.auth_id for u in structure.asym_units}
    polymers = polymer_asym_ids(structure)
    for unit in structure.asym_units:
        if unit.id in polymers and unit.auth_id:
            index.polymers_by_auth.setdefault(unit.auth_id, []).append(unit.id)
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
        index.comp_of[residue] = atom.label_comp_id
        if atom.label_seq_id is not None:
            index.by_label[(atom.label_asym_id, atom.label_seq_id)] = residue
        if atom.auth_asym_id is not None and atom.auth_seq_id is not None:
            key = (atom.auth_asym_id, atom.auth_seq_id, ins)
            index.by_auth.setdefault(key, []).append(residue)
    _index_sequence(index, structure)
    return index


def _index_sequence(index: _Index, structure: Structure) -> None:
    """Add the sequence positions (with or without atoms) to `index`.

    `_entity_poly_seq` survives canonicalisation, which can drop
    `_pdbx_poly_seq_scheme`; it covers Pandora's naming only.
    """

    asyms_by_entity: dict[str, list[str]] = {}
    for unit in structure.asym_units:
        asyms_by_entity.setdefault(unit.entity_id, []).append(unit.id)
    for row in structure.raw.get("_entity_poly_seq", []):
        num = _as_int(row.get("num"))
        for asym_id in asyms_by_entity.get(row.get("entity_id") or "", []):
            if num is not None:
                index.unobserved_label.add((asym_id, num))
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
        f"comp_id{side}": index.comp_of[found[0]],
    }


def _resolve_chain(
    index: _Index, row: Mapping[str, Any], side: str
) -> dict[str, Any] | str:
    """The filled key columns of one chain side, or a reason.

    An author chain id resolves among polymer chains only: ligands and
    waters usually share their polymer's author chain id.
    """

    label = row.get(f"chain_id{side}")
    auth = row.get(f"auth_chain_id{side}")
    if not _present(label) and not _present(auth):
        return "missing_key"
    chains: set[str] = set()
    if _present(label):
        if label not in index.asym_ids:
            return "chain_not_found"
        chains.add(label)
    if _present(auth):
        hits = index.polymers_by_auth.get(auth, [])
        if not hits:
            return "chain_not_found"
        if len(hits) > 1:
            return "ambiguous_author_id"
        chains.add(hits[0])
    if len(chains) > 1:
        return "naming_conflict"
    (chain_id,) = chains
    return {
        f"chain_id{side}": chain_id,
        f"auth_chain_id{side}": index.auth_of.get(chain_id),
    }


def _bad_value(row: Mapping[str, Any], template: LayerTemplate) -> bool:
    """Whether the row's required value has the wrong type."""

    required = _TEMPLATES.get(template, ((), "", None))[2]
    value = row.get(required) if required else None
    if required == "value":
        return (
            isinstance(value, bool)
            or not isinstance(value, int | float)
            or not _present(value)
        )
    if required == "label":
        return not isinstance(value, str) or not value
    if template == "residue_contacts" and _present(row.get("distance")):
        return not _is_number(row["distance"])
    return False


def _is_number(value: Any) -> bool:
    """Whether value is a number or a string holding one (not a bool)."""

    if isinstance(value, bool):
        return False
    if isinstance(value, int | float):
        return True
    try:
        float(str(value))
    except ValueError:
        return False
    return True


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
    seen: set[tuple[Any, ...]] = set()
    sides = _TEMPLATES.get(template, ((), "", None))[0]
    for row in rows:
        result = _resolve_row(row, template, by_entry, indexes)
        if not isinstance(result, str) and template in _UNIQUE:
            key = tuple(
                result.get(f"{stem}{side}")
                for side in sides
                for stem in _KEY_STEMS
            )
            if (result["entry_id"], key) in seen:
                result = "duplicate_key"
            else:
                seen.add((result["entry_id"], key))
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
    present = [column for column, value in row.items() if _present(value)]
    for side, level in _sides(template, present):
        resolve = _resolve_residue if level == "residue" else _resolve_chain
        keys = resolve(index, row, side)
        if isinstance(keys, str):
            return keys
        filled.update(keys)
    if _bad_value(row, template):
        return "bad_value"
    return filled
