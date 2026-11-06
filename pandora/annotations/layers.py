"""Import annotation layers produced by external tools.

A layer is a table in one of the `LayerTemplate` shapes. Each row names
an entry and, per side, a chain or residue, in Pandora's naming
(`chain_id`, `seq_id`) or the author's (`auth_chain_id`, `auth_seq_id`,
`ins_code`).
"""

from __future__ import annotations

from collections.abc import Collection

from pandora.schemas.annotation import LayerTemplate

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
