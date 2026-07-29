from __future__ import annotations

from pathlib import Path

import gemmi

from pandora.schemas.structure import Structure

_ATOM_SITE_FIELDS = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_PDB_model_num",
]

_SS_FIELDS = [
    "beg_label_asym_id",
    "beg_label_seq_id",
    "end_label_asym_id",
    "end_label_seq_id",
    "beg_auth_asym_id",
    "beg_auth_seq_id",
    "end_auth_asym_id",
    "end_auth_seq_id",
]

# Categories rebuilt from Structure's typed fields — canonicalisation can
# rewrite these (chain remap, altloc selection, ligand filtering, ...), so
# they must win over `structure.raw`, which still holds the pre-transform
# values. `_struct` isn't in this set: nothing canonicalises entry title,
# so it's fine to let it flow straight through from raw.
_TYPED_CATEGORIES = frozenset(
    {
        "_entity",
        "_entity_poly",
        "_struct_asym",
        "_atom_site",
        "_struct_conn",
        "_pdbx_struct_assembly",
        "_pdbx_struct_assembly_gen",
        "_pdbx_struct_oper_list",
        "_struct_conf",
        "_struct_sheet_range",
    }
)


def _unwrap(value):
    """Strip a CIF quote/semicolon delimiter `_cs()` leaves embedded.

    Values read via gemmi's low-level accessor keep their literal
    `'...'`/`"..."` or `;...;` delimiters as part of the string (see
    CLAUDE.md's note on `parsing/mmcif.py::_cs`). Writing them back
    through gemmi's own quoting would double-wrap them, so unwrap first
    and let gemmi decide how to re-quote.
    """

    if not isinstance(value, str) or len(value) < 2:
        return value
    if value[0] == ";" and value.endswith(";"):
        return value[1:-1].rstrip("\n")
    if value[0] in "'\"" and value[-1] == value[0]:
        return value[1:-1]
    return value


def _set_category(
    block: gemmi.cif.Block, name: str, columns: dict[str, list]
) -> None:
    cleaned = {
        col: [_unwrap(v) for v in values] for col, values in columns.items()
    }
    block.set_mmcif_category(name, cleaned)


def structure_to_mmcif(structure: Structure, path: str | Path) -> Path:
    """Write a Structure back out as an mmCIF file.

    Rebuilds the standard mmCIF categories from Structure's typed
    fields (entity/entity_poly/struct_asym/atom_site/struct_conn/
    pdbx_struct_assembly*/struct_conf/struct_sheet_range) and re-emits
    every category preserved verbatim in `structure.raw`.

    This is a lossy round-trip: loop ordering, comments, and any field
    mmCIF supports but Pandora's schema doesn't carry (e.g. anisotropic
    B-factors) won't survive.
    """

    doc = gemmi.cif.Document()
    block = doc.add_new_block(structure.entry_id)

    if "_struct" not in structure.raw and structure.entry.title:
        block.set_pair(
            "_struct.title", gemmi.cif.quote(_unwrap(structure.entry.title))
        )

    _write_entities(block, structure)
    _write_struct_asym(block, structure)
    _write_atom_site(block, structure)
    _write_struct_conn(block, structure)
    _write_assemblies(block, structure)
    _write_secondary_structure(block, structure)
    _write_raw(block, structure)

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    doc.write_file(str(out_path))
    return out_path


def _columns_from_fields(records: list, fields: list[str]) -> dict[str, list]:
    return {field: [getattr(r, field) for r in records] for field in fields}


def _columns(rows: list[dict]) -> dict[str, list]:
    if not rows:
        return {}
    return {key: [row.get(key) for row in rows] for key in rows[0]}


def _write_entities(block: gemmi.cif.Block, structure: Structure) -> None:
    if structure.entities:
        _set_category(
            block,
            "_entity",
            _columns_from_fields(
                structure.entities,
                [
                    "id",
                    "type",
                    "pdbx_description",
                    "formula_weight",
                    "src_method",
                ],
            ),
        )

    poly_rows = [
        {
            "entity_id": entity.id,
            "type": entity.poly.type,
            "pdbx_seq_one_letter_code": (entity.poly.pdbx_seq_one_letter_code),
            "pdbx_seq_one_letter_code_can": (
                entity.poly.pdbx_seq_one_letter_code_can
            ),
            "pdbx_strand_id": entity.poly.pdbx_strand_id,
        }
        for entity in structure.entities
        if entity.poly is not None
    ]
    if poly_rows:
        _set_category(block, "_entity_poly", _columns(poly_rows))


def _write_struct_asym(block: gemmi.cif.Block, structure: Structure) -> None:
    if structure.asym_units:
        _set_category(
            block,
            "_struct_asym",
            _columns_from_fields(structure.asym_units, ["id", "entity_id"]),
        )


def _write_atom_site(block: gemmi.cif.Block, structure: Structure) -> None:
    if structure.atoms:
        _set_category(
            block,
            "_atom_site",
            _columns_from_fields(structure.atoms, _ATOM_SITE_FIELDS),
        )


def _write_struct_conn(block: gemmi.cif.Block, structure: Structure) -> None:
    if not structure.connections:
        return
    rows = [
        {
            "id": conn.id,
            "conn_type_id": conn.conn_type_id,
            "ptnr1_label_asym_id": conn.ptnr1.label_asym_id,
            "ptnr1_label_comp_id": conn.ptnr1.label_comp_id,
            "ptnr1_label_seq_id": conn.ptnr1.label_seq_id,
            "ptnr1_label_atom_id": conn.ptnr1.label_atom_id,
            "ptnr1_auth_asym_id": conn.ptnr1.auth_asym_id,
            "ptnr1_auth_seq_id": conn.ptnr1.auth_seq_id,
            "ptnr1_PDB_ins_code": conn.ptnr1.pdbx_PDB_ins_code,
            "ptnr1_symmetry": conn.ptnr1.symmetry,
            "ptnr2_label_asym_id": conn.ptnr2.label_asym_id,
            "ptnr2_label_comp_id": conn.ptnr2.label_comp_id,
            "ptnr2_label_seq_id": conn.ptnr2.label_seq_id,
            "ptnr2_label_atom_id": conn.ptnr2.label_atom_id,
            "ptnr2_auth_asym_id": conn.ptnr2.auth_asym_id,
            "ptnr2_auth_seq_id": conn.ptnr2.auth_seq_id,
            "ptnr2_PDB_ins_code": conn.ptnr2.pdbx_PDB_ins_code,
            "ptnr2_symmetry": conn.ptnr2.symmetry,
            "pdbx_dist_value": conn.pdbx_dist_value,
            "details": conn.details,
        }
        for conn in structure.connections
    ]
    _set_category(block, "_struct_conn", _columns(rows))


def _write_assemblies(block: gemmi.cif.Block, structure: Structure) -> None:
    if not structure.assemblies:
        return

    _set_category(
        block,
        "_pdbx_struct_assembly",
        _columns_from_fields(
            structure.assemblies,
            [
                "id",
                "details",
                "method_details",
                "oligomeric_details",
                "oligomeric_count",
            ],
        ),
    )

    gen_rows = [
        {
            "assembly_id": gen.assembly_id,
            "oper_expression": gen.oper_expression,
            "asym_id_list": ",".join(gen.asym_id_list),
        }
        for asm in structure.assemblies
        for gen in asm.generators
    ]
    if gen_rows:
        _set_category(block, "_pdbx_struct_assembly_gen", _columns(gen_rows))

    seen_ops: set[str] = set()
    op_rows = []
    for asm in structure.assemblies:
        for op in asm.operators:
            if op.id in seen_ops:
                continue
            seen_ops.add(op.id)
            matrix = op.matrix or [[None] * 3 for _ in range(3)]
            vector = op.vector or [None, None, None]
            row = {"id": op.id, "type": op.type}
            for i in range(3):
                for j in range(3):
                    row[f"matrix[{i + 1}][{j + 1}]"] = matrix[i][j]
                row[f"vector[{i + 1}]"] = vector[i]
            op_rows.append(row)
    if op_rows:
        _set_category(block, "_pdbx_struct_oper_list", _columns(op_rows))


def _write_secondary_structure(
    block: gemmi.cif.Block, structure: Structure
) -> None:
    ss = structure.secondary_structure
    if ss.conf_records:
        _set_category(
            block,
            "_struct_conf",
            _columns_from_fields(
                ss.conf_records, ["id", "conf_type_id", *_SS_FIELDS]
            ),
        )
    if ss.sheet_strands:
        _set_category(
            block,
            "_struct_sheet_range",
            _columns_from_fields(
                ss.sheet_strands, ["sheet_id", "id", *_SS_FIELDS]
            ),
        )


def _write_raw(block: gemmi.cif.Block, structure: Structure) -> None:
    for category, rows in structure.raw.items():
        if category in _TYPED_CATEGORIES:
            continue
        _set_category(block, category, _columns(rows))
