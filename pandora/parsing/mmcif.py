from collections.abc import Callable
from functools import partial
from typing import TypeVar

import gemmi

from pandora.schemas.common import Diagnostic, DiagnosticBundle, ResultStatus
from pandora.schemas.structure import (
    AssemblyGenRecord,
    AssemblyOperRecord,
    AssemblyRecord,
    AsymRecord,
    AtomSiteRecord,
    ConfRecord,
    ConnPartner,
    ConnRecord,
    EntityPolyRecord,
    EntityRecord,
    EntryRecord,
    SheetStrandRecord,
    SSRecord,
    Structure,
)

_ENTITY_TYPE_MAP = {
    gemmi.EntityType.Polymer: "polymer",
    gemmi.EntityType.NonPolymer: "non-polymer",
    gemmi.EntityType.Water: "water",
    gemmi.EntityType.Branched: "branched",
}

# mmcif_to_structure helpers --------------------------------------------
_NULL_CIF = frozenset({".", "?"})
_SKIP_RAW = frozenset({"_atom_site", "_atom_site_anisou"})

# Malformed/inconsistent optional mmCIF categories should be skipped, not
# crash the whole parse; but we don't want to swallow unrelated bugs
# (e.g. MemoryError, KeyboardInterrupt), so we name the concrete types.
_MALFORMED_CIF_ERRORS = (
    RuntimeError,
    ValueError,
    KeyError,
    IndexError,
    TypeError,
    AttributeError,
)


def _cs(v: str) -> str | None:
    """mmCIF string value with CIF quoting removed, or None if it's a null
    token ("."/"?")."""

    return None if v in _NULL_CIF else gemmi.cif.as_string(v)


def _cseq(v: str) -> str | None:
    """mmCIF one-letter sequence value with line breaks removed."""

    seq = _cs(v)
    return "".join(seq.split()) if seq is not None else None


def _ci(v: str) -> int | None:
    """mmCIF value as int, or None if null/not parseable."""

    if v in _NULL_CIF:
        return None
    try:
        return int(v)
    except (ValueError, TypeError):
        return None


def _cf(v: str) -> float | None:
    """mmCIF value as float, or None if null/not parseable."""

    if v in _NULL_CIF:
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def _extract_raw(block: object) -> dict[str, list[dict[str, str | None]]]:
    """Every mmCIF category in block, as raw rows, excluding the
    categories in _SKIP_RAW."""

    raw: dict[str, list[dict[str, str | None]]] = {}
    for item in block:  # type: ignore[union-attr]
        loop = getattr(item, "loop", None)
        if loop is not None and loop.tags:
            category = loop.tags[0].rsplit(".", 1)[0]
            if category in _SKIP_RAW:
                continue
            col_names = [t.rsplit(".", 1)[1] for t in loop.tags]
            n_cols = len(col_names)
            raw[category] = [
                {col_names[j]: _cs(loop[i, j]) for j in range(n_cols)}
                for i in range(loop.length())
            ]
            continue

        pair = getattr(item, "pair", None)
        if pair:
            tag, value = pair
            category, col_name = tag.rsplit(".", 1)
            if category in _SKIP_RAW:
                continue
            raw.setdefault(category, [{}])[0][col_name] = _cs(value)
    return raw


def _load(
    path_to_mmcif: str, diag: DiagnosticBundle
) -> tuple[gemmi.Structure, gemmi.cif.Block] | None:
    """gemmi structure and CIF block of the file, or None (with an error
    in diag) if it can't be read or has no model."""

    def _fail(code: str, message: str) -> None:
        diag.errors.append(
            Diagnostic(code=code, severity="error", message=message)
        )

    if not path_to_mmcif.strip():
        return _fail("EMPTY_CONTENT", "mmCIF content is empty")
    try:
        st = gemmi.read_structure(path_to_mmcif, format=gemmi.CoorFormat.Mmcif)
    except (OSError, ValueError, RuntimeError) as exc:
        return _fail("PARSE_ERROR", str(exc))
    try:
        block = gemmi.cif.read(path_to_mmcif).sole_block()
    except (OSError, ValueError, RuntimeError) as exc:
        return _fail("CIF_PARSE_ERROR", str(exc))
    if len(st) == 0:
        return _fail("NO_MODEL", "Structure contains no models")
    return st, block


def _select_model(
    st: gemmi.Structure, model_num: int, diag: DiagnosticBundle
) -> gemmi.Model:
    """Model model_num of st, or its first model with a warning."""

    model = next((m for m in st if m.num == model_num), None)
    if model is None:
        model = st[0]
        diag.warnings.append(
            Diagnostic(
                code="MODEL_NOT_FOUND",
                severity="warning",
                message=(
                    f"Model {model_num} not found; falling back "
                    f"to model {model.num}"
                ),
                entry_id=st.name,
            )
        )
    return model


_RowT = TypeVar("_RowT")


def _read_rows(
    block: gemmi.cif.Block,
    eid: str,
    diag: DiagnosticBundle,
    category: str,
    tags: list[str],
    build: Callable[[gemmi.cif.Table.Row], _RowT],
    code: str,
) -> list[_RowT]:
    """build(row) for each row of an optional category. On malformed data,
    keep the rows read so far and add a `code` warning to diag."""

    out: list[_RowT] = []
    try:
        for row in block.find(category, tags):
            out.append(build(row))
    except _MALFORMED_CIF_ERRORS as exc:
        diag.warnings.append(
            Diagnostic(
                code=code, severity="warning", message=str(exc), entry_id=eid
            )
        )
    return out


# Entities ------------------------------------------------------------
_ENTITY_TAGS = ["id", "pdbx_description", "formula_weight", "src_method"]
_ENTITY_POLY_TAGS = [
    "entity_id",
    "type",
    "pdbx_seq_one_letter_code",
    "pdbx_seq_one_letter_code_can",
    "pdbx_strand_id",
]


def _entity_row(
    row: gemmi.cif.Table.Row,
) -> tuple[str, tuple[str | None, float | None, str | None]]:
    """(entity id, (description, formula weight, source method))."""

    return _cs(row[0]) or "", (_cs(row[1]), _cf(row[2]), _cs(row[3]))


def _entity_poly_row(
    row: gemmi.cif.Table.Row,
) -> tuple[str, EntityPolyRecord]:
    """(entity id, its EntityPolyRecord)."""

    return _cs(row[0]) or "", EntityPolyRecord(
        type=_cs(row[1]),
        pdbx_seq_one_letter_code=_cseq(row[2]),
        pdbx_seq_one_letter_code_can=_cseq(row[3]),
        pdbx_strand_id=_cs(row[4]),
    )


def _entities(
    st: gemmi.Structure,
    desc_map: dict[str, tuple[str | None, float | None, str | None]],
    poly_map: dict[str, EntityPolyRecord],
) -> list[EntityRecord]:
    """One EntityRecord per gemmi entity, with its _entity/_entity_poly
    metadata."""

    entities: list[EntityRecord] = []
    for ent in st.entities:
        desc, fw, src = desc_map.get(ent.name, (None, None, None))
        entities.append(
            EntityRecord(
                id=ent.name,
                type=_ENTITY_TYPE_MAP.get(ent.entity_type, "non-polymer"),
                pdbx_description=desc,
                formula_weight=fw,
                src_method=src,
                poly=poly_map.get(ent.name),
            )
        )
    return entities


# Atoms and asym units ------------------------------------------------
def _residue_seq(res: gemmi.Residue) -> tuple[int | None, str | None]:
    """label_seq_id and insertion code of res, or (None, None) if gemmi
    can't read them."""

    try:
        seq_num = res.label_seq
        ins = res.seqid.icode
    except _MALFORMED_CIF_ERRORS:
        return None, None
    return seq_num, (ins if ins not in (" ", "\x00", "") else None)


def _atoms(
    model: gemmi.Model, subchain_to_entity: dict[str, str]
) -> tuple[list[AtomSiteRecord], dict[str, str]]:
    """Atom records of model, and each label_asym_id's auth_asym_id."""

    asym_auth_map: dict[str, str] = {}
    atoms: list[AtomSiteRecord] = []
    counter = 0

    for chain in model:
        for res in chain:
            sc = res.subchain
            if sc and sc not in asym_auth_map:
                asym_auth_map[sc] = chain.name

            is_polymer = res.entity_type == gemmi.EntityType.Polymer
            seq_num, ins_code = _residue_seq(res)

            for atom in res:
                counter += 1
                atoms.append(
                    AtomSiteRecord(
                        group_PDB="ATOM" if is_polymer else "HETATM",
                        id=atom.serial if atom.serial > 0 else counter,
                        type_symbol=atom.element.name,
                        label_atom_id=atom.name,
                        label_alt_id=(
                            atom.altloc
                            if atom.altloc not in ("\x00", " ", "")
                            else None
                        ),
                        label_comp_id=res.name,
                        label_asym_id=sc or "",
                        label_entity_id=subchain_to_entity.get(sc or "", ""),
                        label_seq_id=seq_num if is_polymer else None,
                        pdbx_PDB_ins_code=ins_code,
                        Cartn_x=atom.pos.x,
                        Cartn_y=atom.pos.y,
                        Cartn_z=atom.pos.z,
                        occupancy=atom.occ,
                        B_iso_or_equiv=atom.b_iso,
                        auth_seq_id=str(res.seqid.num),
                        auth_comp_id=res.name,
                        auth_asym_id=chain.name,
                        auth_atom_id=atom.name,
                        pdbx_PDB_model_num=model.num,
                    )
                )
    return atoms, asym_auth_map


def _asym_units(
    st: gemmi.Structure, asym_auth_map: dict[str, str]
) -> list[AsymRecord]:
    """One AsymRecord per subchain, in entity order."""

    seen: set[str] = set()
    asym_units: list[AsymRecord] = []
    for ent in st.entities:
        for sc in ent.subchains:
            if sc not in seen:
                seen.add(sc)
                asym_units.append(
                    AsymRecord(
                        id=sc,
                        entity_id=ent.name,
                        auth_id=asym_auth_map.get(sc),
                    )
                )
    return asym_units


# Connections ---------------------------------------------------------
_PARTNER_TAGS = [
    "label_asym_id",
    "label_comp_id",
    "label_seq_id",
    "label_atom_id",
    "auth_asym_id",
    "auth_seq_id",
    "PDB_ins_code",
    "symmetry",
]
_CONN_TAGS = [
    "id",
    "conn_type_id",
    *(f"ptnr1_{t}" for t in _PARTNER_TAGS),
    *(f"ptnr2_{t}" for t in _PARTNER_TAGS),
    "pdbx_dist_value",
    "details",
]


def _conn_partner(row: gemmi.cif.Table.Row, i: int) -> ConnPartner:
    """The bond partner whose eight columns start at row[i]."""

    return ConnPartner(
        label_asym_id=_cs(row[i]) or "",
        label_comp_id=_cs(row[i + 1]) or "",
        label_seq_id=_ci(row[i + 2]),
        label_atom_id=_cs(row[i + 3]),
        auth_asym_id=_cs(row[i + 4]),
        auth_seq_id=_cs(row[i + 5]),
        pdbx_PDB_ins_code=_cs(row[i + 6]),
        symmetry=_cs(row[i + 7]),
    )


def _conn_row(row: gemmi.cif.Table.Row) -> ConnRecord:
    """ConnRecord of a _struct_conn row."""

    return ConnRecord(
        id=_cs(row[0]) or "",
        conn_type_id=_cs(row[1]) or "covale",
        ptnr1=_conn_partner(row, 2),
        ptnr2=_conn_partner(row, 10),
        pdbx_dist_value=_cf(row[18]),
        details=_cs(row[19]),
    )


# Assemblies ----------------------------------------------------------
_ASSEMBLY_TAGS = [
    "id",
    "details",
    "method_details",
    "oligomeric_details",
    "oligomeric_count",
]


def _assembly_row(row: gemmi.cif.Table.Row) -> tuple[str, dict]:
    """(assembly id, its _pdbx_struct_assembly metadata)."""

    return _cs(row[0]) or "", {
        "details": _cs(row[1]),
        "method_details": _cs(row[2]),
        "oligomeric_details": _cs(row[3]),
        "oligomeric_count": _ci(row[4]),
    }


def _assembly(asm: gemmi.Assembly, meta: dict) -> AssemblyRecord:
    """AssemblyRecord of a gemmi assembly, each operator listed once."""

    gens: list[AssemblyGenRecord] = []
    operators: dict[str, AssemblyOperRecord] = {}
    for gen in asm.generators:
        gens.append(
            AssemblyGenRecord(
                assembly_id=asm.name,
                oper_expression=",".join(op.name for op in gen.operators),
                asym_id_list=list(gen.subchains),
            )
        )
        for op in gen.operators:
            if op.name not in operators:
                vec = op.transform.vec
                operators[op.name] = AssemblyOperRecord(
                    id=op.name,
                    type=getattr(op, "type", None) or None,
                    matrix=op.transform.mat.tolist(),
                    vector=[vec.x, vec.y, vec.z],
                )

    return AssemblyRecord(
        id=asm.name,
        details=meta.get("details"),
        method_details=meta.get("method_details"),
        oligomeric_details=meta.get("oligomeric_details"),
        oligomeric_count=meta.get("oligomeric_count"),
        generators=gens,
        operators=list(operators.values()),
        author_determined=asm.author_determined,
        software_determined=asm.software_determined,
    )


# Secondary structure -------------------------------------------------
_SS_RANGE_TAGS = [
    "beg_label_asym_id",
    "beg_label_seq_id",
    "end_label_asym_id",
    "end_label_seq_id",
    "beg_auth_asym_id",
    "beg_auth_seq_id",
    "end_auth_asym_id",
    "end_auth_seq_id",
]


def _ss_range(row: gemmi.cif.Table.Row) -> dict:
    """Begin/end fields shared by _struct_conf and _struct_sheet_range,
    from row[2:10]."""

    return {
        "beg_label_asym_id": _cs(row[2]) or "",
        "beg_label_seq_id": _ci(row[3]),
        "end_label_asym_id": _cs(row[4]) or "",
        "end_label_seq_id": _ci(row[5]),
        "beg_auth_asym_id": _cs(row[6]),
        "beg_auth_seq_id": _cs(row[7]),
        "end_auth_asym_id": _cs(row[8]),
        "end_auth_seq_id": _cs(row[9]),
    }


def _conf_row(row: gemmi.cif.Table.Row) -> ConfRecord:
    """ConfRecord of a _struct_conf row."""

    return ConfRecord(
        id=_cs(row[0]) or "", conf_type_id=_cs(row[1]) or "", **_ss_range(row)
    )


def _sheet_row(row: gemmi.cif.Table.Row) -> SheetStrandRecord:
    """SheetStrandRecord of a _struct_sheet_range row."""

    return SheetStrandRecord(
        sheet_id=_cs(row[0]) or "", id=_cs(row[1]) or "", **_ss_range(row)
    )


# mmcif_to_structure ----------------------------------
def mmcif_to_structure(
    path_to_mmcif: str,
    # entry_id: str = "",
    model_num: int = 1,
) -> tuple[Structure | None, DiagnosticBundle, ResultStatus]:
    """Convert raw mmCIF text to a Structure (mmCIF data model).

    Parses the mmCIF file with gemmi, selects the requested model
    (falling back to the first model with a warning if not found),
    and extracts entities, atoms, asym units, connections, assemblies,
    and secondary structure into Pandora's typed `Structure` schema.
    Any other mmCIF categories are preserved verbatim in `Structure.raw`.

    Args:
        path_to_mmcif: Path to the mmCIF file on disk.
        model_num: Which model number to select from a multi-model
            entry.

    Returns:
        A tuple of:
            - The parsed `Structure`, or `None` if parsing failed.
            - A `DiagnosticBundle` of warnings/errors encountered.
            - A `ResultStatus`: "success", "warning", or "failed".
    """
    diag = DiagnosticBundle()
    loaded = _load(path_to_mmcif, diag)
    if loaded is None:
        return None, diag, "failed"
    st, block = loaded
    eid = st.name
    model = _select_model(st, model_num, diag)
    read = partial(_read_rows, block, eid, diag)

    # Read in this order so warnings come out in category order.
    entity_desc = dict(
        read(
            "_entity.", _ENTITY_TAGS, _entity_row, "ENTITY_METADATA_PARSE_ERROR"
        )
    )
    entity_poly = dict(
        read(
            "_entity_poly.",
            _ENTITY_POLY_TAGS,
            _entity_poly_row,
            "ENTITY_POLY_PARSE_ERROR",
        )
    )

    subchain_to_entity = {
        sc: ent.name for ent in st.entities for sc in ent.subchains
    }
    atoms, asym_auth_map = _atoms(model, subchain_to_entity)
    if not atoms:
        diag.warnings.append(
            Diagnostic(
                code="MISSING_ATOM_SITE",
                severity="warning",
                message="No atoms found in structure",
                entry_id=eid,
            )
        )

    connections = read(
        "_struct_conn.", _CONN_TAGS, _conn_row, "CONN_PARSE_ERROR"
    )
    asm_meta = dict(
        read(
            "_pdbx_struct_assembly.",
            _ASSEMBLY_TAGS,
            _assembly_row,
            "ASSEMBLY_METADATA_PARSE_ERROR",
        )
    )
    conf_records = read(
        "_struct_conf.",
        ["id", "conf_type_id", *_SS_RANGE_TAGS],
        _conf_row,
        "STRUCT_CONF_PARSE_ERROR",
    )
    sheet_strands = read(
        "_struct_sheet_range.",
        ["sheet_id", "id", *_SS_RANGE_TAGS],
        _sheet_row,
        "SHEET_RANGE_PARSE_ERROR",
    )

    raw_title = block.find_value("_struct.title")
    return (
        Structure(
            entry_id=eid,
            entry=EntryRecord(
                id=eid, title=_cs(raw_title) if raw_title else None
            ),
            entities=_entities(st, entity_desc, entity_poly),
            asym_units=_asym_units(st, asym_auth_map),
            atoms=atoms,
            connections=connections,
            assemblies=[
                _assembly(asm, asm_meta.get(asm.name, {}))
                for asm in st.assemblies
            ],
            secondary_structure=SSRecord(
                conf_records=conf_records, sheet_strands=sheet_strands
            ),
            raw=_extract_raw(block),
        ),
        diag,
        "warning" if diag.warnings else "success",
    )
