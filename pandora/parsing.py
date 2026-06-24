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
    gemmi.EntityType.Polymer:    "polymer",
    gemmi.EntityType.NonPolymer: "non-polymer",
    gemmi.EntityType.Water:      "water",
    gemmi.EntityType.Branched:   "branched",
}


# def parse_mmcif(
#     raw_content: str,
#     entry_id: str = "",
# ) -> tuple[ParsedStructure | None, DiagnosticBundle, str]:
#     """Parse raw mmCIF text into a ParsedStructure using gemmi."""
#     diag = DiagnosticBundle()

#     if not raw_content.strip():
#         diag.errors.append(Diagnostic(
#             code="EMPTY_CONTENT", 
#             severity="error",
#             message="mmCIF content is empty", 
#             entry_id=entry_id or None,
#         ))
#         return None, diag, "failed"

#     try:
#         st = gemmi.read_structure_string(raw_content, format=gemmi.CoorFormat.Mmcif)
#     except Exception as exc:
#         diag.errors.append(Diagnostic(
#             code="PARSE_ERROR", severity="error",
#             message=str(exc), entry_id=entry_id or None,
#         ))
#         return None, diag, "failed"

#     if len(st) == 0:
#         diag.errors.append(Diagnostic(
#             code="NO_MODEL", severity="error",
#             message="Structure contains no models", entry_id=entry_id or None,
#         ))
#         return None, diag, "failed"

#     model = st[0]
#     eid = entry_id or st.name

#     atoms_out: list[Atom] = []
#     residues_out: list[Residue] = []
#     chains_out: list[Chain] = []
#     ligands_out: list[Ligand] = []
#     _serial = 0

#     # Build a subchain → entity_id lookup once
#     subchain_to_entity: dict[str, str] = {}
#     for ent in st.entities:
#         for sc in ent.subchains:
#             subchain_to_entity[sc] = ent.name

#     for chain in model:
#         residues_in_chain: list[Residue] = []
#         chain_type: str = "non-polymer"
#         entity_id: str = ""

#         # Determine chain type + entity_id from the first subchain
#         for span in chain.subchains():
#             sc_name = span.subchain_id()
#             entity_id = subchain_to_entity.get(sc_name, "")
#             for ent in st.entities:
#                 if sc_name in ent.subchains:
#                     chain_type = _ENTITY_TYPE_MAP.get(ent.entity_type, "non-polymer")
#                     break
#             break  # use the first subchain only for chain-level metadata

#         for res in chain:
#             et = res.entity_type
#             is_polymer = et == gemmi.EntityType.Polymer
#             is_water = et == gemmi.EntityType.Water
#             rid = f"{eid}:{chain.name}:{res.seqid}:{res.name}"

#             atoms_in_res: list[Atom] = []
#             for atom in res:
#                 _serial += 1
#                 altloc = atom.altloc if atom.altloc not in ("\x00", " ", "") else None
#                 atoms_in_res.append(Atom(
#                     atom_id=f"{rid}:{atom.name}:{_serial}",
#                     atom_name=atom.name,
#                     element=atom.element.name,
#                     x=atom.pos.x,
#                     y=atom.pos.y,
#                     z=atom.pos.z,
#                     occupancy=atom.occ,
#                     b_factor=atom.b_iso,
#                     altloc=altloc,
#                     residue_id=rid,
#                     chain_id=chain.name,
#                 ))
#             atoms_out.extend(atoms_in_res)

#             try:
#                 seq_id_num: int | None = res.seqid.num
#                 ins = res.seqid.icode
#                 ins_code: str | None = ins if ins not in (" ", "\x00", "") else None
#             except Exception:
#                 seq_id_num, ins_code = None, None

#             r = Residue(
#                 residue_id=rid,
#                 comp_id=res.name,
#                 seq_id=seq_id_num,
#                 auth_seq_id=str(res.seqid),
#                 insertion_code=ins_code,
#                 chain_id=chain.name,
#                 atoms=atoms_in_res,
#                 is_polymer=is_polymer,
#             )
#             residues_in_chain.append(r)
#             residues_out.append(r)

#             if not is_polymer and not is_water:
#                 atom_names = [a.element.name.upper() for a in res]
#                 ligands_out.append(Ligand(
#                     ligand_id=rid,
#                     chem_comp_id=res.name,
#                     chain_id=chain.name,
#                     residue_id=rid,
#                     is_water=False,
#                     is_ion=(len(atom_names) == 1 and atom_names[0] in _METALS),
#                 ))

#         if not residues_in_chain:
#             diag.warnings.append(Diagnostic(
#                 code="EMPTY_CHAIN", severity="warning",
#                 message=f"Chain {chain.name!r} has no residues",
#                 entry_id=eid, context={"chain_id": chain.name},
#             ))

#         chains_out.append(Chain(
#             chain_id=chain.name,        # auth_asym_id; C02 normalises to label_asym_id
#             auth_chain_id=chain.name,
#             entity_id=entity_id,
#             chain_type=chain_type,      # type: ignore[arg-type]
#             residues=residues_in_chain,
#         ))

#     entities_out: list[Entity] = []
#     for ent in st.entities:
#         entities_out.append(Entity(
#             entity_id=ent.name,
#             entity_type=_ENTITY_TYPE_MAP.get(ent.entity_type, "non-polymer"),
#             description=None,
#             chain_ids=list(ent.subchains),
#             sequence=_entity_sequence(ent),
#         ))

#     assemblies_out: list[Assembly] = []
#     for asm in st.assemblies:
#         gens = [
#             AssemblyGen(
#                 asym_id_list=list(gen.chains),
#                 oper_expression=",".join(str(op) for op in gen.operators),
#             )
#             for gen in asm.generators
#         ]
#         assemblies_out.append(Assembly(
#             assembly_id=asm.name,
#             assembly_gen=gens,
#         ))

#     if not atoms_out:
#         diag.warnings.append(Diagnostic(
#             code="MISSING_ATOM_SITE", severity="warning",
#             message="No atoms found in structure", entry_id=eid,
#         ))

#     return (
#         ParsedStructure(
#             atoms=atoms_out,
#             residues=residues_out,
#             chains=chains_out,
#             entities=entities_out,
#             assemblies=assemblies_out,
#             ligands=ligands_out,
#         ),
#         diag,
#         "warning" if diag.warnings else "success",
#     )

# mmcif_to_structure helpers --------------------------------------------
_NULL_CIF  = frozenset({".", "?"})
_SKIP_RAW  = frozenset({"_atom_site", "_atom_site_anisou"})


def _cs(v: str) -> str | None:
    return None if v in _NULL_CIF else v


def _ci(v: str) -> int | None:
    if v in _NULL_CIF:
        return None
    try:
        return int(v)
    except (ValueError, TypeError):
        return None


def _cf(v: str) -> float | None:
    if v in _NULL_CIF:
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def _extract_raw(block: object) -> dict[str, list[dict[str, str | None]]]:
    raw: dict[str, list[dict[str, str | None]]] = {}
    for item in block:  # type: ignore[union-attr]
        loop = getattr(item, "loop", None)
        if loop is None or not loop.tags:
            continue
        category = loop.tags[0].rsplit(".", 1)[0]
        if category in _SKIP_RAW:
            continue
        col_names = [t.rsplit(".", 1)[1] for t in loop.tags]
        n_cols = len(col_names)
        raw[category] = [
            {col_names[j]: _cs(loop[i, j]) for j in range(n_cols)}
            for i in range(loop.length())
        ]
    return raw


# mmcif_to_structure ----------------------------------
def mmcif_to_structure(
    path_to_mmcif: str,
    # entry_id: str = "",
    model_num: int = 1,
) -> tuple[Structure | None, DiagnosticBundle, ResultStatus]:
    """Convert raw mmCIF text to a Structure (mmCIF data model)."""
    diag = DiagnosticBundle()

    if not path_to_mmcif.strip():
        diag.errors.append(Diagnostic(
            code="EMPTY_CONTENT", severity="error",
            message="mmCIF content is empty"
            # entry_id=entry_id or None,
        ))
        return None, diag, "failed"
    
    try:
        st = gemmi.read_structure(path_to_mmcif, format=gemmi.CoorFormat.Mmcif)
    except Exception as exc:
        diag.errors.append(Diagnostic(
            code="PARSE_ERROR", severity="error",
            message=str(exc)
            # entry_id=entry_id or None,
        ))
        return None, diag, "failed"

    try:
        doc   = gemmi.cif.read(path_to_mmcif)
        block = doc.sole_block()
    except Exception as exc:
        diag.errors.append(Diagnostic(
            code="CIF_PARSE_ERROR", severity="error",
            message=str(exc)
            # entry_id=entry_id or None,
        ))
        return None, diag, "failed"

    if len(st) == 0:
        diag.errors.append(Diagnostic(
            code="NO_MODEL", severity="error",
            message="Structure contains no models"
            # entry_id=entry_id or None,
        ))
        return None, diag, "failed"

    eid = st.name

    # Select model --------------------------
    model = next((m for m in st if m.num == model_num), None)
    if model is None:
        model = st[0]
        diag.warnings.append(Diagnostic(
            code="MODEL_NOT_FOUND", severity="warning",
            message=f"Model {model_num} not found; falling back to model {model.num}",
            entry_id=eid,
        ))

    # Entry ------------------------------------------
    raw_title = block.find_value("_struct.title")
    entry = EntryRecord(
        id=eid,
        title=_cs(raw_title) if raw_title else None,
    )

    # Subchain to entity lookup -----------------------
    subchain_to_entity: dict[str, str] = {
        sc: ent.name
        for ent in st.entities
        for sc in ent.subchains
    }

    # Entity metadata from cif ------------------------------
    entity_desc_map: dict[str, tuple[str | None, float | None, str | None]] = {}
    try:
        for row in block.find("_entity.", ["id", "pdbx_description", "formula_weight", "src_method"]):
            entity_desc_map[row[0]] = (_cs(row[1]), _cf(row[2]), _cs(row[3]))
    except Exception:
        pass

    entity_poly_map: dict[str, EntityPolyRecord] = {}
    try:
        for row in block.find("_entity_poly.", [
            "entity_id", "type",
            "pdbx_seq_one_letter_code", "pdbx_seq_one_letter_code_can",
            "pdbx_strand_id",
        ]):
            entity_poly_map[row[0]] = EntityPolyRecord(
                type=_cs(row[1]),
                pdbx_seq_one_letter_code=_cs(row[2]),
                pdbx_seq_one_letter_code_can=_cs(row[3]),
                pdbx_strand_id=_cs(row[4]),
            )
    except Exception:
        pass

    entities_out: list[EntityRecord] = []
    for ent in st.entities:
        desc, fw, src = entity_desc_map.get(ent.name, (None, None, None))
        entities_out.append(EntityRecord(
            id=ent.name,
            type=_ENTITY_TYPE_MAP.get(ent.entity_type, "non-polymer"),
            pdbx_description=desc,
            formula_weight=fw,
            src_method=src,
            poly=entity_poly_map.get(ent.name),
        ))

    # Atoms + asym tracking ----------------------------------------------
    asym_auth_map: dict[str, str] = {}   # label_asym_id → auth_asym_id
    atoms_out: list[AtomSiteRecord] = []
    _serial = 0

    for chain in model:
        for res in chain:
            sc = res.subchain
            if sc and sc not in asym_auth_map:
                asym_auth_map[sc] = chain.name

            is_polymer = res.entity_type == gemmi.EntityType.Polymer
            group = "ATOM" if is_polymer else "HETATM"

            try:
                seq_num: int | None = res.seqid.num
                ins = res.seqid.icode
                ins_code: str | None = ins if ins not in (" ", "\x00", "") else None
            except Exception:
                seq_num, ins_code = None, None

            for atom in res:
                _serial += 1
                serial = atom.serial if atom.serial > 0 else _serial
                altloc = atom.altloc if atom.altloc not in ("\x00", " ", "") else None
                atoms_out.append(AtomSiteRecord(
                    group_PDB=group,
                    id=serial,
                    type_symbol=atom.element.name,
                    label_atom_id=atom.name,
                    label_alt_id=altloc,
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
                    auth_seq_id=str(res.seqid),
                    auth_comp_id=res.name,
                    auth_asym_id=chain.name,
                    auth_atom_id=atom.name,
                    pdbx_PDB_model_num=model.num,
                ))

    if not atoms_out:
        diag.warnings.append(Diagnostic(
            code="MISSING_ATOM_SITE", severity="warning",
            message="No atoms found in structure", entry_id=eid,
        ))

    # Asym units ------------------------------------------------
    seen_asyms: set[str] = set()
    asym_units: list[AsymRecord] = []
    for ent in st.entities:
        for sc in ent.subchains:
            if sc not in seen_asyms:
                seen_asyms.add(sc)
                asym_units.append(AsymRecord(
                    id=sc,
                    entity_id=ent.name,
                    auth_id=asym_auth_map.get(sc),
                ))

    # Connections -------------------------------------------
    connections_out: list[ConnRecord] = []
    try:
        for row in block.find("_struct_conn.", [
            "id", "conn_type_id",
            "ptnr1_label_asym_id", "ptnr1_label_comp_id", "ptnr1_label_seq_id",
            "ptnr1_label_atom_id", "ptnr1_auth_asym_id", "ptnr1_auth_seq_id",
            "ptnr1_PDB_ins_code", "ptnr1_symmetry",
            "ptnr2_label_asym_id", "ptnr2_label_comp_id", "ptnr2_label_seq_id",
            "ptnr2_label_atom_id", "ptnr2_auth_asym_id", "ptnr2_auth_seq_id",
            "ptnr2_PDB_ins_code", "ptnr2_symmetry",
            "pdbx_dist_value", "details",
        ]):
            connections_out.append(ConnRecord(
                id=row[0],
                conn_type_id=_cs(row[1]) or "covale",
                ptnr1=ConnPartner(
                    label_asym_id=_cs(row[2]) or "",
                    label_comp_id=_cs(row[3]) or "",
                    label_seq_id=_ci(row[4]),
                    label_atom_id=_cs(row[5]),
                    auth_asym_id=_cs(row[6]),
                    auth_seq_id=_cs(row[7]),
                    pdbx_PDB_ins_code=_cs(row[8]),
                    symmetry=_cs(row[9]),
                ),
                ptnr2=ConnPartner(
                    label_asym_id=_cs(row[10]) or "",
                    label_comp_id=_cs(row[11]) or "",
                    label_seq_id=_ci(row[12]),
                    label_atom_id=_cs(row[13]),
                    auth_asym_id=_cs(row[14]),
                    auth_seq_id=_cs(row[15]),
                    pdbx_PDB_ins_code=_cs(row[16]),
                    symmetry=_cs(row[17]),
                ),
                pdbx_dist_value=_cf(row[18]),
                details=_cs(row[19]),
            ))
    except Exception:
        pass

    # Assemblies ------------------------------------
    asm_meta: dict[str, dict] = {}
    try:
        for row in block.find("_pdbx_struct_assembly.", [
            "id", "details", "method_details", "oligomeric_details", "oligomeric_count",
        ]):
            asm_meta[row[0]] = {
                "details":           _cs(row[1]),
                "method_details":    _cs(row[2]),
                "oligomeric_details": _cs(row[3]),
                "oligomeric_count":  _ci(row[4]),
            }
    except Exception:
        pass

    assemblies_out: list[AssemblyRecord] = []
    for asm in st.assemblies:
        meta = asm_meta.get(asm.name, {})
        gens: list[AssemblyGenRecord] = []
        operators: list[AssemblyOperRecord] = []
        seen_ops: set[str] = set()

        for gen in asm.generators:
            gens.append(AssemblyGenRecord(
                assembly_id=asm.name,
                oper_expression=",".join(op.name for op in gen.operators),
                asym_id_list=list(gen.chains),
            ))
            for op in gen.operators:
                if op.name not in seen_ops:
                    seen_ops.add(op.name)
                    mat = op.transform.mat
                    vec = op.transform.vec
                    operators.append(AssemblyOperRecord(
                        id=op.name,
                        type=getattr(op, "type", None) or None,
                        matrix=mat.tolist(),
                        vector=[vec.x, vec.y, vec.z],
                    ))

        assemblies_out.append(AssemblyRecord(
            id=asm.name,
            details=meta.get("details"),
            method_details=meta.get("method_details"),
            oligomeric_details=meta.get("oligomeric_details"),
            oligomeric_count=meta.get("oligomeric_count"),
            generators=gens,
            operators=operators,
        ))

    # Secondary structure --------------------------------
    conf_records: list[ConfRecord] = []
    try:
        for row in block.find("_struct_conf.", [
            "id", "conf_type_id",
            "beg_label_asym_id", "beg_label_seq_id",
            "end_label_asym_id", "end_label_seq_id",
            "beg_auth_asym_id",  "beg_auth_seq_id",
            "end_auth_asym_id",  "end_auth_seq_id",
        ]):
            conf_records.append(ConfRecord(
                id=row[0],
                conf_type_id=_cs(row[1]) or "",
                beg_label_asym_id=_cs(row[2]) or "",
                beg_label_seq_id=_ci(row[3]),
                end_label_asym_id=_cs(row[4]) or "",
                end_label_seq_id=_ci(row[5]),
                beg_auth_asym_id=_cs(row[6]),
                beg_auth_seq_id=_cs(row[7]),
                end_auth_asym_id=_cs(row[8]),
                end_auth_seq_id=_cs(row[9]),
            ))
    except Exception:
        pass

    sheet_strands: list[SheetStrandRecord] = []
    try:
        for row in block.find("_struct_sheet_range.", [
            "sheet_id", "id",
            "beg_label_asym_id", "beg_label_seq_id",
            "end_label_asym_id", "end_label_seq_id",
            "beg_auth_asym_id",  "beg_auth_seq_id",
            "end_auth_asym_id",  "end_auth_seq_id",
        ]):
            sheet_strands.append(SheetStrandRecord(
                sheet_id=row[0],
                id=row[1],
                beg_label_asym_id=_cs(row[2]) or "",
                beg_label_seq_id=_ci(row[3]),
                end_label_asym_id=_cs(row[4]) or "",
                end_label_seq_id=_ci(row[5]),
                beg_auth_asym_id=_cs(row[6]),
                beg_auth_seq_id=_cs(row[7]),
                end_auth_asym_id=_cs(row[8]),
                end_auth_seq_id=_cs(row[9]),
            ))
    except Exception:
        pass

    return (
        Structure(
            entry_id=eid,
            entry=entry,
            entities=entities_out,
            asym_units=asym_units,
            atoms=atoms_out,
            connections=connections_out,
            assemblies=assemblies_out,
            secondary_structure=SSRecord(conf_records=conf_records, sheet_strands=sheet_strands),
            raw=_extract_raw(block),
        ),
        diag,
        "warning" if diag.warnings else "success",
    )
