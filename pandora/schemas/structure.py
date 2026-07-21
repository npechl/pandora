from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

EntityType = Literal["polymer", "non-polymer", "water", "branched"]
ConnType = Literal["disulf", "covale", "hydrog", "metalc"]


class EntryRecord(BaseModel):
    id: str
    title: str | None = None


class EntityPolyRecord(BaseModel):
    type: str | None = None
    pdbx_seq_one_letter_code: str | None = None
    pdbx_seq_one_letter_code_can: str | None = None
    pdbx_strand_id: str | None = None


class EntityRecord(BaseModel):
    id: str
    type: EntityType | str
    pdbx_description: str | None = None
    formula_weight: float | None = None
    src_method: str | None = None
    poly: EntityPolyRecord | None = None


class AsymRecord(BaseModel):
    id: str  # label_asym_id
    entity_id: str
    auth_id: str | None = None  # auth_asym_id


class AtomSiteRecord(BaseModel):
    group_PDB: Literal["ATOM", "HETATM"]
    id: int
    type_symbol: str
    label_atom_id: str
    label_alt_id: str | None = None
    label_comp_id: str
    label_asym_id: str
    label_entity_id: str
    label_seq_id: int | None = None
    pdbx_PDB_ins_code: str | None = None
    Cartn_x: float
    Cartn_y: float
    Cartn_z: float
    occupancy: float
    B_iso_or_equiv: float
    auth_seq_id: str | None = None
    auth_comp_id: str | None = None
    auth_asym_id: str | None = None
    auth_atom_id: str | None = None
    pdbx_PDB_model_num: int = 1


class ConnPartner(BaseModel):
    label_asym_id: str
    label_comp_id: str
    label_seq_id: int | None = None
    label_atom_id: str | None = None
    auth_asym_id: str | None = None
    auth_seq_id: str | None = None
    pdbx_PDB_ins_code: str | None = None
    symmetry: str | None = None


class ConnRecord(BaseModel):
    id: str
    conn_type_id: ConnType | str
    ptnr1: ConnPartner
    ptnr2: ConnPartner
    pdbx_dist_value: float | None = None
    details: str | None = None


class AssemblyOperRecord(BaseModel):
    id: str
    type: str | None = None
    matrix: list[list[float]] | None = None  # 3×3 rotation
    vector: list[float] | None = None  # translation


class AssemblyGenRecord(BaseModel):
    assembly_id: str
    oper_expression: str
    asym_id_list: list[str]


class AssemblyRecord(BaseModel):
    id: str
    details: str | None = None
    method_details: str | None = None
    oligomeric_details: str | None = None
    oligomeric_count: int | None = None
    generators: list[AssemblyGenRecord] = Field(default_factory=list)
    operators: list[AssemblyOperRecord] = Field(default_factory=list)


class ConfRecord(BaseModel):
    """One row from _struct_conf (HELX_P, TURN_P, STRN, etc.)."""

    id: str
    conf_type_id: str
    beg_label_asym_id: str
    beg_label_seq_id: int | None = None
    end_label_asym_id: str
    end_label_seq_id: int | None = None
    beg_auth_asym_id: str | None = None
    beg_auth_seq_id: str | None = None
    end_auth_asym_id: str | None = None
    end_auth_seq_id: str | None = None


class SheetStrandRecord(BaseModel):
    sheet_id: str
    id: str
    beg_label_asym_id: str
    beg_label_seq_id: int | None = None
    end_label_asym_id: str
    end_label_seq_id: int | None = None
    beg_auth_asym_id: str | None = None
    beg_auth_seq_id: str | None = None
    end_auth_asym_id: str | None = None
    end_auth_seq_id: str | None = None


class SSRecord(BaseModel):
    conf_records: list[ConfRecord] = Field(default_factory=list)
    sheet_strands: list[SheetStrandRecord] = Field(default_factory=list)


class Structure(BaseModel):
    entry_id: str
    entry: EntryRecord
    entities: list[EntityRecord] = Field(default_factory=list)
    asym_units: list[AsymRecord] = Field(default_factory=list)
    atoms: list[AtomSiteRecord] = Field(default_factory=list)
    connections: list[ConnRecord] = Field(default_factory=list)
    assemblies: list[AssemblyRecord] = Field(default_factory=list)
    secondary_structure: SSRecord = Field(default_factory=SSRecord)
    # all non-atom_site cif categories not explicitly modelled above
    raw: dict[str, list[dict[str, str | None]]] = Field(default_factory=dict)
