from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

EntityType = Literal["polymer", "non-polymer", "water", "branched"]
ConnType = Literal["disulf", "covale", "hydrog", "metalc"]


class EntryRecord(BaseModel):
    """Entry-level identity and title, from `_struct`.

    Attributes:
        id: The entry's PDB/PDBe identifier (the mmCIF data block name).
        title: The entry's title (`_struct.title`), if present.
    """

    id: str
    title: str | None = None


class EntityPolyRecord(BaseModel):
    """Polymer-specific entity fields, from `_entity_poly`.

    Attributes:
        type: The polymer type (e.g. "polypeptide(L)",
            "polyribonucleotide"), from `_entity_poly.type`.
        pdbx_seq_one_letter_code: The one-letter sequence as deposited,
            including any non-standard residue codes in parentheses.
        pdbx_seq_one_letter_code_can: The canonicalized one-letter
            sequence (non-standard residues mapped to their closest
            standard equivalent).
        pdbx_strand_id: The comma-separated auth chain ids this entity
            spans.
    """

    type: str | None = None
    pdbx_seq_one_letter_code: str | None = None
    pdbx_seq_one_letter_code_can: str | None = None
    pdbx_strand_id: str | None = None


class EntityRecord(BaseModel):
    """One entity (molecule) in the structure, from `_entity`.

    Attributes:
        id: The entity id (`label_entity_id`).
        type: The entity type: "polymer", "non-polymer", "water", or
            "branched".
        pdbx_description: A free-text description of the entity (e.g.
            protein name).
        formula_weight: The entity's molecular formula weight, in
            daltons.
        src_method: How the entity was produced (e.g. "man" for
            genetically manipulated, "nat" for natural, "syn" for
            synthetic).
        poly: Polymer-specific fields, if this entity is a polymer.
    """

    id: str
    type: EntityType | str
    pdbx_description: str | None = None
    formula_weight: float | None = None
    src_method: str | None = None
    poly: EntityPolyRecord | None = None


class AsymRecord(BaseModel):
    """One asym unit (`label_asym_id`), mapping a chain to its entity.

    Attributes:
        id: The `label_asym_id`.
        entity_id: The id of the entity this asym unit belongs to.
        auth_id: The `auth_asym_id` (author-assigned chain id), if
            different from `id`.
    """

    id: str  # label_asym_id
    entity_id: str
    auth_id: str | None = None  # auth_asym_id


class AtomSiteRecord(BaseModel):
    """One row from `_atom_site`: a single atom's identity, position, and
    metadata.

    Attributes:
        group_PDB: "ATOM" for polymer residues, "HETATM" for everything
            else.
        id: The atom's serial number.
        type_symbol: The element symbol.
        label_atom_id: The atom's name within its residue (e.g. "CA").
        label_alt_id: The alternate location (altloc) identifier, if
            the atom has alternate conformers.
        label_comp_id: The residue/component id (e.g. "ALA").
        label_asym_id: The asym unit (chain) this atom belongs to.
        label_entity_id: The entity this atom belongs to.
        label_seq_id: The residue's sequence position within its
            entity; `None` for non-polymer atoms.
        pdbx_PDB_ins_code: The insertion code, if any.
        Cartn_x: The atom's x coordinate, in angstroms.
        Cartn_y: The atom's y coordinate, in angstroms.
        Cartn_z: The atom's z coordinate, in angstroms.
        occupancy: The atom's occupancy (fraction of the site it's
            present at).
        B_iso_or_equiv: The atom's isotropic (or equivalent) B-factor.
        auth_seq_id: The author-assigned residue sequence number.
        auth_comp_id: The author-assigned residue/component id.
        auth_asym_id: The author-assigned chain id.
        auth_atom_id: The author-assigned atom name.
        pdbx_PDB_model_num: Which model this atom belongs to, for
            multi-model entries.
    """

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
    """One side of a `_struct_conn` connection (bond/contact).

    Attributes:
        label_asym_id: The chain of this side of the connection.
        label_comp_id: The residue/component id of this side of the
            connection.
        label_seq_id: The residue sequence position of this side of
            the connection, if polymeric.
        label_atom_id: The atom name of this side of the connection.
        auth_asym_id: The author-assigned chain id of this side of the
            connection.
        auth_seq_id: The author-assigned residue sequence number of
            this side of the connection.
        pdbx_PDB_ins_code: The insertion code of this side of the
            connection, if any.
        symmetry: The symmetry operation applied to this side of the
            connection, if any.
    """

    label_asym_id: str
    label_comp_id: str
    label_seq_id: int | None = None
    label_atom_id: str | None = None
    auth_asym_id: str | None = None
    auth_seq_id: str | None = None
    pdbx_PDB_ins_code: str | None = None
    symmetry: str | None = None


class ConnRecord(BaseModel):
    """One row from `_struct_conn`: a bond or non-bonded connection
    between two atoms.

    Attributes:
        id: The connection's id.
        conn_type_id: The connection type: "disulf", "covale",
            "hydrog", or "metalc".
        ptnr1: The first atom in the connection.
        ptnr2: The second atom in the connection.
        pdbx_dist_value: The measured distance between the two atoms,
            in angstroms.
        details: Free-text details about the connection.
    """

    id: str
    conn_type_id: ConnType | str
    ptnr1: ConnPartner
    ptnr2: ConnPartner
    pdbx_dist_value: float | None = None
    details: str | None = None


class AssemblyOperRecord(BaseModel):
    """One symmetry operator (rotation matrix + translation vector) used
    to generate an assembly.

    Attributes:
        id: The operator's id.
        type: The operator type (e.g. "identity operation", "crystal
            symmetry operation").
        matrix: The 3x3 rotation matrix.
        vector: The 3-element translation vector.
    """

    id: str
    type: str | None = None
    matrix: list[list[float]] | None = None  # 3×3 rotation
    vector: list[float] | None = None  # translation


class AssemblyGenRecord(BaseModel):
    """One assembly-generation instruction: which operators apply to
    which asym units.

    Attributes:
        assembly_id: The assembly this generator instruction belongs
            to.
        oper_expression: The comma-separated operator ids to apply.
        asym_id_list: The asym units (chains) the operators are
            applied to.
    """

    assembly_id: str
    oper_expression: str
    asym_id_list: list[str]


class AssemblyRecord(BaseModel):
    """One biological/crystallographic assembly, from
    `_pdbx_struct_assembly` and its generators/operators.

    Attributes:
        id: The assembly's id.
        details: Free-text details about the assembly.
        method_details: The method used to determine the assembly.
        oligomeric_details: A description of the oligomeric state
            (e.g. "dimeric").
        oligomeric_count: The number of polymer chains in the
            assembly.
        generators: The generator instructions (which operators apply
            to which chains) that build this assembly.
        operators: The symmetry operators referenced by `generators`.
    """

    id: str
    details: str | None = None
    method_details: str | None = None
    oligomeric_details: str | None = None
    oligomeric_count: int | None = None
    generators: list[AssemblyGenRecord] = Field(default_factory=list)
    operators: list[AssemblyOperRecord] = Field(default_factory=list)


class ConfRecord(BaseModel):
    """One row from _struct_conf (HELX_P, TURN_P, STRN, etc.).

    Attributes:
        id: The secondary structure element's id.
        conf_type_id: The secondary structure type (e.g. "HELX_P",
            "TURN_P", "STRN").
        beg_label_asym_id: The chain the element starts in.
        beg_label_seq_id: The residue sequence position the element
            starts at.
        end_label_asym_id: The chain the element ends in.
        end_label_seq_id: The residue sequence position the element
            ends at.
        beg_auth_asym_id: The author-assigned chain id the element
            starts in.
        beg_auth_seq_id: The author-assigned residue number the
            element starts at.
        end_auth_asym_id: The author-assigned chain id the element
            ends in.
        end_auth_seq_id: The author-assigned residue number the
            element ends at.
    """

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
    """One row from `_struct_sheet_range`: one strand of a beta sheet.

    Attributes:
        sheet_id: The beta sheet this strand belongs to.
        id: The strand's id.
        beg_label_asym_id: The chain the strand starts in.
        beg_label_seq_id: The residue sequence position the strand
            starts at.
        end_label_asym_id: The chain the strand ends in.
        end_label_seq_id: The residue sequence position the strand
            ends at.
        beg_auth_asym_id: The author-assigned chain id the strand
            starts in.
        beg_auth_seq_id: The author-assigned residue number the strand
            starts at.
        end_auth_asym_id: The author-assigned chain id the strand ends
            in.
        end_auth_seq_id: The author-assigned residue number the strand
            ends at.
    """

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
    """A structure's secondary structure: helix/turn/strand records and
    sheet strand ranges.

    Attributes:
        conf_records: The helix/turn/strand elements from
            `_struct_conf`.
        sheet_strands: The individual beta-sheet strand ranges from
            `_struct_sheet_range`.
    """

    conf_records: list[ConfRecord] = Field(default_factory=list)
    sheet_strands: list[SheetStrandRecord] = Field(default_factory=list)


class Structure(BaseModel):
    """Pandora's typed mmCIF data model: one parsed structure, with
    everything not promoted to a typed field kept verbatim in `raw`.

    Attributes:
        entry_id: The structure's PDB/PDBe entry identifier.
        entry: Entry-level identity and title.
        entities: Every entity (molecule) in the structure.
        asym_units: Every asym unit (chain), mapping chains to
            entities.
        atoms: Every atom in the structure.
        connections: Every bond/contact connection between atoms.
        assemblies: Every biological/crystallographic assembly.
        secondary_structure: The structure's secondary structure
            elements.
        raw: Every other mmCIF category not promoted to a typed field,
            keyed by category name, preserved verbatim.
    """

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
