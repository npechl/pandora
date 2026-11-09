from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

AnnotationScope = Literal[
    "entry",
    "entity",
    "chain",
    "residue",
    "ligand",
    "interface",
    "pair",
    "cluster",
    "dataset",
]

ContactAtomSet = Literal["all", "heavy", "backbone", "ca"]
"""Atoms that count towards a chain-chain contact: every atom, heavy
(non-H/D) atoms, backbone atoms (protein N/CA/C/O, nucleic-acid
O5'/C5'/C4'/C3'/O3'), or CA atoms only."""


class AnnotationLayer(BaseModel):
    """One derived annotation result (counts, contacts, identity, etc.)
    for an entry, pair, or other scope.

    Attributes:
        layer_name: A human-readable name for the layer.
        layer_type: A machine-readable type identifying the layer
            (e.g. "pairwise_sequence_identity").
        scope: What the layer describes: an entry, a pair, an
            interface, etc.
        method: The identifier of the method/algorithm used to
            compute the layer.
        target_ids: The entry id(s) the layer applies to.
        parameters: The parameters the layer was computed with.
        data: The layer's actual computed output.
        provenance: Additional provenance details about how the layer
            was computed.
    """

    layer_name: str
    layer_type: str
    scope: AnnotationScope
    method: str
    target_ids: list[str] = Field(default_factory=list)
    parameters: dict[str, Any] = Field(default_factory=dict)
    data: dict[str, Any] = Field(default_factory=dict)
    provenance: dict[str, Any] = Field(default_factory=dict)


LayerTemplate = Literal[
    "residue_contacts",
    "residue_values",
    "residue_labels",
    "chain_values",
    "chain_labels",
    "interface_labels",
    "free_form",
]
"""The table shapes an imported annotation layer can take."""

UnmatchedReason = Literal[
    "entry_not_found",
    "missing_key",
    "chain_not_found",
    "residue_not_found",
    "no_coordinates",
    "naming_conflict",
    "ambiguous_author_id",
    "bad_value",
    "duplicate_key",
]
"""Why an imported row could not be linked to a structure."""


class LayerDescription(BaseModel):
    """What an imported annotation layer is and where it came from.

    Written by the user once per layer, as YAML.

    Attributes:
        name: The layer's name, e.g. "cocada".
        template: The table shape the layer's rows follow.
        tool: The tool that produced the data.
        tool_version: The tool's version, if known.
        parameters: The tool's parameters, recorded as given.
        description: Free text about the layer.
    """

    model_config = ConfigDict(extra="forbid")

    name: str
    template: LayerTemplate
    tool: str
    tool_version: str | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)
    description: str = ""


class LayerRow(BaseModel):
    """One row of an imported layer table: any columns, kept as given."""

    model_config = ConfigDict(extra="allow")


_NO_INS_CODE = {"", "?", "."}


class ResidueContactRow(BaseModel):
    """One row of a `residue_contacts` layer: two residues in contact.

    Built-in contact annotations return these, and `annotate import`
    stores them, so any consumer works on either. Each side names its
    residue in both Pandora's naming (`chain_id_s` = `label_asym_id`,
    `seq_id_s` = `label_seq_id`, None for non-polymer residues) and the
    author's (`auth_chain_id_s`, `auth_seq_id_s`, `ins_code_s`).

    Attributes:
        entry_id: The structure's entry id.
        chain_id_1: Side 1 chain (`label_asym_id`).
        seq_id_1: Side 1 residue (`label_seq_id`), None if non-polymer.
        auth_chain_id_1: Side 1 author chain id.
        auth_seq_id_1: Side 1 author residue number.
        ins_code_1: Side 1 insertion code, None for none.
        comp_id_1: Side 1 residue name (e.g. "ALA").
        chain_id_2: Side 2 chain (`label_asym_id`).
        seq_id_2: Side 2 residue (`label_seq_id`), None if non-polymer.
        auth_chain_id_2: Side 2 author chain id.
        auth_seq_id_2: Side 2 author residue number.
        ins_code_2: Side 2 insertion code, None for none.
        comp_id_2: Side 2 residue name.
        contact_type: The kind of contact, if the tool reports one.
        distance: The contact distance in angstroms, if known.
    """

    model_config = ConfigDict(extra="allow")

    entry_id: str
    chain_id_1: str
    seq_id_1: int | None = None
    auth_chain_id_1: str | None = None
    auth_seq_id_1: str | None = None
    ins_code_1: str | None = None
    comp_id_1: str | None = None
    chain_id_2: str
    seq_id_2: int | None = None
    auth_chain_id_2: str | None = None
    auth_seq_id_2: str | None = None
    ins_code_2: str | None = None
    comp_id_2: str | None = None
    contact_type: str | None = None
    distance: float | None = None

    @field_validator("seq_id_1", "seq_id_2", mode="before")
    @classmethod
    def _whole_float_seq_id(cls, value: Any) -> Any:
        """Parquet stores a nullable int column as float (42.0)."""

        if isinstance(value, float) and value.is_integer():
            return int(value)
        return value

    @field_validator("ins_code_1", "ins_code_2", mode="before")
    @classmethod
    def _no_ins_code(cls, value: Any) -> Any:
        """ "", "?" and "." all mean no insertion code."""

        if isinstance(value, str) and value.strip() in _NO_INS_CODE:
            return None
        return value
