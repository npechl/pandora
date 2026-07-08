from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from .common import DiagnosticBundle, ResultStatus


# ── Core structural types ────────────────────────────────────────────────────


class Atom(BaseModel):
    atom_id: str
    atom_name: str
    element: str
    x: float
    y: float
    z: float
    occupancy: float
    b_factor: float
    altloc: str | None = None
    residue_id: str
    chain_id: str
    model_num: int = 1


class Residue(BaseModel):
    residue_id: str
    comp_id: str
    seq_id: int | None = None
    auth_seq_id: str
    insertion_code: str | None = None
    chain_id: str
    atoms: list[Atom] = Field(default_factory=list)
    is_polymer: bool


class Chain(BaseModel):
    chain_id: str
    auth_chain_id: str
    entity_id: str
    chain_type: Literal["polymer", "non-polymer", "water", "branched"]
    residues: list[Residue] = Field(default_factory=list)


class Entity(BaseModel):
    entity_id: str
    entity_type: Literal["polymer", "non-polymer", "water", "branched"]
    description: str | None = None
    chain_ids: list[str] = Field(default_factory=list)
    sequence: str | None = None


class AssemblyGen(BaseModel):
    asym_id_list: list[str]
    oper_expression: str


class Assembly(BaseModel):
    assembly_id: str
    details: str | None = None
    assembly_gen: list[AssemblyGen] = Field(default_factory=list)


class Ligand(BaseModel):
    ligand_id: str
    chem_comp_id: str
    chain_id: str
    residue_id: str
    is_water: bool
    is_ion: bool


class ParsedStructure(BaseModel):
    atoms: list[Atom] = Field(default_factory=list)
    residues: list[Residue] = Field(default_factory=list)
    chains: list[Chain] = Field(default_factory=list)
    entities: list[Entity] = Field(default_factory=list)
    assemblies: list[Assembly] = Field(default_factory=list)
    ligands: list[Ligand] = Field(default_factory=list)


# ── Input schemas ────────────────────────────────────────────────────────────

ProviderType = Literal["pdbe", "pdb", "local", "custom"]
StaleBehavior = Literal["use_stale", "warn", "fail"]


class FetchOptions(BaseModel):
    allow_partial: bool = False
    use_cache: bool = True
    decompress: bool = True
    max_age_seconds: int | None = None
    stale_behavior: StaleBehavior = "use_stale"


class ParallelOptions(BaseModel):
    max_workers: int | None = None
    fail_fast: bool = False


class MmCIFIngestionInput(BaseModel):
    entry_id: str
    provider: ProviderType
    source_uri: str | None = None
    raw_content: str | None = None
    fetch_options: FetchOptions = Field(default_factory=FetchOptions)


class MmCIFBatchInput(BaseModel):
    entries: list[MmCIFIngestionInput]
    mode: Literal["sequential", "parallel"] = "sequential"
    fetch_options: FetchOptions = Field(default_factory=FetchOptions)
    parallel_options: ParallelOptions = Field(default_factory=ParallelOptions)


# ── Output schemas ───────────────────────────────────────────────────────────


class IngestionProvenance(BaseModel):
    provider: str
    source_uri: str | None = None
    retrieved_at: str | None = None
    from_cache: bool = False


class MmCIFIngestionResult(BaseModel):
    entry_id: str
    status: ResultStatus
    parsed_structure: ParsedStructure | None = None
    diagnostics: DiagnosticBundle = Field(default_factory=DiagnosticBundle)
    provenance: IngestionProvenance


class MmCIFBatchSummary(BaseModel):
    total: int
    success: int
    warning: int
    failed: int


class MmCIFBatchResult(BaseModel):
    mode: Literal["sequential", "parallel"]
    summary: MmCIFBatchSummary
    results: list[MmCIFIngestionResult]
