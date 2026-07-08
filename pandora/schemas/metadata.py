from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from pandora.schemas.annotation import AnnotationLayer, AnnotationScope


MetadataSource = Literal["mmcif", "pdbe", "pdb", "uniprot", "sifts", "custom"]


class MetadataProvenance(BaseModel):
    source: MetadataSource | str = "mmcif"
    source_category: str | None = None
    source_record_id: str | None = None


class EntryMetadataRecord(BaseModel):
    entry_id: str
    title: str | None = None
    keywords: str | None = None
    citation_title: str | None = None
    doi: str | None = None
    pubmed_id: str | None = None
    provenance: MetadataProvenance = Field(default_factory=MetadataProvenance)


class TaxonomyRecord(BaseModel):
    entity_id: str | None = None
    ncbi_taxon_id: int | None = None
    organism_scientific: str | None = None
    organism_common: str | None = None
    host_ncbi_taxon_id: int | None = None
    host_scientific: str | None = None
    host_common: str | None = None
    expression_system: str | None = None
    provenance: MetadataProvenance = Field(default_factory=MetadataProvenance)


class QualityRecord(BaseModel):
    experimental_method: str | None = None
    resolution: float | None = None
    r_work: float | None = None
    r_free: float | None = None
    observed_reflections: int | None = None
    percent_possible_observed: float | None = None
    mean_b_factor: float | None = None
    provenance: MetadataProvenance = Field(default_factory=MetadataProvenance)


class EntityMetadataRecord(BaseModel):
    entity_id: str
    entity_type: str | None = None
    description: str | None = None
    formula_weight: float | None = None
    source_method: str | None = None
    ec_number: str | None = None
    mutation: str | None = None
    fragment: str | None = None
    polymer_type: str | None = None
    sequence: str | None = None
    canonical_sequence: str | None = None
    chain_ids: list[str] = Field(default_factory=list)
    provenance: MetadataProvenance = Field(default_factory=MetadataProvenance)


class LigandMetadataRecord(BaseModel):
    entity_id: str | None = None
    comp_id: str
    name: str | None = None
    formula: str | None = None
    formula_weight: float | None = None
    chain_ids: list[str] = Field(default_factory=list)
    provenance: MetadataProvenance = Field(default_factory=MetadataProvenance)


class UniProtMappingRecord(BaseModel):
    entity_id: str | None = None
    asym_id: str | None = None
    accession: str
    db_name: str | None = None
    db_code: str | None = None
    seq_id_start: int | None = None
    seq_id_end: int | None = None
    uniprot_start: int | None = None
    uniprot_end: int | None = None
    identity: float | None = None
    provenance: MetadataProvenance = Field(default_factory=MetadataProvenance)


class MetadataRecord(BaseModel):
    entry_id: str
    entry: EntryMetadataRecord | None = None
    quality: QualityRecord | None = None
    taxonomies: list[TaxonomyRecord] = Field(default_factory=list)
    entities: list[EntityMetadataRecord] = Field(default_factory=list)
    ligands: list[LigandMetadataRecord] = Field(default_factory=list)
    uniprot_mappings: list[UniProtMappingRecord] = Field(default_factory=list)
    raw_categories: list[str] = Field(default_factory=list)


__all__ = [
    "AnnotationLayer",
    "AnnotationScope",
    "EntityMetadataRecord",
    "EntryMetadataRecord",
    "LigandMetadataRecord",
    "MetadataProvenance",
    "MetadataRecord",
    "MetadataSource",
    "QualityRecord",
    "TaxonomyRecord",
    "UniProtMappingRecord",
]
