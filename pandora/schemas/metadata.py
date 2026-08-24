from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field

from pandora.schemas.annotation import AnnotationLayer, AnnotationScope

MetadataSource = Literal["mmcif", "pdbe", "pdb", "uniprot", "sifts", "custom"]


class MetadataProvenance(BaseModel):
    """Record of which source (and, for mmCIF, which category/record)
    one metadata value came from.

    Attributes:
        source: Which source database the value came from.
        source_category: The mmCIF category the value was read from,
            if `source="mmcif"`.
        source_record_id: The id of the specific record within
            `source_category`, if applicable.
    """

    source: MetadataSource | str = "mmcif"
    source_category: str | None = None
    source_record_id: str | None = None


class EntryMetadataRecord(BaseModel):
    """Entry-level source-backed metadata: title, keywords, and
    citation details.

    Attributes:
        entry_id: The structure's entry id.
        title: The entry's title.
        keywords: The entry's keywords.
        citation_title: The primary citation's title.
        doi: The primary citation's DOI.
        pubmed_id: The primary citation's PubMed id.
        provenance: Where this record's data came from.
    """

    entry_id: str
    title: str | None = None
    keywords: str | None = None
    citation_title: str | None = None
    doi: str | None = None
    pubmed_id: str | None = None
    provenance: MetadataProvenance = Field(default_factory=MetadataProvenance)


class TaxonomyRecord(BaseModel):
    """Source organism and expression host metadata for one entity.

    Attributes:
        entity_id: The entity this taxonomy applies to.
        ncbi_taxon_id: The source organism's NCBI taxonomy id.
        organism_scientific: The source organism's scientific name.
        organism_common: The source organism's common name.
        host_ncbi_taxon_id: The expression host's NCBI taxonomy id, if
            genetically manipulated.
        host_scientific: The expression host's scientific name.
        host_common: The expression host's common name.
        expression_system: The expression system or plasmid used.
        provenance: Where this record's data came from.
    """

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
    """Experimental quality metrics (resolution, R-factors, reflection
    counts) for one entry.

    Attributes:
        experimental_method: The experimental method(s) used (e.g.
            "X-RAY DIFFRACTION").
        resolution: The structure's resolution, in angstroms.
        r_work: The working set R-factor.
        r_free: The free set R-factor.
        observed_reflections: The number of observed reflections.
        percent_possible_observed: The percentage of possible
            reflections that were observed.
        mean_b_factor: The mean isotropic B-factor.
        provenance: Where this record's data came from.
    """

    experimental_method: str | None = None
    resolution: float | None = None
    r_work: float | None = None
    r_free: float | None = None
    observed_reflections: int | None = None
    percent_possible_observed: float | None = None
    mean_b_factor: float | None = None
    provenance: MetadataProvenance = Field(default_factory=MetadataProvenance)


class EntityMetadataRecord(BaseModel):
    """Descriptions, sequence, and chain membership for one entity.

    Attributes:
        entity_id: The entity's id.
        entity_type: The entity type (e.g. "polymer").
        description: A free-text description of the entity.
        formula_weight: The entity's molecular formula weight, in
            daltons.
        source_method: How the entity was produced.
        ec_number: The entity's Enzyme Commission number, if
            applicable.
        mutation: A description of any mutation(s), if applicable.
        fragment: A description of the fragment represented, if
            applicable.
        polymer_type: The polymer type, if this entity is a polymer.
        sequence: The entity's deposited one-letter sequence.
        canonical_sequence: The entity's canonicalized one-letter
            sequence.
        chain_ids: The chain (asym) ids belonging to this entity.
        provenance: Where this record's data came from.
    """

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
    """Name, formula, and chain membership for one non-polymer ligand.

    Attributes:
        entity_id: The id of the entity this ligand belongs to, if
            known.
        comp_id: The ligand's component id (e.g. "ZN").
        name: The ligand's full chemical name.
        formula: The ligand's chemical formula.
        formula_weight: The ligand's molecular formula weight, in
            daltons.
        chain_ids: The chain (asym) ids this ligand appears in.
        provenance: Where this record's data came from.
    """

    entity_id: str | None = None
    comp_id: str
    name: str | None = None
    formula: str | None = None
    formula_weight: float | None = None
    chain_ids: list[str] = Field(default_factory=list)
    provenance: MetadataProvenance = Field(default_factory=MetadataProvenance)


class UniProtMappingRecord(BaseModel):
    """One UniProt/SIFTS sequence mapping segment for an entity or chain.

    Attributes:
        entity_id: The entity this mapping applies to, if known.
        asym_id: The chain this mapping applies to, if known.
        accession: The UniProt accession.
        db_name: The reference database name.
        db_code: The reference database entry code.
        seq_id_start: The mapping's start position in the structure's
            sequence numbering.
        seq_id_end: The mapping's end position in the structure's
            sequence numbering.
        uniprot_start: The mapping's start position in the UniProt
            sequence.
        uniprot_end: The mapping's end position in the UniProt
            sequence.
        identity: The sequence identity of the mapped segment, if
            reported (SIFTS mappings only).
        provenance: Where this record's data came from.
    """

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
    """Every source-backed metadata collected for one structure: entry,
    quality, taxonomy, entities, ligands, and UniProt mappings.

    Attributes:
        entry_id: The structure's entry id.
        entry: Entry-level metadata (title, keywords, citation).
        quality: Experimental quality metrics.
        taxonomies: Source organism and expression host metadata.
        entities: Per-entity metadata.
        ligands: Per-ligand metadata.
        uniprot_mappings: UniProt/SIFTS sequence mappings.
        raw_categories: The mmCIF categories that were consulted to
            build this record.
    """

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
