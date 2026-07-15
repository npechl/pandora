from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from typing import Any

from pandora.metadata.utils import (
    RawRow,
    WATER_COMP_IDS,
    as_float,
    as_int,
    clean,
    clean_sequence,
    first_row,
    first_value,
    provenance,
    raw_rows,
)
from pandora.schemas.metadata import (
    EntityMetadataRecord,
    EntryMetadataRecord,
    LigandMetadataRecord,
    QualityRecord,
    TaxonomyRecord,
    UniProtMappingRecord,
)
from pandora.schemas.structure import Structure


def extract_metadata_category(
    structure: Structure,
    category: str,
    columns: Iterable[str] | None = None,
) -> list[RawRow]:
    """Return raw mmCIF rows for one category.

    This is the escape hatch for metadata that Pandora has not promoted to a
    typed schema yet.

    Args:
        structure: The parsed structure to read raw categories from.
        category: The mmCIF category name (e.g. "_entity").
        columns: If given, restrict each row to these column names
            (missing columns come back as `None`); otherwise return
            every column present in the category.

    Returns:
        A list of rows, each a mapping from column name to its
        (already null-normalized) string value.
    """

    rows = raw_rows(structure, category)
    if columns is None:
        return [dict(row) for row in rows]

    selected = list(columns)
    return [{column: row.get(column) for column in selected} for row in rows]


def extract_metadata(
    structure: Structure,
    category: str,
    keys: Iterable[str] | None = None,
) -> list[RawRow]:
    """Compatibility wrapper around `extract_metadata_category()`.

    Args:
        structure: The parsed structure to read raw categories from.
        category: The mmCIF category name (e.g. "_entity").
        keys: If given, restrict each row to these column names;
            otherwise return every column present in the category.

    Returns:
        A list of rows, each a mapping from column name to its
        (already null-normalized) string value.
    """

    return extract_metadata_category(structure, category, keys)


def extract_entry_metadata(structure: Structure) -> EntryMetadataRecord:
    """Extract entry-level source-backed metadata from the parsed structure.

    Args:
        structure: The parsed structure to extract metadata from.

    Returns:
        An `EntryMetadataRecord` with entry id, title, keywords, and
        citation details (title, DOI, PubMed id).
    """

    struct_keywords = first_row(structure, "_struct_keywords")
    citation = first_row(structure, "_citation")

    return EntryMetadataRecord(
        entry_id=structure.entry_id,
        title=structure.entry.title,
        keywords=first_value(
            struct_keywords.get("text"),
            struct_keywords.get("pdbx_keywords"),
        ),
        citation_title=clean(citation.get("title")),
        doi=clean(citation.get("pdbx_database_id_DOI")),
        pubmed_id=clean(citation.get("pdbx_database_id_PubMed")),
        provenance=provenance("_struct"),
    )


def extract_taxonomies(structure: Structure) -> list[TaxonomyRecord]:
    """Extract source organism and expression host metadata.

    Reads the `_entity_src_gen`, `_entity_src_nat`, and
    `_pdbx_entity_src_syn` mmCIF categories, covering genetically
    manipulated, naturally isolated, and synthetic sources
    respectively.

    Args:
        structure: The parsed structure to extract taxonomies from.

    Returns:
        A `TaxonomyRecord` per source row found, in category order
        (generated, then natural, then synthetic sources).
    """

    taxonomies = [
        _taxonomy_from_generated_source(row)
        for row in raw_rows(structure, "_entity_src_gen")
    ]
    taxonomies.extend(
        _taxonomy_from_natural_source(row)
        for row in raw_rows(structure, "_entity_src_nat")
    )
    taxonomies.extend(
        _taxonomy_from_synthetic_source(row)
        for row in raw_rows(structure, "_pdbx_entity_src_syn")
    )
    return taxonomies


def extract_taxonomy(structure: Structure) -> TaxonomyRecord | None:
    """Return the first taxonomy record, if present.

    Args:
        structure: The parsed structure to extract taxonomy from.

    Returns:
        The first `TaxonomyRecord` from `extract_taxonomies()`, or
        `None` if the structure has no source records.
    """

    taxonomies = extract_taxonomies(structure)
    return taxonomies[0] if taxonomies else None


def extract_quality(structure: Structure) -> QualityRecord | None:
    """Extract experimental quality metrics reported in mmCIF.

    Reads the `_exptl`, `_refine`, and `_reflns` categories for the
    experimental method, resolution, R-work/R-free, reflection counts,
    and mean B-factor.

    Args:
        structure: The parsed structure to extract quality data from.

    Returns:
        A `QualityRecord`, or `None` if none of `_exptl`, `_refine`,
        or `_reflns` are present in the structure.
    """

    exptl_rows = raw_rows(structure, "_exptl")
    refine = first_row(structure, "_refine")
    reflns = first_row(structure, "_reflns")

    if not exptl_rows and not refine and not reflns:
        return None

    methods = [
        method
        for method in (clean(row.get("method")) for row in exptl_rows)
        if method is not None
    ]

    return QualityRecord(
        experimental_method="; ".join(methods) if methods else None,
        resolution=as_float(refine.get("ls_d_res_high")),
        r_work=as_float(refine.get("ls_R_factor_R_work")),
        r_free=as_float(refine.get("ls_R_factor_R_free")),
        observed_reflections=as_int(reflns.get("number_obs")),
        percent_possible_observed=as_float(reflns.get("percent_possible_obs")),
        mean_b_factor=as_float(refine.get("B_iso_mean")),
        provenance=provenance("_exptl,_refine,_reflns"),
    )


def extract_entity_metadata(structure: Structure) -> list[EntityMetadataRecord]:
    """Extract entity-level descriptions, sequences, and chain membership.

    Args:
        structure: The parsed structure to extract entity metadata
            from.

    Returns:
        An `EntityMetadataRecord` per entity, combining `_entity`
        fields (EC number, mutation, fragment), polymer sequence data,
        and the chain (asym) ids belonging to that entity.
    """

    entity_rows = {row.get("id"): row for row in raw_rows(structure, "_entity")}
    chain_ids_by_entity: dict[str, list[str]] = defaultdict(list)
    for asym in structure.asym_units:
        chain_ids_by_entity[asym.entity_id].append(asym.id)

    records: list[EntityMetadataRecord] = []
    for entity in structure.entities:
        raw = entity_rows.get(entity.id, {})
        poly = entity.poly
        records.append(
            EntityMetadataRecord(
                entity_id=entity.id,
                entity_type=entity.type,
                description=entity.pdbx_description,
                formula_weight=entity.formula_weight,
                source_method=entity.src_method,
                ec_number=clean(raw.get("pdbx_ec")),
                mutation=clean(raw.get("pdbx_mutation")),
                fragment=clean(raw.get("pdbx_fragment")),
                polymer_type=poly.type if poly else None,
                sequence=(
                    clean_sequence(poly.pdbx_seq_one_letter_code)
                    if poly
                    else None
                ),
                canonical_sequence=(
                    clean_sequence(poly.pdbx_seq_one_letter_code_can)
                    if poly
                    else None
                ),
                chain_ids=chain_ids_by_entity.get(entity.id, []),
                provenance=provenance("_entity", entity.id),
            )
        )
    return records


def extract_ligand_metadata(
    structure: Structure,
    include_waters: bool = False,
) -> list[LigandMetadataRecord]:
    """Extract non-polymer ligand metadata from mmCIF categories.

    Reads `_pdbx_entity_nonpoly` (joined with `_chem_comp` for name and
    formula) for one record per distinct (entity, ligand) pair. Falls
    back to deriving ligand records from HETATM records in
    `structure.atoms` when neither category is present.

    Args:
        structure: The parsed structure to extract ligands from.
        include_waters: If True, include water molecules as ligands.

    Returns:
        A `LigandMetadataRecord` per distinct ligand entity.
    """

    chem_comp = {
        row.get("id"): row for row in raw_rows(structure, "_chem_comp")
    }
    nonpoly_rows = raw_rows(structure, "_pdbx_entity_nonpoly")

    chain_ids_by_entity: dict[str, list[str]] = defaultdict(list)
    for asym in structure.asym_units:
        chain_ids_by_entity[asym.entity_id].append(asym.id)

    records: list[LigandMetadataRecord] = []
    seen: set[tuple[str | None, str]] = set()

    for row in nonpoly_rows:
        comp_id = clean(row.get("comp_id"))
        if comp_id is None:
            continue
        if not include_waters and comp_id.upper() in WATER_COMP_IDS:
            continue

        entity_id = clean(row.get("entity_id"))
        key = (entity_id, comp_id)
        if key in seen:
            continue
        seen.add(key)

        chem = chem_comp.get(comp_id, {})
        records.append(
            LigandMetadataRecord(
                entity_id=entity_id,
                comp_id=comp_id,
                name=first_value(row.get("name"), chem.get("name")),
                formula=clean(chem.get("formula")),
                formula_weight=as_float(chem.get("formula_weight")),
                chain_ids=chain_ids_by_entity.get(entity_id or "", []),
                provenance=provenance("_pdbx_entity_nonpoly", entity_id),
            )
        )

    if records:
        return records

    return _ligands_from_atoms(structure, include_waters)


def extract_uniprot_mappings(
    structure: Structure,
) -> list[UniProtMappingRecord]:
    """Extract UniProt/SIFTS mappings reported in the mmCIF file.

    Reads `_struct_ref`/`_struct_ref_seq` and
    `_pdbx_sifts_unp_segments`, deduplicating mappings that describe
    the same entity/asym/accession/range.

    Args:
        structure: The parsed structure to extract mappings from.

    Returns:
        A `UniProtMappingRecord` per distinct mapping segment found
        across both source categories.
    """

    records: list[UniProtMappingRecord] = []
    seen: set[tuple[Any, ...]] = set()

    struct_ref = {
        row.get("id"): row for row in raw_rows(structure, "_struct_ref")
    }

    for row in raw_rows(structure, "_struct_ref_seq"):
        ref = struct_ref.get(row.get("ref_id"), {})
        accession = first_value(
            row.get("pdbx_db_accession"),
            ref.get("pdbx_db_accession"),
        )
        if accession is None:
            continue
        record = UniProtMappingRecord(
            entity_id=clean(ref.get("entity_id")),
            asym_id=clean(row.get("pdbx_strand_id")),
            accession=accession,
            db_name=clean(ref.get("db_name")),
            db_code=clean(ref.get("db_code")),
            seq_id_start=as_int(row.get("seq_align_beg")),
            seq_id_end=as_int(row.get("seq_align_end")),
            uniprot_start=as_int(row.get("db_align_beg")),
            uniprot_end=as_int(row.get("db_align_end")),
            provenance=provenance("_struct_ref_seq", row.get("align_id")),
        )
        _append_mapping(records, seen, record)

    for row in raw_rows(structure, "_pdbx_sifts_unp_segments"):
        accession = clean(row.get("unp_acc"))
        if accession is None:
            continue
        record = UniProtMappingRecord(
            entity_id=clean(row.get("entity_id")),
            asym_id=clean(row.get("asym_id")),
            accession=accession,
            db_name="UNP",
            seq_id_start=as_int(row.get("seq_id_start")),
            seq_id_end=as_int(row.get("seq_id_end")),
            uniprot_start=as_int(row.get("unp_start")),
            uniprot_end=as_int(row.get("unp_end")),
            identity=as_float(row.get("identity")),
            provenance=provenance(
                "_pdbx_sifts_unp_segments", row.get("segment_id")
            ),
        )
        _append_mapping(records, seen, record)

    return records


def _taxonomy_from_generated_source(row: RawRow) -> TaxonomyRecord:
    return TaxonomyRecord(
        entity_id=clean(row.get("entity_id")),
        ncbi_taxon_id=as_int(row.get("pdbx_gene_src_ncbi_taxonomy_id")),
        organism_scientific=first_value(
            row.get("pdbx_gene_src_scientific_name"),
            row.get("gene_src_scientific_name"),
        ),
        organism_common=clean(row.get("gene_src_common_name")),
        host_ncbi_taxon_id=as_int(row.get("pdbx_host_org_ncbi_taxonomy_id")),
        host_scientific=clean(row.get("pdbx_host_org_scientific_name")),
        host_common=clean(row.get("host_org_common_name")),
        expression_system=first_value(
            row.get("expression_system_id"),
            row.get("plasmid_name"),
        ),
        provenance=provenance("_entity_src_gen", row.get("entity_id")),
    )


def _taxonomy_from_natural_source(row: RawRow) -> TaxonomyRecord:
    return TaxonomyRecord(
        entity_id=clean(row.get("entity_id")),
        ncbi_taxon_id=as_int(
            first_value(
                row.get("pdbx_ncbi_taxonomy_id"),
                row.get("ncbi_taxonomy_id"),
            )
        ),
        organism_scientific=first_value(
            row.get("pdbx_organism_scientific"),
            row.get("organism_scientific"),
        ),
        organism_common=first_value(
            row.get("common_name"),
            row.get("organism_common_name"),
        ),
        provenance=provenance("_entity_src_nat", row.get("entity_id")),
    )


def _taxonomy_from_synthetic_source(row: RawRow) -> TaxonomyRecord:
    return TaxonomyRecord(
        entity_id=clean(row.get("entity_id")),
        ncbi_taxon_id=as_int(row.get("ncbi_taxonomy_id")),
        organism_scientific=clean(row.get("organism_scientific")),
        organism_common=clean(row.get("organism_common_name")),
        provenance=provenance("_pdbx_entity_src_syn", row.get("entity_id")),
    )


def _ligands_from_atoms(
    structure: Structure,
    include_waters: bool,
) -> list[LigandMetadataRecord]:
    by_comp: dict[tuple[str | None, str], set[str]] = defaultdict(set)
    for atom in structure.atoms:
        if atom.group_PDB != "HETATM":
            continue
        comp_id = atom.label_comp_id
        if not include_waters and comp_id.upper() in WATER_COMP_IDS:
            continue
        by_comp[(atom.label_entity_id, comp_id)].add(atom.label_asym_id)

    return [
        LigandMetadataRecord(
            entity_id=entity_id,
            comp_id=comp_id,
            chain_ids=sorted(chain_ids),
            provenance=provenance("_atom_site", entity_id),
        )
        for (entity_id, comp_id), chain_ids in sorted(by_comp.items())
    ]


def _append_mapping(
    records: list[UniProtMappingRecord],
    seen: set[tuple[Any, ...]],
    record: UniProtMappingRecord,
) -> None:
    key = (
        record.entity_id,
        record.asym_id,
        record.accession,
        record.seq_id_start,
        record.seq_id_end,
        record.uniprot_start,
        record.uniprot_end,
    )
    if key not in seen:
        seen.add(key)
        records.append(record)
