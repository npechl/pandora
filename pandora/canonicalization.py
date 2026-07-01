from __future__ import annotations

import string
from collections import defaultdict
from datetime import datetime, timezone
from typing import Iterator

from pandora.schemas.structure import (
    AsymRecord,
    AssemblyRecord,
    AtomSiteRecord,
    EntityRecord,
    Structure,
)
from pandora.schemas.canonicalization import (
    AltlocSelectionMapping,
    AltlocSelectionMappingItem,
    AssemblyMapping,
    AssemblyMappingItem,
    CanonicalizationPolicy,
    CanonicalMappings,
    ChainIdMapping,
    ChainIdMappingItem,
    EntityMapping,
    EntityMappingItem,
    ResidueNumberMapping,
    ResidueNumberMappingItem,
    CanonicalizationProvenance
)

from pandora.schemas.common import Diagnostic, DiagnosticBundle


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()

# dump function for remap ----------------
def _sequential_chain_ids() -> Iterator[str] :
    letters = string.ascii_uppercase
    for c in letters:
        yield c
    for c1 in letters:
        for c2 in letters:
            yield c1 + c2


# normalize_chain_ids ------------------------

def _normalize_chain_ids(
    asym_units: list[AsymRecord],
    strategy: str,
    record: bool,
) -> tuple[dict[str, str], ChainIdMapping]:
    chain_map: dict[str, str] = {}
    mapping = ChainIdMapping()
    id_gen = _sequential_chain_ids()

    for asym in asym_units:
        if strategy == "preserve":
            canonical = asym.id
        elif strategy == "use_auth_chain_id":
            canonical = asym.auth_id or asym.id
        else:  # remap
            canonical = next(id_gen)

        chain_map[asym.id] = canonical
        if record:
            mapping.items.append(ChainIdMappingItem(
                canonical_chain_id=canonical,
                original_chain_id=asym.id,
                original_auth_chain_id=asym.auth_id or asym.id,
            ))

    return chain_map, mapping


def _apply_chain_map(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    structure: Structure,
    chain_map: dict[str, str],
) -> tuple[list[AtomSiteRecord], list[AsymRecord], Structure]:
    
    new_atoms = [
        a.model_copy(update={"label_asym_id": chain_map.get(a.label_asym_id, a.label_asym_id)})
        if chain_map.get(a.label_asym_id, a.label_asym_id) != a.label_asym_id else a
        for a in atoms
    ]
    
    new_asyms = [
        AsymRecord(id=chain_map.get(a.id, a.id), entity_id=a.entity_id, auth_id=a.auth_id)
        for a in asym_units
    ]

    # Update assembly generators
    new_assemblies = []
    for asm in structure.assemblies:
        new_gens = []
        for gen in asm.generators:
            new_asym_list = [chain_map.get(aid, aid) for aid in gen.asym_id_list]
            new_gens.append(gen.model_copy(update={"asym_id_list": new_asym_list}))
        new_assemblies.append(asm.model_copy(update={"generators": new_gens}))

    # Update connections
    new_conns = []
    for conn in structure.connections:
        p1 = conn.ptnr1.model_copy(update={
            "label_asym_id": chain_map.get(conn.ptnr1.label_asym_id, conn.ptnr1.label_asym_id)
        })
        p2 = conn.ptnr2.model_copy(update={
            "label_asym_id": chain_map.get(conn.ptnr2.label_asym_id, conn.ptnr2.label_asym_id)
        })
        new_conns.append(conn.model_copy(update={"ptnr1": p1, "ptnr2": p2}))

    # Update secondary structure
    new_conf = [
        r.model_copy(update={
            "beg_label_asym_id": chain_map.get(r.beg_label_asym_id, r.beg_label_asym_id),
            "end_label_asym_id": chain_map.get(r.end_label_asym_id, r.end_label_asym_id),
        })
        for r in structure.secondary_structure.conf_records
    ]

    new_strands = [
        r.model_copy(update={
            "beg_label_asym_id": chain_map.get(r.beg_label_asym_id, r.beg_label_asym_id),
            "end_label_asym_id": chain_map.get(r.end_label_asym_id, r.end_label_asym_id),
        })
        for r in structure.secondary_structure.sheet_strands
    ]

    new_ss = structure.secondary_structure.model_copy(
        update={"conf_records": new_conf, "sheet_strands": new_strands}
    )

    new_structure = structure.model_copy(update={
        "assemblies": new_assemblies,
        "connections": new_conns,
        "secondary_structure": new_ss,
    })

    return new_atoms, new_asyms, new_structure


# normalize_residue_numbering ------------------------

def _normalize_residue_numbering(
    atoms: list[AtomSiteRecord],
    chain_map: dict[str, str],
    strategy: str,
    preserve_insertion_codes: bool,
    record: bool,
) -> tuple[list[AtomSiteRecord], ResidueNumberMapping]:
    mapping = ResidueNumberMapping()

    if strategy == "preserve":
        if preserve_insertion_codes:
            return atoms, mapping
        result = []
        for a in atoms:
            ins = a.pdbx_PDB_ins_code
            if ins and ins not in (".", "?"):
                a = a.model_copy(update={"pdbx_PDB_ins_code": None})
            result.append(a)
        return result, mapping

    seen: set[tuple] = set()

    if strategy == "use_auth_seq":
        result = []
        for a in atoms:
            raw = a.auth_seq_id or ""
            digits = "".join(c for c in raw if c.isdigit() or c == "-")
            try:
                new_seq = int(digits)
            except ValueError:
                new_seq = a.label_seq_id

            ins = a.pdbx_PDB_ins_code if preserve_insertion_codes else None

            key = (a.label_asym_id, a.label_seq_id, a.auth_seq_id)
            if record and key not in seen:
                seen.add(key)
                mapping.items.append(ResidueNumberMappingItem(
                    canonical_chain_id=chain_map.get(a.label_asym_id, a.label_asym_id),
                    canonical_seq_id=new_seq or 0,
                    original_chain_id=a.label_asym_id,
                    original_seq_id=a.label_seq_id,
                    original_auth_seq_id=a.auth_seq_id or "",
                    original_insertion_code=a.pdbx_PDB_ins_code,
                ))
            result.append(a.model_copy(update={"label_seq_id": new_seq, "pdbx_PDB_ins_code": ins}))
        return result, mapping

    # strategy == "renumber"
    chain_res_order: dict[str, list[tuple]] = defaultdict(list)
    chain_res_seen: dict[str, set] = defaultdict(set)
    for a in atoms:
        rk = (a.label_seq_id, a.label_comp_id, a.pdbx_PDB_ins_code)
        if rk not in chain_res_seen[a.label_asym_id]:
            chain_res_seen[a.label_asym_id].add(rk)
            chain_res_order[a.label_asym_id].append(rk)

    renumber: dict[str, dict[tuple, int]] = {
        ch: {rk: i + 1 for i, rk in enumerate(rks)}
        for ch, rks in chain_res_order.items()
    }

    result = []
    for a in atoms:
        rk = (a.label_seq_id, a.label_comp_id, a.pdbx_PDB_ins_code)
        new_seq = renumber[a.label_asym_id].get(rk, a.label_seq_id)

        key = (a.label_asym_id, a.label_seq_id, a.auth_seq_id)
        if record and key not in seen:
            seen.add(key)
            mapping.items.append(ResidueNumberMappingItem(
                canonical_chain_id=chain_map.get(a.label_asym_id, a.label_asym_id),
                canonical_seq_id=new_seq or 0,
                original_chain_id=a.label_asym_id,
                original_seq_id=a.label_seq_id,
                original_auth_seq_id=a.auth_seq_id or "",
                original_insertion_code=a.pdbx_PDB_ins_code,
            ))
        result.append(a.model_copy(update={"label_seq_id": new_seq, "pdbx_PDB_ins_code": None}))

    return result, mapping


# normalize_assemblies ------------------------------

def _normalize_assemblies(
    assemblies: list[AssemblyRecord],
    assembly_rules,
    id_strategy: str,
    record: bool,
) -> tuple[list[AssemblyRecord], AssemblyMapping]:
    mapping = AssemblyMapping()
    result = list(assemblies)

    if assembly_rules.strategy == "select_first_assembly" and result:
        result = [result[0]]

    # standardize_biological_assembly: without external metadata we keep order as-is;
    # the preferred assembly (if annotated) would be identified from pdbx flags in raw data.

    if id_strategy == "preserve":
        if record:
            for asm in result:
                mapping.items.append(AssemblyMappingItem(
                    canonical_assembly_id=asm.id,
                    original_assembly_id=asm.id,
                ))
        return result, mapping

    # remap or standardize → sequential integers
    new_result = []
    for i, asm in enumerate(result, 1):
        new_id = str(i)
        if record:
            mapping.items.append(AssemblyMappingItem(
                canonical_assembly_id=new_id,
                original_assembly_id=asm.id,
            ))
        new_result.append(asm.model_copy(update={"id": new_id}))

    return new_result, mapping


# handle_missing_data ---------------------------

_BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})


def _handle_missing_atoms(
    atoms: list[AtomSiteRecord],
    rules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> list[AtomSiteRecord]:
    strategy = rules.strategy
    if strategy in ("preserve", "impute"):
        return atoms

    polymer_by_residue: dict[tuple, list[AtomSiteRecord]] = defaultdict(list)
    hetatm: list[AtomSiteRecord] = []
    for a in atoms:
        if a.group_PDB == "ATOM":
            polymer_by_residue[(a.label_asym_id, a.label_seq_id, a.label_comp_id)].append(a)
        else:
            hetatm.append(a)

    result: list[AtomSiteRecord] = []
    for key, res_atoms in polymer_by_residue.items():
        present = {a.label_atom_id for a in res_atoms}
        missing = _BACKBONE_ATOMS - present
        if missing:
            asym_id, seq_id, comp_id = key
            if rules.record_missingness:
                diagnostics.warnings.append(Diagnostic(
                    code="MISSING_ATOMS",
                    severity="warning",
                    message=f"Residue {comp_id} {seq_id} in chain {asym_id} missing backbone atoms",
                    entry_id=entry_id,
                    context={"chain": asym_id, "seq_id": seq_id, "comp_id": comp_id,
                             "missing": sorted(missing)},
                ))
            if strategy == "drop_partial_residue":
                continue
        result.extend(res_atoms)

    result.extend(hetatm)
    return result


def _handle_missing_residues(
    atoms: list[AtomSiteRecord],
    rules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> list[AtomSiteRecord]:
    strategy = rules.strategy
    record = rules.record_gaps

    if strategy == "preserve" and not record:
        return atoms

    chain_seq_ids: dict[str, list[int]] = defaultdict(list)
    chain_seen: dict[str, set] = defaultdict(set)
    for a in atoms:
        if a.group_PDB == "ATOM" and a.label_seq_id is not None:
            if a.label_seq_id not in chain_seen[a.label_asym_id]:
                chain_seen[a.label_asym_id].add(a.label_seq_id)
                chain_seq_ids[a.label_asym_id].append(a.label_seq_id)

    gap_chains: set[str] = set()
    for chain, seq_ids in chain_seq_ids.items():
        sorted_ids = sorted(seq_ids)
        for i in range(len(sorted_ids) - 1):
            if sorted_ids[i + 1] - sorted_ids[i] > 1:
                gap_chains.add(chain)
                if record:
                    diagnostics.warnings.append(Diagnostic(
                        code="SEQUENCE_GAP",
                        severity="warning",
                        message=(
                            f"Sequence gap in chain {chain} between "
                            f"residues {sorted_ids[i]} and {sorted_ids[i + 1]}"
                        ),
                        entry_id=entry_id,
                        context={
                            "chain": chain,
                            "gap_start": sorted_ids[i],
                            "gap_end": sorted_ids[i + 1],
                        },
                    ))

    if strategy == "drop_chain_segment" and gap_chains:
        return [a for a in atoms if a.label_asym_id not in gap_chains]

    return atoms


def _handle_incomplete_chains(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    rules,
) -> tuple[list[AtomSiteRecord], list[AsymRecord]]:
    strategy = rules.strategy
    if strategy == "preserve":
        return atoms, asym_units

    chain_seq_ids: dict[str, list[int]] = defaultdict(list)
    chain_seen: dict[str, set] = defaultdict(set)
    for a in atoms:
        if a.group_PDB == "ATOM" and a.label_seq_id is not None:
            if a.label_seq_id not in chain_seen[a.label_asym_id]:
                chain_seen[a.label_asym_id].add(a.label_seq_id)
                chain_seq_ids[a.label_asym_id].append(a.label_seq_id)

    incomplete: set[str] = set()
    for chain, seq_ids in chain_seq_ids.items():
        s = sorted(seq_ids)
        if s[-1] - s[0] + 1 > len(s):
            incomplete.add(chain)

    if strategy == "exclude":
        keep = {a.id for a in asym_units} - incomplete
        return (
            [a for a in atoms if a.label_asym_id in keep],
            [a for a in asym_units if a.id in keep],
        )

    # truncate_to_complete_regions: keep longest contiguous run per incomplete chain
    keep_ids: dict[str, set[int]] = {}
    for chain in incomplete:
        s = sorted(chain_seq_ids[chain])
        best: list[int] = []
        run: list[int] = [s[0]]
        for i in range(1, len(s)):
            if s[i] == s[i - 1] + 1:
                run.append(s[i])
            else:
                if len(run) > len(best):
                    best = run
                run = [s[i]]
        if len(run) > len(best):
            best = run
        keep_ids[chain] = set(best)

    result = []
    for a in atoms:
        if a.label_asym_id in incomplete:
            if a.label_seq_id in keep_ids.get(a.label_asym_id, set()):
                result.append(a)
        else:
            result.append(a)

    return result, asym_units


# resolve_altlocs ------------------------------

def _resolve_altlocs(
    atoms: list[AtomSiteRecord],
    rules,
) -> tuple[list[AtomSiteRecord], AltlocSelectionMapping]:
    mapping = AltlocSelectionMapping()
    strategy = rules.strategy

    if strategy == "preserve":
        return atoms, mapping

    no_alt: list[AtomSiteRecord] = []
    by_residue: dict[tuple, dict[str, list[AtomSiteRecord]]] = defaultdict(lambda: defaultdict(list))

    for a in atoms:
        alt = a.label_alt_id
        if not alt or alt in (".", "?"):
            no_alt.append(a)
        else:
            key = (a.label_asym_id, a.label_seq_id, a.label_comp_id)
            by_residue[key][alt].append(a)

    result: list[AtomSiteRecord] = list(no_alt)

    for (asym_id, seq_id, comp_id), alt_groups in by_residue.items():
        altlocs = sorted(alt_groups.keys())

        if len(altlocs) == 1:
            selected = altlocs[0]
        elif strategy == "select_first":
            selected = altlocs[0]
        elif strategy == "select_user_defined":
            udl = rules.user_defined_altloc
            selected = udl if udl in altlocs else altlocs[0]
        else:  # select_best_occupancy
            def _mean_occ(alt: str) -> float:
                grp = alt_groups[alt]
                return sum(a.occupancy for a in grp) / len(grp)

            occ = {alt: _mean_occ(alt) for alt in altlocs}
            max_occ = max(occ.values())
            tied = [alt for alt, o in occ.items() if o == max_occ]

            if len(tied) == 1:
                selected = tied[0]
            else:
                tb = rules.tie_breaker
                if tb == "alphabetical_last":
                    selected = max(tied)
                elif tb == "lowest_b_factor":
                    def _mean_b(alt: str) -> float:
                        grp = alt_groups[alt]
                        return sum(a.B_iso_or_equiv for a in grp) / len(grp)
                    selected = min(tied, key=_mean_b)
                elif tb == "highest_b_factor":
                    def _mean_b(alt: str) -> float:
                        grp = alt_groups[alt]
                        return sum(a.B_iso_or_equiv for a in grp) / len(grp)
                    selected = max(tied, key=_mean_b)
                else:  # alphabetical_first
                    selected = min(tied)

        for a in alt_groups[selected]:
            result.append(a.model_copy(update={"label_alt_id": None}))

        if rules.record_selection and len(altlocs) > 1:
            reason_map = {
                "select_first": "first_alphabetical",
                "select_user_defined": "user_defined",
                "select_best_occupancy": "best_occupancy",
            }
            mapping.items.append(AltlocSelectionMappingItem(
                canonical_chain_id=asym_id,
                residue_id=f"{seq_id}:{comp_id}",
                selected_altloc=selected,
                available_altlocs=altlocs,
                selection_reason=reason_map.get(strategy, "first_alphabetical"),
            ))

    return result, mapping


# normalize_entities -----------------------------------

def _normalize_entities(
    entities: list[EntityRecord],
    asym_units: list[AsymRecord],
    atoms: list[AtomSiteRecord],
    rules,
    record: bool,
) -> tuple[list[EntityRecord], list[AsymRecord], list[AtomSiteRecord], EntityMapping]:
    mapping = EntityMapping()

    if rules.strategy == "preserve":
        if record:
            for ent in entities:
                mapping.items.append(EntityMappingItem(
                    canonical_entity_id=ent.id,
                    original_entity_ids=[ent.id],
                ))
        return entities, asym_units, atoms, mapping

    entity_id_map: dict[str, str] = {}

    if rules.strategy == "standardize":
        new_entities = []
        for i, ent in enumerate(entities, 1):
            new_id = str(i)
            entity_id_map[ent.id] = new_id
            new_entities.append(ent.model_copy(update={"id": new_id}))
            if record:
                orig = ent.id if rules.preserve_original_entity_ids else new_id
                mapping.items.append(EntityMappingItem(
                    canonical_entity_id=new_id,
                    original_entity_ids=[orig],
                ))

    elif rules.strategy == "merge_equivalent_entities":
        seq_to_id: dict[str, str] = {}
        merged_orig: dict[str, list[str]] = defaultdict(list)
        new_entities = []
        counter = 1
        for ent in entities:
            seq_key = (
                ent.poly.pdbx_seq_one_letter_code_can
                if ent.poly and ent.poly.pdbx_seq_one_letter_code_can
                else None
            )
            if seq_key and seq_key in seq_to_id:
                canonical_id = seq_to_id[seq_key]
                entity_id_map[ent.id] = canonical_id
                merged_orig[canonical_id].append(ent.id)
            else:
                canonical_id = str(counter)
                counter += 1
                entity_id_map[ent.id] = canonical_id
                if seq_key:
                    seq_to_id[seq_key] = canonical_id
                merged_orig[canonical_id] = [ent.id]
                new_entities.append(ent.model_copy(update={"id": canonical_id}))

        if record:
            for cid, orig_ids in merged_orig.items():
                mapping.items.append(EntityMappingItem(
                    canonical_entity_id=cid,
                    original_entity_ids=orig_ids,
                ))
    else:
        return entities, asym_units, atoms, mapping

    new_asyms = [
        a.model_copy(update={"entity_id": entity_id_map.get(a.entity_id, a.entity_id)})
        for a in asym_units
    ]
    new_atoms = [
        a.model_copy(update={"label_entity_id": entity_id_map.get(a.label_entity_id, a.label_entity_id)})
        if entity_id_map.get(a.label_entity_id, a.label_entity_id) != a.label_entity_id
        else a
        for a in atoms
    ]
    return new_entities, new_asyms, new_atoms, mapping


# filter_ligands ---------------------------

def _filter_ligands(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    entities: list[EntityRecord],
    rules,
) -> tuple[list[AtomSiteRecord], list[AsymRecord]]:
    if rules.strategy != "filter":
        return atoms, asym_units

    entity_type: dict[str, str] = {e.id: e.type for e in entities}
    entity_desc: dict[str, str] = {
        e.id: (e.pdbx_description or "").upper() for e in entities
    }
    asym_entity: dict[str, str] = {a.id: a.entity_id for a in asym_units}

    _ION_KEYWORDS = frozenset({
        "ION", "ZINC", "CALCIUM", "MAGNESIUM", "SODIUM", "POTASSIUM",
        "IRON", "COPPER", "MANGANESE", "COBALT", "NICKEL", "CHLORIDE",
        "SULFATE", "PHOSPHATE",
    })

    def _keep(asym_id: str) -> bool:
        eid = asym_entity.get(asym_id)
        if eid is None:
            return True
        etype = entity_type.get(eid, "polymer")
        if etype == "polymer":
            return True
        if etype == "water":
            return rules.keep_waters
        desc = entity_desc.get(eid, "")
        is_ion = any(kw in desc for kw in _ION_KEYWORDS)
        if is_ion:
            return rules.keep_ions
        return rules.keep_nonpolymer_ligands

    keep = {a.id for a in asym_units if _keep(a.id)}
    return (
        [a for a in atoms if a.label_asym_id in keep],
        [a for a in asym_units if a.id in keep],
    )


# validate_canonical_structure ------------------------------

def _validate(
    atoms: list[AtomSiteRecord],
    asym_units: list[AsymRecord],
    rules,
    diagnostics: DiagnosticBundle,
    entry_id: str,
) -> str:
    chain_ids = [a.id for a in asym_units]
    if len(chain_ids) != len(set(chain_ids)):
        diagnostics.errors.append(Diagnostic(
            code="CANONICAL_CHAIN_ID_COLLISION",
            severity="error",
            message="Duplicate canonical chain IDs",
            entry_id=entry_id,
        ))

    chain_res_keys: dict[str, set] = defaultdict(set)
    for a in atoms:
        if a.group_PDB == "ATOM" and a.label_seq_id is not None:
            rk = (a.label_seq_id, a.pdbx_PDB_ins_code)
            if rk in chain_res_keys[a.label_asym_id]:
                diagnostics.errors.append(Diagnostic(
                    code="RESIDUE_NUMBER_COLLISION",
                    severity="error",
                    message=(
                        f"Residue number collision in chain {a.label_asym_id} "
                        f"at seq_id {a.label_seq_id}"
                    ),
                    entry_id=entry_id,
                ))
            chain_res_keys[a.label_asym_id].add(rk)

    if rules.strictness == "strict":
        for w in list(diagnostics.warnings):
            diagnostics.errors.append(
                w.model_copy(update={"severity": "error"})
            )

    has_errors = bool(diagnostics.errors)
    has_warnings = bool(diagnostics.warnings)

    if has_errors and rules.fail_on_unresolved_issues:
        return "failed"
    if rules.warnings_as_errors and has_warnings:
        return "failed"
    if has_errors or has_warnings:
        return "warning"
    return "success"


# Public function ----------------------------------

def canonicalize_structure(
    structure: Structure,
    policy: CanonicalizationPolicy,
) -> tuple[Structure, CanonicalMappings, CanonicalizationProvenance]:
    """Convert a parsed structure into a canonical structure."""
    ir = policy.identifier_rules
    mdr = policy.missing_data_rules
    ar = policy.altloc_rules
    asmr = policy.assembly_rules
    er = policy.entity_rules
    lr = policy.ligand_rules
    vr = policy.validation_rules
    pr = policy.provenance_rules

    record = pr.record_original_mappings
    transforms: list[str] = []
    diagnostics = DiagnosticBundle()

    atoms: list[AtomSiteRecord] = list(structure.atoms)
    asym_units: list[AsymRecord] = list(structure.asym_units)
    entities: list[EntityRecord] = list(structure.entities)
    assemblies: list[AssemblyRecord] = list(structure.assemblies)

    # normalize_chain_ids --------------------------------
    chain_map, chain_id_mapping = _normalize_chain_ids(
        asym_units, ir.chain_id.strategy, record
    )
    atoms, asym_units, structure = _apply_chain_map(atoms, asym_units, structure, chain_map)
    assemblies = list(structure.assemblies)
    if ir.chain_id.strategy != "preserve":
        transforms.append(f"chain_id:{ir.chain_id.strategy}")

    # normalize_residue_numbering ---------------------------------
    atoms, residue_number_mapping = _normalize_residue_numbering(
        atoms,
        chain_map,
        ir.residue_numbering.strategy,
        ir.residue_numbering.preserve_insertion_codes,
        record,
    )
    if ir.residue_numbering.strategy != "preserve":
        transforms.append(f"residue_numbering:{ir.residue_numbering.strategy}")

    # normalize_assemblies --------------------------------
    assemblies, assembly_mapping = _normalize_assemblies(
        assemblies, asmr, ir.assembly_id.strategy, record
    )
    if ir.assembly_id.strategy != "preserve" or asmr.strategy != "preserve_as_reported":
        transforms.append(
            f"assembly:{asmr.strategy}/id:{ir.assembly_id.strategy}"
        )

    # handle_missing_atoms ---------------------------------------
    atoms = _handle_missing_atoms(atoms, mdr.missing_atoms, diagnostics, structure.entry_id)
    if mdr.missing_atoms.strategy not in ("preserve",):
        transforms.append(f"missing_atoms:{mdr.missing_atoms.strategy}")

    # handle_missing_residues ---------------------------------------
    atoms = _handle_missing_residues(atoms, mdr.missing_residues, diagnostics, structure.entry_id)
    if mdr.missing_residues.strategy not in ("preserve",):
        transforms.append(f"missing_residues:{mdr.missing_residues.strategy}")

    # handle_incomplete_chains -------------------------------------
    atoms, asym_units = _handle_incomplete_chains(atoms, asym_units, mdr.incomplete_chains)
    if mdr.incomplete_chains.strategy != "preserve":
        transforms.append(f"incomplete_chains:{mdr.incomplete_chains.strategy}")

    # esolve_altlocs ------------------------------
    atoms, altloc_selection_mapping = _resolve_altlocs(atoms, ar)
    if ar.strategy != "preserve":
        transforms.append(f"altloc:{ar.strategy}")

    # normalize_entities ---------------------------------------
    entities, asym_units, atoms, entity_mapping = _normalize_entities(
        entities, asym_units, atoms, er, record
    )
    if er.strategy != "preserve":
        transforms.append(f"entity:{er.strategy}")

    # filter_ligands ---------------------------------------
    atoms, asym_units = _filter_ligands(atoms, asym_units, entities, lr)
    if lr.strategy != "preserve":
        transforms.append(f"ligands:{lr.strategy}")

    # validate_canonical_structure ---------------------------------------
    _validate(atoms, asym_units, vr, diagnostics, structure.entry_id)

    canonical = structure.model_copy(update={
        "atoms": atoms,
        "asym_units": asym_units,
        "entities": entities,
        "assemblies": assemblies,
    })

    mappings = CanonicalMappings(
        chain_id_mapping=chain_id_mapping,
        residue_number_mapping=residue_number_mapping,
        assembly_mapping=assembly_mapping,
        entity_mapping=entity_mapping,
        altloc_selection_mapping=altloc_selection_mapping,
    )

    provenance = CanonicalizationProvenance(
        canonicalized_at=_now_iso(),
        policy_id=policy.policy_id,
        policy_name=policy.policy_name,
        policy_version=policy.policy_version,
        transforms=transforms,
        report={
            "warnings": len(diagnostics.warnings),
            "errors": len(diagnostics.errors),
        } if pr.emit_canonicalization_report else {},
    )

    return canonical, mappings, provenance
