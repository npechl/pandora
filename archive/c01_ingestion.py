"""Component 01 — mmCIF Ingestion: public functions."""
from __future__ import annotations

import gzip
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

import gemmi
import httpx

from archive.schemas.c01_ingestion import (
    Assembly,
    AssemblyGen,
    Atom,
    Chain,
    Entity,
    FetchOptions,
    IngestionProvenance,
    Ligand,
    MmCIFBatchInput,
    MmCIFBatchResult,
    MmCIFBatchSummary,
    MmCIFIngestionInput,
    MmCIFIngestionResult,
    ParsedStructure,
    Residue,
)
from pandora.schemas.common import Diagnostic, DiagnosticBundle


# Metals for simple is_ion detection (single-atom residue whose element is a metal)
_METALS: frozenset[str] = frozenset({
    "LI", "NA", "K", "RB", "CS", "BE", "MG", "CA", "SR", "BA",
    "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
    "Y", "ZR", "NB", "MO", "RU", "RH", "PD", "AG", "CD",
    "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",
    "AL", "GA", "IN", "SN", "TL", "PB", "BI",
    "LA", "CE", "PR", "ND", "SM", "EU", "GD", "TB", "DY",
    "HO", "ER", "TM", "YB", "LU",
})

_PROVIDER_URLS: dict[str, str] = {
    "pdbe": "https://www.ebi.ac.uk/pdbe/entry-files/download/{id}_updated.cif",
    "pdb":  "https://files.rcsb.org/download/{id}.cif",
}

_CACHE_DIR = Path.home() / ".pandora" / "cache" / "mmcif"

_ENTITY_TYPE_MAP = {
    gemmi.EntityType.Polymer:    "polymer",
    gemmi.EntityType.NonPolymer: "non-polymer",
    gemmi.EntityType.Water:      "water",
    gemmi.EntityType.Branched:   "branched",
}


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _cache_fresh(path: Path, max_age_seconds: int | None) -> bool:
    if max_age_seconds is None:
        return True
    return (time.time() - path.stat().st_mtime) < max_age_seconds


def _entity_sequence(ent: gemmi.Entity) -> str | None:
    """One-letter sequence for a polymer entity from ent.full_sequence (3-letter codes)."""
    if ent.entity_type != gemmi.EntityType.Polymer or not ent.full_sequence:
        return None
    return "".join(
        gemmi.find_tabulated_residue(code).one_letter_code.strip() or "X"
        for code in ent.full_sequence
    )


# ── Public API ────────────────────────────────────────────────────────────────

def fetch_mmcif(
    entry_id: str,
    provider: str,
    source_uri: str | None,
    fetch_options: FetchOptions,
) -> tuple[str, IngestionProvenance]:
    """Fetch a raw mmCIF file from a provider URL or local path, with optional disk cache."""
    cache = _CACHE_DIR / f"{entry_id.lower()}.cif"

    # Cache hit
    if (
        provider != "local"
        and fetch_options.use_cache
        and cache.exists()
        and _cache_fresh(cache, fetch_options.max_age_seconds)
    ):
        return cache.read_text(encoding="utf-8"), IngestionProvenance(
            provider=provider,
            source_uri=str(cache),
            retrieved_at=_now_iso(),
            from_cache=True,
        )

    # Local filesystem read
    if provider == "local":
        path = Path(source_uri or f"{entry_id.lower()}.cif")
        content = path.read_text(encoding="utf-8")
        return content, IngestionProvenance(
            provider=provider,
            source_uri=str(path),
            retrieved_at=_now_iso(),
            from_cache=False,
        )

    # Build remote URL
    if source_uri:
        url = source_uri
    elif provider in _PROVIDER_URLS:
        fmt_id = entry_id.lower() if provider == "pdbe" else entry_id.upper()
        url = _PROVIDER_URLS[provider].format(id=fmt_id)
    else:
        raise ValueError(f"provider={provider!r} requires an explicit source_uri")

    resp = httpx.get(url, follow_redirects=True, timeout=60.0)
    resp.raise_for_status()
    raw_bytes = resp.content

    if fetch_options.decompress and (url.endswith(".gz") or raw_bytes[:2] == b"\x1f\x8b"):
        raw_bytes = gzip.decompress(raw_bytes)
    content = raw_bytes.decode("utf-8")

    if fetch_options.use_cache:
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache.write_text(content, encoding="utf-8")

    return content, IngestionProvenance(
        provider=provider,
        source_uri=url,
        retrieved_at=_now_iso(),
        from_cache=False,
    )


def parse_mmcif(
    raw_content: str,
    entry_id: str = "",
) -> tuple[ParsedStructure | None, DiagnosticBundle, str]:
    """Parse raw mmCIF text into a ParsedStructure using gemmi."""
    diag = DiagnosticBundle()

    if not raw_content.strip():
        diag.errors.append(Diagnostic(
            code="EMPTY_CONTENT", severity="error",
            message="mmCIF content is empty", entry_id=entry_id or None,
        ))
        return None, diag, "failed"

    try:
        st = gemmi.read_structure_string(raw_content, format=gemmi.CoorFormat.Mmcif)
    except Exception as exc:
        diag.errors.append(Diagnostic(
            code="PARSE_ERROR", severity="error",
            message=str(exc), entry_id=entry_id or None,
        ))
        return None, diag, "failed"

    if len(st) == 0:
        diag.errors.append(Diagnostic(
            code="NO_MODEL", severity="error",
            message="Structure contains no models", entry_id=entry_id or None,
        ))
        return None, diag, "failed"

    model = st[0]
    eid = entry_id or st.name

    atoms_out: list[Atom] = []
    residues_out: list[Residue] = []
    chains_out: list[Chain] = []
    ligands_out: list[Ligand] = []
    _serial = 0

    # Build a subchain → entity_id lookup once
    subchain_to_entity: dict[str, str] = {}
    for ent in st.entities:
        for sc in ent.subchains:
            subchain_to_entity[sc] = ent.name

    for chain in model:
        residues_in_chain: list[Residue] = []
        chain_type: str = "non-polymer"
        entity_id: str = ""

        # Determine chain type + entity_id from the first subchain
        for span in chain.subchains():
            sc_name = span.subchain_id()
            entity_id = subchain_to_entity.get(sc_name, "")
            for ent in st.entities:
                if sc_name in ent.subchains:
                    chain_type = _ENTITY_TYPE_MAP.get(ent.entity_type, "non-polymer")
                    break
            break  # use the first subchain only for chain-level metadata

        for res in chain:
            et = res.entity_type
            is_polymer = et == gemmi.EntityType.Polymer
            is_water = et == gemmi.EntityType.Water
            rid = f"{eid}:{chain.name}:{res.seqid}:{res.name}"

            atoms_in_res: list[Atom] = []
            for atom in res:
                _serial += 1
                altloc = atom.altloc if atom.altloc not in ("\x00", " ", "") else None
                atoms_in_res.append(Atom(
                    atom_id=f"{rid}:{atom.name}:{_serial}",
                    atom_name=atom.name,
                    element=atom.element.name,
                    x=atom.pos.x,
                    y=atom.pos.y,
                    z=atom.pos.z,
                    occupancy=atom.occ,
                    b_factor=atom.b_iso,
                    altloc=altloc,
                    residue_id=rid,
                    chain_id=chain.name,
                ))
            atoms_out.extend(atoms_in_res)

            try:
                seq_id_num: int | None = res.seqid.num
                ins = res.seqid.icode
                ins_code: str | None = ins if ins not in (" ", "\x00", "") else None
            except Exception:
                seq_id_num, ins_code = None, None

            r = Residue(
                residue_id=rid,
                comp_id=res.name,
                seq_id=seq_id_num,
                auth_seq_id=str(res.seqid),
                insertion_code=ins_code,
                chain_id=chain.name,
                atoms=atoms_in_res,
                is_polymer=is_polymer,
            )
            residues_in_chain.append(r)
            residues_out.append(r)

            if not is_polymer and not is_water:
                atom_names = [a.element.name.upper() for a in res]
                ligands_out.append(Ligand(
                    ligand_id=rid,
                    chem_comp_id=res.name,
                    chain_id=chain.name,
                    residue_id=rid,
                    is_water=False,
                    is_ion=(len(atom_names) == 1 and atom_names[0] in _METALS),
                ))

        if not residues_in_chain:
            diag.warnings.append(Diagnostic(
                code="EMPTY_CHAIN", severity="warning",
                message=f"Chain {chain.name!r} has no residues",
                entry_id=eid, context={"chain_id": chain.name},
            ))

        chains_out.append(Chain(
            chain_id=chain.name,        # auth_asym_id; C02 normalises to label_asym_id
            auth_chain_id=chain.name,
            entity_id=entity_id,
            chain_type=chain_type,      # type: ignore[arg-type]
            residues=residues_in_chain,
        ))

    entities_out: list[Entity] = []
    for ent in st.entities:
        entities_out.append(Entity(
            entity_id=ent.name,
            entity_type=_ENTITY_TYPE_MAP.get(ent.entity_type, "non-polymer"),
            description=None,
            chain_ids=list(ent.subchains),
            sequence=_entity_sequence(ent),
        ))

    assemblies_out: list[Assembly] = []
    for asm in st.assemblies:
        gens = [
            AssemblyGen(
                asym_id_list=list(gen.chains),
                oper_expression=",".join(str(op) for op in gen.operators),
            )
            for gen in asm.generators
        ]
        assemblies_out.append(Assembly(
            assembly_id=asm.name,
            assembly_gen=gens,
        ))

    if not atoms_out:
        diag.warnings.append(Diagnostic(
            code="MISSING_ATOM_SITE", severity="warning",
            message="No atoms found in structure", entry_id=eid,
        ))

    return (
        ParsedStructure(
            atoms=atoms_out,
            residues=residues_out,
            chains=chains_out,
            entities=entities_out,
            assemblies=assemblies_out,
            ligands=ligands_out,
        ),
        diag,
        "warning" if diag.warnings else "success",
    )


def validate_mmcif(
    parsed_structure: ParsedStructure,
    entry_id: str = "",
) -> tuple[str, DiagnosticBundle]:
    # TODO: implement — V1 error/warning rules per spec Section 6.3
    #   (MISSING_ATOM_SITE, EMPTY_CHAIN, DISCONTINUOUS_SEQID, ...)
    return "valid", DiagnosticBundle()


def ingest_mmcif(inp: MmCIFIngestionInput) -> MmCIFIngestionResult:
    """Run the full single-entry ingestion workflow (fetch → parse → validate)."""
    if inp.raw_content is not None:
        raw = inp.raw_content
        prov = IngestionProvenance(
            provider=inp.provider,
            source_uri=inp.source_uri,
            retrieved_at=None,
            from_cache=False,
        )
    else:
        raw, prov = fetch_mmcif(
            inp.entry_id, inp.provider, inp.source_uri, inp.fetch_options
        )

    parsed, parse_diag, parse_status = parse_mmcif(raw, inp.entry_id)

    if parse_status == "failed":
        return MmCIFIngestionResult(
            entry_id=inp.entry_id,
            status="failed",
            parsed_structure=None,
            diagnostics=parse_diag,
            provenance=prov,
        )

    val_status, val_diag = validate_mmcif(parsed, inp.entry_id)  # type: ignore[arg-type]

    if val_status == "invalid":
        status = "failed"
    elif parse_status == "warning" or val_status == "warning":
        status = "warning"
    else:
        status = "success"

    return MmCIFIngestionResult(
        entry_id=inp.entry_id,
        status=status,
        parsed_structure=parsed,
        diagnostics=DiagnosticBundle(
            warnings=parse_diag.warnings + val_diag.warnings,
            errors=parse_diag.errors + val_diag.errors,
        ),
        provenance=prov,
    )


def ingest_list_mmcif(batch: MmCIFBatchInput) -> MmCIFBatchResult:
    """Run ingestion for a batch of entries (sequential or parallel)."""
    # TODO: implement parallel mode via concurrent.futures
    results = [ingest_mmcif(inp) for inp in batch.entries]
    summary = MmCIFBatchSummary(
        total=len(results),
        success=sum(1 for r in results if r.status == "success"),
        warning=sum(1 for r in results if r.status == "warning"),
        failed=sum(1 for r in results if r.status == "failed"),
    )
    return MmCIFBatchResult(mode=batch.mode, summary=summary, results=results)
