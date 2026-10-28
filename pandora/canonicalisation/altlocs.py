from collections import defaultdict

from pandora.schemas.canonicalisation import (
    AltlocSelectionMapping,
    AltlocSelectionMappingItem,
)
from pandora.schemas.structure import AtomSiteRecord

_SELECTION_REASON = {
    "select_first": "first_alphabetical",
    "select_user_defined": "user_defined",
    "select_best_occupancy": "best_occupancy",
}


def _has_altloc(a: AtomSiteRecord) -> bool:
    """Whether a carries a real altloc code."""

    return bool(a.label_alt_id) and a.label_alt_id not in (".", "?")


def _site(a: AtomSiteRecord) -> tuple:
    """Physical site of a: chain plus residue numbers."""

    return (a.label_asym_id, a.auth_seq_id, a.label_seq_id)


def _mean(atoms: list[AtomSiteRecord], field: str) -> float:
    """Mean of field over atoms."""

    return sum(getattr(a, field) for a in atoms) / len(atoms)


def _select_altloc(alt_groups: dict[str, list[AtomSiteRecord]], rules) -> str:
    """The altloc to keep at one site, per the altloc rules."""

    altlocs = sorted(alt_groups)
    if len(altlocs) == 1 or rules.strategy == "select_first":
        return altlocs[0]
    if rules.strategy == "select_user_defined":
        udl = rules.user_defined_altloc
        return udl if udl in altlocs else altlocs[0]

    # select_best_occupancy
    occ = {alt: _mean(alt_groups[alt], "occupancy") for alt in altlocs}
    max_occ = max(occ.values())
    tied = [alt for alt, o in occ.items() if o == max_occ]
    if len(tied) == 1:
        return tied[0]

    def _mean_b(alt: str) -> float:
        """Mean B-factor of alt's atoms."""

        return _mean(alt_groups[alt], "B_iso_or_equiv")

    tb = rules.tie_breaker
    if tb == "alphabetical_last":
        return max(tied)
    if tb == "lowest_b_factor":
        return min(tied, key=_mean_b)
    if tb == "highest_b_factor":
        return max(tied, key=_mean_b)
    return min(tied)  # alphabetical_first


def _resolve_altlocs(
    atoms: list[AtomSiteRecord],
    rules,
) -> tuple[list[AtomSiteRecord], AltlocSelectionMapping]:
    """Resolve altloc-disordered atoms to a single conformer per residue,
    per the altloc rules."""

    mapping = AltlocSelectionMapping()

    if rules.strategy == "preserve":
        return atoms, mapping

    # Keyed by physical site (chain + residue number), not residue
    # identity: mmCIF also uses label_alt_id for *compositional* disorder
    # (different comp_id occupying the same slot, e.g. two natural-product
    # variants modeled at one position), not just conformational disorder
    # of a single residue. Grouping by comp_id too would treat those as
    # unrelated single-altloc residues, leave both in the output, and
    # produce a residue-number collision downstream (validation.py) that
    # was never actually ambiguous — the file already resolved it via
    # label_alt_id, we just weren't reading it that way.
    by_site: dict[tuple, dict[str, list[AtomSiteRecord]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for a in atoms:
        if _has_altloc(a):
            by_site[_site(a)][a.label_alt_id].append(a)

    selected_altloc: dict[tuple, str] = {}
    for site, alt_groups in by_site.items():
        selected = _select_altloc(alt_groups, rules)
        selected_altloc[site] = selected

        if rules.record_selection and len(alt_groups) > 1:
            asym_id, _auth_seq_id, seq_id = site
            selected_comp_id = alt_groups[selected][0].label_comp_id
            mapping.items.append(
                AltlocSelectionMappingItem(
                    canonical_chain_id=asym_id,
                    residue_id=f"{seq_id}:{selected_comp_id}",
                    selected_altloc=selected,
                    available_altlocs=sorted(alt_groups),
                    selection_reason=_SELECTION_REASON.get(
                        rules.strategy, "first_alphabetical"
                    ),
                )
            )

    result: list[AtomSiteRecord] = []
    for a in atoms:
        if not _has_altloc(a):
            result.append(a)
        elif a.label_alt_id == selected_altloc.get(_site(a)):
            result.append(a.model_copy(update={"label_alt_id": None}))

    return result, mapping
