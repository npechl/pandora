from collections import defaultdict

from pandora.schemas.structure import AtomSiteRecord

from pandora.schemas.canonicalisation import (
    AltlocSelectionMapping,
    AltlocSelectionMappingItem,
)


def _resolve_altlocs(
    atoms: list[AtomSiteRecord],
    rules,
) -> tuple[list[AtomSiteRecord], AltlocSelectionMapping]:
    mapping = AltlocSelectionMapping()
    strategy = rules.strategy

    if strategy == "preserve":
        return atoms, mapping

    by_residue: dict[tuple, dict[str, list[AtomSiteRecord]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for a in atoms:
        alt = a.label_alt_id
        if alt and alt not in (".", "?"):
            key = (
                a.label_asym_id,
                a.auth_seq_id,
                a.label_seq_id,
                a.label_comp_id,
            )
            by_residue[key][alt].append(a)

    selected_altloc: dict[tuple, str] = {}

    for (
        asym_id,
        _auth_seq_id,
        seq_id,
        comp_id,
    ), alt_groups in by_residue.items():
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

        selected_altloc[(asym_id, _auth_seq_id, seq_id, comp_id)] = selected

        if rules.record_selection and len(altlocs) > 1:
            reason_map = {
                "select_first": "first_alphabetical",
                "select_user_defined": "user_defined",
                "select_best_occupancy": "best_occupancy",
            }
            mapping.items.append(
                AltlocSelectionMappingItem(
                    canonical_chain_id=asym_id,
                    residue_id=f"{seq_id}:{comp_id}",
                    selected_altloc=selected,
                    available_altlocs=altlocs,
                    selection_reason=reason_map.get(
                        strategy, "first_alphabetical"
                    ),
                )
            )

    result: list[AtomSiteRecord] = []
    for a in atoms:
        alt = a.label_alt_id
        if not alt or alt in (".", "?"):
            result.append(a)
            continue
        key = (a.label_asym_id, a.auth_seq_id, a.label_seq_id, a.label_comp_id)
        if alt == selected_altloc.get(key):
            result.append(a.model_copy(update={"label_alt_id": None}))

    return result, mapping
