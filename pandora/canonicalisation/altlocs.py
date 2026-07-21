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

    no_alt: list[AtomSiteRecord] = []
    by_residue: dict[tuple, dict[str, list[AtomSiteRecord]]] = defaultdict(
        lambda: defaultdict(list)
    )

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

    return result, mapping
