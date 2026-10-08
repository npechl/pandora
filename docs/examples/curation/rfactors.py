from pandora.datasets import curate_structure
from pandora.metadata import collect_metadata
from pandora.parsing import mmcif_to_structure
from pandora.schemas.dataset import DatasetCurationPolicy, QualityRules


def check(entry_id: str, quality_rules: QualityRules) -> str:
    structure, _, _ = mmcif_to_structure(f"datasets/dev/mmcif/{entry_id}.cif")
    policy = DatasetCurationPolicy(
        policy_id="c",
        policy_name="c",
        policy_version="1.0.0",
        quality_rules=quality_rules,
    )
    curated, exclusions, _ = curate_structure(
        structure, collect_metadata(structure), policy
    )
    if curated is not None:
        return "kept"
    return f"{exclusions[-1].reason_code} ({exclusions[-1].message})"


# 1a08 is X-ray. It reports Rmerge but neither Rsym nor Rfree.
structure, _, _ = mmcif_to_structure("datasets/dev/mmcif/1a08.cif")
quality = collect_metadata(structure).quality
print(
    f"1a08: r_free={quality.r_free} r_work={quality.r_work} "
    f"r_sym={quality.r_sym} r_merge={quality.r_merge}"
)

# Rsym falls back to Rmerge.
print("max_r_sym=0.05:", check("1a08", QualityRules(max_r_sym=0.05)))
# A missing Rfree is kept by default...
print("max_r_free=0.25:", check("1a08", QualityRules(max_r_free=0.25)))
# ...or excluded if you ask for it.
print(
    "max_r_free=0.25, exclude nulls:",
    check(
        "1a08",
        QualityRules(max_r_free=0.25, null_rfactor_behavior="exclude"),
    ),
)
# R-factor rules only check X-ray entries, so cryo-EM 22jy passes.
print(
    "22jy, same rule:",
    check(
        "22jy",
        QualityRules(max_r_free=0.25, null_rfactor_behavior="exclude"),
    ),
)
