"""Mechanical scans for review-pandora. Prints leads, not verdicts.

Run from the repo root:
    uv run python .claude/skills/review-pandora/scripts/scan.py [--limit N]
"""

from __future__ import annotations

import argparse
import importlib
import json
import re
import tempfile
import typing
from collections import defaultdict
from pathlib import Path

import gemmi
import yaml
from pydantic import BaseModel

from pandora.canonicalisation import canonicalise_structure
from pandora.export import structure_to_mmcif
from pandora.parsing import mmcif_to_structure
from pandora.schemas.canonicalisation import canonicalisationPolicy
from pandora.schemas.dataset import DatasetCurationPolicy
from pandora.schemas.similarity import HitFilter

ROOT = Path.cwd()
FIXTURES = sorted((ROOT / "datasets/dev/mmcif").glob("*/*.cif"))
POLICIES = [canonicalisationPolicy, DatasetCurationPolicy, HitFilter]
POLICY_YAMLS = {
    "canonicalisation.yaml": canonicalisationPolicy,
    "dataset.yaml": DatasetCurationPolicy,
}
# mmCIF category -> how many typed records the parser made from it.
TYPED = {
    "_struct_conn.": lambda s: len(s.connections),
    "_entity.": lambda s: len(s.entities),
    "_entity_poly.": lambda s: sum(e.poly is not None for e in s.entities),
    "_pdbx_struct_assembly.": lambda s: len(s.assemblies),
    "_struct_conf.": lambda s: len(s.secondary_structure.conf_records),
    "_struct_sheet_range.": lambda s: len(s.secondary_structure.sheet_strands),
}
_BASE = {"policy_id": "scan", "policy_name": "scan", "policy_version": "0"}
_METADATA = {*_BASE, "description"}


def _submodel(annotation) -> type[BaseModel] | None:
    """The pydantic model inside a field annotation, if any."""

    for t in (annotation, *typing.get_args(annotation)):
        if isinstance(t, type) and issubclass(t, BaseModel):
            return t
    return None


def _fields(model: type[BaseModel], pre: str = ""):
    """(dotted path, FieldInfo) for every field, nested models included."""

    for name, f in model.model_fields.items():
        yield pre + name, f
        sub = _submodel(f.annotation)
        if sub is not None:
            yield from _fields(sub, f"{pre}{name}.")


def _section(title: str, lines: list[str]) -> None:
    """Print one report section."""

    print(f"\n## {title} ({len(lines)})")
    for line in lines or ["none"]:
        print(f"- {line}")


def unread_policy_fields() -> list[str]:
    """Policy fields no module outside pandora/schemas reads."""

    code = "\n".join(
        p.read_text()
        for p in (ROOT / "pandora").rglob("*.py")
        if "schemas" not in p.parts
    )
    out = []
    for model in POLICIES:
        for path, f in _fields(model):
            if _submodel(f.annotation) or path in _METADATA:
                continue
            name = path.split(".")[-1]
            if not re.search(rf"\b{name}\b", code):
                out.append(f"{model.__name__}.{path}")
    return out


def undocumented_policy_fields() -> list[str]:
    """Policy fields docs/reference/policies.md never mentions."""

    doc = (ROOT / "docs/reference/policies.md").read_text()
    return [
        f"{model.__name__}.{path}"
        for model in POLICIES
        for path, f in _fields(model)
        if not _submodel(f.annotation)
        and path not in _METADATA
        and not re.search(rf"`[\w.]*\b{path.split('.')[-1]}`", doc)
    ]


def policy_yaml_drift() -> list[str]:
    """Field paths that differ between docs/policies/*.yaml and models."""

    def paths(d: dict, pre: str = "") -> set[str]:
        out = set()
        for k, v in d.items():
            out.add(pre + k)
            if isinstance(v, dict):
                out |= paths(v, f"{pre}{k}.")
        return out

    out = []
    docs = ROOT / "docs/policies"
    for p in sorted(docs.glob("*.yaml")):
        model = POLICY_YAMLS.get(p.name)
        if model is None:
            out.append(f"{p.name}: no matching policy model")
            continue
        y = paths(next(iter(yaml.safe_load(p.read_text()).values())))
        m = {path for path, _ in _fields(model)}
        out += [f"{p.name}: `{x}` not in model" for x in sorted(y - m)]
        out += [f"{p.name}: `{x}` missing" for x in sorted(m - y)]
    return out


def unreferenced_packages() -> list[str]:
    """Packages with a public API missing from docs/reference/functions.md."""

    doc = (ROOT / "docs/reference/functions.md").read_text()
    out = []
    for init in sorted((ROOT / "pandora").glob("*/__init__.py")):
        pkg = f"pandora.{init.parent.name}"
        if pkg == "pandora.schemas" or f"::: {pkg}\n" in doc + "\n":
            continue
        if getattr(importlib.import_module(pkg), "__all__", None):
            out.append(f"{pkg} exports names but has no `::: {pkg}`")
    return out


def parser_coverage(fixtures: list[Path]) -> list[str]:
    """Typed categories whose parsed record count differs from the file."""

    bad: dict[str, list[str]] = defaultdict(list)
    for p in fixtures:
        block = gemmi.cif.read(str(p)).sole_block()
        s, _, _ = mmcif_to_structure(str(p))
        for cat, count in TYPED.items():
            rows = len(block.find_mmcif_category(cat))
            if rows and count(s) == 0:
                bad[cat].append(f"{p.stem} ({rows} rows -> 0)")
    return [
        f"`{cat}`: {len(v)} fixtures, e.g. {', '.join(v[:3])}"
        for cat, v in bad.items()
    ]


def _variants() -> list[tuple[str, object]]:
    """(field path, value) for every non-default Literal or bool value of
    the canonicalisation policy."""

    out = []
    for path, f in _fields(canonicalisationPolicy):
        if path in _METADATA:
            continue
        args = typing.get_args(f.annotation)
        if f.annotation is bool:
            out.append((path, not f.default))
        elif args and all(isinstance(a, str) for a in args):
            out += [(path, a) for a in args if a != f.default]
    return out


def _policy(path: str, value: object) -> canonicalisationPolicy:
    """Default policy with one field set."""

    data = canonicalisationPolicy(**_BASE).model_dump()
    *parents, leaf = path.split(".")
    node = data
    for key in parents:
        node = node[key]
    node[leaf] = value
    return canonicalisationPolicy.model_validate(data)


def _invariants(c) -> list[str]:
    """Names of the referential checks the canonical structure fails."""

    asyms = {u.id for u in c.asym_units}
    atom_asyms = {a.label_asym_id for a in c.atoms}
    ss = c.secondary_structure
    spans = [*ss.conf_records, *ss.sheet_strands]
    checks = {
        "no atoms left": not c.atoms,
        "asym unit without atoms": bool(asyms - atom_asyms),
        "atoms in unlisted asym": bool(atom_asyms - asyms),
        "entity without atoms": bool(
            {e.id for e in c.entities} - {a.label_entity_id for a in c.atoms}
        ),
        "assembly names missing asym": any(
            x not in asyms
            for a in c.assemblies
            for g in a.generators
            for x in g.asym_id_list
        ),
        "bond names missing asym": any(
            {b.ptnr1.label_asym_id, b.ptnr2.label_asym_id} - asyms
            for b in c.connections
        ),
        "SS span names missing asym": any(
            {s.beg_label_asym_id, s.end_label_asym_id} - asyms for s in spans
        ),
    }
    return [name for name, failed in checks.items() if failed]


def _fingerprint(canonical, mappings) -> str:
    """Comparable dump of a canonicalisation result."""

    return json.dumps(
        [canonical.model_dump(mode="json"), mappings.model_dump(mode="json")],
        sort_keys=True,
    )


def _export_failures(c, path: Path) -> list[str]:
    """Names of the export round-trip checks c fails."""

    if not c.atoms:
        return []
    r, _, status = mmcif_to_structure(str(structure_to_mmcif(c, path)))
    if status == "failed" or len(r.atoms) != len(c.atoms):
        return ["export does not re-parse"]
    if len(r.connections) != len(c.connections):
        return ["export loses bonds"]
    return []


def canonicalisation_sweep(fixtures: list[Path]) -> tuple[list, list]:
    """Run each one-field policy variant over the fixtures; collect
    invariant failures and variants that never change the output."""

    tmp = Path(tempfile.mkdtemp())
    parsed = {p.stem: mmcif_to_structure(str(p))[0] for p in fixtures}
    default = canonicalisationPolicy(**_BASE)
    baseline = {}
    for eid, s in parsed.items():
        try:
            baseline[eid] = _fingerprint(
                *canonicalise_structure(s, default)[:2]
            )
        except ValueError:
            baseline[eid] = None

    hits: dict[tuple[str, str], list[str]] = defaultdict(list)
    no_effect = []
    for path, value in _variants():
        label = f"{path}={value}"
        policy = _policy(path, value)
        changed = False
        for eid, s in parsed.items():
            try:
                c, m, _ = canonicalise_structure(s, policy)
            except ValueError as exc:
                hits[(label, f"raises {type(exc).__name__}")].append(eid)
                changed = True
                continue
            changed |= _fingerprint(c, m) != baseline[eid]
            failed = _invariants(c) + _export_failures(c, tmp / f"{eid}.cif")
            for name in failed:
                hits[(label, name)].append(eid)
        if not changed:
            no_effect.append(label)

    failures = [
        f"{label}: {name} in {len(eids)} entries, e.g. {', '.join(eids[:3])}"
        for (label, name), eids in sorted(hits.items())
    ]
    return failures, no_effect


def main() -> None:
    """Run every scan and print a Markdown report."""

    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, help="only the first N fixtures")
    ap.add_argument("--skip-sweep", action="store_true")
    args = ap.parse_args()
    fixtures = FIXTURES[: args.limit] if args.limit else FIXTURES

    print(f"# review-pandora scan ({len(fixtures)} fixtures)")
    _section("Policy fields no logic reads", unread_policy_fields())
    _section(
        "Policy fields missing from policies.md", undocumented_policy_fields()
    )
    _section("docs/policies YAML drift", policy_yaml_drift())
    _section("Packages missing from functions.md", unreferenced_packages())
    _section("Typed categories parsed as empty", parser_coverage(fixtures))
    if not args.skip_sweep:
        failures, no_effect = canonicalisation_sweep(fixtures)
        _section("Canonicalisation invariant failures", failures)
        _section("Policy values that never change the output", no_effect)


if __name__ == "__main__":
    main()
