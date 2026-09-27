"""검증 세트 실행기 — 결과가 알려진 조합으로 시험 전체를 점검한다 (계획서 7·8절).

kb/validation/assays_v0.yaml 의 각 사례는 레시피 하나와 시험 코드 하나, 기대 조건(checks)을 갖는다.
기준: 안전 게이트(S0) 위험 사례 재현율 100%, 나머지 시험 통과율 90% 이상.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml
from pymatgen.core import Composition

from msl.env import KB_DIR
from msl.registry.load import Registry
from msl.runtime.runner import run_recipe
from msl.schema.recipe import Recipe
from msl.schema.result import AssayResult

DEFAULT_SET = KB_DIR / "validation" / "assays_v0.yaml"
HAZARD = {"부적합", "주의"}


@dataclass
class CaseResult:
    id: str
    assay: str
    passed: bool
    hazard_expected: bool = False
    detail: list[str] = field(default_factory=list)
    summary: str = ""


def _value(res: AssayResult, name: str):
    for v in res.values:
        if v.name == name or v.name.startswith(name):
            return v
    return None


def _formulas(side: str) -> list[Composition]:
    """'0.5 Zn(FeO2)2 + 0.5 CO2' → 각 화학종의 축약 조성 (표기 순서와 무관하게 비교하려고)."""
    out = []
    for term in side.split(" + "):
        term = re.sub(r"^\s*[\d.]+\s+", "", term.strip())
        try:
            out.append(Composition(term).reduced_composition)
        except Exception:
            continue
    return out


def _check(res: AssayResult, check: dict[str, Any]) -> tuple[bool, str]:
    (kind, arg), = check.items()
    data = res.data or {}
    if kind == "status_in":
        return res.status.value in arg, f"status {res.status.value} ∈ {arg}"
    if kind == "verdict_in":
        return data.get("verdict") in arg, f"판정 {data.get('verdict')} ∈ {arg}"
    if kind == "value_between":
        v = _value(res, arg["name"])
        ok = v is not None and isinstance(v.value, (int, float)) and arg["range"][0] <= v.value <= arg["range"][1]
        return ok, f"{arg['name']} = {v.value if v else '없음'} ∈ {arg['range']}"
    if kind == "value_contains":
        v = _value(res, arg["name"])
        ok = v is not None and arg["text"] in str(v.value)
        return ok, f"{arg['name']} = {v.value if v else '없음'} ⊇ '{arg['text']}'"
    if kind == "summary_contains":
        return arg in data.get("summary", ""), f"요약 ⊇ '{arg}'"
    if kind == "table_contains":
        cells = [str(c) for row in (data.get("table") or {}).get("rows", []) for c in row]
        return any(arg in c for c in cells), f"표 ⊇ '{arg}'"
    if kind in ("curve_min_below", "curve_min_above"):
        curves = data.get("curves") or []
        if not curves:
            return False, "곡선 없음"
        low = min(k["e"] for k in curves[-1]["kinks"])
        ok = low < arg if kind == "curve_min_below" else low >= arg
        return ok, f"{curves[-1]['label']} 최저 {low} meV/atom {'<' if kind == 'curve_min_below' else '≥'} {arg}"
    if kind == "product_contains":
        curves = data.get("curves") or []
        if not curves:
            return False, "곡선 없음"
        low = min(curves[-1]["kinks"], key=lambda k: k["e"])
        products = low["rxn"].split("->")[-1]
        options = arg if isinstance(arg, list) else [arg]
        found = _formulas(products)
        want = {Composition(o).reduced_composition for o in options}
        ok = any(any(f.almost_equals(w) for w in want) for f in found)
        return ok, f"생성물 '{products.strip()}' ⊇ {options} (조성 비교)"
    if kind == "recipe_products_contains":  # A2: 레시피 혼합비에서의 평형 생성물
        found = [Composition(f).reduced_composition for f in data.get("products_at_recipe") or []]
        options = arg if isinstance(arg, list) else [arg]
        want = [Composition(o).reduced_composition for o in options]
        ok = any(f.almost_equals(w) for f in found for w in want)
        return ok, f"레시피 혼합비 생성물 {data.get('products_at_recipe')} ⊇ {options}"
    return False, f"알 수 없는 검사: {kind}"


def run_set(path: Path = DEFAULT_SET, registry: Registry | None = None, use_cache: bool = True) -> list[CaseResult]:
    cases = yaml.safe_load(path.read_text(encoding="utf-8"))
    out: list[CaseResult] = []
    for case in cases:
        recipe = Recipe.model_validate({**case["recipe"], "assays": [case["assay"]]})
        report = run_recipe(recipe, use_cache=use_cache, registry=registry)
        res = next((r for r in report.results if r.assay == case["assay"]), None)
        cr = CaseResult(case["id"], case["assay"], False, hazard_expected=bool(case.get("hazard")))
        if res is None:
            cr.detail.append("결과 없음")
        else:
            cr.summary = (res.data or {}).get("summary", "")
            checks = [_check(res, c) for c in case["checks"]]
            cr.passed = all(ok for ok, _ in checks)
            cr.detail = [("✓ " if ok else "✗ ") + msg for ok, msg in checks]
        out.append(cr)
    return out


def score(results: list[CaseResult]) -> dict[str, Any]:
    by_assay: dict[str, list[CaseResult]] = {}
    for r in results:
        by_assay.setdefault(r.assay, []).append(r)
    hazards = [r for r in results if r.assay == "S0" and r.hazard_expected]
    others = [r for r in results if r.assay != "S0" or not r.hazard_expected]
    recall = sum(r.passed for r in hazards) / len(hazards) if hazards else 1.0
    other_rate = sum(r.passed for r in others) / len(others) if others else 1.0
    return {
        "by_assay": {a: (sum(r.passed for r in rs), len(rs)) for a, rs in sorted(by_assay.items())},
        "s0_recall": recall, "other_rate": other_rate,
        "passed": recall == 1.0 and other_rate >= 0.9,
    }
