"""레시피 처리 — 도메인 라우터 (계획서 3.3절 '라우팅 규칙 v0')."""

from __future__ import annotations

from msl.resolve import Resolved
from msl.schema.recipe import Mode, Recipe, State
from msl.schema.refs import Namespace

METALS_EXCLUDED = {"H", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar", "As", "Se", "Br", "Kr", "Te", "I", "Xe", "At", "Rn"}


def route(recipe: Recipe, comps: list[Resolved]) -> list[tuple[str, str]]:
    """적용할 시험과 그 이유. assays 를 직접 지정했으면 그대로 따른다."""
    if recipe.assays != "auto":
        return [(code, "레시피에 직접 지정") for code in recipe.assays]
    if recipe.mode is Mode.SYSTEM:
        return [("A1", "원소 계 탐색"), ("A7", "안정상의 물성 조회")]

    out: dict[str, str] = {"S0": "항상 — 모든 혼합은 안전 판정부터"}
    materials = [c for c in comps if c.props]
    if materials and len(materials) == len(comps):
        out["A8"] = "사용자 정의 소재끼리의 블렌드"
    water = any(c.is_water for c in comps) or any(c.state is State.AQUEOUS for c in comps)
    gas = any(c.state is State.GAS for c in comps)
    solids = [c for c in comps if not c.props and c.state in (None, State.SOLID) and not c.is_water]
    if water:
        out["A4"] = "물 또는 수용액 성분이 있음"
    if gas:
        out["A3"] = "기체 성분이 있음"
    if solids and len(solids) == len(comps):
        out["A1"] = "고체 성분만 있음 → 화학계 안정성"
        if len(comps) == 2:
            out["A2"] = "고체 두 성분 → 계면 반응"
        if all(c.ref.namespace is Namespace.ELEMENT and c.ref.key not in METALS_EXCLUDED for c in comps):
            out["A6"] = "금속 원소만 있음 → 합금 상평형"
            out["A7"] = "금속 원소의 물성"
    if any(c.cas for c in comps):
        out["A10"] = "CAS 번호가 있는 성분 → 국내 규제 확인"
    out["A9"] = "항상 — 원료 공급 위험"
    order = ["S0", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10"]
    return [(c, out[c]) for c in order if c in out]
