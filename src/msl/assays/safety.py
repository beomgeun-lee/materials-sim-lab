"""S0 혼합 안전 게이트 — CAMEO 반응성 그룹(PubChem) × 자체 규칙 v0 (kb/safety_rules.yaml)."""

from __future__ import annotations

import re
from functools import cache
from itertools import combinations
from typing import Any

import yaml

from msl.assays.base import Context, Outcome, not_applicable, val
from msl.env import KB_DIR
from msl.resolve import Resolved
from msl.schema.recipe import Mode
from msl.schema.result import Fidelity, Status, ValueKind

VERDICT = {"incompatible": ("부적합", "bad", 2), "caution": ("주의", "warn", 1)}


@cache
def rules() -> list[dict[str, Any]]:
    return yaml.safe_load((KB_DIR / "safety_rules.yaml").read_text(encoding="utf-8")) or []


def _match(cond: dict[str, Any], comp: Resolved) -> bool:
    def group_ok(pattern: str) -> bool:
        if pattern.endswith("*"):
            return any(g.startswith(pattern[:-1]) for g in comp.groups)
        return pattern in comp.groups

    if not any(group_ok(p) for p in cond.get("groups", [])):
        return False
    formula = cond.get("formula")
    return formula is None or bool(comp.formula and re.search(formula, comp.formula))


def run(ctx: Context) -> Outcome:
    if ctx.recipe.mode is Mode.SYSTEM:
        return not_applicable("원소 계 탐색 — 실제로 섞는 레시피가 아니라 혼합 판정 대상이 아님", "cameo-reactivity")

    comps = [c for c in ctx.comps if c.error is None]
    hits: list[dict[str, Any]] = []
    for a, b in combinations(comps, 2):
        for rule in rules():
            for x, y in ((a, b), (b, a)):
                if _match(rule["a"], x) and _match(rule["b"], y):
                    hits.append({
                        "rule": rule["id"], "name": rule["name"], "a": x.name, "b": y.name,
                        "verdict": VERDICT[rule["verdict"]][0], "code": VERDICT[rule["verdict"]][1],
                        "hazards": rule["hazards"], "gases": rule["gases"], "note": rule["note"],
                    })
                    break

    known = [c for c in comps if c.groups]
    if hits:
        worst = max(hits, key=lambda h: 2 if h["code"] == "bad" else 1)
        verdict, code = worst["verdict"], worst["code"]
        title = worst["name"] + (" — " + ", ".join(worst["gases"]) + " 발생" if worst["gases"] else "")
    elif len(known) < 2:
        verdict, code = "판정 불가", "na"
        title = "반응성 그룹 정보가 있는 성분이 2개 미만"
    else:
        verdict, code = "규칙 해당 없음", "ok"
        title = f"자체 규칙 v0 {len(rules())}개 중 해당하는 위험 조합 없음"

    caveats = [f"NOAA 쌍별 판정표 확보 전의 자체 규칙 v0({len(rules())}개) — 목록에 없는 조합은 '해당 없음'일 뿐 안전 보장이 아님"]
    if len(comps) >= 3:
        caveats.append("2성분 쌍만 평가 — 3성분 이상의 상호작용은 평가하지 않음")
    missing = [c.name for c in ctx.comps if not c.groups]
    if missing:
        caveats.append("반응성 그룹 정보 없음: " + ", ".join(missing))

    status = {"bad": Status.WARNING, "warn": Status.WARNING}.get(code, Status.OK)
    return Outcome(
        status=status, fidelity=Fidelity.L0, engine="cameo-reactivity", engine_version="자체 규칙 v0",
        conditions_basis="상온·상압, 두 성분을 비슷한 양으로 섞는다고 가정 (CAMEO 전제)",
        values=[val("판정", verdict, kind=ValueKind.FLAG), val("해당 규칙 수", float(len(hits)))],
        sources=[("cameo", None), ("pubchem", None)],
        caveats=caveats,
        data={
            "verdict": verdict, "verdict_code": code, "title": title, "pairs": hits,
            "groups": [{"name": c.name, "groups": c.groups} for c in ctx.comps],
            "rules_total": len(rules()),
        },
        summary=title,
    )
