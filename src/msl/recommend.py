"""목표 기반 추천 — 원소와 조건을 주면 조성 후보를 만들어 평가하고 순위를 매긴다.

1. 후보: 필수 원소 + 선택 원소로 만들 수 있는 정수 조성(화학식당 원자 수 상한)을 모두 나열
2. 거르기: 산화수로 전하 균형이 가능한 조성만 (혼합 원자가 허용, 금속끼리는 검사 안 함)
3. 평가: DB 에 있는 조성은 DFT 값(MP), 없는 조성은 L1 예측값(80% 구간)
4. 판정: 조건을 구간까지 만족하면 '충족', 구간이 경계에 걸치면 '가능성 있음', 구간이 조건 밖이면 '불충족'
5. 순위: 판정 → 목표(가장 안정·밴드갭 목표값에 가까움 등)
실행 결과는 data/recommendations/ 에 JSON 으로 남긴다.
"""

from __future__ import annotations

import datetime as dt
import itertools
import json
import math
import re
from functools import reduce
from pathlib import Path
from typing import Any, Callable, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from pymatgen.core import Composition, Element

from msl.env import DATA_DIR
from msl.schema.refs import ELEMENTS

RUNS = DATA_DIR / "recommendations"
Prop = Literal["ehull", "gap", "ef", "density", "p_metal"]
PROP_KO = {"ehull": "hull 거리", "gap": "밴드갭", "ef": "형성에너지", "density": "밀도", "p_metal": "금속일 확률"}
PROP_UNIT = {"ehull": "eV/atom", "gap": "eV", "ef": "eV/atom", "density": "g/cm³", "p_metal": ""}
TIER = {"충족": 0, "가능성 있음": 1, "불충족": 2}
SUSPICIOUS_BELOW_HULL = -0.1  # eV/atom — 알려진 hull 보다 이만큼 낮게 예측되면 모델 과대평가일 가능성이 크다


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Constraint(_Model):
    prop: Prop
    min: float | None = None
    max: float | None = None

    @model_validator(mode="after")
    def _one(self) -> Constraint:
        if self.min is None and self.max is None:
            raise ValueError("min 과 max 중 하나는 있어야 함")
        return self

    def text(self) -> str:
        u = PROP_UNIT[self.prop]
        if self.min is not None and self.max is not None:
            return f"{PROP_KO[self.prop]} {self.min:g}~{self.max:g} {u}".strip()
        return f"{PROP_KO[self.prop]} {'≥ ' + format(self.min, 'g') if self.min is not None else '≤ ' + format(self.max, 'g')} {u}".strip()


class Objective(_Model):
    prop: Prop = "ehull"
    direction: Literal["min", "max", "target"] = "min"
    value: float | None = None

    @model_validator(mode="after")
    def _target(self) -> Objective:
        if self.direction == "target" and self.value is None:
            raise ValueError("target 에는 value 가 필요함")
        return self

    def text(self) -> str:
        if self.direction == "target":
            return f"{PROP_KO[self.prop]} {self.value:g} {PROP_UNIT[self.prop]}에 가까운 순".replace("  ", " ")
        return f"{PROP_KO[self.prop]} {'낮은' if self.direction == 'min' else '높은'} 순"


class Goal(_Model):
    id: str = Field(pattern=r"^goal-[a-z0-9][a-z0-9\-]*$")
    name: str
    description: str | None = None
    required: list[str] = Field(default_factory=list)  # 모든 후보에 들어갈 원소
    optional: list[str] = Field(default_factory=list)  # 들어갈 수도 있는 원소
    max_elements: int = Field(default=4, ge=1, le=6)
    max_atoms: int = Field(default=10, ge=2, le=24)  # 화학식 단위당 원자 수 상한 (약분 후)
    charge_balanced: bool = True
    include_known: bool = True  # DB 에 있는 조성도 후보로
    include_new: bool = True  # DB 에 없는 조성 (L1 예측)
    constraints: list[Constraint] = Field(default_factory=list)
    objective: Objective = Objective()
    top_k: int = Field(default=30, ge=1, le=500)
    max_candidates: int = Field(default=6000, ge=10, le=50000)

    @field_validator("required", "optional")
    @classmethod
    def _elements(cls, v: list[str]) -> list[str]:
        bad = [e for e in v if e not in ELEMENTS]
        if bad:
            raise ValueError(f"원소 기호가 아님: {', '.join(bad)}")
        return list(dict.fromkeys(v))

    @model_validator(mode="after")
    def _pool(self) -> Goal:
        if not self.required and not self.optional:
            raise ValueError("원소를 하나 이상 적어야 함 (required 또는 optional)")
        if set(self.required) & set(self.optional):
            raise ValueError(f"필수와 선택에 같은 원소: {sorted(set(self.required) & set(self.optional))}")
        if len(self.required) > self.max_elements:
            raise ValueError(f"필수 원소({len(self.required)}개)가 max_elements({self.max_elements})보다 많음")
        if not (self.include_known or self.include_new):
            raise ValueError("include_known 과 include_new 중 하나는 참이어야 함")
        return self


# ── 후보 생성 ─────────────────────────────────────────────────────────────


def _states(el: Element) -> tuple[int, ...]:
    return tuple(sorted(set(el.icsd_oxidation_states) | set(el.common_oxidation_states)))


def charge_balanced(comp: Composition) -> bool:
    """산화수로 전하 균형이 가능한가 (혼합 원자가 허용).

    가장 전기음성도가 큰 원소는 음의 산화수(여럿이면 그 범위), 나머지 원소는 가진 산화수 전체 범위에서
    원자마다 달라도 된다고 보고 0 이 합의 범위 안에 드는지 본다. 음의 산화수가 없는 조합(금속끼리)은 검사하지 않는다.
    """
    els = sorted(comp.elements, key=lambda e: e.X if not math.isnan(e.X) else 0.0)
    anion = els[-1]
    neg = [s for s in _states(anion) if s < 0]
    if not neg:
        return True
    lo = comp[anion] * min(neg)
    hi = comp[anion] * max(neg)
    for e in els[:-1]:
        st = _states(e)
        if not st:
            return False
        lo += comp[e] * min(st)
        hi += comp[e] * max(st)
    return lo <= 0 <= hi


def candidates(goal: Goal) -> list[Composition]:
    """필수 원소를 모두 포함하는 원소 부분집합마다, 원자 수 상한 안의 모든 약분 정수 조성."""
    req, opt = goal.required, goal.optional
    seen: set[str] = set()
    out: list[Composition] = []
    for k in range(max(0, 1 - len(req)) if req else 1, len(opt) + 1):
        for extra in itertools.combinations(opt, k):
            els = [*req, *extra]
            if not els or len(els) > goal.max_elements:
                continue
            n = len(els)
            for total in range(n, goal.max_atoms + 1):
                for cut in itertools.combinations(range(1, total), n - 1):
                    counts = [b - a for a, b in zip((0, *cut), (*cut, total))]
                    if reduce(math.gcd, counts) != 1:
                        continue  # 약분되는 조성은 더 작은 total 에서 이미 나옴
                    comp = Composition(dict(zip(els, counts)))
                    key = comp.reduced_formula
                    if key in seen:
                        continue
                    seen.add(key)
                    out.append(comp)
    return out


# ── 평가 · 판정 ───────────────────────────────────────────────────────────


def _values(rec: dict[str, Any]) -> tuple[dict[str, float | None], dict[str, tuple[float, float] | None], str]:
    """(점 값, 구간, 출처). DB 에 있으면 DFT 값(구간 없음), 없으면 L1 예측."""
    k = rec["known"]
    if k is not None:
        v = {"ehull": k["ehull"], "gap": k["gap"], "ef": k["ef"], "density": k["density"], "p_metal": 1.0 if k["gap"] <= 0 else 0.0}
        return v, {p: None for p in v}, "DFT"
    v = {"ehull": rec["ehull"], "gap": rec["gap"], "ef": rec["ef"], "density": rec["density"], "p_metal": rec["p_metal"]}
    iv = {"ehull": tuple(rec["ehull_interval"]) if rec["ehull_interval"] else None, "gap": tuple(rec["gap_interval"]),
          "ef": tuple(rec["ef_interval"]), "density": None, "p_metal": None}
    return v, iv, "L1"


def _judge(c: Constraint, v: float | None, iv: tuple[float, float] | None) -> str:
    if v is None:
        return "불충족"
    lo, hi = iv if iv else (v, v)
    mn = -math.inf if c.min is None else c.min
    mx = math.inf if c.max is None else c.max
    if mn <= lo and hi <= mx:
        return "충족"
    if hi < mn or lo > mx:
        return "불충족"
    return "가능성 있음"


def _objective(o: Objective, v: dict[str, float | None]) -> float:
    x = v.get(o.prop)
    if x is None:
        return math.inf
    return {"min": x, "max": -x, "target": abs(x - (o.value or 0.0))}[o.direction]


def recommend(goal: Goal, progress: Callable[[str], None] | None = None) -> dict[str, Any]:
    from msl.ml.predict import info, predict_batch

    say = progress or (lambda s: None)
    comps = candidates(goal)
    n_all = len(comps)
    if goal.charge_balanced:
        comps = [c for c in comps if len(c.elements) == 1 or charge_balanced(c)]
    n_balanced = len(comps)
    if len(comps) > goal.max_candidates:
        raise ValueError(f"후보 {len(comps):,}개가 max_candidates {goal.max_candidates:,} 를 넘음 — 원소·원자 수를 줄이거나 한도를 올릴 것")
    say(f"후보 {n_all:,}개 → 전하 균형 {n_balanced:,}개 · 평가 중")
    recs = predict_batch([c.formula.replace(" ", "") for c in comps])
    rows = []
    for rec in recs:
        is_known = rec["known"] is not None
        if (is_known and not goal.include_known) or (not is_known and not goal.include_new):
            continue
        v, iv, src = _values(rec)
        checks = [(c, _judge(c, v[c.prop], iv[c.prop])) for c in goal.constraints]
        status = max((j for _, j in checks), key=lambda j: TIER[j], default="충족")
        flags = []
        if not is_known and v["ehull"] is not None and v["ehull"] < SUSPICIOUS_BELOW_HULL:
            flags.append(f"알려진 상보다 {-v['ehull']:.2f} eV/atom 낮게 예측 — 과대평가 의심 (새 바닥 상태가 이만큼 낮은 경우는 드묾)")
        if not is_known and not rec["in_domain"]:
            flags.append("학습 범위 밖 — 신뢰 낮음")
        rows.append({
            "formula": rec["formula"], "source": src, "material_id": rec["known"]["material_id"] if is_known else None,
            "status": status, "flags": flags, "failed": [c.text() for c, j in checks if j == "불충족"],
            "uncertain": [c.text() for c, j in checks if j == "가능성 있음"],
            "values": {p: None if x is None else round(float(x), 4) for p, x in v.items()},
            "intervals": {p: None if x is None else [round(float(x[0]), 4), round(float(x[1]), 4)] for p, x in iv.items()},
            "in_domain": rec["in_domain"] or is_known, "objective": _objective(goal.objective, v),
            "theoretical": rec["known"]["theoretical"] if is_known else None,
        })
    # 같은 판정 안에서는 DB(DFT)로 확인된 후보 먼저, 그다음 의심 표시 없는 L1 예측 — 예측값이 확인된 물질보다 위에 오지 않게
    rows.sort(key=lambda r: (TIER[r["status"]], r["source"] != "DFT", bool(r["flags"]), r["objective"]))
    top = [r for r in rows if r["status"] != "불충족"][: goal.top_k]
    _decompose(top)
    counts = {"candidates": n_all, "charge_balanced": n_balanced, "evaluated": len(rows),
              "known": sum(r["source"] == "DFT" for r in rows), "new": sum(r["source"] == "L1" for r in rows),
              "flagged": sum(bool(r["flags"]) for r in rows),
              "out_of_domain": sum(any("범위 밖" in f for f in r["flags"]) for r in rows),
              "below_hull_suspicious": sum(any("과대평가" in f for f in r["flags"]) for r in rows),
              **{k: sum(r["status"] == k for r in rows) for k in TIER}}
    m = info()
    return {"goal": goal.model_dump(mode="json"), "created": dt.datetime.now().isoformat(timespec="seconds"),
            "constraints_text": [c.text() for c in goal.constraints], "objective_text": goal.objective.text(),
            "counts": counts, "results": top,
            "scatter": [{"f": r["formula"], "e": r["values"]["ehull"], "g": r["values"]["gap"], "s": r["status"], "src": r["source"]}
                        for r in rows if r["values"]["ehull"] is not None][:3000],
            "model": {"version": m["version"], "data": m["data"]["version"], "ef_mae": m["metrics"]["chemsys"]["ef_mae"],
                      "stability_accuracy": m["metrics"]["random"]["stability_accuracy"]}}


def _decompose(rows: list[dict[str, Any]]) -> None:
    """상위 후보 중 hull 위인 새 조성에 분해 생성물을 붙인다."""
    from msl.ml.predict import _pd_for

    for r in rows:
        if r["source"] != "L1" or (r["values"]["ehull"] or 0) <= 0:
            continue
        c = Composition(r["formula"])
        try:
            dec = _pd_for(tuple(sorted(str(e) for e in c.elements))).get_decomposition(c)
            r["decomposition"] = {e.name: round(float(x), 3) for e, x in sorted(dec.items(), key=lambda p: -p[1])}
        except Exception:
            pass


def save_run(result: dict[str, Any]) -> Path:
    RUNS.mkdir(parents=True, exist_ok=True)
    stamp = re.sub(r"[^0-9]", "", result["created"])[:14]
    path = RUNS / f"{stamp}-{result['goal']['id']}.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=1, default=float), encoding="utf-8")
    return path


def list_runs(limit: int = 20) -> list[dict[str, Any]]:
    if not RUNS.exists():
        return []
    out = []
    for p in sorted(RUNS.glob("*.json"), reverse=True)[:limit]:
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
            out.append({"file": p.name, "name": d["goal"]["name"], "created": d["created"], "counts": d["counts"]})
        except Exception:
            continue
    return out


def to_csv(result: dict[str, Any]) -> str:
    import csv
    import io

    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["순위", "화학식", "출처", "MP ID", "판정", "hull 거리 (eV/atom)", "hull 구간", "밴드갭 (eV)", "밴드갭 구간",
                "형성에너지 (eV/atom)", "밀도 (g/cm³)", "금속일 확률", "학습 범위 안", "불확실 조건", "주의", "분해 생성물"])
    for i, r in enumerate(result["results"], 1):
        v, iv = r["values"], r["intervals"]
        fmt = lambda x: "" if x is None else f"{x[0]:.3f}~{x[1]:.3f}"
        w.writerow([i, r["formula"], r["source"], r["material_id"] or "", r["status"], v["ehull"], fmt(iv["ehull"]), v["gap"],
                    fmt(iv["gap"]), v["ef"], v["density"], v["p_metal"], "예" if r["in_domain"] else "아니오",
                    "; ".join(r["uncertain"]), "; ".join(r["flags"]), " + ".join(r.get("decomposition", {}) or {})])
    return buf.getvalue()
