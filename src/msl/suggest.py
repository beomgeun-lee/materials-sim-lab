"""다음에 시험할 조합 — 베이지안 최적화 (BayBE, 계획서 4단계, D25).

목표 기반 추천(D18)은 모든 후보를 L1 로 한 번에 순위 매긴다. 여기서는 평가 비용이 큰 단계(L2 uMLIP·DFT·실험)를 위해
'지금까지 잰 값'을 보고 다음에 잴 조성 몇 개를 고른다.

- 후보 풀: 목표(goal) 추천 상위 pool 개, 또는 화학식 목록.
- 조성 표현: L1 특성 112개 + 목표 물성의 사전값(DB 는 DFT, 새 조성은 L1 예측). BayBE 가 상관 높은 열을 줄인다.
- 측정값: data/campaigns/{id}.csv (formula, value, source, note, added). 저장된 L2 결과는 import_l2 로 가져온다.
- 측정값이 없으면 BayBE 초기 추천(특성 공간에서 고르게), 있으면 가우시안 과정 대리모델 + 획득 함수로 고른다.
- 서로 다른 충실도(L2·DFT·실험)의 값을 한 목표에 섞으면 경고한다 — 같은 척도가 아니다 (D19).
"""

from __future__ import annotations

import datetime as dt
import csv
import warnings
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import Field, model_validator
from pymatgen.core import Composition

from msl.env import DATA_DIR, ROOT
from msl.recommend import PROP_KO, PROP_UNIT, Goal, Objective, _Model

STORE = DATA_DIR / "campaigns"
GOALS = ROOT / "examples" / "goals"
Source = Literal["L2", "DFT", "실험", "기타"]


class CampaignSpec(_Model):
    id: str = Field(pattern=r"^camp-[a-z0-9][a-z0-9\-]*$")
    name: str
    description: str | None = None
    goal: str | None = None  # examples/goals 의 목표 id (goal-…) — 이 목표의 추천 상위 pool 개가 후보
    formulas: list[str] = Field(default_factory=list)  # 또는 후보를 직접
    pool: int = Field(default=300, ge=10, le=500)
    target: Objective = Objective()
    sigma: float | None = Field(default=None, gt=0)  # target 방향일 때 허용 폭 (기본: |값|의 10%, 최소 0.05)
    batch: int = Field(default=5, ge=1, le=50)
    seed: int = 0

    @model_validator(mode="after")
    def _one(self) -> CampaignSpec:
        if bool(self.goal) == bool(self.formulas):
            raise ValueError("goal 과 formulas 중 하나만 적어야 함")
        return self


def load_campaign(path: Path) -> CampaignSpec:
    return CampaignSpec.model_validate(yaml.safe_load(path.read_text(encoding="utf-8")))


def _key(formula: str) -> str:
    from msl.ml.predict import formula_key

    return formula_key(Composition(formula))


# ── 후보 풀 ───────────────────────────────────────────────────────────────


def _goal(goal_id: str) -> Goal:
    for p in sorted(GOALS.glob("*.yaml")):
        d = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
        if d.get("id") == goal_id:
            return Goal.model_validate(d)
    raise ValueError(f"목표를 찾지 못함: {goal_id} (examples/goals)")


def pool(spec: CampaignSpec) -> list[dict[str, Any]]:
    """후보 — [{formula, source(DFT|L1), prior(목표 물성 사전값)}]."""
    prop = spec.target.prop
    if spec.goal:
        from msl.recommend import recommend

        g = _goal(spec.goal)
        res = recommend(g.model_copy(update={"top_k": spec.pool}))
        return [{"formula": r["formula"], "source": r["source"], "prior": r["values"][prop]} for r in res["results"]]
    from msl.ml.predict import predict_batch
    from msl.recommend import _values

    out = []
    for rec in predict_batch(list(dict.fromkeys(_key(f) for f in spec.formulas))):
        v, _, src = _values(rec)
        out.append({"formula": rec["formula"], "source": src, "prior": v[prop]})
    return out


# ── 측정값 ────────────────────────────────────────────────────────────────

FIELDS = ["formula", "value", "source", "note", "added"]


def _path(spec: CampaignSpec) -> Path:
    return STORE / f"{spec.id}.csv"


def measurements(spec: CampaignSpec) -> list[dict[str, Any]]:
    p = _path(spec)
    if not p.exists():
        return []
    with p.open(encoding="utf-8") as f:
        return [dict(r, value=float(r["value"])) for r in csv.DictReader(f)]


def add_measurement(spec: CampaignSpec, formula: str, value: float, source: Source, note: str = "") -> dict[str, Any]:
    row = {"formula": _key(formula), "value": float(value), "source": source, "note": note,
           "added": dt.datetime.now().isoformat(timespec="seconds")}
    STORE.mkdir(parents=True, exist_ok=True)
    new = not _path(spec).exists()
    with _path(spec).open("a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new:
            w.writeheader()
        w.writerow(row)
    return row


def import_l2(spec: CampaignSpec, formulas: list[str]) -> list[dict[str, Any]]:
    """저장된 L2 결과(두 모델 평균 hull 거리)를 측정값으로 — 목표가 ehull 일 때만. 이미 있는 조성은 건너뛴다."""
    from msl import l2

    if spec.target.prop != "ehull":
        raise ValueError("L2 결과는 hull 거리만 있어 목표 물성이 ehull 일 때만 가져올 수 있음")
    have = {m["formula"] for m in measurements(spec)}
    added = []
    for f in formulas:
        k = _key(f)
        if k in have:
            continue
        res = l2.saved(k)
        if res and res.get("ehull_mean") is not None:
            added.append(add_measurement(spec, k, res["ehull_mean"], "L2", f"두 모델 평균, 차이 {res.get('ehull_spread')}"))
    return added


# ── 제안 ──────────────────────────────────────────────────────────────────


def _target(spec: CampaignSpec, name: str):
    from baybe.targets import NumericalTarget

    o = spec.target
    if o.direction == "target":
        sigma = spec.sigma or max(abs(o.value) * 0.1, 0.05)
        return NumericalTarget.match_bell(name, match_value=o.value, sigma=sigma)
    return NumericalTarget(name=name, minimize=o.direction == "min")


def suggest(spec: CampaignSpec, batch: int | None = None, candidates: list[dict[str, Any]] | None = None,
            meas: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """다음에 잴 조성 batch 개. candidates·meas 를 주면 파일 대신 그것을 쓴다 (벤치마크용)."""
    import numpy as np
    import pandas as pd
    from baybe import Campaign
    from baybe.parameters import CustomDiscreteParameter
    from baybe.searchspace import SearchSpace
    from baybe.utils.random import set_random_seed

    from msl.ml.features import feature_names, featurize_many

    set_random_seed(spec.seed)
    cands = candidates if candidates is not None else pool(spec)
    if len(cands) < 2:
        raise ValueError("후보가 2개보다 적음")
    prop = spec.target.prop
    formulas = [c["formula"] for c in cands]
    X = pd.DataFrame(featurize_many(formulas), index=formulas, columns=feature_names())
    priors = np.array([np.nan if c["prior"] is None else c["prior"] for c in cands], dtype=float)
    X[f"prior_{prop}"] = np.where(np.isnan(priors), np.nanmedian(priors) if np.isfinite(priors).any() else 0.0, priors)
    X = X.loc[~X.index.duplicated(), X.std() > 1e-9]

    meas = measurements(spec) if meas is None else meas
    in_pool = [m for m in meas if m["formula"] in X.index]
    outside = sorted({m["formula"] for m in meas} - set(X.index))
    sources = sorted({m["source"] for m in in_pool})
    name = f"{prop}"
    param = CustomDiscreteParameter(name="조성", data=X, decorrelate=True)
    camp = Campaign(searchspace=SearchSpace.from_product([param]), objective=_target(spec, name))
    camp.allow_recommending_already_measured = False
    if in_pool:
        camp.add_measurements(pd.DataFrame({"조성": [m["formula"] for m in in_pool], name: [m["value"] for m in in_pool]}))
    n = min(batch or spec.batch, len(X) - len({m["formula"] for m in in_pool}))
    if n < 1:
        raise ValueError("모든 후보를 이미 쟀음")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rec = camp.recommend(n)
        stats = camp.posterior_stats(rec) if in_pool else None
    by = {c["formula"]: c for c in cands}
    rows = []
    for i, f in enumerate(rec["조성"].tolist()):
        row = {"formula": f, "prior": by[f]["prior"], "prior_source": by[f]["source"]}
        if stats is not None:
            row["mean"] = round(float(stats[f"{name}_mean"].iloc[i]), 4)
            row["std"] = round(float(stats[f"{name}_std"].iloc[i]), 4)
        rows.append(row)
    notes = []
    if not in_pool:
        notes.append("측정값이 없어 특성 공간에서 서로 멀리 떨어진 조성을 고름 (초기 탐색)")
    elif len(in_pool) < 5:
        notes.append(f"측정값 {len(in_pool)}개 — 대리모델이 아직 거칠다. 몇 번 더 재고 나면 제안이 좋아진다")
    if len(sources) > 1:
        notes.append(f"측정값 출처가 섞여 있음 ({', '.join(sources)}) — L2·DFT·실험 값은 같은 척도가 아니다 (D19)")
    if outside:
        notes.append(f"후보 풀에 없는 측정 조성 {len(outside)}개는 쓰지 않음: {', '.join(outside[:5])}")
    return {"campaign": spec.id, "target": f"{PROP_KO[prop]} ({PROP_UNIT[prop]}) · {spec.target.text()}",
            "n_candidates": len(X), "n_measured": len(in_pool), "sources": sources, "suggestions": rows, "notes": notes,
            "best_measured": (min if spec.target.direction == "min" else max)(in_pool, key=lambda m: m["value"])
            if in_pool and spec.target.direction != "target" else None}
