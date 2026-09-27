"""L1 조성 예측 — 학습된 대리모델로 화학식 하나의 물성을 추정한다 (충실도 L1, D17).

각 예측에는 신뢰도 정보를 함께 붙인다.
- 80% 예측 구간 (분위수 모델)
- 가장 비슷한 알려진 물질 3개 (특성 공간 최근접, DFT 값과 함께)
- 학습 범위 밖 여부 (최근접 거리 > 시험 세트 95% 분위)
- DB 에 이미 있는 조성이면 DFT 값 (예측보다 그 값을 쓰라고 안내)
"""

from __future__ import annotations

from functools import cache
from typing import Any

import numpy as np
from pymatgen.core import Composition

from msl.ml.features import featurize
from msl.ml.train import MODEL_DIR, QUANTILES, VERSION, HullIndex

AMU_PER_A3 = 1.66053907  # g/cm³ per (amu/Å³)
NEAR_HULL = 0.05  # eV/atom


class ModelMissing(RuntimeError):
    """학습된 모델 파일이 없을 때 (msl ml train 필요)."""


@cache
def model() -> dict[str, Any]:
    import joblib

    path = MODEL_DIR / f"{VERSION}.joblib"
    if not path.exists():
        raise ModelMissing(f"모델 파일 없음: {path} — `uv run msl ml train` 으로 먼저 학습")
    return joblib.load(path)


@cache
def _hull() -> HullIndex:
    return HullIndex(model()["nn_table"])


@cache
def _known() -> dict[str, dict[str, Any]]:
    t = model()["nn_table"]
    return {formula_key(Composition(r.formula)): r._asdict() for r in t.itertuples(index=False)}


def formula_key(comp: Composition) -> str:
    """DB(MP formula_pretty)와 맞춘 약분 정수 화학식: Li1.2Ni0.6Mn0.2O2 → Li6MnNi3O10, Li0.5Co0.5O → LiCoO2."""
    return Composition(comp.get_integer_formula_and_factor()[0]).reduced_formula


def info() -> dict[str, Any]:
    """모델 판·학습 데이터·평가 지표 (화면·리포트 표시용)."""
    m = model()
    return {"version": m["version"], "created": m["created"], "data": m["data"], "metrics": m["metrics"],
            "nn_threshold": m["nn_threshold"]}


def stability_label(ehull: float, known: bool = False) -> str:
    if ehull <= 0:
        return "안정 예측 — 다른 알려진 상보다 낮음" + ("" if known else " (새 바닥 상태 후보)")
    if ehull <= NEAR_HULL:
        return "준안정 — hull 가까움 (합성 가능성 있음)"
    return "불안정 예측 — 다른 상들로 분해되는 쪽이 유리"


def predict(formula: str) -> dict[str, Any]:
    comp = Composition(formula)
    if comp.num_atoms <= 0:
        raise ValueError(f"빈 조성: {formula}")
    m = model()
    md = m["models"]
    x = featurize(comp)[None, :]
    reduced = formula_key(comp)
    q_lo, q_hi = QUANTILES

    ef = float(md["ef"].predict(x)[0])
    cq = m.get("conformal", {})  # 분할 등각 보정량 (구간을 넓혀 80% 포함률을 맞춘다)
    ef_lo = float(md[f"ef_q{q_lo}"].predict(x)[0]) - cq.get("ef", 0.0)
    ef_hi = float(md[f"ef_q{q_hi}"].predict(x)[0]) + cq.get("ef", 0.0)
    ef_lo, ef_hi = min(ef_lo, ef), max(ef_hi, ef)  # 분위수·점 모델을 따로 학습해 어긋날 때 구간이 점을 포함하도록
    p_metal = float(md["metal"].predict_proba(x)[0, 1])
    gap_raw = float(md["gap"].predict(x)[0])
    gap = 0.0 if p_metal > 0.5 else max(gap_raw, 0.0)
    gap_lo = max(float(md[f"gap_q{q_lo}"].predict(x)[0]) - cq.get("gap", 0.0), 0.0)
    gap_hi = max(float(md[f"gap_q{q_hi}"].predict(x)[0]) + cq.get("gap", 0.0), 0.0)
    gap_lo, gap_hi = min(gap_lo, gap), max(gap_hi, gap)
    vpa = float(md["vpa"].predict(x)[0])
    density = comp.weight / comp.num_atoms / vpa * AMU_PER_A3 if vpa > 0 else None

    xs = m["scaler"].transform(x).astype(np.float32)[0]
    d = np.sqrt(((m["nn_X"] - xs) ** 2).sum(axis=1))
    order = np.argsort(d)[:4]
    table = m["nn_table"]
    nearest = [{"formula": table["formula"].iloc[i], "material_id": table["material_id"].iloc[i],
                "distance": round(float(d[i]), 2), "ef": round(float(table["ef"].iloc[i]), 3),
                "gap": round(float(table["gap"].iloc[i]), 2), "ehull": round(float(table["ehull"].iloc[i]), 3)}
               for i in order if table["formula"].iloc[i] != reduced][:3]
    dmin = float(d[order[0]]) if table["formula"].iloc[order[0]] != reduced else float(d[order[1]])

    out: dict[str, Any] = {
        "input": formula, "formula": reduced, "elements": sorted(str(e) for e in comp.elements),
        "ef": round(ef, 3), "ef_interval": [round(ef_lo, 3), round(ef_hi, 3)],
        "gap": round(gap, 2), "gap_interval": [round(gap_lo, 2), round(gap_hi, 2)], "p_metal": round(p_metal, 3),
        "vpa": round(vpa, 2), "density": None if density is None else round(density, 3),
        "nearest": nearest, "nn_distance": round(dmin, 2), "in_domain": dmin <= m["nn_threshold"],
        "known": None, "ehull": None, "stability": None, "decomposition": None,
    }
    if (k := _known().get(reduced)) is not None:
        out["known"] = {"material_id": k["material_id"], "ef": round(float(k["ef"]), 3), "ehull": round(float(k["ehull"]), 3),
                        "gap": round(float(k["gap"]), 2), "density": None if k["density_g_cm3"] is None else round(float(k["density_g_cm3"]), 3),
                        "theoretical": bool(k["theoretical"])}
    if len(comp.elements) >= 2:
        res = _hull().analyze(comp, exclude=reduced)
        if res is not None:
            hull_e, decomp = res
            out["ehull"] = round(ef - hull_e, 3)
            out["stability"] = stability_label(ef - hull_e, known=out["known"] is not None)
            out["decomposition"] = decomp
    return out


@cache
def _pd_for(elements: tuple[str, ...]):
    """원소 집합의 로컬 상태도 (알려진 바닥 다형 전부). 일괄 예측에서 계 마다 한 번만 만든다."""
    import warnings
    from itertools import combinations

    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

    h = _hull()
    entries = [PDEntry(Composition(e), 0.0, name=e) for e in elements]
    for k in range(2, len(elements) + 1):
        for sub in combinations(elements, k):
            for f, amt, ef in h.by_set.get(frozenset(sub), []):
                c = Composition(amt)
                entries.append(PDEntry(c, ef * c.num_atoms, name=f))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return PhaseDiagram(entries)


def predict_batch(formulas: list[str], chunk: int = 256) -> list[dict[str, Any]]:
    """여러 조성을 한 번에 — 모델 예측·최근접 거리를 묶어서 계산 (추천에서 수천 개).

    predict() 와 같은 값을 내되, 분해 생성물·비슷한 물질 목록은 빼고 최근접 거리만 준다 (필요하면 predict() 로 자세히).
    DB 에 있는 조성은 known 에 DFT 값을 붙인다. hull 거리는 DB 에 없는 조성만 계산한다 (있는 조성은 DFT 값 사용).
    """
    from msl.ml.features import featurize_many

    if not formulas:
        return []
    m = model()
    md, cq = m["models"], m.get("conformal", {})
    comps = [Composition(f) for f in formulas]
    X = featurize_many(formulas)
    q_lo, q_hi = QUANTILES
    ef = md["ef"].predict(X)
    ef_lo, ef_hi = md[f"ef_q{q_lo}"].predict(X) - cq.get("ef", 0.0), md[f"ef_q{q_hi}"].predict(X) + cq.get("ef", 0.0)
    p_metal = md["metal"].predict_proba(X)[:, 1]
    gap = np.where(p_metal > 0.5, 0.0, np.clip(md["gap"].predict(X), 0, None))
    gap_lo = np.clip(md[f"gap_q{q_lo}"].predict(X) - cq.get("gap", 0.0), 0, None)
    gap_hi = np.clip(md[f"gap_q{q_hi}"].predict(X) + cq.get("gap", 0.0), 0, None)
    vpa = md["vpa"].predict(X)
    ef_lo, ef_hi = np.minimum(ef_lo, ef), np.maximum(ef_hi, ef)
    gap_lo, gap_hi = np.minimum(gap_lo, gap), np.maximum(gap_hi, gap)

    Xs = m["scaler"].transform(X).astype(np.float32)
    B = m["nn_X"]
    b2 = (B * B).sum(axis=1)
    dmin = np.empty(len(Xs), dtype=np.float32)
    for i in range(0, len(Xs), chunk):
        a = Xs[i:i + chunk]
        d2 = (a * a).sum(axis=1)[:, None] + b2[None, :] - 2 * a @ B.T
        dmin[i:i + chunk] = np.sqrt(np.clip(d2.min(axis=1), 0, None))

    known = _known()
    out = []
    for i, c in enumerate(comps):
        key = formula_key(c)
        k = known.get(key)
        density = c.weight / c.num_atoms / vpa[i] * AMU_PER_A3 if vpa[i] > 0 else None
        rec: dict[str, Any] = {
            "formula": key, "elements": sorted(str(e) for e in c.elements),
            "ef": float(ef[i]), "ef_interval": [float(ef_lo[i]), float(ef_hi[i])],
            "gap": float(gap[i]), "gap_interval": [float(gap_lo[i]), float(gap_hi[i])], "p_metal": float(p_metal[i]),
            "density": None if density is None else float(density), "nn_distance": float(dmin[i]),
            "in_domain": bool(dmin[i] <= m["nn_threshold"]), "known": None, "ehull": None, "ehull_interval": None,
        }
        if k is not None:
            rec["known"] = {"material_id": k["material_id"], "ef": float(k["ef"]), "ehull": float(k["ehull"]), "gap": float(k["gap"]),
                            "density": None if k["density_g_cm3"] is None else float(k["density_g_cm3"]), "theoretical": bool(k["theoretical"])}
        elif len(c.elements) >= 2:
            h = float(_pd_for(tuple(rec["elements"])).get_hull_energy_per_atom(c))
            rec["ehull"] = float(ef[i]) - h
            rec["ehull_interval"] = [float(ef_lo[i]) - h, float(ef_hi[i]) - h]
        out.append(rec)
    return out
