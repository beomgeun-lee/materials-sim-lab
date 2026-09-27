"""조성 특성 — 화학식만으로 만드는 수치 벡터 (Magpie 방식, Ward et al. npj Comput Mater 2016).

원소 물성 20종마다 조성 가중 평균·평균 편차·최솟값·최댓값·범위(5 × 20 = 100),
화학량론 노름 6개, 원자가 전자 궤도 분율 4개, 이온성 2개 → 112개.
원소 데이터는 pymatgen 내장값만 쓴다 (오프라인, 추가 패키지 없음). 빠진 값은 원소 전체의 중앙값으로 채운다.
"""

from __future__ import annotations

from functools import cache

import numpy as np
from pymatgen.core import Composition, Element

PROPS = [
    "Z", "atomic_mass", "row", "group", "mendeleev_no", "X", "atomic_radius", "atomic_radius_calculated",
    "melting_point", "molar_volume", "electron_affinity", "ionization_energy", "density_of_solid",
    "average_ionic_radius", "n_s", "n_p", "n_d", "n_f", "n_valence", "n_unfilled",
]
STATS = ["mean", "avg_dev", "min", "max", "range"]
NORMS = [0, 2, 3, 5, 7, 10]  # p=0 은 원소 수


def _valence(el: Element) -> tuple[int, int, int, int]:
    """원자가 s·p·d·f 전자 수 (Magpie 정의: 최외각 n 의 s·p, n−1 의 d, n−2 의 f)."""
    shells = el.full_electronic_structure
    n = max(s[0] for s in shells)
    count = {"s": 0, "p": 0, "d": 0, "f": 0}
    for shell, orb, e in shells:
        if (orb in "sp" and shell == n) or (orb == "d" and shell == n - 1) or (orb == "f" and shell == n - 2):
            count[orb] += e
    return count["s"], count["p"], count["d"], count["f"]


def _raw(el: Element, prop: str) -> float | None:
    if prop in ("n_s", "n_p", "n_d", "n_f", "n_valence", "n_unfilled"):
        s, p, d, f = _valence(el)
        if prop == "n_valence":
            return float(s + p + d + f)
        if prop == "n_unfilled":
            return float(sum(cap - x for x, cap in ((s, 2), (p, 6), (d, 10), (f, 14)) if x > 0))
        return float({"n_s": s, "n_p": p, "n_d": d, "n_f": f}[prop])
    try:
        v = getattr(el, prop)
    except Exception:
        return None
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


@cache
def element_table() -> dict[str, np.ndarray]:
    """원소 → 물성 벡터 (빠진 값은 그 물성의 원소 중앙값)."""
    import warnings

    els = [Element.from_Z(z) for z in range(1, 104)]
    with warnings.catch_warnings():  # pymatgen 이 빠진 원소 데이터마다 경고를 낸다
        warnings.simplefilter("ignore")
        raw = np.array([[v if (v := _raw(e, p)) is not None else np.nan for p in PROPS] for e in els], dtype=float)
    med = np.nanmedian(raw, axis=0)
    filled = np.where(np.isnan(raw), med, raw)
    return {e.symbol: filled[i] for i, e in enumerate(els)}


def feature_names() -> list[str]:
    names = [f"{s}_{p}" for p in PROPS for s in STATS]
    names += [f"norm_{p}" for p in NORMS]
    names += ["frac_s", "frac_p", "frac_d", "frac_f", "ionic_max", "ionic_mean"]
    return names


def featurize(comp: Composition | str) -> np.ndarray:
    """조성 하나 → 특성 벡터 (길이 len(feature_names()))."""
    comp = Composition(comp) if isinstance(comp, str) else comp
    fr = comp.fractional_composition
    syms = [str(e) for e in fr.elements]
    x = np.array([fr[s] for s in syms], dtype=float)
    tab = element_table()
    P = np.array([tab[s] for s in syms])  # (원소 수, 물성 수)
    mean = x @ P
    avg_dev = x @ np.abs(P - mean)
    mn, mx = P.min(axis=0), P.max(axis=0)
    stats = np.stack([mean, avg_dev, mn, mx, mx - mn], axis=1).reshape(-1)  # 물성마다 5개
    norms = [float(len(x))] + [float(np.sum(x ** p) ** (1 / p)) for p in NORMS[1:]]
    idx = {p: PROPS.index(p) for p in ("n_s", "n_p", "n_d", "n_f", "n_valence", "X")}
    total_v = float(x @ P[:, idx["n_valence"]]) or 1.0
    orb = [float(x @ P[:, idx[k]]) / total_v for k in ("n_s", "n_p", "n_d", "n_f")]
    chi = P[:, idx["X"]]
    if len(x) > 1:
        ic = 1 - np.exp(-0.25 * (chi[:, None] - chi[None, :]) ** 2)
        ionic = [float(ic.max()), float(x @ ic @ x)]
    else:
        ionic = [0.0, 0.0]
    return np.concatenate([stats, norms, orb, ionic]).astype(np.float32)


def featurize_many(formulas: list[str]) -> np.ndarray:
    return np.vstack([featurize(f) for f in formulas])
