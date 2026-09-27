"""L1 조성 대리모델 학습 — MP(CC BY 4.0) 화학식별 바닥 다형으로 형성에너지·밴드갭·금속 여부·원자 부피를 배운다.

평가 두 가지 (계획서 4단계 L1, D17):
- 무작위 분할(조성 80/20): 이미 아는 화학계 안의 다른 조성
- 화학계 분할(chemsys 통째로 20% 빼기): 학습에 없던 원소 조합 — 새 조성의 현실적 기대치
그리고 형성에너지 예측으로 hull 거리를 계산해 '안정성 판정'이 얼마나 맞는지 잰다 (Bartel 2020 이 지적한 약점).
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import time
import warnings
from itertools import combinations
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from msl.env import DATA_DIR
from msl.ml.features import feature_names, featurize

MODEL_DIR = DATA_DIR / "models"
VERSION = "l1-composition-v1"
TARGETS = ("ef", "gap", "vpa")  # 회귀: 형성에너지(eV/atom), 밴드갭(eV), 원자당 부피(Å³)
QUANTILES = (0.1, 0.9)  # 80% 예측 구간 (분위수 모델 + 분할 등각 보정)
MAX_EHULL = 0.5  # eV/atom — 바닥 다형이 이보다 hull 위면 비물리적 가상 구조뿐인 조성 → 학습에서 뺀다
MAX_VPA = 150.0  # Å³/atom — 원자 하나에 이보다 큰 부피는 진공 상자 속 분자·실패한 이완


def ground_states() -> pd.DataFrame:
    """화학식별 바닥 다형 (hull 에 가장 가까운 것). deprecated 제외."""
    from msl.db import connect

    df = connect().sql("""
        select material_id, formula_pretty as formula, chemsys, nsites, volume_A3, density_g_cm3,
               formation_energy_per_atom_eV as ef, energy_above_hull_eV as ehull, band_gap_eV as gap,
               is_metal, theoretical, source_version
        from phases_calc where not deprecated and formation_energy_per_atom_eV is not null
              and band_gap_eV is not null and volume_A3 is not null and nsites > 0
        qualify row_number() over (partition by formula_pretty order by energy_above_hull_eV, formation_energy_per_atom_eV) = 1
    """).df()
    df["vpa"] = df["volume_A3"] / df["nsites"]
    df["metal"] = df["is_metal"].fillna(df["gap"] <= 0).astype(bool)  # is_metal 이 빈 행은 밴드갭 0 으로 판단
    keep = (df["ehull"] <= MAX_EHULL) & (df["vpa"] <= MAX_VPA)
    df.attrs["dropped"] = int((~keep).sum())
    return df[keep].reset_index(drop=True)


def _hgb(kind: str, **kw):
    from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor

    common = dict(max_iter=1500, learning_rate=0.08, max_leaf_nodes=63, min_samples_leaf=20, l2_regularization=1e-3,
                  early_stopping=True, validation_fraction=0.1, n_iter_no_change=30, random_state=0)
    return (HistGradientBoostingClassifier if kind == "clf" else HistGradientBoostingRegressor)(**common, **kw)


def _fit_all(X: np.ndarray, df: pd.DataFrame, quantiles: bool, log: Callable[[str], None]) -> dict[str, Any]:
    models: dict[str, Any] = {}
    for t in TARGETS:
        t0 = time.time()
        models[t] = _hgb("reg").fit(X, df[t].to_numpy())
        log(f"  {t} 학습 {time.time() - t0:.0f}초 (반복 {models[t].n_iter_})")
        if quantiles and t in ("ef", "gap"):
            for q in QUANTILES:
                models[f"{t}_q{q}"] = _hgb("reg", loss="quantile", quantile=q).fit(X, df[t].to_numpy())
    models["metal"] = _hgb("clf").fit(X, df["metal"].to_numpy())
    return models


def _score(models: dict[str, Any], X: np.ndarray, df: pd.DataFrame) -> dict[str, float]:
    from sklearn.metrics import accuracy_score, mean_absolute_error, r2_score, roc_auc_score

    out: dict[str, float] = {}
    for t in TARGETS:
        pred = models[t].predict(X)
        out[f"{t}_mae"] = float(mean_absolute_error(df[t], pred))
        out[f"{t}_r2"] = float(r2_score(df[t], pred))
        out[f"{t}_baseline_mae"] = float(np.mean(np.abs(df[t] - df[t].mean())))  # 평균값만 찍을 때
        lo, hi = f"{t}_q{QUANTILES[0]}", f"{t}_q{QUANTILES[1]}"
        if lo in models:
            inside = (df[t] >= models[lo].predict(X)) & (df[t] <= models[hi].predict(X))
            out[f"{t}_interval80_coverage"] = float(inside.mean())
    p = models["metal"].predict_proba(X)[:, 1]
    out["metal_accuracy"] = float(accuracy_score(df["metal"], p > 0.5))
    out["metal_auc"] = float(roc_auc_score(df["metal"], p))
    ins = ~df["metal"].to_numpy()
    out["gap_mae_nonmetal"] = float(np.mean(np.abs(df["gap"].to_numpy()[ins] - np.clip(models["gap"].predict(X[ins]), 0, None))))
    return out


def _conformal(models, X: np.ndarray, df: pd.DataFrame, metrics: dict[str, float]) -> dict[str, float]:
    """분할 등각 보정 (CQR, Romano 2019): 분위수 구간을 넓혀 80% 가 실제로 80% 를 덮게 한다.

    시험 세트 앞 반으로 보정량 q 를 정하고, 뒤 반에서 보정 후 포함률을 잰다 (보정에 쓴 자료로 재지 않는다).
    """
    rng = np.random.default_rng(1)
    order = rng.permutation(len(df))
    cal, ev = order[: len(order) // 2], order[len(order) // 2:]
    out: dict[str, float] = {}
    for t in ("ef", "gap"):
        y = df[t].to_numpy()
        lo, hi = models[f"{t}_q{QUANTILES[0]}"].predict(X), models[f"{t}_q{QUANTILES[1]}"].predict(X)
        score = np.maximum(lo[cal] - y[cal], y[cal] - hi[cal])
        level = min(1.0, (QUANTILES[1] - QUANTILES[0]) * (1 + 1 / len(cal)))
        q = float(np.quantile(score, level))
        out[t] = q
        inside = (y[ev] >= lo[ev] - q) & (y[ev] <= hi[ev] + q)
        metrics[f"{t}_interval80_coverage_calibrated"] = float(inside.mean())
        metrics[f"{t}_interval80_width_calibrated"] = float(np.mean(hi[ev] - lo[ev] + 2 * q))
    return out


class HullIndex:
    """원소 집합 → 그 원소들로만 된 바닥 다형 (형성에너지). 조성의 hull 에너지를 로컬 DB 로 계산한다."""

    def __init__(self, df: pd.DataFrame):
        from pymatgen.core import Composition

        self.by_set: dict[frozenset, list[tuple[str, dict[str, float], float]]] = {}
        for f, ef in zip(df["formula"], df["ef"]):
            c = Composition(f)
            key = frozenset(str(e) for e in c.elements)
            self.by_set.setdefault(key, []).append((f, c.get_el_amt_dict(), float(ef)))

    def analyze(self, comp, exclude: str | None = None) -> tuple[float, dict[str, float]] | None:
        """조성 comp 의 hull 에너지 (eV/atom) 와 그 조성의 분해 생성물 {화학식: 원자 분율}.

        exclude 화학식은 빼고 계산한다 (자기 자신을 뺀 '다른 알려진 상 대비').
        """
        from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
        from pymatgen.core import Composition

        els = sorted(str(e) for e in comp.elements)
        entries = [PDEntry(Composition(e), 0.0, name=e) for e in els]
        for k in range(2, len(els) + 1):
            for sub in combinations(els, k):
                for f, amt, ef in self.by_set.get(frozenset(sub), []):
                    if f != exclude:
                        c = Composition(amt)
                        entries.append(PDEntry(c, ef * c.num_atoms, name=f))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                pd = PhaseDiagram(entries)
                decomp = pd.get_decomposition(comp)
                energy = sum(e.energy_per_atom * x for e, x in decomp.items())
                return float(energy), {e.name: round(float(x), 4) for e, x in sorted(decomp.items(), key=lambda p: -p[1])}
            except Exception:
                return None

    def hull_energy(self, comp, exclude: str | None = None) -> float | None:
        res = self.analyze(comp, exclude)
        return None if res is None else res[0]


def _stability_check(models, X_test, test: pd.DataFrame, hull: HullIndex, n: int, log) -> dict[str, float]:
    """형성에너지 예측 → hull 거리. DFT 로 계산한 같은 양과 비교 (자기 자신은 hull 에서 뺀다)."""
    from pymatgen.core import Composition

    rng = np.random.default_rng(0)
    idx = rng.choice(len(test), size=min(n, len(test)), replace=False)
    idx = [i for i in idx if len(Composition(test["formula"].iloc[i]).elements) >= 2]
    pred = models["ef"].predict(X_test[idx])
    rows = []
    for j, i in enumerate(idx):
        f = test["formula"].iloc[i]
        h = hull.hull_energy(Composition(f), exclude=f)
        if h is not None:
            rows.append((pred[j] - h, float(test["ef"].iloc[i]) - h))
    d = np.array(rows)
    stable_true, stable_pred = d[:, 1] <= 0.0, d[:, 0] <= 0.0
    tp = float(np.sum(stable_true & stable_pred))
    log(f"  안정성 점검 {len(d)}개 조성")
    return {"ehull_mae": float(np.mean(np.abs(d[:, 0] - d[:, 1]))), "stable_fraction": float(stable_true.mean()),
            "stable_precision": tp / max(float(stable_pred.sum()), 1.0), "stable_recall": tp / max(float(stable_true.sum()), 1.0),
            "stability_accuracy": float(np.mean(stable_true == stable_pred)), "n_checked": int(len(d))}


def train(log: Callable[[str], None] = print, stability_n: int = 1500, sample: int | None = None,
          out_dir: Path | None = None) -> Path:
    """sample: 조성 수를 줄여 빠르게 점검 (파이프라인 시험용). out_dir: 저장 폴더 (기본 data/models)."""
    from sklearn.model_selection import GroupShuffleSplit, train_test_split
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import StandardScaler

    t0 = time.time()
    df = ground_states()
    if sample:
        df = df.sample(n=min(sample, len(df)), random_state=0).reset_index(drop=True)
    log(f"학습 데이터: MP 바닥 다형 {len(df):,}개 조성 (원본 {df['source_version'].iloc[0]}) · "
        f"비물리적 조성 {df.attrs.get('dropped', 0):,}개 제외 (hull 위 {MAX_EHULL} eV/atom 초과 또는 원자 부피 {MAX_VPA:g} Å³ 초과)")
    X = np.vstack([featurize(f) for f in df["formula"]])
    log(f"특성 {X.shape[1]}개 · {time.time() - t0:.0f}초")

    metrics: dict[str, Any] = {}
    tr, te = train_test_split(np.arange(len(df)), test_size=0.2, random_state=0)
    log("평가 1 — 무작위 분할 (조성 80/20)")
    m = _fit_all(X[tr], df.iloc[tr], quantiles=True, log=log)
    metrics["random"] = _score(m, X[te], df.iloc[te])
    conformal = _conformal(m, X[te], df.iloc[te], metrics["random"])
    metrics["random"].update(_stability_check(m, X[te], df.iloc[te], HullIndex(df.iloc[tr]), stability_n, log))
    scaler_tr = StandardScaler().fit(X[tr])
    nn_tr = NearestNeighbors(n_neighbors=1).fit(scaler_tr.transform(X[tr]))
    dist_te = nn_tr.kneighbors(scaler_tr.transform(X[te]))[0][:, 0]
    err = np.abs(m["ef"].predict(X[te]) - df["ef"].to_numpy()[te])
    q95 = float(np.quantile(dist_te, 0.95))
    metrics["domain"] = {"nn_distance_q95": q95, "ef_mae_inside": float(err[dist_te <= q95].mean()),
                         "ef_mae_outside": float(err[dist_te > q95].mean())}

    log("평가 2 — 화학계 분할 (원소 조합 통째로 20% 제외)")
    gtr, gte = next(GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=0).split(X, groups=df["chemsys"]))
    m2 = _fit_all(X[gtr], df.iloc[gtr], quantiles=False, log=log)
    metrics["chemsys"] = _score(m2, X[gte], df.iloc[gte])

    log("최종 모델 — 전체 데이터")
    final = _fit_all(X, df, quantiles=True, log=log)
    scaler = StandardScaler().fit(X)
    Xs = scaler.transform(X).astype(np.float32)
    fp = hashlib.sha256("".join(sorted(df["material_id"])).encode()).hexdigest()[:16]
    artifact = {
        "version": VERSION, "created": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "feature_names": feature_names(), "models": final, "scaler": scaler,
        "nn_X": Xs, "nn_table": df[["material_id", "formula", "ef", "ehull", "gap", "density_g_cm3", "theoretical"]].reset_index(drop=True),
        "nn_threshold": q95, "metrics": metrics, "conformal": conformal,
        "data": {"source": "materials-project", "version": str(df["source_version"].iloc[0]), "n": int(len(df)),
                 "license": "CC-BY-4.0", "material_ids_sha256_16": fp,
                 "target_scheme": "MP summary 기본 thermo (GGA/GGA+U/r2SCAN 혼합) 형성에너지"},
    }
    import joblib

    folder = out_dir or MODEL_DIR
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{VERSION}.joblib"
    joblib.dump(artifact, path, compress=3)
    meta = {k: v for k, v in artifact.items() if k in ("version", "created", "feature_names", "nn_threshold", "metrics", "conformal", "data")}
    (folder / f"{VERSION}.json").write_text(json.dumps(meta, ensure_ascii=False, indent=1), encoding="utf-8")
    log(f"저장: {path} · 전체 {time.time() - t0:.0f}초")
    return path
