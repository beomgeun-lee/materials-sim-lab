"""L1 조성 대리모델 — 특성, hull 계산, 예측·평가 지표 (모델 파일이 없으면 예측 시험은 skip)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from pymatgen.core import Composition, Element

from msl.ml.features import feature_names, featurize, _valence
from msl.ml.train import MODEL_DIR, VERSION, HullIndex

HAVE_MODEL = (MODEL_DIR / f"{VERSION}.joblib").exists()
needs_model = pytest.mark.skipif(not HAVE_MODEL, reason="L1 모델 없음 — uv run msl ml train")


# ── 특성 ─────────────────────────────────────────────────────────────────


def test_feature_vector() -> None:
    f = featurize("Fe2O3")
    assert f.shape == (len(feature_names()),) == (112,) and not np.isnan(f).any()
    assert np.allclose(featurize("Fe4O6"), f) and np.allclose(featurize("O3Fe2"), f)  # 분율만 본다, 순서 무관
    names = feature_names()
    cu = featurize("Cu")
    assert cu[names.index("ionic_max")] == 0 and cu[names.index("norm_0")] == 1
    nacl = featurize("NaCl")
    assert nacl[names.index("ionic_max")] > 0.5 and nacl[names.index("range_X")] == pytest.approx(3.16 - 0.93, abs=1e-3)


@pytest.mark.parametrize(("el", "spdf"), [("Fe", (2, 0, 6, 0)), ("O", (2, 4, 0, 0)), ("Pb", (2, 2, 10, 14)), ("Au", (1, 0, 10, 14))])
def test_valence_counts(el: str, spdf: tuple) -> None:
    assert _valence(Element(el)) == spdf


# ── hull (로컬 DB 기반) ───────────────────────────────────────────────────


def test_hull_index_toy() -> None:
    df = pd.DataFrame({"formula": ["LiF", "Li2O", "LiF3"], "ef": [-3.0, -2.0, -0.5]})
    h = HullIndex(df)
    e, dec = h.analyze(Composition("LiF"))
    assert e == pytest.approx(-3.0) and dec == {"LiF": 1.0}
    e, dec = h.analyze(Composition("LiF"), exclude="LiF")  # 자기 자신을 빼면 Li + LiF3 로 (2/3 × −0.5)
    assert e == pytest.approx(-1 / 3) and set(dec) == {"Li", "LiF3"}
    e, dec = h.analyze(Composition("LiF2"))  # LiF + F 사이 → hull 위
    assert e == pytest.approx(-3.0 * 2 / 3) and set(dec) == {"LiF", "F"}


# ── 예측 · 지표 ───────────────────────────────────────────────────────────


@needs_model
def test_metrics_meet_floor() -> None:
    """학습 지표 하한 — 문헌의 조성 모델(Magpie + 앙상블) 수준. 떨어지면 데이터·특성 회귀를 의심한다."""
    from msl.ml.predict import info

    m = info()["metrics"]
    assert m["random"]["ef_mae"] < 0.09 and m["chemsys"]["ef_mae"] < 0.10  # 학습 시 0.069 · 0.072 eV/atom
    assert m["random"]["metal_accuracy"] > 0.87 and m["random"]["gap_mae"] < 0.5 and m["random"]["vpa_mae"] < 1.2
    assert m["random"]["stability_accuracy"] > 0.65  # 조성 모델의 약점 — 0.73 (Bartel 2020 과 같은 경향)
    for t in ("ef", "gap"):  # 보정한 80% 구간이 (보정에 안 쓴 절반에서) 실제로 약 80% 를 덮는다
        assert 0.75 < m["random"][f"{t}_interval80_coverage_calibrated"] < 0.9
    for split in ("random", "chemsys"):
        assert m[split]["ef_mae"] < 0.5 * m[split]["ef_baseline_mae"]  # 평균만 찍는 것보다 훨씬 낫다


@needs_model
def test_predict_known_and_new() -> None:
    from msl.ml.predict import predict

    nacl = predict("NaCl")
    assert nacl["known"] and nacl["p_metal"] < 0.2 and nacl["gap"] > 3 and nacl["in_domain"]
    assert abs(nacl["ef"] - nacl["known"]["ef"]) < 0.3
    cu3au = predict("Cu3Au")
    assert cu3au["p_metal"] > 0.8 and cu3au["gap"] == 0
    new = predict("Li1.2Ni0.6Mn0.2O2")  # 정수로 바꾸면 Li6MnNi3O10 — DB 에 없는 조성
    assert new["formula"] == "Li6MnNi3O10" and new["known"] is None and new["ehull"] is not None
    assert predict("Li0.5Co0.5O")["known"]["material_id"] == predict("LiCoO2")["known"]["material_id"]  # 분수 표기도 DB 와 맞춘다
    assert new["ef_interval"][0] <= new["ef"] <= new["ef_interval"][1]
    assert set(new["decomposition"]) and len(new["nearest"]) == 3
