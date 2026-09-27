"""다음에 시험할 조합 (BayBE, D25)."""

from __future__ import annotations

import pytest

pytest.importorskip("baybe")

from msl import suggest as sg  # noqa: E402


@pytest.fixture(autouse=True)
def _store(tmp_path, monkeypatch):
    monkeypatch.setattr(sg, "STORE", tmp_path)


def _spec(**kw):
    return sg.CampaignSpec(**({"id": "camp-t", "name": "t", "formulas": ["LiCoO2"], "batch": 3} | kw))


def _cands(n: int = 12):
    fs = ["LiCoO2", "LiMnO2", "LiNiO2", "LiFeO2", "Li2MnO3", "LiMn2O4", "LiCo2O4", "Li2NiO3", "NaCoO2", "NaMnO2", "KCoO2", "MgTiO3"]
    return [{"formula": f, "source": "L1", "prior": 0.01 * i} for i, f in enumerate(fs[:n])]


def test_spec_needs_goal_or_formulas() -> None:
    with pytest.raises(ValueError):
        sg.CampaignSpec(id="camp-t", name="t")
    with pytest.raises(ValueError):
        sg.CampaignSpec(id="camp-t", name="t", goal="goal-x", formulas=["LiCoO2"])


def test_cold_start_then_uses_measurements() -> None:
    spec = _spec()
    r0 = sg.suggest(spec, candidates=_cands())
    assert len(r0["suggestions"]) == 3 and r0["n_measured"] == 0 and "초기 탐색" in r0["notes"][0]
    for f, v in [("LiCoO2", 0.0), ("NaMnO2", 0.12), ("MgTiO3", 0.3)]:
        sg.add_measurement(spec, f, v, "L2")
    r1 = sg.suggest(spec, candidates=_cands())
    assert r1["n_measured"] == 3 and all("mean" in s for s in r1["suggestions"])
    assert not {s["formula"] for s in r1["suggestions"]} & {"LiCoO2", "NaMnO2", "MgTiO3"}  # 잰 것은 다시 고르지 않음
    assert r1["best_measured"]["formula"] == "LiCoO2"


def test_mixed_sources_warn() -> None:
    spec = _spec()
    sg.add_measurement(spec, "LiCoO2", 0.0, "L2")
    sg.add_measurement(spec, "MgTiO3", 0.1, "실험")
    assert any("섞여" in n for n in sg.suggest(spec, candidates=_cands())["notes"])


def test_measurement_formula_is_normalized() -> None:
    spec = _spec()
    sg.add_measurement(spec, "Li0.5Co0.5O", 0.0, "DFT")
    assert sg.measurements(spec)[0]["formula"] == "LiCoO2"


def test_target_direction_needs_value() -> None:
    spec = _spec(target={"prop": "gap", "direction": "target", "value": 1.5})
    assert sg._target(spec, "gap") is not None
