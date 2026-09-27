"""목표 기반 추천 — 전하 균형, 후보 나열, 목표 검증, 판정 규칙, 실제 추천·저장·API (모델이 없으면 뒤쪽은 skip)."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError
from pymatgen.core import Composition

from msl.ml.train import MODEL_DIR, VERSION
from msl.recommend import Constraint, Goal, _judge, candidates, charge_balanced

needs_model = pytest.mark.skipif(not (MODEL_DIR / f"{VERSION}.joblib").exists(), reason="L1 모델 없음 — uv run msl ml train")


@pytest.mark.parametrize(("formula", "ok"), [
    ("Li2O", True), ("LiO", False), ("LiO2", False), ("LiMn2O4", True), ("Li6MnNi3O10", True),  # 혼합 원자가 허용
    ("Fe3O4", True), ("NaCl", True), ("NaCl2", False), ("Li2SO4", True), ("Cu2ZnSnS4", True),
    ("CuAu", True), ("Al3Mg2", True),  # 금속끼리는 검사하지 않음
])
def test_charge_balanced(formula: str, ok: bool) -> None:
    assert charge_balanced(Composition(formula)) is ok


def test_candidates_enumeration() -> None:
    g = Goal(id="goal-t", name="t", required=["Li", "O"], max_atoms=4)
    got = sorted(c.reduced_formula for c in candidates(g))
    want = sorted(Composition(f).reduced_formula for f in ["LiO", "Li2O", "LiO2", "Li3O", "LiO3"])  # pymatgen 은 LiO 를 Li2O2 로 표기
    assert got == want and len(got) == len(set(got))  # Li2O2 와 LiO 는 한 번만
    g2 = Goal(id="goal-t", name="t", required=["Li"], optional=["O", "S"], max_elements=2, max_atoms=3)
    names = {c.reduced_formula for c in candidates(g2)}
    assert "Li" in names and {"Li2O", "Li2S"} <= names and not any("O" in n and "S" in n for n in names)


def test_goal_validation() -> None:
    with pytest.raises(ValidationError):
        Goal(id="goal-t", name="t")  # 원소 없음
    with pytest.raises(ValidationError):
        Goal(id="goal-t", name="t", required=["Li"], optional=["Li"])  # 겹침
    with pytest.raises(ValidationError):
        Goal(id="goal-t", name="t", required=["Xx"])
    with pytest.raises(ValidationError):
        Goal(id="goal-t", name="t", required=["Li"], objective={"prop": "gap", "direction": "target"})  # 목표값 없음
    with pytest.raises(ValidationError):
        Goal(id="bad", name="t", required=["Li"])


def test_judge_with_intervals() -> None:
    c = Constraint(prop="ehull", max=0.05)
    assert _judge(c, 0.01, (-0.02, 0.04)) == "충족"
    assert _judge(c, 0.01, (-0.02, 0.10)) == "가능성 있음"  # 점은 충족이지만 구간이 경계를 넘음
    assert _judge(c, 0.08, (0.02, 0.15)) == "가능성 있음"  # 점은 불충족이지만 구간이 조건 안에 걸침
    assert _judge(c, 0.20, (0.10, 0.30)) == "불충족"
    assert _judge(Constraint(prop="gap", min=1, max=2), 1.5, None) == "충족"  # DFT 값은 구간 없음
    assert _judge(c, None, None) == "불충족"


@needs_model
def test_recommend_li_co_o(tmp_path, monkeypatch) -> None:
    import msl.recommend as rc

    monkeypatch.setattr(rc, "RUNS", tmp_path)
    g = Goal(id="goal-t", name="Li-Co-O", required=["Li", "Co", "O"], max_atoms=8,
             constraints=[{"prop": "ehull", "max": 0.02}, {"prop": "p_metal", "max": 0.5}], top_k=50)
    res = rc.recommend(g)
    top = res["results"]
    assert top and all(r["status"] != "불충족" for r in top)
    tiers = [rc.TIER[r["status"]] for r in top]
    assert tiers == sorted(tiers)  # 판정 순
    for t in set(tiers):  # 같은 판정 안에서는 DB(DFT) 후보가 L1 예측보다 먼저 (웹 사용 점검 11번)
        src = [r["source"] for r in top if rc.TIER[r["status"]] == t]
        assert src == sorted(src, key=lambda x: x != "DFT")
    lco = next(r for r in top if r["formula"] == "LiCoO2")  # 실제 양극재 LiCoO2 는 DB 값으로 충족
    assert lco["source"] == "DFT" and lco["status"] == "충족"
    for r in top:
        if r["source"] == "L1" and r["intervals"]["ehull"]:
            lo, hi = r["intervals"]["ehull"]
            assert lo <= r["values"]["ehull"] <= hi
    c = res["counts"]
    assert c["candidates"] >= c["charge_balanced"] >= c["evaluated"] == c["known"] + c["new"]
    path = rc.save_run(res)
    assert path.parent == tmp_path and rc.list_runs()[0]["file"] == path.name
    assert rc.to_csv(res).splitlines()[0].startswith("순위,화학식,출처")


@needs_model
def test_recommend_api(tmp_path, monkeypatch) -> None:
    from fastapi import HTTPException

    import msl.recommend as rc
    from msl.web import app as web

    monkeypatch.setattr(rc, "RUNS", tmp_path)
    out = web.recommend_api(web.RecommendRequest(goal={"id": "goal-web-t", "name": "t", "required": ["Al"], "optional": ["Mg"],
                                                       "max_atoms": 5, "charge_balanced": False, "constraints": [{"prop": "ehull", "max": 0.03}]}))
    assert out["result"]["results"] and (tmp_path / out["file"]).exists()
    assert web.recommend_run(out["file"])["goal"]["id"] == "goal-web-t"
    with pytest.raises(HTTPException):
        web.recommend_api(web.RecommendRequest(goal={"id": "goal-t", "name": "t"}))
    with pytest.raises(HTTPException):
        web.recommend_run("../../etc/passwd")
    assert len(web.goals()) >= 3
