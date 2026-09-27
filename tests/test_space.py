"""조합 공간 생성기 — 격자·양 스윕·부분집합 열거, 값 수집, 리포트 (네트워크 없이)."""

from __future__ import annotations

import datetime as dt
from pathlib import Path

import pytest

from msl.recipe.space import Collect, Space, SpaceResult, expand, load_space, pick, run_space
from msl.registry.load import load_registry
from msl.report.space import to_csv, to_html, to_markdown
from msl.schema.provenance import RunProvenance
from msl.schema.recipe import RecipeError
from msl.schema.result import AssayResult, Fidelity, ResultValue, Status, ValueKind

EXAMPLES = Path(__file__).resolve().parents[1] / "examples" / "spaces"


def space(generator: dict, **kw) -> Space:
    return Space.model_validate({"id": "spc-t", "name": "t", "generator": generator,
                                 "collect": kw.pop("collect", [{"assay": "S0", "field": "verdict"}]), **kw})


# ── 생성 ─────────────────────────────────────────────────────────────────


def test_ratio_grid_counts() -> None:
    two = {"kind": "ratio", "components": ["element:Cu", "element:Ni"], "step": 0.1, "unit": "at%"}
    assert len(expand(space(two))) == 9  # 끝점(순수 성분) 제외
    assert len(expand(space({**two, "include_pure": True}))) == 11
    three = {"kind": "ratio", "components": ["element:Cu", "element:Ni", "element:Zn"], "step": 0.25}
    assert len(expand(space(three))) == 3  # 4 를 양의 정수 3개로: (2,1,1) 순열
    assert len(expand(space({**three, "include_pure": True}))) == 15
    v = expand(space(two))[2]
    assert v.params == {"Cu": 0.3, "Ni": 0.7} and [str(c.amount) for c in v.recipe.components] == ["30 at%", "70 at%"]


def test_ratio_total_and_bad_step() -> None:
    v = expand(space({"kind": "ratio", "components": ["formula:MgO", "formula:Al2O3"], "step": 0.5, "total": "2 mol"}))
    assert [c.amount.value for c in v[0].recipe.components] == [1.0, 1.0]
    with pytest.raises(ValueError):
        space({"kind": "ratio", "components": ["element:Cu", "element:Ni"], "step": 0.3}).generator.fractions()


def test_amount_sweep_drops_zero_and_keeps_fixed() -> None:
    sp = space({"kind": "amount", "component": {"ref": "formula:NaOH", "state": "aqueous", "label": "NaOH"},
                "start": "0 mol", "stop": "0.2 mol", "steps": 5},
               base={"components": [{"ref": "formula:HCl", "amount": "0.1 mol", "state": "aqueous"}]})
    vs = expand(sp)
    assert len(vs) == 5 and [str(c.ref) for c in vs[0].recipe.components] == ["formula:HCl"]
    assert [str(c.ref) for c in vs[1].recipe.components] == ["formula:NaOH", "formula:HCl"]
    assert vs[-1].params == {"NaOH": 0.2} and vs[-1].recipe.id == "rcp-t-004"


def test_subsets_and_limit() -> None:
    pool = [f"element:{e}" for e in ("Fe", "Cu", "Zn", "Al")]
    vs = expand(space({"kind": "subsets", "pool": pool, "size": [1, 2], "require": ["element:O"]}, base={"mode": "system"}))
    assert len(vs) == 4 + 6 and all(str(v.recipe.components[-1].ref) == "element:O" for v in vs)
    with pytest.raises(RecipeError, match="max_variants"):
        expand(space({"kind": "subsets", "pool": pool, "size": 2}, max_variants=3))


def test_examples_expand() -> None:
    counts = {p.stem: len(expand(load_space(p))) for p in sorted(EXAMPLES.glob("*.yaml"))}
    assert counts == {"cu-ni-composition": 9, "hcl-naoh-titration": 21, "li-tm-oxides": 4,
                      "metals-in-water": 8, "reagent-shelf-s0": 21}


# ── 값 수집 ──────────────────────────────────────────────────────────────


def _res(assay: str, values: list[tuple[str, object, str | None]], **data: object) -> AssayResult:
    now = dt.datetime(2026, 9, 27, tzinfo=dt.UTC)
    return AssayResult.model_validate(dict(
        assay=assay, recipe_id="rcp-t", status=Status.OK, fidelity=Fidelity.T, engine="x", engine_version="x",
        conditions_basis="-", data=data,
        values=[ResultValue(name=n, kind=ValueKind.PROPERTY, value=v, unit=u) for n, v, u in values],
        provenance=RunProvenance(input_hash="0" * 64, msl_version="t", engine_versions={}, started_at=now, finished_at=now)))


def test_pick() -> None:
    results = [_res("A6", [("고상선 (평형, 액상 첫 출현)", 1467.0, "K"), ("고상선 (평형, 액상 첫 출현)", 1193.9, "°C")], summary="s"),
               _res("S0", [("판정", "부적합", None)], verdict="부적합")]
    assert pick(results, Collect(assay="A6", value="고상선", unit="°C")) == 1193.9  # 앞부분 일치 + 단위
    assert pick(results, Collect(assay="A6", value="고상선")) == 1467.0
    assert pick(results, Collect(assay="S0", field="verdict")) == "부적합"
    assert pick(results, Collect(assay="S0", field="status")) == "ok"
    assert pick(results, Collect(assay="A9", field="summary")) is None
    with pytest.raises(ValueError):
        Collect(assay="A6", value="x", field="y")


# ── 실행·리포트 ───────────────────────────────────────────────────────────


def test_titration_space_runs_and_reports() -> None:
    """강산 0.1 mol 을 강염기로 적정 — 당량점(0.1 mol)에서 pH 7 근처로 급변 (A4, PHREEQC 로컬)."""
    sp = space({"kind": "amount", "component": {"ref": "formula:NaOH", "state": "aqueous", "label": "NaOH"},
                "start": "0 mol", "stop": "0.2 mol", "steps": 5},
               base={"components": [{"ref": "formula:HCl", "amount": "0.1 mol", "state": "aqueous"},
                                    {"ref": "formula:H2O", "amount": "1 L", "state": "liquid"}], "assays": ["A4"]},
               collect=[{"assay": "A4", "value": "pH", "label": "pH"}])
    res = run_space(sp, load_registry(), use_cache=False)
    ph = [r["values"]["pH"] for r in res.rows]
    assert ph == sorted(ph) and ph[0] == pytest.approx(1.08, abs=0.1) and ph[2] == pytest.approx(7.0, abs=0.3)
    assert res.axis() == ("NaOH (mol)", [0.0, 0.05, 0.1, 0.15, 0.2])
    assert to_csv(res).splitlines()[0] == "변형,NaOH,pH,recipe_id"
    assert "| NaOH 0.1 mol |" in to_markdown(res)
    assert "<polyline" in to_html(res)


def test_matrix_for_pairs() -> None:
    sp = space({"kind": "subsets", "pool": ["element:Fe", "element:Cu", "element:Zn"], "size": 2})
    res = SpaceResult(sp, rows=[{"label": "·", "params": {"members": list(m)}, "recipe_id": "r", "values": {"S0 verdict": v}}
                                for m, v in [(("Fe", "Cu"), "부적합"), (("Fe", "Zn"), "주의"), (("Cu", "Zn"), "규칙 해당 없음")]])
    m = res.matrix()
    assert m["cells"][1][0] == m["cells"][0][1] == "부적합" and m["cells"][0][0] is None
    page = to_html(res)
    assert 'class="mx bad"' in page and "안전하다는 뜻이 아니다" in page
