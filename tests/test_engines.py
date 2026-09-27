"""계산 엔진·시험 v0 — 네트워크 없이 도는 것만 (MP·PubChem 호출 없음)."""

from __future__ import annotations

from pathlib import Path

import pytest

from msl.assays import blend, equilibrium, safety
from msl.assays.base import Context
from msl.engines import nasa
from msl.recipe import route
from msl.registry.load import load_registry
from msl.resolve import Resolved
from msl.schema.quantity import Quantity
from msl.schema.recipe import Component, Recipe, State
from msl.schema.refs import SubstanceRef
from msl.schema.result import Status

REG = load_registry()


def comp(ref: str, formula: str | None, amount: str | None, state: State | None = None,
         groups: list[str] | None = None, props: dict | None = None) -> Resolved:
    return Resolved(ref=SubstanceRef.model_validate(ref), name=ref, formula=formula, groups=groups or [],
                    state=state, amount=Quantity.model_validate(amount) if amount else None, props=props or {})


def ctx(comps: list[Resolved], **conditions: str) -> Context:
    recipe = Recipe.model_validate({
        "id": "rcp-test", "name": "t",
        "components": [{"ref": str(c.ref), "amount": str(c.amount) if c.amount else None,
                        "state": c.state.value if c.state else None} for c in comps],
        "conditions": conditions,
    })
    return Context(recipe=recipe, comps=comps, registry=REG)


# ── NASA 평형 (Cantera) ───────────────────────────────────────────────────


def test_gas_lookup_keeps_molecules() -> None:
    """회귀: 약분 조성으로 비교하면 O₂ 가 산소 원자(O)로 잘못 연결됐다 (2026-09-26 발견)."""
    assert nasa.find_gas("O2").name == "O2"
    assert nasa.find_gas("N2").name == "N2"
    assert nasa.find_gas("CH4").name == "CH4"


def test_methane_air_adiabatic_flame_temperature() -> None:
    gas = nasa.gas_equilibrium({"CH4": 1, "O2": 2, "N2": 7.52}, 298.15, 101325, "HP")
    assert gas.T == pytest.approx(2225, abs=15)  # 교과서 값 약 2,225 K


@pytest.mark.parametrize(("T", "decomposed"), [(1073.15, False), (1173.15, True)])
def test_calcite_decomposition_temperature(T: float, decomposed: bool) -> None:
    """방해석 분해온도는 1 atm CO₂ 에서 약 898 °C — 800 °C 안정, 900 °C 분해."""
    eq = nasa.multiphase_equilibrium({"CaCO3": 1.0}, T, 101325)
    assert ("CaCO3(caL)" not in eq.condensed) is decomposed


def test_nasa_out_of_range_is_explicit() -> None:
    with pytest.raises(nasa.OutOfRange):
        nasa.multiphase_equilibrium({"CaCO3": 1.0}, 1273.15, 101325)


# ── 시험 ──────────────────────────────────────────────────────────────────


def test_a3_methane_air() -> None:
    c = ctx([comp("cas:74-82-8", "CH4", "1 mol", State.GAS), comp("cas:7782-44-7", "O2", "2 mol", State.GAS),
             comp("cas:7727-37-9", "N2", "7.52 mol", State.GAS)], T="25 °C")
    out = equilibrium.a3(c)
    assert out.status is Status.OK
    t_ad = next(v.value for v in out.values if v.unit == "K")
    assert t_ad == pytest.approx(2225, abs=15)


def test_a4_calcite_water_equilibrium() -> None:
    c = ctx([comp("mineral:calcite", "CaCO3", "1 g", State.SOLID), comp("cas:7732-18-5", "H2O", "1 L", State.LIQUID)],
            T="25 °C", atmosphere="air")
    out = equilibrium.a4(c)
    assert out.status is Status.OK
    ph = next(v.value for v in out.values if v.name == "pH")
    assert ph == pytest.approx(8.2, abs=0.1)  # 대기 CO₂ 와 평형인 방해석 포화 용액


def _reaktoro() -> bool:
    from msl.engines import reaktoro_bridge

    return reaktoro_bridge.available()


@pytest.mark.skipif(not _reaktoro(), reason="Reaktoro 환경 없음 (~/micromamba/envs/reaktoro)")
@pytest.mark.parametrize("T, ph", [("25 °C", 9.91), ("75 °C", 8.87)])
def test_a4_reaktoro_cross_check_closed_calcite(T: str, ph: float) -> None:
    """닫힌계 방해석 + 물 — 같은 phreeqc.dat 로 푼 PHREEQC 와 Reaktoro 의 pH 가 0.1 안에서 맞아야 한다 (D22)."""
    c = ctx([comp("mineral:calcite", "CaCO3", "1 g", State.SOLID), comp("cas:7732-18-5", "H2O", "1 L", State.LIQUID)], T=T)
    out = equilibrium.a4(c)
    v = {x.name: x.value for x in out.values}
    assert v["pH"] == pytest.approx(ph, abs=0.05)
    assert v["pH (Reaktoro 교차검증)"] == pytest.approx(v["pH"], abs=equilibrium.PH_CROSS_TOL)
    assert not any("차이" in cv for cv in out.caveats)


@pytest.mark.skipif(not _reaktoro(), reason="Reaktoro 환경 없음 (~/micromamba/envs/reaktoro)")
def test_a4_reaktoro_cross_check_dissolved_salts() -> None:
    c = ctx([comp("cas:144-55-8", "NaHCO3", "0.05 mol", State.AQUEOUS), comp("cas:7647-14-5", "NaCl", "0.1 mol", State.AQUEOUS),
             comp("cas:7732-18-5", "H2O", "1 L", State.LIQUID)])
    out = equilibrium.a4(c)
    v = {x.name: x.value for x in out.values}
    assert v["pH (Reaktoro 교차검증)"] == pytest.approx(v["pH"], abs=equilibrium.PH_CROSS_TOL)


def test_a4_skips_reaktoro_under_air_co2() -> None:
    """대기 CO₂ 고정 조건은 Reaktoro 쪽에 같은 제약을 걸지 않으므로 교차검증을 건너뛴다."""
    c = ctx([comp("mineral:calcite", "CaCO3", "1 g", State.SOLID), comp("cas:7732-18-5", "H2O", "1 L", State.LIQUID)],
            atmosphere="air")
    out = equilibrium.a4(c)
    assert "reaktoro" not in out.data


def test_reaktoro_input_is_charge_balanced() -> None:
    import phreeqpython

    from msl.engines import reaktoro_bridge as rb

    db = str(Path(phreeqpython.__file__).parent / "database" / "phreeqc.dat")
    spec = rb.from_a4({"Na": 0.1, "C": 0.05}, ["Calcite"], db, 1.0, 298.15, 1.0)
    assert sum(rb._charge(s) * n for s, n in spec["species"].items()) == pytest.approx(0, abs=1e-12)
    assert {"Ca", "C", "Na"} <= set(spec["elements"]) and spec["minerals"] == {"Calcite": 10.0}


def test_a4_refuses_unknown_oxidation_state() -> None:
    """차아염소산(Cl +1)은 phreeqc.dat 에 없음 — 조용히 염화물로 계산하지 않고 거부해야 한다."""
    c = ctx([comp("cas:7681-52-9", "NaClO", "0.05 mol", State.AQUEOUS), comp("cas:7732-18-5", "H2O", "1 L", State.LIQUID)])
    out = equilibrium.a4(c)
    assert out.status is Status.NOT_APPLICABLE
    assert "NaClO" in out.summary


def test_a8_blend_bounds() -> None:
    epoxy = {"density": 1.2, "youngs_modulus": 3.5, "poisson_ratio": 0.35}
    glass = {"density": 2.55, "youngs_modulus": 72.0, "poisson_ratio": 0.22}
    c = ctx([comp("material:epoxy-generic", None, "70 vol%", props=epoxy),
             comp("material:e-glass-fiber", None, "30 vol%", props=glass)])
    out = blend.a8(c)
    v = {x.name: x.value for x in out.values}
    assert v["밀도 (혼합법칙)"] == pytest.approx(0.7 * 1.2 + 0.3 * 2.55)
    assert v["영률 상한 · Voigt (섬유 방향)"] == pytest.approx(0.7 * 3.5 + 0.3 * 72, abs=0.01)
    assert v["영률 하한 · Reuss (섬유 수직)"] == pytest.approx(1 / (0.7 / 3.5 + 0.3 / 72), abs=0.01)
    lo, hi = v["영률 · Hashin–Shtrikman 하한"], v["영률 · Hashin–Shtrikman 상한"]
    assert v["영률 하한 · Reuss (섬유 수직)"] <= lo <= hi <= v["영률 상한 · Voigt (섬유 방향)"]


def test_s0_bleach_acid_is_incompatible() -> None:
    c = ctx([comp("cas:7681-52-9", "NaClO", "0.05 mol", State.AQUEOUS,
                  ["Oxidizing Agents, Strong", "Salts, Basic", "Water and Aqueous Solutions"]),
             comp("cas:7647-01-0", "HCl", "0.1 mol", State.AQUEOUS, ["Acids, Strong Non-oxidizing", "Water and Aqueous Solutions"])])
    out = safety.run(c)
    assert out.data["verdict"] == "부적합"
    assert any("Cl₂" in g for p in out.data["pairs"] for g in p["gases"])


def test_s0_no_rule_is_not_called_safe() -> None:
    c = ctx([comp("mineral:calcite", "CaCO3", "1 mol", State.SOLID, ["Carbonate Salts"]),
             comp("mineral:quartz", "SiO2", "1 mol", State.SOLID, ["Siloxanes"])])
    out = safety.run(c)
    assert out.data["verdict"] == "규칙 해당 없음"
    assert any("안전 보장이 아님" in cv for cv in out.caveats)


# ── 라우터 ────────────────────────────────────────────────────────────────


def test_router_auto() -> None:
    solids = [comp("formula:MgO", "MgO", "1 mol", State.SOLID), comp("formula:Al2O3", "Al2O3", "1 mol", State.SOLID)]
    codes = [c for c, _ in route(ctx(solids).recipe, solids)]
    assert codes == ["S0", "A1", "A2", "A9"]
    water = [comp("mineral:calcite", "CaCO3", "1 g", State.SOLID), comp("cas:7732-18-5", "H2O", "1 L", State.LIQUID)]
    assert "A4" in [c for c, _ in route(ctx(water).recipe, water)]


# ── A10 국내 규제 (네트워크 대신 가짜 응답) ───────────────────────────────


def test_a10_flags_key_regulations(monkeypatch: pytest.MonkeyPatch) -> None:
    from msl.assays import regulation
    from msl.engines import datagokr as dg

    monkeypatch.setattr(dg, "keco_substance", lambda cas: {"sbstnNmKor": "염화 수소", "korexst": "KE-20189", "typeList": [
        {"sbstnClsfTypeNm": "기존화학물질", "unqNo": "KE-20189"},
        {"sbstnClsfTypeNm": "사고대비물질", "unqNo": "42", "contInfo": "10% 이상 함유 혼합물"}]})
    monkeypatch.setattr(dg, "keco_ghs", lambda cas: {"sfsgwd": "위험", "pctgrmCd": "GHS05^GHS06",
                                                    "hrmflnList": [{"hrmDngrCd": "H314"}, {"hrmDngrCd": "H331"}]})
    monkeypatch.setattr(dg, "kosha_chem", lambda cas: None)
    c = ctx([comp("cas:7647-01-0", "HCl", "0.1 mol", State.AQUEOUS)])
    c.comps[0].cas = "7647-01-0"
    out = regulation.a10(c)
    assert out.status is Status.WARNING
    row = out.data["table"]["rows"][0]
    assert "사고대비물질(42)" in row[3] and row[4] == "위험" and row[6] == "H314 H331"


def test_a10_pending_without_key(monkeypatch: pytest.MonkeyPatch) -> None:
    from msl.assays import regulation
    from msl.engines import datagokr as dg

    def boom(cas: str) -> None:
        raise dg.DataGoKrUnavailable("DATA_GO_KR_SERVICE_KEY 가 .env 에 없음")

    monkeypatch.setattr(dg, "keco_substance", boom)
    c = ctx([comp("cas:7647-01-0", "HCl", "0.1 mol", State.AQUEOUS)])
    c.comps[0].cas = "7647-01-0"
    assert regulation.a10(c).status is Status.PENDING
