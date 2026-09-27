"""A5 Pourbaix — 표기, 대상 선택, 전위 시나리오, 판정 분류, 하이브리드 에너지 (네트워크 없이).

Pourbaix 엔트리가 필요한 시험은 data/cache/mp 캐시와 이온 스냅샷(data/raw)이 있을 때만 (없으면 skip).
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from msl.assays.aqueous import _scenarios, _verdict, pretty, target
from msl.assays.base import Context
from msl.connectors import mp_ion_ref_data
from msl.db import raw_dir
from msl.engines import mp
from msl.engines import pourbaix as pb
from msl.registry.load import load_registry
from msl.resolve import Resolved
from msl.schema.quantity import Quantity
from msl.schema.recipe import Recipe, State
from msl.schema.refs import SubstanceRef

REG = load_registry()


def res(ref: str, formula: str | None, state: State | None = State.SOLID) -> Resolved:
    return Resolved(ref=SubstanceRef.model_validate(ref), name=ref, formula=formula, state=state,
                    amount=Quantity.model_validate("1 g"))


def needs_pbx(*elements: str):
    path = mp.CACHE / f"pourbaix_{'-'.join(sorted(elements))}_{mp.THERMO_TYPE}_ion{mp_ion_ref_data.VERSION}.pkl"
    seen = mp.CACHE / f"pourbaix_observed_{'-'.join(sorted(elements))}_{mp.THERMO_TYPE}.json"
    ions = raw_dir(mp_ion_ref_data.SOURCE, mp_ion_ref_data.VERSION) / mp_ion_ref_data.FILE
    return pytest.mark.skipif(not (path.exists() and seen.exists() and ions.exists()), reason=f"Pourbaix 캐시·이온 스냅샷 없음: {elements}")


def ctx(conditions: dict) -> Context:
    recipe = Recipe.model_validate({"id": "rcp-t", "name": "t", "conditions": conditions,
                                    "components": [{"ref": "element:Fe", "amount": "1 g"}]})
    return Context(recipe=recipe, comps=[], registry=REG)


# ── 표기·대상 ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize(("raw", "shown"), [
    ("Fe[+2]", "Fe²⁺"), ("FeOH[+1]", "FeOH⁺"), ("Al(OH)4[-1]", "Al(OH)₄⁻"), ("CuO2[-2]", "CuO₂²⁻"),
    ("Fe2O3(s)", "Fe₂O₃(s)"), ("Mg(HO)2(s)", "Mg(OH)₂(s)"), ("H3AuO3(aq)", "H₃AuO₃(aq)"),
])
def test_pretty(raw: str, shown: str) -> None:
    assert pretty(raw) == shown


def test_target_selection() -> None:
    assert target(res("element:Fe", "Fe")) == ["Fe"]
    assert target(res("formula:ZnFe2O4", "ZnFe2O4")) == ["Fe", "Zn"]
    assert target(res("formula:Al(OH)3", "Al(OH)3")) == ["Al"]
    assert target(res("mineral:calcite", "CaCO3")) is None  # 탄소 — 금속·O·H 만이 아님
    assert target(res("cas:7647-14-5", "NaCl")) is None
    assert target(res("cas:7732-18-5", "H2O", State.LIQUID)) is None
    assert target(res("cas:1310-73-2", "NaHO", State.AQUEOUS)) is None  # 녹아 있는 성분은 대상이 아님


# ── 전위 시나리오 ─────────────────────────────────────────────────────────


def test_water_lines() -> None:
    h2, o2 = pb.water_lines(7)
    assert h2 == pytest.approx(-0.414, abs=1e-3) and o2 == pytest.approx(0.815, abs=1e-3)


def test_scenarios_from_conditions() -> None:
    assert _scenarios(ctx({"Eh": -0.5}), 7) == [("지정 전위", -0.5)]
    air = _scenarios(ctx({"atmosphere": "air"}), 7)
    assert len(air) == 1 and air[0][1] == pytest.approx(0.815 - 0.0100, abs=1e-3)  # pO₂ 0.21
    assert _scenarios(ctx({"atmosphere": "N2"}), 7)[0][1] == pytest.approx(-0.414, abs=1e-3)
    assert _scenarios(ctx({"atmosphere": "CO2"}), 7) == _scenarios(ctx({}), 7)  # 'co2' 를 산소로 오인하지 않음
    assert len(_scenarios(ctx({}), 7)) == 2


# ── 판정 분류 ─────────────────────────────────────────────────────────────


def part(name: str, kind: str) -> SimpleNamespace:
    return SimpleNamespace(name=name, phase_type=kind)


def test_verdicts() -> None:
    assert _verdict(True, True, [part("Fe(s)", "Solid")])[0] == "immune"
    assert _verdict(True, False, [part("Fe2O3(s)", "Solid")]) == ("passive", "부동태 — 표면에 Fe₂O₃(s) 생성")
    assert _verdict(True, False, [part("Fe[+2]", "Ion")]) == ("corrosion", "부식 — Fe²⁺(으)로 녹음")
    assert _verdict(False, False, [part("Fe[+3]", "Ion")])[0] == "dissolve"
    assert _verdict(False, True, [part("Fe2O3(s)", "Solid")])[0] == "stable"
    assert _verdict(True, False, [part("Zn[+2]", "Ion"), part("Fe2O3(s)", "Solid")])[0] == "partial"


# ── 하이브리드 에너지 (D15) ───────────────────────────────────────────────


def _upper_immunity(pd, metal: str, pH: float = 0.0) -> float:
    lo, hi = -3.0, 2.0  # 이분 탐색: 금속이 안정한 최고 전위
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if pd.get_stable_entry(pH, mid).name == f"{metal}(s)" else (lo, mid)
    return lo


@needs_pbx("Cu")
@needs_pbx("Ag")
@needs_pbx("Zn")
@pytest.mark.parametrize(("el", "e0", "n"), [("Cu", 0.3419, 2), ("Ag", 0.7996, 1), ("Zn", -0.7618, 2)])
def test_hybrid_immunity_matches_standard_potential(el: str, e0: float, n: int) -> None:
    """면역 상한 = E°(Mⁿ⁺/M) + (0.05916/n)·log10(10⁻⁶) — 실험 이온 에너지를 그대로 쓰므로 ±10 mV."""
    _, pd, _ = pb.diagram([el])
    assert _upper_immunity(pd, el) == pytest.approx(e0 + 0.05916 / n * -6, abs=0.01)


@needs_pbx("Cu")
def test_hybrid_sources_and_mp_comparison() -> None:
    _, pd, src = pb.diagram(["Cu"])
    assert src["Cu(s)"] == pb.ELEMENT and src["CuO(s)"].startswith(pb.EXP_NASA) and src["Cu[+2]"] == pb.EXP_ION
    _, pd_mp, _ = pb.diagram(["Cu"], method="mp")
    assert _upper_immunity(pd_mp, "Cu") < 0.05  # MP 원본은 Cu 를 0 V 근처에서 녹는 것으로 본다 (도감 +0.16 V)


@needs_pbx("Mg")
def test_hybrid_magnesium_hydroxide_passivation() -> None:
    """Mg(OH)₂ 석출: Ksp 5.6e-12, Mg²⁺ 10⁻⁶ M → pH 11.4. MP 원본 에너지는 수산화물을 크게 불안정하게 본다."""
    _, pd, _ = pb.diagram(["Mg"])
    assert pd.get_stable_entry(11.0, 0.0).name == "Mg[+2]"
    assert pd.get_stable_entry(11.8, 0.0).name == "Mg(HO)2(s)"
    _, pd_mp, _ = pb.diagram(["Mg"], method="mp")
    assert pd_mp.get_stable_entry(12.0, 0.0).phase_type == "Ion"


@needs_pbx("Cu")
def test_observed_only_drops_theoretical_phases() -> None:
    """MP 의 이론 상 Cu₂O₃ (ICSD 없음) 는 빼고, 실재하는 상만으로 도감과 같은 Cu 도표가 된다."""
    _, pd, _ = pb.diagram(["Cu"])
    assert {e.name for e in pd.stable_entries if e.phase_type == "Solid"} == {"Cu(s)", "Cu2O(s)", "CuO(s)"}
    ents, _ = pb.hybrid_entries(["Cu"], observed_only=False)
    assert "Cu2O3(s)" in {e.name for e in ents}


@needs_pbx("Fe")
def test_domains_cover_chart() -> None:
    _, pd, _ = pb.diagram(["Fe"])
    doms = pb.domains(pd, [[0, 14], [-2, 2]])
    names = {d["name"] for d in doms}
    assert {"Fe(s)", "Fe[+2]", "Fe2O3(s)"} <= names
    assert all(0 <= x <= 14 and -2 <= y <= 2 for d in doms for x, y in d["vertices"])


# ── 라우터·A4 연계 ────────────────────────────────────────────────────────


def _recipe(components: list[dict]) -> Recipe:
    return Recipe.model_validate({"id": "rcp-t", "name": "t", "components": components})


def test_router_adds_a5_for_metal_in_water() -> None:
    from msl.recipe import route

    water = res("cas:7732-18-5", "H2O", State.LIQUID)
    iron = [res("element:Fe", "Fe"), water]
    calcite = [res("mineral:calcite", "CaCO3"), water]
    comps = lambda cs: [{"ref": str(c.ref), "amount": "1 g"} for c in cs]
    assert "A5" in dict(route(_recipe(comps(iron)), iron))
    assert "A5" not in dict(route(_recipe(comps(calcite)), calcite))


def test_a4_accepts_molar_concentration() -> None:
    """A5 가 용액 pH 를 얻으려고 부르는 A4 — mol/L 을 물 부피로 환산한다 (0.1 M NaOH → pH ≈ 12.9)."""
    from msl.assays.equilibrium import a4

    comps = [Resolved(ref=SubstanceRef.model_validate("cas:1310-73-2"), name="NaOH", formula="NaHO", state=State.AQUEOUS,
                      amount=Quantity.model_validate("0.1 mol/L")),
             Resolved(ref=SubstanceRef.model_validate("cas:7732-18-5"), name="water", formula="H2O", state=State.LIQUID,
                      amount=Quantity.model_validate("1 L"))]
    c = Context(recipe=_recipe([{"ref": "cas:1310-73-2", "amount": "0.1 mol/L", "state": "aqueous"},
                                {"ref": "cas:7732-18-5", "amount": "1 L", "state": "liquid"}]), comps=comps, registry=REG)
    out = a4(c)
    assert out.status.value == "ok" and c.shared["pH"] == pytest.approx(12.9, abs=0.15)
