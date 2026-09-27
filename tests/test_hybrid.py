"""A2 v1 하이브리드 열역학 — NASA ΔGf(T) 계산, 원소 기준상, 하이브리드 엔트리 구성, 교차검증.

네트워크 없이 돈다. MP 엔트리가 필요한 시험은 data/cache/mp 캐시가 있을 때만 (없으면 skip).
thermo.inp 원본(data/raw, 저장소 밖)에만 있는 화학종·상 이름을 보는 시험도 원본이 없으면 skip.
"""

from __future__ import annotations

import pytest
from pymatgen.analysis.interface_reactions import InterfacialReactivity
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Composition

from msl.assays import mp_assays
from msl.assays.base import Context
from msl.engines import mp, nasa
from msl.engines import thermo_hybrid as th
from msl.registry.load import load_registry
from msl.resolve import Resolved
from msl.schema.quantity import Quantity
from msl.schema.recipe import Recipe, State
from msl.schema.refs import SubstanceRef
from msl.schema.result import Status

KJ = nasa.EV / 1000  # kJ/mol per eV
RAW = pytest.mark.skipif(not nasa.THERMO_INP.exists(), reason="thermo.inp 원본 없음 (data/raw)")
REG = load_registry()


def cached(*elements: str) -> bool:
    return (mp.CACHE / f"entries_{'-'.join(sorted(elements))}_{mp.THERMO_TYPE}.pkl").exists()


def needs_mp(*elements: str):
    return pytest.mark.skipif(not cached(*elements), reason=f"MP 캐시 없음: {'-'.join(sorted(elements))}")


def dgf(formula: str, T: float) -> float:
    """ΔGf(T), kJ/mol (화학식 단위)."""
    return nasa.formation_gibbs(formula, T).dgf_per_atom * Composition(formula).num_atoms * KJ


def comp(ref: str, formula: str, amount: str) -> Resolved:
    return Resolved(ref=SubstanceRef.model_validate(ref), name=ref, formula=formula, groups=[],
                    state=State.SOLID, amount=Quantity.model_validate(amount), props={})


def ctx(comps: list[Resolved], T: str) -> Context:
    recipe = Recipe.model_validate({
        "id": "rcp-test", "name": "t", "conditions": {"T": T},
        "components": [{"ref": str(c.ref), "amount": str(c.amount), "state": "solid"} for c in comps],
    })
    return Context(recipe=recipe, comps=comps, registry=REG)


# ── thermo.inp 파서 ───────────────────────────────────────────────────────

SNIPPET = """\
! 주석
thermo
    200.00   1000.00   6000.00  20000.   9/8/2021
CaCO3(cr)         Hexagonal Gurvich,1996a pt1 p483 pt2 p376.
 2 tpis96 CA  1.00C   1.00O   3.00    0.00    0.00 1  100.0869000   -1206600.000
    300.000    500.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        14480.499
-1.329862425D+07 2.500517168D+05-1.907177167D+03 7.616666270D+00-1.655870860D-02
 1.879382277D-05-8.720713270D-09                -1.271058215D+06 9.957497560D+03
    500.000   1603.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        14480.499
-2.583555736D+05 0.000000000D+00 1.197256363D+01 3.263812299D-03 0.000000000D+00
 0.000000000D+00 0.000000000D+00                -1.497009803D+05-5.961133653D+01
CaCO3+            가짜 이온 — 전자(E)가 있으면 빠져야 한다
 1 tpis96 CA  1.00C   1.00O   3.00E  -1.00    0.00 0  100.0869000   -1206600.000
    300.000    500.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0        14480.499
-1.329862425D+07 2.500517168D+05-1.907177167D+03 7.616666270D+00-1.655870860D-02
 1.879382277D-05-8.720713270D-09                -1.271058215D+06 9.957497560D+03
END PRODUCTS
"""


def test_parse_thermo_inp_record() -> None:
    (sp,) = nasa._parse_thermo_inp(SNIPPET)  # 이온은 빠진다
    assert sp.name == "CaCO3(cr)" and sp.condensed
    assert sp.comp.almost_equals(Composition("CaCO3"))
    assert (sp.tmin, sp.tmax) == (300.0, 1603.0)
    assert sp.h(298.15) / 1000 == pytest.approx(-1206.6, abs=0.5)  # H(298.15) = ΔHf (원소 기준 0)
    assert sp.s(298.15) == pytest.approx(91.7, abs=0.5)  # 방해석 S°298 (JANAF 91.71 J/mol·K)
    assert sp.valid(298.15) and not sp.valid(1700)  # 하한 300 K 는 2 K 허용


def test_nasa9_matches_cantera_bundle() -> None:
    """같은 화학종(CO₂·O₂)의 H·S 가 Cantera 번들과 0.1% 안에서 같다 — 계수 해석이 맞는지."""
    bundle = {s.name: s for s in nasa.gas_species()}
    for name in ("CO2", "O2"):
        mine = next(s for s in nasa.nasa9_species() if s.name == name and not s.condensed)
        for T in (500.0, 1273.15, 1800.0):
            assert mine.h(T) == pytest.approx(bundle[name].thermo.h(T) / 1000, rel=1e-3, abs=100)
            assert mine.s(T) == pytest.approx(bundle[name].thermo.s(T) / 1000, rel=1e-3)


# ── ΔGf(T) 와 원소 기준상 ────────────────────────────────────────────────


@pytest.mark.parametrize(("formula", "T", "janaf"), [
    ("CO2", 298.15, -394.4), ("CO2", 1300.0, -396.1),
    ("H2O", 298.15, -237.1), ("H2O", 1000.0, -192.6),  # 298 K 는 액체가 안정형, 1,000 K 는 기체
])
def test_formation_gibbs_matches_janaf(formula: str, T: float, janaf: float) -> None:
    """ΔGf(T) = G(화합물) − Σ nᵢ G(원소 기준상) 이 JANAF 표와 같다 (H 기준상이 D₂ 로 잡히면 틀린다)."""
    assert dgf(formula, T) == pytest.approx(janaf, abs=1.0)


def test_element_reference_is_zero_formation() -> None:
    for el in ("O", "C", "Si"):
        g = nasa.formation_gibbs(el, 1273.15)
        assert g.dgf_per_atom == pytest.approx(0.0, abs=1e-12)


@RAW
@pytest.mark.parametrize(("el", "T", "name"), [
    ("Ca", 1073.15, "Ca(b)"), ("Ca", 1273.15, "Ca(L)"), ("C", 1273.15, "C(gr)"), ("O", 1273.15, "O2"),
    ("H", 1000.0, "H2"), ("Mg", 1673.15, "Mg"),  # Mg 끓는점 1,363 K 위 → 1 bar 단원자 기체 (JANAF 관례)
])
def test_element_reference_phase(el: str, T: float, name: str) -> None:
    assert nasa.element_reference(el, T)[0].name == name


@RAW
def test_stable_form_follows_melting() -> None:
    """Li₂CO₃ 녹는점 1,005 K — 아래는 결정, 위는 액체를 고른다."""
    assert nasa.stable_form("Li2CO3", 900.0)[0].name == "Li2CO3(cr)"
    assert nasa.stable_form("Li2CO3", 1173.15)[0].name == "Li2CO3(L)"


@RAW
def test_calcite_decomposition_from_thermo_inp() -> None:
    """thermo.inp 는 CaCO₃(cr) 를 1,603 K 까지 준다 (번들은 1,200 K). ΔG(CaCO₃ → CaO + CO₂) 부호가 1 bar CO₂
    분해온도(실측 약 1,171 K) 앞뒤로 바뀐다."""
    def rxn(T: float) -> float:
        return dgf("CaO", T) + dgf("CO2", T) - dgf("CaCO3", T)

    assert rxn(1073.15) > 5 and rxn(1273.15) < -5
    assert rxn(1120.0) > 0 > rxn(1200.0)


@RAW
def test_element_shift_pmg_vs_nasa() -> None:
    """두 원소 G 표는 같은 관례 — O·C·Si 는 수 meV/atom 안. Ba 는 크게 다르다 (그래서 기준을 맞춘다)."""
    for el in ("O", "C", "Si", "Ca"):
        assert abs(th.element_shift(el, 1300.0)) < 0.01
    assert th.element_shift("Ba", 1300.0) < -0.03
    assert th.element_shift("Au", 1300.0) is None  # NASA 에 Au 없음 → 보정 0


# ── 하이브리드 엔트리 구성 (가짜 에너지) ──────────────────────────────────


def fg(formula: str, dgf_per_atom: float, condensed: bool = True) -> nasa.FormationG:
    return nasa.FormationG(Composition(formula).reduced_formula, formula + ("(cr)" if condensed else ""),
                           condensed, 1273.15, dgf_per_atom, {})


MP_FAKE = {  # eV/atom — g: SISSO(원소 기준 NASA), e0: 0 K 형성에너지. MP 는 NASA 보다 과안정 (실제 경향)
    "CaO": th.MPPhase(g=-2.65, e0=-3.30, entry_id="mp-cao"),
    "SiO2": th.MPPhase(g=-2.65, e0=-3.27, entry_id="mp-sio2"),
    "CaSiO3": th.MPPhase(g=-2.80, e0=-3.45, entry_id="mp-casio3"),
}
NASA_FAKE = {"CaO": fg("CaO", -2.59), "SiO2": fg("SiO2", -2.36)}


def test_nasa_replaces_mp_on_same_composition() -> None:
    h = th.combine(MP_FAKE, NASA_FAKE, ["Ca", "O", "Si"], 1273.15)
    same = [e for e in h.pd.all_entries if e.composition.reduced_formula == "CaO"]
    assert len(same) == 1 and same[0].attribute["source"] == th.NASA
    assert h.sources == {"Ca": th.ELEMENT, "O2": th.ELEMENT, "Si": th.ELEMENT,  # pymatgen 약분식: O → O2
                         "CaO": th.NASA, "SiO2": th.NASA, "CaSiO3": th.MP_0K}
    assert h.compare["CaO"]["diff"] == pytest.approx(60.0)  # NASA − MP, meV/atom
    assert h.compare["SiO2"]["diff"] == pytest.approx(290.0)


@pytest.mark.parametrize(("method", "expected"), [
    ("nk", -3.45 + 0.4 * (-2.59 + 3.30) + 0.6 * (-2.36 + 3.27)),  # 0 K 반응에너지 + NASA 앵커
    ("sisso", -2.80 + 0.4 * (-2.59 + 2.65) + 0.6 * (-2.36 + 2.65)),  # SISSO 반응에너지 + NASA 앵커
    ("raw", -2.80),  # 원시 혼합 (비교용)
])
def test_anchor_methods(method: str, expected: float) -> None:
    h = th.combine(MP_FAKE, NASA_FAKE, ["Ca", "O", "Si"], 1273.15, method=method)
    assert h.best("CaSiO3").energy_per_atom == pytest.approx(expected)


def test_nk_anchor_keeps_0k_reaction_energy() -> None:
    """0 K 앵커는 CaO + SiO₂ → CaSiO₃ 반응에너지를 0 K DFT 값 그대로 둔다 (MP 생성에너지 오차는 상쇄)."""
    h = th.combine(MP_FAKE, NASA_FAKE, ["Ca", "O", "Si"], 1273.15)
    e = {f: h.best(f).energy_per_atom for f in ("CaO", "SiO2", "CaSiO3")}
    assert e["CaSiO3"] - (0.4 * e["CaO"] + 0.6 * e["SiO2"]) == pytest.approx(-3.45 - (0.4 * -3.30 + 0.6 * -3.27))
    assert h.anchoring["CaSiO3"]["via"] == {"CaO": 0.4, "SiO2": 0.6}


def test_gas_product_falls_back_to_sisso() -> None:
    """분해 생성물에 O₂ 가 끼는 상(과산화물)은 기체 엔트로피가 커서 0 K 반응에너지를 쓰지 않는다."""
    mp_ = {**MP_FAKE, "CaO2": th.MPPhase(g=-1.90, e0=-2.40, entry_id="mp-cao2")}
    h = th.combine(mp_, NASA_FAKE, ["Ca", "O", "Si"], 1273.15, condensed_elements={"O": False})
    assert h.anchoring["CaO2"]["method"] == th.MP_GIBBS
    assert h.anchoring["CaO2"]["via"] == {"CaO": pytest.approx(2 / 3, abs=1e-3)}


def test_gas_and_solid_of_same_formula_both_kept() -> None:
    """NASA 의 T 안정형이 기체(SiO)이고 MP 엔트리가 고체면 서로 다른 상 — 둘 다 두고 앵커로 쓰지 않는다.
    MP 엔트리가 실험 기체값(CO₂)이면 같은 상이라 NASA 로 바꾼다."""
    mp_ = {**MP_FAKE, "SiO": th.MPPhase(g=-1.50, e0=-1.90, entry_id="mp-sio"),
           "CO2": th.MPPhase(g=-1.37, e0=-1.40, entry_id="mp-co2", experimental=True)}
    nasa_ = {**NASA_FAKE, "SiO": fg("SiO", -1.08, condensed=False), "CO2": fg("CO2", -1.3686, condensed=False)}
    h = th.combine(mp_, nasa_, ["C", "Ca", "O", "Si"], 1273.15)
    sio = {e.attribute["source"] for e in h.pd.all_entries if e.composition.reduced_formula == "SiO"}
    co2 = {e.attribute["source"] for e in h.pd.all_entries if e.composition.reduced_formula == "CO2"}
    assert sio == {th.NASA, th.MP_0K} and "SiO" not in h.compare
    assert co2 == {th.NASA} and "CO2" in h.compare


def test_nasa_only_hybrid_decomposes_calcite_at_1273() -> None:
    """MP 상 없이 NASA 만으로도 (T 에서 유효하면) 방해석은 1,273 K 에서 CaO + CO₂ 쪽."""
    T = 1273.15
    nasa_g, _ = th.nasa_gf(["C", "Ca", "O"], T)
    if "CaCO3" not in nasa_g:
        pytest.skip("CaCO3 가 이 온도에서 유효한 NASA 데이터 없음 (번들 대체)")
    h = th.combine({}, nasa_g, ["C", "Ca", "O"], T)
    assert h.e_hull("CaCO3") > 0.01
    assert {e.composition.reduced_formula for e in h.pd.get_decomposition(Composition("CaCO3"))} == {"CaO", "CO2"}


# ── 교차검증 ──────────────────────────────────────────────────────────────


def test_cross_check_flags_only_real_mismatch() -> None:
    """NASA 평형은 방해석이 분해된다는데 하이브리드가 안정으로 보면 불일치(WARNING 근거), 같은 판정이면 일치."""
    T = 1273.15
    c = ctx([comp("mineral:calcite", "CaCO3", "1 mol"), comp("mineral:quartz", "SiO2", "1 mol")], "1000 °C")
    overstable = th.combine({"CaCO3": th.MPPhase(g=-2.10, e0=-2.69, entry_id="mp-cal")},
                            {"CaO": fg("CaO", -2.59), "CO2": fg("CO2", -1.3686, condensed=False)},
                            ["C", "Ca", "O"], T)
    bad = mp_assays._cross_check(c, overstable, None, None)
    assert bad["mismatch"] and bad["checks"][0]["mismatch"]

    nasa_g, _ = th.nasa_gf(["C", "Ca", "O", "Si"], T)
    if "CaCO3" not in nasa_g:
        pytest.skip("CaCO3 가 이 온도에서 유효한 NASA 데이터 없음 (번들 대체)")
    good = mp_assays._cross_check(c, th.combine({}, nasa_g, ["C", "Ca", "O", "Si"], T), None, None)
    assert not good["mismatch"]
    assert "일치" in good["checks"][0]["note"]


# ── MP 캐시가 있을 때: 실제 계 ────────────────────────────────────────────


def hybrid(elements: list[str], T: float, method: str = "nk") -> th.HybridPD:
    return th.build(PhaseDiagram(mp.entries_in_chemsys(elements)), T, method=method)


def lowest(h: th.HybridPD, r1: str, r2: str) -> tuple[float, str]:
    ir = InterfacialReactivity(Composition(r1), Composition(r2), h.pd, norm=True, use_hull_energy=False)
    _, _, e, rxn, _ = min(ir.get_kinks(), key=lambda k: k[2])
    return float(e), str(rxn)


@RAW
@needs_mp("C", "Ca", "O", "Si")
def test_a2_limestone_quartz_is_favorable() -> None:
    """회귀: v0(MP+Gibbs 근사)은 1,273 K 에서 '반응 없음'. 하이브리드는 CaSiO₃/Ca₂SiO₄ + CO₂ 가 유리."""
    out = mp_assays.a2(ctx([comp("mineral:calcite", "CaCO3", "1 mol"), comp("mineral:quartz", "SiO2", "1 mol")],
                           "1000 °C"))
    assert out.status is Status.OK and not out.warnings
    main = out.data["curves"][-1]
    assert "하이브리드" in main["label"]
    assert min(k["e"] for k in main["kinks"]) < -50
    assert set(out.data["products_at_recipe"]) == {"CaSiO3", "CO2"}
    assert out.data["sources_by_phase"]["CaCO3"] == th.NASA
    assert out.data["sources_by_phase"]["CaSiO3"] in (th.MP_0K, th.MP_GIBBS)
    assert out.energy_reference and "혼합" in out.energy_reference
    assert any("출처 혼합" in c for c in out.caveats)


@RAW
@needs_mp("C", "Ca", "O", "Si")
def test_calcite_e_hull_sign_in_hybrid() -> None:
    assert hybrid(["C", "Ca", "O", "Si"], 1073.15).e_hull("CaCO3") < 1e-3
    assert hybrid(["C", "Ca", "O", "Si"], 1273.15).e_hull("CaCO3") > 0.01


@RAW
@needs_mp("Ba", "C", "O", "Ti")
def test_baco3_tio2_gives_batio3() -> None:
    h = hybrid(["Ba", "C", "O", "Ti"], 1373.15)
    assert h.e_hull("BaTiO3") < 1e-3  # SISSO 앵커는 BaTiO₃ 를 hull 위 11 meV/atom 로 밀어낸다 — 0 K 앵커로 해결
    mix = Composition("BaCO3").fractional_composition * 0.625 + Composition("TiO2").fractional_composition * 0.375
    assert {e.composition.reduced_formula for e in h.pd.get_decomposition(mix)} == {"BaTiO3", "CO2"}
    assert lowest(h, "BaCO3", "TiO2")[0] < -0.05


@RAW
@needs_mp("Al", "Mg", "O")
def test_mgo_alumina_spinel() -> None:
    e, rxn = lowest(hybrid(["Al", "Mg", "O"], 1673.15), "MgO", "Al2O3")
    assert e < -0.03 and "MgAl2O4" in rxn


@needs_mp("Au", "O", "Si")
def test_gold_quartz_no_reaction() -> None:
    assert lowest(hybrid(["Au", "O", "Si"], 1273.15), "Au", "SiO2")[0] > -0.001


BENCH = [  # NASA 에 있는 삼원 산화물 ← 이원 산화물
    (["Mg", "O", "Si"], "MgSiO3", "MgO", "SiO2"), (["Mg", "O", "Si"], "Mg2SiO4", "MgO", "SiO2"),
    (["Mg", "O", "Ti"], "MgTiO3", "MgO", "TiO2"), (["Al", "O", "Si"], "Al2SiO5", "Al2O3", "SiO2"),
    (["Al", "Li", "O"], "LiAlO2", "Li2O", "Al2O3"), (["Al", "Mg", "O"], "MgAl2O4", "MgO", "Al2O3"),
]


@RAW
def test_anchor_benchmark_against_nasa_ternaries() -> None:
    """방법 선택의 근거: NASA 삼원 산화물 생성 반응에너지(이원 산화물 기준)를 NASA 값 없이 MP 로 예측해 비교.
    NASA 삼원 산화물 엔트리를 빼고(=MP 전용 상으로 취급) 하이브리드를 만든다."""
    T = 1300.0
    errs: dict[str, list[float]] = {m: [] for m in th.METHODS}
    for els, x, a, b in BENCH:
        if not cached(*els):
            continue
        pd0 = PhaseDiagram(mp.entries_in_chemsys(els))
        mp_, _ = th.mp_phases(pd0, T)
        nasa_g, _ = th.nasa_gf(els, T)
        truth = nasa_g.pop(x).dgf_per_atom
        for m in th.METHODS:
            h = th.combine(mp_, nasa_g, els, T, method=m)
            ca, cb = Composition(a), Composition(b)
            dec = PhaseDiagram([th._entry(a, 0), th._entry(b, 0), *(th._entry(e, 0) for e in els)]).get_decomposition(
                Composition(x))
            ref = sum(amt * nasa_g[e.name].dgf_per_atom for e, amt in dec.items() if e.name in (ca.reduced_formula, cb.reduced_formula))
            pred = h.best(x).energy_per_atom - ref
            errs[m].append(abs(pred - (truth - ref)))
    if not errs["nk"]:
        pytest.skip("벤치마크 계 MP 캐시 없음")
    mae = {m: sum(v) / len(v) for m, v in errs.items()}
    assert mae["nk"] < 0.02 and mae["nk"] < mae["sisso"] < mae["raw"]
