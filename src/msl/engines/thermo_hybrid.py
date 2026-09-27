"""유한온도 하이브리드 상태도 — NASA 실험 ΔGf(T) 우선 + MP 보충 (A2 v1, 계획서 3.2·3.5절).

좌표계: 모든 엔트리 에너지 = T 에서 NASA 원소 기준상(기체는 1 bar)으로부터의 생성 깁스 에너지 (원소 = 0).

1. NASA 조성: ΔGf(T) = G(화합물,T) − Σ nᵢ G(원소 기준상,T)                        (nasa.formation_gibbs)
   같은 조성에 둘 다 있으면 NASA 값만 쓰고 MP 다형은 뺀다. 단 NASA 의 T 안정형이 기체이고 MP 엔트리가
   실험 기체(G_GASES)가 아닌 고체면 서로 다른 상이므로 둘 다 둔다 (예: SiO 기체 vs 고체).
2. MP 조성:   GibbsComputedStructureEntry = N·(ΔHf,0K + G^δ_SISSO(T)) − Σ nᵢ Gᵢ^pmg(T)  (Bartel 2018)
   - 원소 기준 맞춤: Σ nᵢ [Gᵢ^pmg(T) − Gᵢ^NASA(T)] 를 더한다. 두 원소 표는 같은 관례(G − H_SER(298 K))라
     대부분 수 meV/atom 안에서 같지만 Ba 는 1,300 K 에서 51 meV/atom 다르다.
   - 다형: 조성마다 0 K 바닥 다형만 쓴다. SISSO 는 부피 항 때문에 저밀도 가상 구조를 과안정화한다
     (예: SiO₂ V=58 Å³/atom 구조는 0 K +199 meV/atom 인데 1,273 K 에서 바닥 다형보다 160 meV/atom 낮아짐).
3. 앵커 보정: 원소 기준이 같아도 MP 의 생성에너지는 화합물마다 크게 틀린다 (1,273 K NASA 대비
   CaO −63, CaCO₃ −230, SiO₂ −290 meV/atom — MP2020 ΔHf 오차 + SISSO G^δ 오차). 그래서 MP 에만 있는 조성 X 는
   'NASA·MP 양쪽에 있는 앵커 상 p' 로의 분해 반응에너지만 MP 에서 가져오고 앵커는 NASA 값을 쓴다.
       G(X) = Σ_p c_p G_NASA(p) + ΔG_MP(Σ c_p p → X)
   Pourbaix 가 실험 이온 에너지를 기준 고체로 DFT 에 잇는 방식(Persson 2012)과 같다. 앵커 없는 원소 몫은 보정 0.
   ΔG_MP 는 두 가지:
   - MP-0K (노이만–코프): 분해 생성물이 모두 T 에서 응축상이면 0 K DFT 반응에너지 ΔE₀ 를 그대로 쓴다
     (고체↔고체 반응은 ΔCp·ΔS 가 작다). 온도 의존은 NASA 앵커가 전부 가진다.
   - MP-Gibbs (SISSO): 분해 생성물에 기체(O₂·CO₂ 등)가 있으면 기체 엔트로피가 커서 SISSO ΔG(T) 차이를 쓴다.
   검증 (NASA 에 있는 삼원 산화물 11종 × 1,000/1,300/1,600 K, 이원 산화물로부터의 반응에너지):
   원시 혼합 MAE 94 · SISSO 앵커 29 · 0 K 앵커 11 meV/atom. 원시 혼합은 MP 상을 체계적으로 과안정화한다.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache

import numpy as np
from pymatgen.analysis.compatibility.computed_entries import G_ELEMS
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import GibbsComputedStructureEntry

from msl.engines import nasa

T_MIN, T_MAX = 300.0, 2000.0  # GibbsComputedStructureEntry(SISSO) 적용 범위
NASA, MP_0K, MP_GIBBS, ELEMENT = "NASA", "MP-0K", "MP-Gibbs", "원소 기준"
METHODS = ("nk", "sisso", "raw")  # 앵커 0 K(기본) / 앵커 SISSO / 원시 혼합(비교용)
ENERGY_REFERENCE = ("혼합: NASA 실험 ΔGf(T) 우선 + MP(MP2020 GGA/GGA+U) 반응에너지 — 원소 기준 NASA, "
                    "MP 전용 상은 NASA 앵커 상으로부터의 0 K 반응에너지(고체 생성물) 또는 SISSO ΔG(T)(기체 생성물)")
ANCHOR_MAE = {"nk": 11, "sisso": 29, "raw": 94}  # meV/atom — 위 검증 (삼원 산화물 11종 × 3 온도)
SISSO_MAE = 60  # meV/atom — Bartel 2018 보고 G^δ MAE


@dataclass(frozen=True)
class MPPhase:
    """MP 0 K 바닥 다형 하나. 에너지는 eV/atom."""

    g: float  # SISSO ΔGf(T), 원소 기준 NASA 로 옮긴 값
    e0: float  # 0 K 형성에너지 (MP2020)
    entry_id: str
    experimental: bool = False  # pymatgen G_GASES 실험 기체값(CO₂ 등)을 쓴 엔트리


@dataclass
class HybridPD:
    T: float
    pd: PhaseDiagram
    sources: dict[str, str]  # 약분 화학식 → NASA | MP-0K | MP-Gibbs | 원소 기준
    nasa_species: dict[str, str]  # 약분 화학식 → NASA 화학종 이름 (상 포함, 예: 'CaCO3(cr)', 'Li2CO3(L)')
    compare: dict[str, dict[str, float]] = field(default_factory=dict)  # 앵커 상: NASA vs MP Gibbs, meV/atom
    anchoring: dict[str, dict] = field(default_factory=dict)  # MP 전용 상 → 방법·앵커 분해·에너지
    element_shift: dict[str, float] = field(default_factory=dict)  # 원소 → Gᵢ^pmg − Gᵢ^NASA, meV/atom
    nasa_out_of_range: list[str] = field(default_factory=list)  # NASA 에 있으나 T 에서 유효 구간 밖 → MP 로 채움
    method: str = "nk"

    def best(self, formula: str) -> PDEntry | None:
        """같은 조성 엔트리 중 hull 에 가장 가까운 것."""
        target = Composition(formula).reduced_formula
        same = [e for e in self.pd.all_entries if e.composition.reduced_formula == target]
        return min(same, key=self.pd.get_e_above_hull) if same else None

    def e_hull(self, formula: str) -> float | None:
        """E_hull (eV/atom). 조성이 상태도에 없으면 None."""
        e = self.best(formula)
        return None if e is None else float(self.pd.get_e_above_hull(e))

    def source_of(self, formula: str) -> str | None:
        return self.sources.get(Composition(formula).reduced_formula)


def g_elem_pmg(el: str, T: float) -> float:
    """pymatgen(Bartel) 원소 G(T), eV/atom — GibbsComputedStructureEntry 와 같은 선형 보간."""
    temps = sorted(int(t) for t in G_ELEMS)
    return float(np.interp(T, temps, [G_ELEMS[str(t)][el] for t in temps]))


def element_shift(el: str, T: float) -> float | None:
    """Gᵢ^pmg(T) − Gᵢ^NASA(T), eV/atom. 어느 한쪽 표에 원소가 없으면 None (보정 0)."""
    ref = nasa.element_reference(el, T)
    if ref is None or el not in G_ELEMS[str(int(T_MIN))]:
        return None
    return g_elem_pmg(el, T) - ref[1] / nasa.EV


def element_condensed(el: str, T: float) -> bool:
    """원소 기준상이 T 에서 응축상인가. NASA 에 없는 원소(Au 등 중금속)는 응축상으로 본다."""
    ref = nasa.element_reference(el, T)
    return True if ref is None else ref[0].condensed


@cache
def _nasa_comps() -> tuple[tuple[str, frozenset[str]], ...]:
    return tuple((f, frozenset(str(e) for e in Composition(f).elements)) for f in nasa._by_reduced())


def nasa_gf(elements: list[str], T: float) -> tuple[dict[str, nasa.FormationG], list[str]]:
    """화학계 안 NASA 화합물(2원소 이상)의 ΔGf(T)와, 조성은 있으나 T 에서 유효하지 않은 화학식."""
    els = set(elements)
    out: dict[str, nasa.FormationG] = {}
    missing = []
    for f, fe in _nasa_comps():
        if len(fe) < 2 or not fe <= els:
            continue
        g = nasa.formation_gibbs(f, T)
        if g is None:
            missing.append(f)
        else:
            out[f] = g
    return out, sorted(missing)


def mp_phases(pd0: PhaseDiagram, T: float) -> tuple[dict[str, MPPhase], dict[str, float]]:
    """0 K MP 상태도 → 조성마다 바닥 다형의 0 K 형성에너지와 SISSO ΔGf(T) (원소 기준 NASA). 원소는 제외."""
    shifts = {str(el): element_shift(str(el), T) for el in pd0.elements}
    gs: dict[str, object] = {}
    for e in pd0.all_entries:
        f = e.composition.reduced_formula
        if not e.composition.is_element and (f not in gs or e.energy_per_atom < gs[f].energy_per_atom):
            gs[f] = e
    out: dict[str, MPPhase] = {}
    for f, e in gs.items():
        e0 = pd0.get_form_energy_per_atom(e)
        g = GibbsComputedStructureEntry(e.structure, formation_enthalpy_per_atom=e0, temp=T, correction=0,
                                        data=e.data, entry_id=e.entry_id)
        val = g.energy_per_atom
        if not g.experimental:  # 실험 기체(G_GASES)는 원소 표를 쓰지 않는다
            val += sum(x * (shifts[str(el)] or 0.0) for el, x in e.composition.fractional_composition.items())
        out[f] = MPPhase(float(val), float(e0), str(e.entry_id), g.experimental)
    return out, {k: round(v * 1000, 1) for k, v in shifts.items() if v is not None}


def _entry(f: str, g_per_atom: float, **attr) -> PDEntry:
    c = Composition(f)
    return PDEntry(c, g_per_atom * c.num_atoms, name=f, attribute=attr)


def combine(mp: dict[str, MPPhase], nasa_g: dict[str, nasa.FormationG], elements: list[str], T: float, *,
            method: str = "nk", condensed_elements: dict[str, bool] | None = None,
            missing: list[str] | None = None, element_shift_meV: dict[str, float] | None = None) -> HybridPD:
    """MP·NASA 에너지(eV/atom, 같은 원소 기준)를 합쳐 상태도를 만든다. 네트워크·결정 구조 없이 동작."""
    if method not in METHODS:
        raise ValueError(f"method 는 {METHODS} 중 하나")
    cond_el = condensed_elements or {}
    anchors = {f: m for f, m in mp.items() if f in nasa_g and (nasa_g[f].condensed or m.experimental)}
    compare = {f: {"nasa": round(nasa_g[f].dgf_per_atom * 1000, 1), "mp_gibbs": round(m.g * 1000, 1),
                   "diff": round((nasa_g[f].dgf_per_atom - m.g) * 1000, 1)} for f, m in anchors.items()}
    base = [_entry(el, 0.0, source=ELEMENT) for el in elements]
    ref0 = PhaseDiagram(base + [_entry(f, m.e0) for f, m in anchors.items()])
    refg = PhaseDiagram(base + [_entry(f, m.g) for f, m in anchors.items()])

    def condensed(name: str) -> bool:
        return nasa_g[name].condensed and not mp[name].experimental if name in anchors else cond_el.get(name, True)

    entries = list(base)
    anchoring: dict[str, dict] = {}
    for f, m in mp.items():
        if f in anchors:
            continue
        g, src, via = m.g, MP_GIBBS, {}
        if method != "raw" and not m.experimental:  # NASA 에 없는 실험 기체(G_GASES)는 JANAF 값 그대로
            dec0 = {p.name: a for p, a in ref0.get_decomposition(Composition(f)).items()}
            if method == "nk" and all(condensed(p) for p in dec0):
                via = {p: a for p, a in dec0.items() if p in anchors}
                g = m.e0 + sum(a * (nasa_g[p].dgf_per_atom - anchors[p].e0) for p, a in via.items())
                src = MP_0K
            else:
                via = {p.name: a for p, a in refg.get_decomposition(Composition(f)).items() if p.name in anchors}
                g = m.g + sum(a * (nasa_g[p].dgf_per_atom - anchors[p].g) for p, a in via.items())
            anchoring[f] = {"method": src, "via": {p: round(float(a), 3) for p, a in via.items()},
                            "g_meV": round(float(g) * 1000, 1), "mp_gibbs_meV": round(m.g * 1000, 1)}
        entries.append(_entry(f, g, source=src, entry_id=m.entry_id))
    for f, ng in nasa_g.items():
        entries.append(_entry(f, ng.dgf_per_atom, source=NASA, species=ng.species))

    pd = PhaseDiagram(entries)
    sources: dict[str, str] = {}
    for e in sorted(entries, key=pd.get_e_above_hull, reverse=True):  # NASA 기체·MP 고체가 겹치면 hull 에 가까운 쪽
        sources[e.composition.reduced_formula] = e.attribute["source"]
    return HybridPD(
        T=T, pd=pd, sources=sources, nasa_species={f: g.species for f, g in nasa_g.items()},
        compare=compare, anchoring=anchoring, element_shift=element_shift_meV or {},
        nasa_out_of_range=[f for f in (missing or []) if f in mp], method=method,
    )


def build(pd0: PhaseDiagram, T: float, *, method: str = "nk") -> HybridPD:
    """0 K MP 상태도에서 T(300–2,000 K)의 하이브리드 상태도를 만든다."""
    if not T_MIN <= T <= T_MAX:
        raise ValueError(f"하이브리드 상태도는 {T_MIN:.0f}–{T_MAX:.0f} K 만 (요청 {T:.0f} K)")
    elements = sorted(str(e) for e in pd0.elements)
    mp, shifts = mp_phases(pd0, T)
    nasa_g, missing = nasa_gf(elements, T)
    return combine(mp, nasa_g, elements, T, method=method, missing=missing, element_shift_meV=shifts,
                   condensed_elements={el: element_condensed(el, T) for el in elements})


def nasa_only_pd(elements: list[str], T: float) -> PhaseDiagram:
    """NASA thermo.inp 만으로 만든 T 상태도 (평형 트랙 교차검증용). 원소는 기준상 = 0, 기체는 1 bar."""
    nasa_g, _ = nasa_gf(elements, T)
    return PhaseDiagram([_entry(el, 0.0, source=ELEMENT) for el in elements]
                        + [_entry(f, g.dgf_per_atom, source=NASA, species=g.species) for f, g in nasa_g.items()])
