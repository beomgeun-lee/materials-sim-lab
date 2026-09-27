"""NASA 열역학 데이터 평형 계산 — Cantera 번들 nasa_gas.yaml(748종)·nasa_condensed.yaml(382종).

평형 트랙(T)의 v0 엔진이다. 레시피 원소로만 이루어진 모든 화학종을 넣어 평형을 푼다.
DB 에 없는 화학종은 결과에 나올 수 없고, 응축상은 화학종마다 유효 온도 구간이 있다.

생성 깁스 에너지 ΔGf(T)(A2 하이브리드용)는 원본 thermo.inp 를 직접 읽어 계산한다. Cantera 번들보다
화학종·온도 구간이 넓다(예: CaCO3(cr) 1,603 K 까지, BaCO3·TiO2(cr) 포함). 원본이 없으면 번들로 대신한다.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from functools import cache
from importlib.metadata import version

import cantera as ct
from pymatgen.core import Composition

from msl.env import DATA_DIR

ENERGY_REFERENCE = "NASA 열화학 (298.15 K 원소 기준)"
THERMO_INP = DATA_DIR / "raw" / "nasa-cea-thermo" / "cantera-3.2.0" / "thermo.inp"
R = 8.314462618  # J/(mol·K) — thermo.inp 계수가 이 값 기준 (NASA TP-2002-211556)
EV = 96_485.33212  # J/mol per eV
T_TOL = 2.0  # K — 유효 구간 경계 허용치 (2021 개정에서 하한을 298.15 → 300 K 로 올린 화학종이 많음)


class OutOfRange(ValueError):
    """화학종은 있으나 요청 온도가 NASA 데이터 유효 구간 밖일 때."""

    def __init__(self, formula: str, name: str, tmin: float, tmax: float):
        super().__init__(f"{formula}({name}) 데이터 유효 구간 {tmin:.0f}–{tmax:.0f} K 밖")
        self.formula, self.name, self.tmin, self.tmax = formula, name, tmin, tmax


def engine_version() -> str:
    return version("cantera")


def _comp(sp: ct.Species) -> Composition:
    return Composition({el: n for el, n in sp.composition.items() if el != "E"})


def _elements(comp: Composition) -> set[str]:
    return {str(e) for e in comp.elements}


@cache
def gas_species() -> tuple[ct.Species, ...]:
    return tuple(s for s in ct.Species.list_from_file("nasa_gas.yaml") if "E" not in s.composition)


@cache
def condensed_species() -> tuple[ct.Species, ...]:
    return tuple(ct.Species.list_from_file("nasa_condensed.yaml"))


def _same(sp: ct.Species, formula: str) -> bool:
    """조성이 정확히 같은가. 약분하지 않는다 (O₂ ≠ O)."""
    return _comp(sp).almost_equals(Composition(formula))


def find_gas(formula: str) -> ct.Species | None:
    hits = sorted((s for s in gas_species() if _same(s, formula)), key=lambda s: ("," in s.name, len(s.name)))
    return hits[0] if hits else None


def find_condensed(formula: str, T: float | None = None) -> list[ct.Species]:
    """같은 조성의 응축상(고체 다형·액체). T 를 주면 그 온도에서 유효한 것만."""
    hits = [s for s in condensed_species() if _same(s, formula)]
    if T is not None:
        hits = [s for s in hits if s.thermo.min_temp <= T <= s.thermo.max_temp]
    return hits


def gas_equilibrium(moles: dict[str, float], T: float, P: float, mode: str = "TP") -> ct.Solution:
    """이상기체 혼합물 평형. moles 는 {화학식: mol}. mode='HP' 면 단열(엔탈피 일정) 평형."""
    elements = set().union(*(_elements(Composition(f)) for f in moles))
    species = [s for s in gas_species() if _elements(_comp(s)) <= elements]
    gas = ct.Solution(thermo="ideal-gas", species=species)
    x: dict[str, float] = {}
    missing = []
    for formula, n in moles.items():
        sp = find_gas(formula)
        if sp is None:
            missing.append(formula)
        else:
            x[sp.name] = x.get(sp.name, 0.0) + n
    if missing:
        raise KeyError(", ".join(missing))
    gas.TPX = T, P, x
    gas.equilibrate(mode)
    return gas


@dataclass
class Equilibrium:
    T: float
    P: float
    condensed: dict[str, float]  # 응축상 이름 → mol
    gas: dict[str, float]  # 기체 화학종 → 몰분율 (상위)
    gas_moles: float
    n_species: int


def multiphase_equilibrium(moles: dict[str, float], T: float, P: float) -> Equilibrium:
    """고체·기체 혼합물 평형 (Cantera Mixture). moles 는 {화학식: mol}.

    각 응축상은 순수 고정 조성상으로 넣는다(고용체 없음). 반응물이 T 에서 유효한 응축상
    데이터가 없으면 OutOfRange 를 낸다.
    """
    elements = set().union(*(_elements(Composition(f)) for f in moles))
    gases = [s for s in gas_species() if _elements(_comp(s)) <= elements]
    solids = [
        s for s in condensed_species()
        if _elements(_comp(s)) <= elements and s.thermo.min_temp <= T <= s.thermo.max_temp
    ]
    start_gas: dict[str, float] = {}
    start_solid: dict[str, float] = {}
    for formula, n in moles.items():
        valid = [s for s in solids if _same(s, formula)]
        if valid:
            start_solid[valid[0].name] = start_solid.get(valid[0].name, 0.0) + n
            continue
        anywhere = find_condensed(formula)
        if anywhere:
            s = anywhere[0]
            raise OutOfRange(formula, s.name, s.thermo.min_temp, s.thermo.max_temp)
        g = find_gas(formula)
        if g is None:
            raise KeyError(formula)
        start_gas[g.name] = start_gas.get(g.name, 0.0) + n

    gas = ct.Solution(thermo="ideal-gas", species=gases)
    gas.TPX = T, P, (start_gas or {gas.species_names[0]: 1.0})
    phases: list[tuple[ct.Solution, float]] = [(gas, sum(start_gas.values()) or 1e-10)]
    for sp in solids:
        ph = ct.Solution(thermo="fixed-stoichiometry", species=[sp])
        ph.TP = T, P
        phases.append((ph, start_solid.get(sp.name, 0.0)))
    mix = ct.Mixture(phases)
    mix.T, mix.P = T, P
    try:
        mix.equilibrate("TP", solver="gibbs", max_steps=5000)
    except ct.CanteraError:
        mix.equilibrate("TP", solver="vcs", max_steps=5000)
    condensed = {
        sp.name: float(mix.phase_moles(i + 1)) for i, sp in enumerate(solids) if mix.phase_moles(i + 1) > 1e-7
    }
    g_moles = float(mix.phase_moles(0))
    top = sorted(zip(gas.species_names, gas.X), key=lambda p: -p[1])[:6] if g_moles > 1e-7 else []
    return Equilibrium(
        T=T, P=P, condensed=condensed, gas={k: float(v) for k, v in top if v > 1e-4},
        gas_moles=g_moles, n_species=len(gases) + len(solids),
    )


# ── 생성 깁스 에너지 ΔGf(T) — thermo.inp NASA 9계수 ───────────────────────


@dataclass(frozen=True)
class Nasa9:
    """NASA 다항식 화학종 하나 (thermo.inp 9계수; 번들 대체 시 7계수도 같은 표현).
    H 는 298.15 K 원소 기준상 = 0, S 는 제3법칙 절대값 (1 bar)."""

    name: str
    comp: Composition
    condensed: bool
    ranges: tuple[tuple[float, float, tuple[float, ...], tuple[float, ...], float, float], ...]
    # 구간마다 (Tmin, Tmax, T 지수들, 계수 a, b1, b2):  Cp/R = Σ aₖ T^eₖ

    @property
    def tmin(self) -> float:
        return self.ranges[0][0]

    @property
    def tmax(self) -> float:
        return self.ranges[-1][1]

    def valid(self, T: float) -> bool:
        return self.tmin - T_TOL <= T <= self.tmax + T_TOL

    def _range(self, T: float):
        for r in self.ranges:
            if T <= r[1]:
                return r
        return self.ranges[-1]

    def h(self, T: float) -> float:
        """몰 엔탈피 J/mol."""
        _, _, ex, a, b1, _ = self._range(T)
        s = sum(ak * (math.log(T) if e == -1 else T ** (e + 1) / (e + 1)) for e, ak in zip(ex, a))
        return R * (s + b1)

    def s(self, T: float) -> float:
        """몰 엔트로피 J/(mol·K), 1 bar."""
        _, _, ex, a, _, b2 = self._range(T)
        s = sum(ak * (math.log(T) if e == 0 else T ** e / e) for e, ak in zip(ex, a))
        return R * (s + b2)

    def g(self, T: float) -> float:
        """몰 깁스 에너지 G = H − TS, J/mol (1 bar)."""
        return self.h(T) - T * self.s(T)


def _d(x: str) -> float:
    return float(x.replace("D", "E"))


def _parse_thermo_inp(text: str) -> list[Nasa9]:
    """thermo.inp 의 생성물(products) 구역을 읽는다. 이온(전자 포함)·반응물 전용 화학종은 제외.

    형식: NASA TP-2002-211556 부록 — 레코드 1 이름, 2 구간 수·조성·상 플래그, 구간마다 3(온도·지수)·4·5(계수).
    """
    lines = text.splitlines()
    i = next(k for k, ln in enumerate(lines) if ln.strip().lower() == "thermo") + 2
    out: list[Nasa9] = []
    while i < len(lines) and not lines[i].startswith("END"):
        name, rec2 = lines[i].split()[0], lines[i + 1]
        n = int(rec2[0:2])
        els: dict[str, float] = {}
        for k in range(5):
            f = rec2[10 + 8 * k: 18 + 8 * k]
            sym, cnt = f[:2].strip(), float(f[2:] or 0)
            if sym and cnt:
                els[sym[0] + sym[1:].lower()] = cnt
        ranges = []
        for j in range(n):
            r3, r4, r5 = lines[i + 2 + 3 * j: i + 5 + 3 * j]
            ex = tuple(float(r3[23 + 5 * k: 28 + 5 * k]) for k in range(int(r3[22])))
            a = [_d(r4[16 * k: 16 * k + 16]) for k in range(5)] + [_d(r5[0:16]), _d(r5[16:32])]
            ranges.append((float(r3[0:11]), float(r3[11:22]), ex, tuple(a[:len(ex)]), _d(r5[48:64]), _d(r5[64:80])))
        i += 2 + 3 * n if n else 3
        if not ranges or {"E", "D", "T"} & els.keys():  # 이온(전자)·중수소·삼중수소 제외
            continue
        try:
            comp = Composition(els)
        except ValueError:
            continue
        out.append(Nasa9(name, comp, rec2[50:52].strip() not in ("", "0"), tuple(ranges)))
    return out


def _from_cantera(sp: ct.Species, condensed: bool) -> Nasa9 | None:
    """Cantera 번들 화학종 → 같은 표현. NASA7 은 지수 0–4 + (a6, a7) 로, NASA9 은 −2–4 + (b1, b2) 로 옮긴다."""
    th = sp.input_data.get("thermo", {})
    n = {"NASA7": 5, "NASA9": 7}.get(th.get("model"))
    if n is None or {"E", "D", "T"} & sp.composition.keys():
        return None
    t = th["temperature-ranges"]
    ex = (0.0, 1.0, 2.0, 3.0, 4.0) if n == 5 else (-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0)
    ranges = tuple((t[k], t[k + 1], ex, tuple(d[:n]), d[n], d[n + 1]) for k, d in enumerate(th["data"]))
    return Nasa9(sp.name, _comp(sp), condensed, ranges)


@cache
def nasa9_species() -> tuple[Nasa9, ...]:
    """ΔGf 계산용 화학종 목록. 원본 thermo.inp 우선, 없으면 Cantera 번들."""
    if THERMO_INP.exists():
        return tuple(_parse_thermo_inp(THERMO_INP.read_text(encoding="latin-1")))
    conv = [_from_cantera(s, False) for s in gas_species()] + [_from_cantera(s, True) for s in condensed_species()]
    return tuple(s for s in conv if s is not None)


def nasa9_source() -> str:
    return "NASA CEA thermo.inp (2021-09 개정, Cantera 3.2.0 동봉 원본)" if THERMO_INP.exists() else "Cantera nasa_*.yaml"


@cache
def _by_reduced() -> dict[str, tuple[Nasa9, ...]]:
    idx: dict[str, list[Nasa9]] = {}
    with warnings.catch_warnings():  # 비활성 기체(He·Ne·Ar)는 전기음성도가 없다는 pymatgen 경고 — 무해
        warnings.simplefilter("ignore", UserWarning)
        for s in nasa9_species():
            idx.setdefault(s.comp.reduced_formula, []).append(s)
    return {k: tuple(v) for k, v in idx.items()}


def stable_form(formula: str, T: float) -> tuple[Nasa9, float] | None:
    """조성(약분 기준)이 같은 화학종 중 T 에서 유효하고 원자당 G 가 가장 낮은 것 → (화학종, G J/mol-atom).

    고체 다형·액체·기체(1 bar)를 모두 비교하므로 T 에서의 안정 형태를 고른다 (예: 1,173 K Li2CO3 → 액체).
    """
    cands = [s for s in _by_reduced().get(Composition(formula).reduced_formula, ()) if s.valid(T)]
    if not cands:
        return None
    best = min(cands, key=lambda s: s.g(T) / s.comp.num_atoms)
    return best, best.g(T) / best.comp.num_atoms


def element_reference(el: str, T: float) -> tuple[Nasa9, float] | None:
    """원소 기준상 — T 에서 원자당 G 가 가장 낮은 순수 원소 화학종 (예: Ca(b)·Ca(L), C(gr), O2).

    JANAF/NASA 관례(끓는점 위는 1 bar 단원자 기체)와 같다. 2,000 K 이하에서 O₂·H₂·N₂ 는 해리되지 않는다.
    """
    return stable_form(el, T)


@dataclass(frozen=True)
class FormationG:
    formula: str  # 약분 화학식
    species: str  # NASA 화학종 이름 (상 포함, 예: 'CaCO3(cr)')
    condensed: bool
    T: float
    dgf_per_atom: float  # eV/atom — T 의 원소 기준상으로부터의 생성 깁스 에너지
    refs: dict[str, str]  # 원소 → 기준상 이름


def formation_gibbs(formula: str, T: float) -> FormationG | None:
    """ΔGf(T) = G(화합물, T) − Σ nᵢ G(원소 기준상, T), 원자당 eV. 데이터·유효 구간이 없으면 None."""
    comp = Composition(formula)
    hit = stable_form(formula, T)
    if hit is None:
        return None
    refs: dict[str, str] = {}
    g_el = 0.0
    for el, x in comp.fractional_composition.items():
        ref = element_reference(str(el), T)
        if ref is None:
            return None
        refs[str(el)] = ref[0].name
        g_el += x * ref[1]
    sp, g = hit
    return FormationG(comp.reduced_formula, sp.name, sp.condensed, T, (g - g_el) / EV, refs)
