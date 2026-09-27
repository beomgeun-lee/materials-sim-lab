"""NASA 열역학 데이터 평형 계산 — Cantera 번들 nasa_gas.yaml(748종)·nasa_condensed.yaml(382종).

평형 트랙(T)의 v0 엔진이다. 레시피 원소로만 이루어진 모든 화학종을 넣어 평형을 푼다.
DB 에 없는 화학종은 결과에 나올 수 없고, 응축상은 화학종마다 유효 온도 구간이 있다.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from importlib.metadata import version

import cantera as ct
from pymatgen.core import Composition

ENERGY_REFERENCE = "NASA 열화학 (298.15 K 원소 기준)"


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
