"""평형 트랙(T) 시험 — A3 기체 평형(Cantera·NASA), A4 수용액 평형(PHREEQC)."""

from __future__ import annotations

import re
from functools import cache
from importlib.metadata import version
from pathlib import Path

from pymatgen.core import Composition

from msl.assays.base import Context, Outcome, not_applicable, val
from msl.engines import nasa
from msl.schema.quantity import Dimension
from msl.schema.recipe import State
from msl.schema.result import Fidelity, Status

# ── A3 기체 평형·연소 ─────────────────────────────────────────────────────


def a3(ctx: Context) -> Outcome:
    comps = ctx.comps
    if any(c.formula is None for c in comps):
        return not_applicable("화학식이 없는 성분이 있음", "cantera", Fidelity.T)
    no_gas = [c.formula for c in comps if nasa.find_gas(c.formula) is None]
    if no_gas:
        return not_applicable(f"NASA 기체 DB 에 없는 성분: {', '.join(no_gas)}", "cantera", Fidelity.T)
    moles: dict[str, float] = {}
    for c in comps:
        a = c.amount
        if a is None:
            return not_applicable("기체 평형은 모든 성분의 양이 필요함", "cantera", Fidelity.T)
        if a.dimension is Dimension.AMOUNT:
            moles[c.formula] = a.to_si()
        elif a.dimension in (Dimension.MOLE_FRACTION, Dimension.VOLUME_FRACTION):  # 이상기체: vol% = mol%
            moles[c.formula] = a.to_si()
        else:
            return not_applicable(f"기체 성분 양은 mol 또는 mol%/vol% 로 적어야 함 ({c.formula}: {a})", "cantera", Fidelity.T)

    T0, P = ctx.T or 298.15, ctx.P
    tp = nasa.gas_equilibrium(moles, T0, P, "TP")
    hp = nasa.gas_equilibrium(moles, T0, P, "HP")
    xt, xh = dict(zip(tp.species_names, tp.X)), dict(zip(hp.species_names, hp.X))
    names = sorted(xt, key=lambda n: -max(xt[n], xh[n]))[:8]
    rows = [[n, round(float(xt[n]), 4), round(float(xh[n]), 4)] for n in names]
    return Outcome(
        status=Status.OK, fidelity=Fidelity.T, engine="cantera", engine_version=nasa.engine_version(),
        conditions_basis=f"시작 {T0:.0f} K, {P / 101325:.2f} atm, 이상기체",
        energy_reference=nasa.ENERGY_REFERENCE,
        values=[val("단열 평형 온도", round(hp.T, 0), "K"), val("단열 평형 온도", round(hp.T - 273.15, 0), "°C"),
                val("평형 계산에 쓴 화학종 수", float(len(tp.species_names)))],
        sources=[("nasa-cea-thermo", "Cantera nasa_gas.yaml")],
        summary=f"완전히 반응해 열을 잃지 않으면 {hp.T - 273.15:,.0f} °C 까지 오름 (단열 평형)",
        caveats=["평형 계산 — 점화 여부·반응 속도는 모름 (반응 메커니즘 필요)", "열손실·해리 이외의 비평형 효과 미반영"],
        data={"table": {"columns": ["화학종", f"몰분율 · 등온 {T0:.0f} K", f"몰분율 · 단열 {hp.T:.0f} K"], "rows": rows}},
    )


# ── A4 수용액 평형 ────────────────────────────────────────────────────────


@cache
def _db() -> tuple[dict[str, Composition], dict[str, dict[int, str]]]:
    """phreeqc.dat 의 광물(PHASES) 화학식과 원소별 마스터 화학종(산화수 → 이름)."""
    import phreeqpython

    txt = (Path(phreeqpython.__file__).parent / "database" / "phreeqc.dat").read_text(encoding="latin-1")
    phases_txt = txt.split("\nPHASES", 1)[1].split("\nEXCHANGE_MASTER_SPECIES", 1)[0]
    phases: dict[str, Composition] = {}
    for name, formula in re.findall(r"^([A-Z][A-Za-z0-9()_\-]*)\s*\n\s+([^=\n]+?)\s*=", phases_txt, re.M):
        if "(g)" in name:
            continue
        try:
            phases[name] = Composition(formula.split("+")[0].strip().replace(":", "."))
        except Exception:
            continue
    masters: dict[str, dict[int, str]] = {}
    msm = txt.split("SOLUTION_MASTER_SPECIES", 1)[1].split("SOLUTION_SPECIES", 1)[0]
    for line in msm.splitlines():
        parts = line.split()
        if len(parts) < 2 or parts[0].startswith("#") or parts[0] in ("H", "H(0)", "H(1)", "E", "O", "O(0)", "O(-2)"):
            continue
        name, species = parts[0], parts[1]
        m = re.match(r"^([A-Z][a-z]?)(?:\((-?\d+)\))?$", name)
        if not m:
            continue
        el, ox = m.group(1), m.group(2)
        oxi = int(ox) if ox is not None else _oxi_from_species(el, species)
        if oxi is not None:
            masters.setdefault(el, {})[oxi] = name
    return phases, masters


def _oxi_from_species(el: str, species: str) -> int | None:
    """마스터 화학종 식(예: CO3-2)에서 중심 원소의 산화수 (O=-2, H=+1 가정)."""
    m = re.match(r"^(.+?)([+-]\d*)?$", species)
    if not m:
        return None
    body, chg = m.group(1), m.group(2)
    charge = 0 if not chg else (int(chg[1:] or 1) * (1 if chg[0] == "+" else -1))
    try:
        comp = Composition(body)
    except Exception:
        return None
    n = comp[el]
    if n == 0:
        return None
    rest = sum(comp[e] * {"O": -2, "H": 1}.get(str(e), 0) for e in comp.elements if str(e) != el)
    return round((charge - rest) / n)


def a4(ctx: Context) -> Outcome:
    comps = ctx.comps
    water = [c for c in comps if c.is_water]
    if not water and not any(c.state is State.AQUEOUS for c in comps):
        return not_applicable("물이나 수용액 성분이 없음", "phreeqpython", Fidelity.T)
    phases, masters = _db()
    liters = 1.0
    if water and water[0].amount is not None and water[0].amount.dimension is Dimension.VOLUME:
        liters = water[0].amount.to_si() * 1000

    totals: dict[str, float] = {}  # 마스터 화학종 → mol/kgw
    minerals: dict[str, float | None] = {}  # PHREEQC 광물 이름 → 넣은 양(mol)
    for c in comps:
        if c.is_water:
            continue
        if c.formula is None:
            return not_applicable(f"화학식이 없는 성분: {c.name}", "phreeqpython", Fidelity.T)
        comp = Composition(c.formula)
        phase = next((n for n, pc in phases.items() if pc.almost_equals(comp)), None)
        if phase and c.state in (None, State.SOLID):
            minerals[phase] = c.moles()
            continue
        n = c.moles()
        if n is None:
            return not_applicable(f"{c.formula}: 양을 mol 또는 질량으로 적어야 함", "phreeqpython", Fidelity.T)
        guesses = comp.oxi_state_guesses(max_sites=-1)
        if not guesses:
            have = {el: sorted(masters.get(str(el), {})) for el in comp.elements if str(el) not in ("H", "O")}
            return not_applicable(
                f"{c.formula}: 산화수를 판정할 수 없어 PHREEQC 입력으로 바꿀 수 없음 "
                f"(phreeqc.dat 의 산화수: {', '.join(f'{k} {v}' for k, v in have.items())})", "phreeqpython", Fidelity.T)
        for el, amt in comp.get_el_amt_dict().items():
            if el in ("H", "O"):
                continue
            ox = round(guesses[0][el])
            name = masters.get(el, {}).get(ox)
            if name is None:
                return not_applicable(f"{c.formula}: phreeqc.dat 에 {el}({ox:+d}) 화학종이 없음", "phreeqpython", Fidelity.T)
            totals[name] = totals.get(name, 0.0) + n * amt / liters

    from phreeqpython import PhreeqPython

    T_c = (ctx.T or 298.15) - 273.15
    pp = PhreeqPython()
    sol = pp.add_solution({"temp": round(T_c, 2), "units": "mol/kgw", "pH": "7 charge",
                           **{k: f"{v:.8g}" for k, v in totals.items()}})
    eq_phases, targets = list(minerals), [0.0] * len(minerals)
    atm = (ctx.recipe.conditions.atmosphere or "").lower()
    if atm in ("air", "공기"):
        eq_phases.append("CO2(g)")
        targets.append(-3.4)
    if eq_phases:
        sol.equalize(eq_phases, targets)

    species = sorted(((k, v) for k, v in sol.species.items() if k != "H2O"), key=lambda p: -p[1])[:8]
    si = {k: v for k, v in sol.phases.items() if "(g)" not in k and k != "Fix_pH"}
    supersat = sorted(((k, v) for k, v in si.items() if v > 0.05), key=lambda p: -p[1])[:5]
    values = [val("pH", round(sol.pH, 2)), val("이온 세기", round(sol.I, 4), "mol/kgw")]
    caveats = ["phreeqc.dat 기본 DB — 고농도(대략 1 mol/kgw 이상)는 Pitzer DB 가 필요", "평형 계산 — 용해·침전 속도는 모름"]
    for m, given in minerals.items():
        el = next(str(e) for e in phases[m].elements if str(e) not in ("C", "O", "H"))
        dissolved = sol.total_element(el, "mol") if hasattr(sol, "total_element") else None
        if dissolved is not None:
            values.append(val(f"{m} 용해량", round(dissolved * 1000, 4), "mmol"))
            if given is not None and dissolved > given:
                caveats.append(f"{m}: 넣은 양({given:.4g} mol)보다 평형 용해량이 커서 과잉 고체를 가정한 결과와 다름")
    if supersat:
        values.append(val("과포화 광물 (SI>0, 침전 가능)", ", ".join(f"{k} {v:+.2f}" for k, v in supersat)))
    summary = f"pH {sol.pH:.2f}" + (f", {', '.join(minerals)} 와 평형" if minerals else "") + (
        ", 대기 CO₂ 평형" if "CO2(g)" in eq_phases else "")
    return Outcome(
        status=Status.OK, fidelity=Fidelity.T, engine="phreeqpython", engine_version=version("phreeqpython"),
        conditions_basis=f"{T_c:.0f} °C, 물 {liters:g} L(≈kgw)" + (", pCO₂ 10^-3.4 atm" if "CO2(g)" in eq_phases else ""),
        values=values, sources=[("phreeqc-db", "phreeqc.dat (phreeqpython 1.6.2 동봉)")], caveats=caveats,
        summary=summary,
        data={"table": {"columns": ["화학종", "몰랄농도 (mol/kgw)"], "rows": [[k, f"{v:.3e}"] for k, v in species]},
              "si": {k: round(v, 3) for k, v in sorted(si.items(), key=lambda p: -p[1])[:10]}},
    )
