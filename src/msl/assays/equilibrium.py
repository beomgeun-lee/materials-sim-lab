"""평형 트랙(T) 시험 — A3 기체 평형(Cantera·NASA), A4 수용액 평형(PHREEQC)."""

from __future__ import annotations

import math
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


LLNL = Path(__file__).resolve().parents[1] / "data" / "phreeqc" / "llnl.dat"  # USGS PHREEQC 배포본 동봉 (LLNL thermo.com.V8.R6)


def db_path(name: str = "phreeqc.dat") -> Path:
    """A4 가 쓰는 DB 파일 — phreeqc.dat 은 phreeqpython 동봉본, llnl.dat 은 패키지에 넣은 PHREEQC 배포본 파일 (D27)."""
    if name == "llnl.dat":
        return LLNL
    import phreeqpython

    return Path(phreeqpython.__file__).parent / "database" / name


@cache
def _db(name: str = "phreeqc.dat") -> tuple[dict[str, Composition], dict[str, dict[int, str]]]:
    """DB 의 광물(PHASES) 화학식과 원소별 마스터 화학종(산화수 → 이름)."""
    txt = db_path(name).read_text(encoding="latin-1")
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


class _Unsupported:
    """DB 가 모르는 산화 상태 — 다른 DB 로 다시 해 볼 수 있다."""

    def __init__(self, reason: str):
        self.reason = reason


def _oxi_guess(comp: Composition, masters: dict[str, dict[int, str]]) -> dict[str, float] | None:
    """산화수 배정 — pymatgen 기본 추정이 DB 에 있는 상태면 그대로, 아니면 DB 가 아는 산화 상태 안에서 다시 찾는다
    (NaClO: pymatgen 은 Cl⁺¹ 을 흔한 상태로 보지 않아 추정이 비지만, llnl.dat 에는 Cl(1)=ClO⁻ 가 있다)."""
    els = [str(e) for e in comp.elements if str(e) not in ("H", "O")]
    guesses = comp.oxi_state_guesses(max_sites=-1)
    for g in guesses:
        if all(round(g[el]) in masters.get(el, {}) for el in els):
            return g
    override = {el: sorted(masters[el]) for el in els if masters.get(el)}
    if len(override) < len(els):
        return None
    guesses = comp.oxi_state_guesses(oxi_states_override=override, max_sites=-1)
    return guesses[0] if guesses else None


def _totals(comps, liters: float, db: str):
    """성분 → (마스터 화학종별 mol/kgw, 평형 광물 {이름: 넣은 mol}). 못 하면 Outcome(해당 없음) 또는 _Unsupported."""
    phases, masters = _db(db)
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
        if n is None and c.amount is not None and c.amount.dimension is Dimension.CONCENTRATION:
            n = c.amount.to_si() / 1000 * liters  # mol/m³ → mol/L × 물 부피 (질량수지와 같은 환산)
        elif n is None and c.amount is not None and c.amount.dimension is Dimension.MOLALITY:
            n = c.amount.to_si() * liters  # mol/kg × kgw(≈ L)
        if n is None:
            return not_applicable(f"{c.formula}: 양을 mol·질량·몰농도로 적어야 함", "phreeqpython", Fidelity.T)
        guess = _oxi_guess(comp, masters)
        if guess is None:
            have = {el: sorted(masters.get(str(el), {})) for el in comp.elements if str(el) not in ("H", "O")}
            states = "; ".join(f"{k} {', '.join(f'{x:+d}' for x in v) or '없음'}" for k, v in have.items())
            return _Unsupported(
                f"{c.formula}: 수용액 계산 DB 가 다루는 이온으로 바꿀 수 없어 수용액 평형은 계산하지 않음 "
                f"— DB 가 아는 산화 상태는 {states} 뿐")
        for el, amt in comp.get_el_amt_dict().items():
            if el in ("H", "O"):
                continue
            name = masters[el][round(guess[el])]
            totals[name] = totals.get(name, 0.0) + n * amt / liters
    return totals, minerals


def a4(ctx: Context) -> Outcome:
    comps = ctx.comps
    water = [c for c in comps if c.is_water]
    if not water and not any(c.state is State.AQUEOUS for c in comps):
        return not_applicable("물이나 수용액 성분이 없음", "phreeqpython", Fidelity.T)
    liters = 1.0
    if water and water[0].amount is not None and water[0].amount.dimension is Dimension.VOLUME:
        liters = water[0].amount.to_si() * 1000

    # phreeqc.dat 으로 먼저, 그 DB 가 모르는 산화 상태(락스의 Cl⁺¹ 등)가 있으면 llnl.dat 으로 전체를 계산한다 (D27)
    db = "phreeqc.dat"
    got = _totals(comps, liters, db)
    if isinstance(got, _Unsupported) and LLNL.exists():
        alt = _totals(comps, liters, "llnl.dat")
        if not isinstance(alt, (_Unsupported, Outcome)):
            db, got = "llnl.dat", alt
    if isinstance(got, _Unsupported):
        return not_applicable(got.reason, "phreeqpython", Fidelity.T)
    if isinstance(got, Outcome):
        return got
    totals, minerals = got
    phases, _ = _db(db)

    from phreeqpython import PhreeqPython

    T_c = (ctx.T or 298.15) - 273.15
    P_atm = ctx.P / 101_325
    pp = PhreeqPython(database=db, database_directory=db_path(db).parent)  # 기본값은 vitens.dat(Stimela 파생) — 출처 표기와 맞춘다
    sol = pp.add_solution({"temp": round(T_c, 2), "pressure": f"{P_atm:.6g}", "units": "mol/kgw", "pH": "7 charge",
                           **{k: f"{v:.8g}" for k, v in totals.items()}})
    eq_phases, targets = list(minerals), [0.0] * len(minerals)
    atm = (ctx.recipe.conditions.atmosphere or "").lower()
    if atm in ("air", "공기"):
        eq_phases.append("CO2(g)")
        targets.append(-3.4)
    if eq_phases:
        sol.equalize(eq_phases, targets)

    ctx.shared["pH"] = float(sol.pH)  # A5 Pourbaix 가 쓴다
    species = sorted(((k, v) for k, v in sol.species.items() if k != "H2O"), key=lambda p: -p[1])[:8]
    si = {k: v for k, v in sol.phases.items() if "(g)" not in k and k != "Fix_pH"}
    supersat = sorted(((k, v) for k, v in si.items() if v > 0.05), key=lambda p: -p[1])[:5]
    values = [val("pH", round(sol.pH, 2)), val("이온 세기", round(sol.I, 4), "mol/kgw")]
    caveats = ["phreeqc.dat 기본 DB — 고농도(대략 1 mol/kgw 이상)는 Pitzer DB 가 필요" if db == "phreeqc.dat" else
               "llnl.dat 사용 — phreeqc.dat 에 없는 산화 상태가 있어 전체를 LLNL DB(B-dot 활동도, 0–300 °C)로 계산. 대략 1 mol/kgw 이상은 부정확",
               "평형 계산 — 용해·침전 속도는 모름"]
    if db == "llnl.dat":
        caveats.append("산화 상태는 섞은 그대로 둔다 — 차아염소산이 염화물·산소로 분해되는 느린 반응은 반영하지 않음")
    cl2 = _chlorine(sol, totals, db, T_c + 273.15)
    if cl2 is not None:
        values.append(val("염소 기체 평형 분압 (Cl₂)", f"{10 ** cl2:.2e}", "atm"))  # 10⁻⁸ atm 도 0 으로 반올림되지 않게 지수 표기
        values.append(val("염소 기체 평형 분압 log₁₀", round(cl2, 2), "log atm"))
        if cl2 > 0:
            caveats.append("염소 기체가 빠져나가며 용액 농도가 줄어드는 과정은 반영하지 않음 (평형 분압은 섞은 직후 기준)")
    for m, given in minerals.items():
        el = next(str(e) for e in phases[m].elements if str(e) not in ("C", "O", "H"))
        dissolved = sol.total_element(el, "mol") if hasattr(sol, "total_element") else None
        if dissolved is not None:
            values.append(val(f"{m} 용해량", round(dissolved * 1000, 4), "mmol"))
            if given is not None and dissolved > given:
                caveats.append(f"{m}: 넣은 양({given:.4g} mol)보다 평형 용해량이 커서 과잉 고체를 가정한 결과와 다름")
    if supersat:
        values.append(val("과포화 광물 (SI>0, 침전 가능)", ", ".join(f"{k} {v:+.2f}" for k, v in supersat)))
    if T_c > 100 or P_atm > 1.5:
        caveats.append("phreeqc.dat 의 온도 의존성은 대략 0–100 °C 에서 검증됨 — 고온·고압은 경향만 볼 것 (D22)")
    # Reaktoro 는 산화 상태까지 모두 평형으로 풀어(차아염소산 → 염화물+산소) 섞은 그대로 둔 llnl 계산과 비교할 수 없다
    cross = (_reaktoro_cross(totals, list(minerals), liters, ctx.T or 298.15, P_atm * 1.01325)
             if "CO2(g)" not in eq_phases and db == "phreeqc.dat" else None)
    if cross:
        if "pH" in cross:
            values.append(val("pH (Reaktoro 교차검증)", round(cross["pH"], 2)))
            if abs(cross["pH"] - sol.pH) > PH_CROSS_TOL:
                caveats.append(f"PHREEQC 와 Reaktoro 의 pH 가 {abs(cross['pH'] - sol.pH):.2f} 차이 — 두 엔진 해석이 갈리는 조건 (결과 신중히)")
        else:
            caveats.append(f"Reaktoro 교차검증 실패: {cross['error']}")
    summary = f"pH {sol.pH:.2f}" + (f", {', '.join(minerals)} 와 평형" if minerals else "") + (
        ", 대기 CO₂ 평형" if "CO2(g)" in eq_phases else "")
    if cl2 is not None and cl2 > 0:
        summary = f"염소(Cl₂) 기체 발생 — 평형 분압 {10 ** cl2:.2g} atm (대기압 초과) · " + summary
    elif cl2 is not None and cl2 > -6:
        summary += f" · 염소(Cl₂) 기체 조금 발생 — 평형 분압 {10 ** cl2:.2g} atm (약 {10 ** cl2 * 1e6:.3g} ppm)"
    return Outcome(
        status=Status.OK, fidelity=Fidelity.T, engine="phreeqpython", engine_version=version("phreeqpython"),
        conditions_basis=f"{T_c:.0f} °C, 물 {liters:g} L(≈kgw)" + (", pCO₂ 10^-3.4 atm" if "CO2(g)" in eq_phases else ""),
        values=values, caveats=caveats,
        sources=[("phreeqc-db", "phreeqc.dat (phreeqpython 1.6.2 동봉)" if db == "phreeqc.dat" else "llnl.dat (PHREEQC 동봉, LLNL thermo.com.V8.R6)")],
        summary=summary,
        data={"table": {"columns": ["화학종", "몰랄농도 (mol/kgw)"], "rows": [[k, f"{v:.3e}"] for k, v in species]},
              "si": {k: round(v, 3) for k, v in sorted(si.items(), key=lambda p: -p[1])[:10]},
              **({"reaktoro": cross} if cross else {})},
    )


PH_CROSS_TOL = 0.1


@cache
def _thermo(db: str, block: str, name: str) -> tuple[float, float | None, list[float] | None]:
    """DB 의 반응 상수 (log_k, delta_H kJ/mol, analytic 계수). block: 'species' = 반응 결과가 name 인 화학종, 'phase' = 상 이름."""
    txt = db_path(db).read_text(encoding="latin-1")
    if block == "species":
        body = txt.split("\nSOLUTION_SPECIES", 1)[1].split("\nPHASES", 1)[0]
        m = re.search(rf"^[^#\n]*= *{re.escape(name)}\s*\n(.*?)(?=^[^\s#]|^ *[^#\s-][^\n]*=)", body, re.M | re.S)
    else:
        body = txt.split("\nPHASES", 1)[1]
        m = re.search(rf"^{re.escape(name)}\s*\n\s+[^\n]*=[^\n]*\n(.*?)(?=^[^\s#])", body, re.M | re.S)
    if m is None:
        raise KeyError(f"{db}: {name} 없음")
    blk = m.group(1)
    lk = float(re.search(r"log_k\s+(\S+)", blk).group(1))
    dh = re.search(r"-delta_H\s+(\S+)\s+kJ", blk)
    an = re.search(r"^\s*-analytic\s+([^\n#]+)", blk, re.M)
    return lk, float(dh.group(1)) if dh else None, [float(x) for x in an.group(1).split()] if an else None


def _logk(db: str, block: str, name: str, T: float) -> float:
    """온도 T(K)의 log K — analytic 식이 있으면 그것, 없으면 van 't Hoff (PHREEQC 와 같은 우선순위)."""
    lk, dh, an = _thermo(db, block, name)
    if an:
        a = an + [0.0] * (6 - len(an))
        return a[0] + a[1] * T + a[2] / T + a[3] * math.log10(T) + a[4] / T ** 2 + a[5] * T ** 2
    if dh is not None:
        return lk - dh / (8.314462618e-3 * math.log(10)) * (1 / T - 1 / 298.15)
    return lk


def _chlorine(sol, totals: dict[str, float], db: str, T: float) -> float | None:
    """차아염소산(Cl⁺¹)과 염화 이온이 함께 있으면 Cl₂(g) 평형 분압의 log (atm), 아니면 None (D27).

    llnl.dat: Cl⁻ + ½O₂ = ClO⁻ (k1), Cl₂(g) + H₂O = ½O₂ + 2Cl⁻ + 2H⁺ (K2) → O₂ 를 소거하면
    log P(Cl₂) = log a(ClO⁻) + log a(Cl⁻) + 2 log a(H⁺) − log a(H₂O) − k1 − K2.
    25 °C 에서 HOCl + H⁺ + Cl⁻ ⇌ Cl₂(g) + H₂O 의 K = 10^4.53 (문헌 약 10^4.55). 물 활동도는 1 로 둔다 (묽은 용액)."""
    if db != "llnl.dat" or "Cl(1)" not in totals or not any(k in totals for k in ("Cl", "Cl(-1)")):
        return None
    try:
        la = {sp: math.log10(sol.activity(sp, "mol")) for sp in ("ClO-", "Cl-", "H+")}
        return la["ClO-"] + la["Cl-"] + 2 * la["H+"] - _logk(db, "species", "ClO-", T) - _logk(db, "phase", "Cl2(g)", T)
    except (ValueError, KeyError, AttributeError):
        return None


def _reaktoro_cross(totals: dict[str, float], minerals: list[str], liters: float, T_K: float, P_bar: float) -> dict | None:
    """같은 phreeqc.dat 로 Reaktoro 가 푼 pH (D22). Reaktoro 환경이 없으면 None — 교차검증을 건너뛴다."""
    from msl.engines import reaktoro_bridge as rb

    if not rb.available():
        return None
    import phreeqpython

    db = str(Path(phreeqpython.__file__).parent / "database" / "phreeqc.dat")
    try:
        out = rb.solve(rb.from_a4(totals, minerals, db, liters, T_K, P_bar))
    except Exception as e:  # 교차검증 실패는 A4 결과를 막지 않는다
        return {"error": str(e)[:200]}
    if not out.get("ok"):
        return {"error": "Reaktoro 평형 계산이 수렴하지 않음"}
    return {"pH": out["pH"], "I": out["I"], "engine": f"reaktoro {rb.version()}", "si": out["si"]}
