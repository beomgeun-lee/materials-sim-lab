"""A6 합금 상평형 — pycalphad + 공개 TDB (평형 트랙 T, 조사 문서 04).

레시피 원소를 모두 덮는 TDB 를 tdb_systems 에서 고르고, 레시피 조성에서 온도를 훑어 안정 상과
상분율을 구한다. 고상선 = 액상이 처음 나타나는 온도, 액상선 = 고상이 모두 사라지는 온도 (평형 기준,
구간을 격자 세분해 0.5 K 안으로 좁힌다). 2원계 크기면 scheil 로 비평형 응고 종료 온도도 낸다.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from typing import Any

import numpy as np

import msl.db as mdb
from msl.assays.base import Context, Outcome, not_applicable, pending, val
from msl.connectors._tdb import read_database
from msl.registry.models import Partition
from msl.schema.quantity import Dimension
from msl.schema.refs import Namespace
from msl.schema.result import Fidelity, Status, ValueKind

ENGINE = "pycalphad"
ENERGY_REFERENCE = "CALPHAD (TDB 기준 상태: SER)"
QUALITY = {"evaluated": 2, "compiled": 1}
EPS = 1e-6
ROOM_T = 298.15
MAX_POINTS = 121  # 스윕 점이 이보다 많으면 솎는다 (다원계 TDB 는 점당 0.1 s 안팎)
DEFAULT_SWEEP = (300.0, 2000.0, 35)
AUX_GRID = np.arange(300.0, 4001.0, 50.0)  # 스윕 안에서 고상선·액상선을 못 잡았을 때
SCHEIL_MAX_PHASES = 15


@dataclass
class Tdb:
    source: str
    version: str
    file: str
    path: Path
    reference: str
    quality: str
    license: str
    restricted: bool


# ── TDB 선택 ──────────────────────────────────────────────────────────────


def _pairs(elements: list[str]) -> set[str]:
    return {f"{a}-{b}" for i, a in enumerate(elements) for b in elements[i + 1:]}


def rank(row: dict[str, Any], elements: list[str], partition: Partition) -> tuple:
    """클수록 좋다: 모든 쌍 평가됨 > 평가 DB > 계가 딱 맞음 > 대표 평가 > 재배포 가능 > 원소 수 적음."""
    have = set(filter(None, str(row.get("binaries") or "").split(",")))
    return (_pairs(elements) <= have, QUALITY.get(row.get("quality"), 0), row["n_elements"] == len(elements),
            bool(row.get("preferred")), partition is Partition.OPEN, -row["n_elements"])


def choose(elements: list[str], registry: Any) -> tuple[Tdb | None, str]:
    """레시피 원소를 덮고 모든 원소 쌍이 평가된 TDB 중 가장 나은 것. 없으면 (None, 이유)."""
    if "tdb_systems" not in mdb.tables():
        return None, "pending"
    rows = mdb.connect().sql("SELECT * FROM tdb_systems WHERE parse_ok").df().to_dict("records")
    need = set(elements)
    cover = [r for r in rows if need <= set(str(r["system"]).split("-"))]
    if not cover:
        return None, "no-cover"

    def part(r: dict[str, Any]) -> Partition:
        lic = registry.licenses.get(r["license"])
        return lic.partition if lic else Partition.RESTRICTED

    best = max(cover, key=lambda r: rank(r, elements, part(r)))
    if not rank(best, elements, part(best))[0]:
        return None, "no-pairs"
    path = mdb.DATA_DIR / "raw" / best["source"] / str(best["source_version"]).replace("/", "-") / best["file"]
    return Tdb(best["source"], str(best["source_version"]), best["file"], path, best["reference"] or "",
               best.get("quality") or "", best["license"], part(best) is Partition.RESTRICTED), ""


# ── 조성·온도 ─────────────────────────────────────────────────────────────


def mole_fractions(comps: list[Any]) -> dict[str, float]:
    """레시피 양(at%·mol%·mol·wt%·g) → 원소 몰분율."""
    from pymatgen.core import Element

    moles: dict[str, float] = {}
    for c in comps:
        a, el = c.amount, c.ref.key
        if a is None:
            raise ValueError("합금 조성(at%·wt%·mol·g)이 필요함")
        if a.dimension in (Dimension.MOLE_FRACTION, Dimension.AMOUNT):
            moles[el] = a.to_si()
        elif a.dimension is Dimension.MASS_FRACTION:
            moles[el] = a.to_si() / float(Element(el).atomic_mass)
        elif a.dimension is Dimension.MASS:  # kg → g → mol (mol 로 적은 성분과 섞여도 단위가 맞게)
            moles[el] = a.to_si() * 1000 / float(Element(el).atomic_mass)
        else:
            raise ValueError(f"합금 조성은 at%·wt%·mol·g 로 적어야 함 ({el}: {a})")
    total = sum(moles.values())
    return {el: n / total for el, n in moles.items()}


def temperatures(ctx: Context) -> np.ndarray:
    sw = ctx.recipe.sweep
    if sw is not None and sw.param == "T":
        T = np.linspace(sw.start.to_si(), sw.stop.to_si(), sw.steps)
        return T if len(T) <= MAX_POINTS else np.linspace(T[0], T[-1], MAX_POINTS)
    if ctx.T is not None:
        return np.array([ctx.T])
    return np.linspace(*DEFAULT_SWEEP)


# ── pycalphad 계산 ────────────────────────────────────────────────────────


class Calc:
    """TDB 하나·조성 하나에 대한 평형 계산 묶음. 모델은 한 번만 만든다."""

    def __init__(self, dbf: Any, elements: list[str], x: dict[str, float], P: float) -> None:
        from pycalphad import Model
        from pycalphad.core.utils import filter_phases

        self.dbf, self.P = dbf, P
        self.comps = [e.upper() for e in elements] + ["VA"]
        self.models: dict[str, Any] = {}
        self.dropped: list[str] = []
        for ph in filter_phases(dbf, self.comps):
            try:
                self.models[ph] = Model(dbf, self.comps, ph)
            except Exception as exc:  # pycalphad 가 못 만드는 상 모델(규칙-불규칙 비율 불일치 등)은 빼고 계산
                self.dropped.append(f"{ph} ({str(exc)[:80]})")
        self.phases = sorted(self.models)
        major = max(x, key=x.get)  # 주원소를 종속 변수로 두고 나머지 몰분율을 조건으로
        self.xconds = {e.upper(): v for e, v in x.items() if e != major}

    def _label(self, phase: str, y: np.ndarray) -> str:
        """규칙상 이름으로 나온 불규칙 고용체(FCC_L12 인데 부격자 점유가 같음)를 불규칙상 이름으로 바꾼다."""
        dis = self.dbf.phases[phase].model_hints.get("disordered_phase")
        if not dis or dis == phase or dis not in self.dbf.phases:
            return phase
        n = len(self.dbf.phases[phase].sublattices) - len(self.dbf.phases[dis].sublattices) + 1
        sub: dict[int, dict[str, float]] = {}
        for sf, value in zip(self.models[phase].site_fractions, y):
            if sf.sublattice_index < n:
                sub.setdefault(sf.sublattice_index, {})[sf.species.name] = float(value)
        first = sub.get(0, {})
        same = all(s.keys() == first.keys() and all(abs(s[k] - first[k]) < 1e-4 for k in s) for s in sub.values())
        return dis if same else phase

    def at(self, T: np.ndarray) -> list[dict[str, float] | None]:
        """온도마다 {상 이름: 몰분율}. 같은 상의 조성 세트는 이름에 #2 를 붙여 따로 센다. 수렴 실패는 None."""
        from pycalphad import equilibrium, variables as v

        T = np.atleast_1d(np.asarray(T, dtype=float))
        conds = {v.T: T, v.P: self.P, v.N: 1, **{v.X(e): xe for e, xe in self.xconds.items()}}
        with warnings.catch_warnings():  # pycalphad 0.11.2 의 NumPy 2.5 DeprecationWarning (점마다 수백 줄)
            warnings.simplefilter("ignore", DeprecationWarning)
            eq = equilibrium(self.dbf, self.comps, self.phases, conds, model=self.models)
        names = eq.Phase.values.reshape(len(T), -1)
        amounts = eq.NP.values.reshape(len(T), -1)
        ys = eq.Y.values.reshape(len(T), names.shape[1], -1)
        out: list[dict[str, float] | None] = []
        for i in range(len(T)):
            fr: dict[str, float] = {}
            for j, ph in enumerate(names[i]):
                amt = amounts[i, j]
                if not ph or np.isnan(amt) or amt < EPS:
                    continue
                label = self._label(str(ph), ys[i, j])
                key, k = label, 2
                while key in fr:
                    key, k = f"{label}#{k}", k + 1
                fr[key] = float(amt)
            out.append(fr or None)
        return out


def is_liquid(phase: str) -> bool:
    return "LIQ" in phase.upper()


def is_solid(phase: str) -> bool:
    return not is_liquid(phase) and not phase.upper().startswith("GAS")


def _liquid(fr: dict[str, float] | None) -> float:
    return sum(v for k, v in (fr or {}).items() if is_liquid(k))


def _solid(fr: dict[str, float] | None) -> float:
    return sum(v for k, v in (fr or {}).items() if is_solid(k))


def _first(T: np.ndarray, fracs: list[dict[str, float] | None], pred: Any, start: int = 0) -> int | None:
    for i in range(start, len(T)):
        if fracs[i] is not None and pred(fracs[i]):
            return i
    return None


def refine(calc: Calc, lo: float, hi: float, pred: Any, below: dict[str, float] | None = None,
           tol: float = 0.5) -> tuple[float, dict[str, float] | None]:
    """pred 가 lo 에서 거짓, hi 에서 참일 때 처음 참이 되는 온도를 tol 안으로 좁힌다. (온도, lo 쪽 상분율)."""
    for _ in range(6):
        if hi - lo <= tol:
            break
        grid = np.linspace(lo, hi, 11)
        fr = calc.at(grid)
        k = _first(grid, fr, pred)
        if k is None:  # 세분 격자에서 못 잡으면 바깥 구간 유지
            break
        if k > 0:
            lo, below = float(grid[k - 1]), fr[k - 1]
        hi = float(grid[k])
    return (lo + hi) / 2, below


def transitions(calc: Calc, T: np.ndarray, fracs: list[dict[str, float] | None]) -> dict[str, Any]:
    """고상선·액상선과 고상선 바로 아래 상분율. 스윕에서 못 잡으면 보조 격자(300–4000 K)로."""
    has_liq = lambda f: _liquid(f) > EPS  # noqa: E731
    no_sol = lambda f: _solid(f) < EPS  # noqa: E731
    out: dict[str, Any] = {"solidus": None, "liquidus": None, "below": None, "grid": "sweep"}
    i = _first(T, fracs, has_liq)
    j = _first(T, fracs, no_sol, i or 0) if i is not None else None
    if i is None or j is None or i == 0:
        T, fracs, out["grid"] = AUX_GRID, calc.at(AUX_GRID), "aux"
        i = _first(T, fracs, has_liq)
        j = _first(T, fracs, no_sol, i or 0) if i is not None else None
    if i is not None and i > 0:
        out["solidus"], out["below"] = refine(calc, T[i - 1], T[i], has_liq, fracs[i - 1])
    if j is not None and j > 0:
        out["liquidus"], _ = refine(calc, T[j - 1], T[j], no_sol)
    return out


def scheil(calc: Calc, start: float) -> dict[str, Any] | None:
    """Scheil–Gulliver 응고 (액상 완전 혼합·고상 확산 없음). 상이 많은 다원계 TDB 는 느려서 건너뛴다."""
    if len(calc.phases) > SCHEIL_MAX_PHASES or "LIQUID" not in calc.phases:
        return None
    from pycalphad import variables as v
    from scheil import simulate_scheil_solidification

    comp = {v.X(e): xe for e, xe in calc.xconds.items()}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        r = simulate_scheil_solidification(calc.dbf, calc.comps, calc.phases, comp, start, step_temperature=1.0)
    formed = sorted(p for p, amt in r.cum_phase_amounts.items() if amt and amt[-1] > 1e-3)
    return {"end": float(r.temperatures[-1]), "liquid_left": float(r.fraction_liquid[-1]), "phases": formed}


# ── 표시 ──────────────────────────────────────────────────────────────────


def _base(name: str) -> str:
    return name.split("#")[0]


def describe(fr: dict[str, float] | None) -> str:
    """'FCC_A1', 'FCC_A1 ×2 (상분리)', 'LIQUID 62% + FCC_A1 38%' 같은 한 줄 설명."""
    if not fr:
        return "수렴 실패"
    if len(fr) == 1:
        return next(iter(fr))
    if len({_base(k) for k in fr}) == 1:
        return f"{_base(next(iter(fr)))} ×{len(fr)} (상분리)"
    return " + ".join(f"{k} {v:.0%}" for k, v in sorted(fr.items(), key=lambda kv: -kv[1]))


def _series(T: np.ndarray, fracs: list[dict[str, float] | None]) -> dict[str, Any]:
    names = sorted({_base(k) for f in fracs if f for k in f}, key=lambda n: (is_liquid(n), n))
    phases = {n: [round(float(sum(v for k, v in f.items() if _base(k) == n)), 6) if f else None for f in fracs]
              for n in names}
    return {"T": [round(float(t), 2) for t in T], "unit": "K", "phases": phases}


def _table(series: dict[str, Any], fracs: list[dict[str, float] | None]) -> dict[str, Any]:
    names = list(series["phases"])
    rows = [[t, round(t - 273.15, 1), describe(f), *[None if series["phases"][n][i] is None
                                                       else round(series["phases"][n][i], 4) for n in names]]
            for i, (t, f) in enumerate(zip(series["T"], fracs))]
    return {"columns": ["T (K)", "T (°C)", "안정 상", *names], "rows": rows}


# ── 시험 ──────────────────────────────────────────────────────────────────


def a6(ctx: Context) -> Outcome:
    comps = ctx.comps
    if any(c.ref.namespace is not Namespace.ELEMENT for c in comps):
        return not_applicable("합금 상평형은 원소 성분만 받음", ENGINE, Fidelity.T)
    elements = sorted(c.ref.key for c in comps)
    system = "–".join(elements)
    try:
        x = mole_fractions(comps)
    except ValueError as exc:
        return not_applicable(str(exc), ENGINE, Fidelity.T)

    tdb, why = choose(elements, ctx.registry)
    if tdb is None:
        if why == "pending":
            return pending("tdb_systems 미적재 — msl db load sgte-binary-collection (조사 문서 04)", ENGINE, Fidelity.T)
        detail = "원소를 모두 덮는 TDB 없음" if why == "no-cover" else "덮는 TDB 는 있으나 평가되지 않은 원소 쌍이 있음"
        return not_applicable(f"공개 TDB 에 {system} 없음 ({detail}) — 조사 문서 04 참조", ENGINE, Fidelity.T)
    if not tdb.path.exists():
        return pending(f"TDB 원본 스냅샷 없음 — msl db fetch {tdb.source}", ENGINE, Fidelity.T)

    calc = Calc(read_database(str(tdb.path)), elements, x, ctx.P)
    T = temperatures(ctx)
    fracs = calc.at(T)
    tr = transitions(calc, T, fracs)
    room = calc.at(np.array([ROOM_T]))[0]
    sch = None
    if tr["liquidus"] is not None:
        try:
            sch = scheil(calc, tr["liquidus"] + 10.0)
        except Exception as exc:  # scheil 은 선택 — 실패해도 평형 결과는 낸다
            calc.dropped.append(f"scheil 실패: {type(exc).__name__}: {str(exc)[:80]}")

    comp_txt = ", ".join(f"{el} {x[el] * 100:.1f}" for el in elements) + " at%"
    values = []
    sol, liq = tr["solidus"], tr["liquidus"]
    if sol is not None:
        values += [val("고상선 (평형, 액상 첫 출현)", round(sol, 1), "K"), val("고상선 (평형, 액상 첫 출현)", round(sol - 273.15, 1), "°C")]
    else:
        values.append(val("고상선", "300–4000 K 안에서 못 찾음", kind=ValueKind.TEXT))
    if liq is not None:
        values += [val("액상선 (평형, 고상 소멸)", round(liq, 1), "K"), val("액상선 (평형, 고상 소멸)", round(liq - 273.15, 1), "°C")]
    else:
        values.append(val("액상선", "300–4000 K 안에서 못 찾음", kind=ValueKind.TEXT))
    if sol is not None and liq is not None:
        values.append(val("응고 구간 (액상선 − 고상선)", round(liq - sol, 1), "K"))
    main = None
    if tr["below"]:
        solids = {k: v for k, v in tr["below"].items() if is_solid(k)}
        main = max(solids, key=solids.get) if solids else None
        values.append(val("고상선 바로 아래 고체상", describe(tr["below"]), kind=ValueKind.TEXT))
    values.append(val(f"상온({ROOM_T:.0f} K) 평형 상", describe(room), kind=ValueKind.TEXT))
    if sch and sch["liquid_left"] < 0.01:
        values += [val("Scheil 응고 종료 온도 (비평형)", round(sch["end"], 1), "K"),
                   val("Scheil 응고 종료 온도 (비평형)", round(sch["end"] - 273.15, 1), "°C")]

    if sol is not None and liq is not None:
        summary = (f"{system} ({comp_txt}): {sol - 273.15:,.0f} °C 에서 녹기 시작해 {liq - 273.15:,.0f} °C 에서 모두 녹음"
                   + (f" — 굳으면 {main}" if main else ""))
    else:
        summary = f"{system} ({comp_txt}): 고상선·액상선을 300–4000 K 안에서 모두 잡지 못함 — 표의 상분율 참조"

    caveats = [
        "평형 계산 — 실제 주조 응고는 편석 때문에 고상선보다 낮은 온도까지 액상이 남을 수 있음. Scheil 값은 고상 확산이 "
        "없다고 본 반대쪽 극단 (전율고용계에서는 저융점 순원소 융점까지 내려감) — 실제는 둘 사이",
        "상온 상은 평형 기준 — 확산이 느려 실제 재료에는 고온 상이 그대로 남는 경우가 많음",
        f"정확도는 TDB 평가 품질을 따름: {tdb.reference}",
    ]
    if tdb.restricted:
        caveats.insert(0, f"연구용 TDB — 라이선스 {tdb.license} 소스({tdb.source})라 restricted 파티션, 배포 빌드에서는 쓰지 않음 (D2)")
    if tdb.quality == "compiled":
        caveats.append("문헌 파라미터를 모은 진행 중 DB — 평가 DB 보다 오차가 클 수 있음")
    series = _series(T, fracs)
    engine_version = version("pycalphad") + (f" + scheil {version('scheil')}" if sch else "")
    return Outcome(
        status=Status.OK, fidelity=Fidelity.T, engine=ENGINE, engine_version=engine_version,
        conditions_basis=(f"평형 계산, P = {ctx.P / 101325:.2f} atm, 조성 {comp_txt}, "
                          f"T {T[0]:.0f}–{T[-1]:.0f} K {len(T)}점" + (" (고상선·액상선은 300–4000 K 보조 격자)" if tr["grid"] == "aux" else "")),
        energy_reference=ENERGY_REFERENCE, values=values, summary=summary, caveats=caveats,
        warnings=[f"계산에서 뺀 상: {d}" for d in calc.dropped],
        sources=[(tdb.source, tdb.version)],
        data={"table": _table(series, fracs), "phase_fractions": series,
              "tdb": {"source": tdb.source, "file": tdb.file, "license": tdb.license, "reference": tdb.reference}},
    )
