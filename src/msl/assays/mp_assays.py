"""Materials Project 기반 시험 — A1 조성 안정성, A2 고상·계면 반응, A7 물성 (L0)."""

from __future__ import annotations

from typing import Any

from pymatgen.analysis.interface_reactions import InterfacialReactivity
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import GibbsComputedStructureEntry

from msl.assays.base import Context, Outcome, not_applicable, pending, val
from msl.engines import mp, nasa
from msl.schema.recipe import Mode, State
from msl.schema.result import Fidelity, Status, ValueKind

MP_SRC = ("materials-project", mp.DB_VERSION)
NEAR_HULL = 0.05  # eV/atom — 준안정 후보로 보여 줄 범위


def _elements(ctx: Context) -> list[str] | None:
    if ctx.recipe.mode is Mode.SYSTEM:
        return [c.ref.key for c in ctx.comps]
    if any(c.formula is None for c in ctx.comps):
        return None
    return sorted({str(e) for c in ctx.comps for e in Composition(c.formula).elements})


def _phase_diagram(ctx: Context, elements: list[str]) -> tuple[list, PhaseDiagram]:
    key = "pd:" + "-".join(elements)
    if key not in ctx.shared:
        entries = mp.entries_in_chemsys(elements)
        ctx.shared[key] = (entries, PhaseDiagram(entries))
    return ctx.shared[key]


def _best(pd: PhaseDiagram, entries: list, formula: str):
    """같은 조성 엔트리 중 hull 에 가장 가까운 것 (다형 중 최안정)."""
    target = Composition(formula).reduced_formula
    same = [e for e in entries if e.composition.reduced_formula == target]
    return min(same, key=pd.get_e_above_hull) if same else None


# ── A1 ────────────────────────────────────────────────────────────────────


def a1(ctx: Context) -> Outcome:
    elements = _elements(ctx)
    if elements is None:
        return not_applicable("화학식이 없는 성분(사용자 정의 소재 등)이 있어 상태도를 만들 수 없음", "pymatgen")
    if len(elements) > 6:
        return not_applicable(f"원소 {len(elements)}개 — v0 은 6원소 이하 화학계만", "pymatgen")
    try:
        entries, pd = _phase_diagram(ctx, elements)
    except mp.MPUnavailable as exc:
        return pending(f"Materials Project 연결 불가: {exc}", "pymatgen")

    rows: list[list[Any]] = []
    stable = sorted(pd.stable_entries, key=lambda e: pd.get_form_energy_per_atom(e))
    for e in stable:
        if len(e.composition.elements) == 1:
            continue
        rows.append([e.composition.reduced_formula, str(e.entry_id).split("-GGA")[0],
                     round(pd.get_form_energy_per_atom(e), 3), 0.0, "안정"])
    n_stable = len(rows)
    values = [val("안정 화합물 수", float(n_stable))]

    if ctx.recipe.mode is Mode.SYSTEM:
        best: dict[str, tuple[float, Any]] = {}
        for e in entries:
            eh = pd.get_e_above_hull(e)
            f = e.composition.reduced_formula
            if 1e-6 < eh <= NEAR_HULL and len(e.composition.elements) > 1 and f not in {r[0] for r in rows}:
                if f not in best or eh < best[f][0]:
                    best[f] = (eh, e)
        near = sorted(best.values(), key=lambda p: p[0])[:8]
        for eh, e in near:
            rows.append([e.composition.reduced_formula, str(e.entry_id).split("-GGA")[0],
                         round(pd.get_form_energy_per_atom(e), 3), round(eh * 1000, 1), "준안정 → L2 재확인 대상"])
        values.append(val("hull 근처 준안정 후보 (≤50 meV/atom)", float(len(near))))
        summary = f"{'-'.join(elements)} 계에서 안정 화합물 {n_stable}개, L2 재확인 후보 {len(near)}개"
    else:
        parts = []
        for c in ctx.comps:
            e = _best(pd, entries, c.formula)
            if e is None:
                parts.append(f"{c.formula}: MP 에 없음")
                continue
            eh = pd.get_e_above_hull(e) * 1000
            values.append(val(f"{c.formula} E_hull", round(eh, 1), "meV/atom", kind=ValueKind.ENERGY))
            parts.append(f"{c.formula} {eh:.0f} meV/atom")
        summary = "반응물 안정성 — " + ", ".join(parts)

    return Outcome(
        status=Status.OK, fidelity=Fidelity.L0, engine="pymatgen", engine_version=_pmg(),
        conditions_basis="0 K DFT", energy_reference=mp.ENERGY_REFERENCE, values=values,
        sources=[MP_SRC], summary=summary,
        caveats=["0 K 계산 — 고온에서만 안정한 다형은 hull 위로 나올 수 있음",
                 "구조 탐색은 MP 에 이미 있는 구조로 한정 — 새 구조는 L2(uMLIP) 단계에서"],
        data={"table": {"columns": ["상", "MP ID", "형성에너지 (eV/atom)", "E_hull (meV/atom)", "판정"], "rows": rows}},
    )


# ── A2 ────────────────────────────────────────────────────────────────────


def _kinks(ir: InterfacialReactivity) -> list[dict[str, Any]]:
    return [{"x": round(float(x), 4), "e": round(float(e) * 1000, 1), "rxn": str(rxn)} for _, x, e, rxn, _ in ir.get_kinks()]


def _interp(kinks: list[dict[str, Any]], x: float) -> float:
    for a, b in zip(kinks, kinks[1:]):
        if a["x"] <= x <= b["x"]:
            t = 0.0 if b["x"] == a["x"] else (x - a["x"]) / (b["x"] - a["x"])
            return a["e"] + t * (b["e"] - a["e"])
    return kinks[-1]["e"]


def _recipe_x(ctx: Context) -> float | None:
    """레시피 양으로 본 혼합비 x (InterfacialReactivity 와 같은 원자 기준)."""
    c1, c2 = ctx.comps
    n1, n2 = c1.moles(), c2.moles()
    if n1 is None or n2 is None:
        return None
    a1_, a2_ = n1 * Composition(c1.formula).num_atoms, n2 * Composition(c2.formula).num_atoms
    return a1_ / (a1_ + a2_)


def _cross_check(ctx: Context, gibbs_pd: PhaseDiagram | None, gibbs_entries: list | None) -> dict[str, Any]:
    """평형 트랙 교차검증 — NASA 실험 열화학으로 각 반응물이 T 에서 분해되는지 본다."""
    T = ctx.T
    out: dict[str, Any] = {"T": T, "checks": [], "mismatch": False}
    if T is None:
        return out
    for c in ctx.comps:
        item: dict[str, Any] = {"formula": c.formula}
        T_used = T
        try:
            eq = nasa.multiphase_equilibrium({c.formula: 1.0}, T, ctx.P)
        except nasa.OutOfRange as exc:
            if exc.tmax < T:
                T_used = exc.tmax - 1
                eq = nasa.multiphase_equilibrium({c.formula: 1.0}, T_used, ctx.P)
                item["note"] = f"NASA 데이터가 {exc.tmax:.0f} K 까지라 {T_used:.0f} K 에서 계산"
            else:
                item["note"] = str(exc)
                out["checks"].append(item)
                continue
        except KeyError:
            item["note"] = "NASA DB 에 없음"
            out["checks"].append(item)
            continue
        same = {s.name for s in nasa.find_condensed(c.formula)}
        left = sum(v for n, v in eq.condensed.items() if n in same)
        item.update({"T_used": T_used, "decomposed": left < 0.5,
                     "products": {**{k: round(v, 3) for k, v in eq.condensed.items()},
                                  **({"기체": round(eq.gas_moles, 3)} if eq.gas_moles > 1e-6 else {})}})
        if gibbs_pd is not None and gibbs_entries is not None:
            e = _best(gibbs_pd, gibbs_entries, c.formula)
            if e is not None:
                eh = gibbs_pd.get_e_above_hull(e)
                item["mp_gibbs_e_hull_meV"] = round(eh * 1000, 1)
                if item["decomposed"] and eh < 0.005:
                    item["mismatch"] = True
                    out["mismatch"] = True
        out["checks"].append(item)
    return out


def a2(ctx: Context) -> Outcome:
    comps = ctx.comps
    if ctx.recipe.mode is not Mode.MIXTURE or len(comps) != 2:
        return not_applicable("고상·계면 반응은 두 성분 혼합에만 적용", "pymatgen")
    if any(c.formula is None for c in comps) or any(c.state not in (None, State.SOLID) for c in comps):
        return not_applicable("두 성분 모두 화학식이 있는 고체여야 함", "pymatgen")
    elements = _elements(ctx)
    try:
        entries, pd = _phase_diagram(ctx, elements)
    except mp.MPUnavailable as exc:
        return pending(f"Materials Project 연결 불가: {exc}", "pymatgen")

    c1, c2 = (Composition(c.formula) for c in comps)
    use_hull = any(_best(pd, entries, c.formula) is None for c in comps)
    curves = [{"label": "0 K (DFT)", "kinks": _kinks(InterfacialReactivity(c1, c2, pd, norm=True, use_hull_energy=use_hull))}]
    gibbs_pd = gibbs_entries = None
    T = ctx.T
    if T and T >= 300:
        Tg = int(min(2000, max(300, round(T / 100) * 100)))
        gibbs_entries = GibbsComputedStructureEntry.from_entries(entries, temp=Tg)
        gibbs_pd = PhaseDiagram(gibbs_entries)
        curves.append({"label": f"{Tg} K (Gibbs 근사)", "kinks": _kinks(
            InterfacialReactivity(c1, c2, gibbs_pd, norm=True, use_hull_energy=use_hull))})

    x = _recipe_x(ctx)
    values = []
    for cv in curves:
        low = min(cv["kinks"], key=lambda k: k["e"])
        values.append(val(f"최저 반응에너지 · {cv['label']}", low["e"], "meV/atom", kind=ValueKind.ENERGY))
        if x is not None:
            values.append(val(f"레시피 혼합비에서 · {cv['label']}", round(_interp(cv["kinks"], x), 1),
                              "meV/atom", kind=ValueKind.ENERGY))
    main = curves[-1]
    low = min(main["kinks"], key=lambda k: k["e"])
    reacts = low["e"] < -1.0
    summary = (f"{main['label']}: {low['rxn']} — 반응에너지 {low['e']:.0f} meV/atom" if reacts
               else f"{main['label']}: 두 성분 사이에 열역학적으로 유리한 반응이 없음")

    check = _cross_check(ctx, gibbs_pd, gibbs_entries)
    caveats = ["열역학 판정만 — 반응 속도·핵생성은 반영하지 않음",
               "고체 유한온도는 Bartel 기술자 근사, 기체(CO₂ 등)는 실험 G(T)"]
    warnings = []
    status = Status.OK
    if check["mismatch"]:
        status = Status.WARNING
        bad = [c for c in check["checks"] if c.get("mismatch")]
        for c in bad:
            warnings.append(
                f"교차검증 불일치: NASA 실험 열화학은 {c['formula']} 가 {c['T_used']:.0f} K 에서 분해된다고 하지만 "
                f"MP+Gibbs 근사는 안정(E_hull {c['mp_gibbs_e_hull_meV']} meV/atom)으로 봄 → 이 온도의 L0 반응 판정은 신뢰 낮음")
    stable_compounds = {e.composition.reduced_formula for e in pd.stable_entries if len(e.composition.elements) > 1}
    missing = sorted(f for f in stable_compounds if not (nasa.find_condensed(f) or nasa.find_gas(f)))
    if missing:
        caveats.append("NASA 실험 DB 에 없는 안정상(평형 트랙으로 확인 불가): " + ", ".join(missing[:8]))

    return Outcome(
        status=status, fidelity=Fidelity.L0, engine="pymatgen", engine_version=_pmg(),
        conditions_basis=f"0 K DFT + {curves[-1]['label']}" if len(curves) > 1 else "0 K DFT",
        energy_reference=mp.ENERGY_REFERENCE, values=values,
        sources=[MP_SRC, ("nasa-cea-thermo", "Cantera nasa_*.yaml")],
        caveats=caveats, warnings=warnings, summary=summary,
        data={"curves": curves, "recipe_x": x, "labels": [comps[0].formula, comps[1].formula], "crosscheck": check},
    )


# ── A7 ────────────────────────────────────────────────────────────────────

FIELDS = ["formula_pretty", "density", "band_gap", "formation_energy_per_atom", "energy_above_hull"]


def a7(ctx: Context) -> Outcome:
    if all(c.props for c in ctx.comps):
        rows = [[c.name, c.props.get("density"), c.props.get("youngs_modulus"), c.props.get("poisson_ratio")] for c in ctx.comps]
        return Outcome(status=Status.OK, fidelity=Fidelity.L0, engine="msl-mixing-rules", engine_version="v0",
                       conditions_basis="사용자 입력 물성", summary="사용자 정의 소재 물성 (kb/materials.yaml)",
                       data={"table": {"columns": ["소재", "밀도 (g/cm³)", "영률 (GPa)", "푸아송비"], "rows": rows}},
                       values=[val("소재 수", float(len(rows)))])
    elements = _elements(ctx)
    if elements is None:
        return not_applicable("화학식이 없는 성분이 있음", "mp-api")
    try:
        entries, pd = _phase_diagram(ctx, elements)
    except mp.MPUnavailable as exc:
        return pending(f"Materials Project 연결 불가: {exc}", "mp-api")

    if ctx.recipe.mode is Mode.SYSTEM:
        targets = [e for e in pd.stable_entries if len(e.composition.elements) > 1]
    else:
        targets = [e for e in (_best(pd, entries, c.formula) for c in ctx.comps) if e is not None]
    ids = [str(e.entry_id).split("-GGA")[0] for e in targets]
    docs = mp.summaries(ids, FIELDS)
    rows = []
    for mid in ids:
        d = docs.get(mid, {})
        rows.append([d.get("formula_pretty"), mid, _r(d.get("density"), 3), _r(d.get("band_gap"), 2),
                     _r(d.get("formation_energy_per_atom"), 3)])
    values = [val("조회한 상 수", float(len(rows)))]

    comps = ctx.comps
    if ctx.recipe.mode is Mode.MIXTURE and all(c.ref.namespace.value == "element" for c in comps):
        fr = [c.amount.to_si() if c.amount and c.amount.dimension.value == "mole_fraction" else None for c in comps]
        dens = [docs.get(i, {}).get("density") for i in ids]
        if None not in fr and None not in dens and len(dens) == len(comps):
            masses = [Composition(c.formula).weight for c in comps]
            rho = sum(f * m for f, m in zip(fr, masses)) / sum(f * m / d for f, m, d in zip(fr, masses, dens))
            values.append(val("합금 밀도 추정 (몰부피 가산)", round(rho, 3), "g/cm³"))

    return Outcome(
        status=Status.OK, fidelity=Fidelity.L0, engine="mp-api", engine_version=mp.engine_version(),
        conditions_basis="0 K DFT (PBE/PBE+U)", energy_reference=mp.ENERGY_REFERENCE, values=values,
        sources=[MP_SRC], summary=f"MP 요약 물성 {len(rows)}개 상 조회",
        caveats=["PBE 밴드갭은 실측보다 약 40% 작게 나오는 경향",
                 "MP 요약 구조에는 r²SCAN 재계산분이 섞여 있어 밀도가 실측과 몇 % 다를 수 있음 (예: Cu 9.22 vs 실측 8.96 g/cm³)"],
        data={"table": {"columns": ["상", "MP ID", "밀도 (g/cm³)", "밴드갭 (eV)", "형성에너지 (eV/atom)"], "rows": rows}},
    )


def _r(v: Any, n: int) -> Any:
    return round(float(v), n) if isinstance(v, (int, float)) else v


def _pmg() -> str:
    from importlib.metadata import version

    return version("pymatgen")
