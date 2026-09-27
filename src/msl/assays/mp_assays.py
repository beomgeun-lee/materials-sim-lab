"""Materials Project 기반 시험 — A1 조성 안정성, A2 고상·계면 반응, A7 물성 (L0)."""

from __future__ import annotations

from typing import Any

from pymatgen.analysis.interface_reactions import InterfacialReactivity
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import GibbsComputedStructureEntry

from msl.assays.base import Context, Outcome, not_applicable, pending, val
from msl.engines import mp, nasa
from msl.engines.mp import mp_int
from msl.engines import thermo_hybrid as th
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


def _polymorph(entries: list, mp_id: str | None):
    """광물이 가리키는 MP 다형 엔트리 (GGA/GGA+U 엔트리 중 같은 material_id)."""
    target = mp_int(mp_id)
    if target is None:
        return None
    return next((e for e in entries if mp_int(e.entry_id) == target), None)


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
        near = sorted(best.values(), key=lambda p: p[0])
        for eh, e in near[:8]:  # 표에는 hull 에 가까운 8개만
            rows.append([e.composition.reduced_formula, str(e.entry_id).split("-GGA")[0],
                         round(pd.get_form_energy_per_atom(e), 3), round(eh * 1000, 1), "준안정 → L2 재확인 대상"])
        values.append(val("hull 근처 준안정 후보 (≤50 meV/atom)", float(len(near))))
        summary = f"{'-'.join(elements)} 계에서 안정 화합물 {n_stable}개, L2 재확인 후보 {len(near)}개" + (
            " (표에는 가까운 8개)" if len(near) > 8 else "")
    else:
        parts, poly_notes = [], []
        for c in ctx.comps:
            poly = _polymorph(entries, getattr(c, "mp_id", None))
            e = poly or _best(pd, entries, c.formula)
            if e is None:
                parts.append(f"{c.formula}: MP 에 없음")
                continue
            eh = pd.get_e_above_hull(e) * 1000
            if poly is not None:  # 광물: 그 광물의 다형 구조로 평가 (예: 석영 ≠ SiO2 바닥 다형)
                label = f"{c.name} ({c.formula}) E_hull"
                parts.append(f"{c.name} {eh:.0f} meV/atom")
                poly_notes.append(f"{c.name}: MP {c.mp_id} 다형으로 평가")
            else:
                label = f"{c.formula} E_hull"
                parts.append(f"{c.formula} {eh:.0f} meV/atom")
                if getattr(c, "mp_id", None) or c.ref.namespace.value == "mineral":
                    poly_notes.append(f"{c.name}: 광물 다형을 GGA/GGA+U 엔트리에서 찾지 못해 같은 조성의 최저 에너지 다형으로 평가")
            values.append(val(label, round(eh, 1), "meV/atom", kind=ValueKind.ENERGY))
        summary = "반응물 안정성 — " + ", ".join(parts)

    return Outcome(
        status=Status.OK, fidelity=Fidelity.L0, engine="pymatgen", engine_version=_pmg(),
        conditions_basis="0 K DFT", energy_reference=mp.ENERGY_REFERENCE, values=values,
        sources=[MP_SRC], summary=summary,
        caveats=["0 K 계산 — 고온에서만 안정한 다형은 hull 위로 나올 수 있음",
                 "구조 탐색은 MP 에 이미 있는 구조로 한정 — 새 구조는 L2(uMLIP) 단계에서",
                 *([] if ctx.recipe.mode is Mode.SYSTEM else poly_notes)],
        data={"table": {"columns": ["상", "MP ID", "형성에너지 (eV/atom)", "E_hull (meV/atom)", "판정"], "rows": rows}},
    )


# ── A2 ────────────────────────────────────────────────────────────────────

MISMATCH_HULL = 0.025  # eV/atom — NASA 평형은 안정인데 하이브리드가 이보다 높게 hull 위로 보면 불일치 (앵커 오차 폭)


def _kinks(ir: InterfacialReactivity) -> list[dict[str, Any]]:
    return [{"x": round(float(x), 4), "e": round(float(e) * 1000, 1), "rxn": str(rxn)} for _, x, e, rxn, _ in ir.get_kinks()]


def _rxn_phases(ir: InterfacialReactivity) -> set[str]:
    """곡선의 반응식에 나오는 상 (약분 화학식)."""
    return {c.reduced_formula for _, _, _, rxn, _ in ir.get_kinks() for c in [*rxn.reactants, *rxn.products]}


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


def _products_at(pd: PhaseDiagram, c1: Composition, c2: Composition, x: float) -> list[str]:
    """혼합비 x(원자 기준)의 평형 생성물 — 상태도 위 분해, 많은 순."""
    f1, f2 = c1.fractional_composition, c2.fractional_composition
    mix = Composition({el: x * f1[el] + (1 - x) * f2[el] for el in {*f1.elements, *f2.elements}})
    dec = sorted(pd.get_decomposition(mix).items(), key=lambda p: -p[1])
    return [e.composition.reduced_formula for e, amt in dec if amt > 1e-6]


def _nasa_eq(formula: str, T: float, P: float) -> dict[str, Any]:
    """반응물 하나가 T 에서 분해되는지 — NASA 실험 열화학. 1순위 Cantera 번들 평형(하이브리드와 다른 판본의
    데이터라 독립적), 번들에 없으면 thermo.inp 단독 상태도. 판정 불가면 note 만 돌려준다."""
    try:
        eq = nasa.multiphase_equilibrium({formula: 1.0}, T, P)
        T_used, note = T, None
    except nasa.OutOfRange as exc:
        if exc.tmax >= T:
            return {"note": str(exc)}
        T_used = exc.tmax - 1
        eq = nasa.multiphase_equilibrium({formula: 1.0}, T_used, P)
        note = f"NASA 데이터(Cantera 번들)가 {exc.tmax:.0f} K 까지라 {T_used:.0f} K 에서 계산"
    except KeyError:
        if nasa.formation_gibbs(formula, T) is None:
            return {"note": "NASA DB(Cantera 번들·thermo.inp)에 없거나 이 온도의 데이터가 없음"}
        comp = Composition(formula)
        npd = th.nasa_only_pd([str(e) for e in comp.elements], T)
        target = next(e for e in npd.all_entries if e.composition.reduced_formula == comp.reduced_formula)
        dec = npd.get_decomposition(comp)
        return {"T_used": T, "decomposed": npd.get_e_above_hull(target) > 5e-4,
                "products": {e.attribute.get("species", e.name): round(a * comp.num_atoms / e.composition.num_atoms, 3)
                             for e, a in dec.items() if a > 1e-6},
                "note": "Cantera 번들에 없어 thermo.inp 단독 상태도로 확인"}
    same = {s.name for s in nasa.find_condensed(formula)}
    left = sum(v for n, v in eq.condensed.items() if n in same)
    out = {"T_used": T_used, "decomposed": left < 0.5,
           "products": {**{k: round(v, 3) for k, v in eq.condensed.items()},
                        **({"기체": round(eq.gas_moles, 3)} if eq.gas_moles > 1e-6 else {})}}
    return {**out, "note": note} if note else out


def _cross_check(ctx: Context, hyb: th.HybridPD | None, gibbs_pd: PhaseDiagram | None,
                 gibbs_entries: list | None) -> dict[str, Any]:
    """평형 트랙 교차검증 — NASA 실험 열화학으로 각 반응물이 T 에서 분해되는지 보고 하이브리드 상태도의
    E_hull 과 비교한다. 둘이 어긋날 때만 불일치로 표시한다."""
    T = ctx.T
    out: dict[str, Any] = {"T": T, "checks": [], "mismatch": False}
    if T is None:
        return out
    for c in ctx.comps:
        item: dict[str, Any] = {"formula": c.formula, **_nasa_eq(c.formula, T, ctx.P)}
        if "decomposed" not in item:
            out["checks"].append(item)
            continue
        notes = [item.pop("note")] if item.get("note") else []
        T_used = item["T_used"]
        if gibbs_pd is not None and gibbs_entries is not None:
            e = _best(gibbs_pd, gibbs_entries, c.formula)
            if e is not None:
                item["mp_gibbs_e_hull_meV"] = round(gibbs_pd.get_e_above_hull(e) * 1000, 1)
        eh = hyb.e_hull(c.formula) if hyb is not None else None
        if eh is not None:
            src = hyb.source_of(c.formula)
            item.update({"hybrid_e_hull_meV": round(eh * 1000, 1), "hybrid_source": src})
            if item["decomposed"] and eh < 0.005 and T_used <= hyb.T + 1:
                verdict, bad = "불일치", True  # 기체가 나오는 분해는 온도가 높을수록 더 유리 — 낮은 T 의 '분해'는 T 에서도 분해
            elif not item["decomposed"] and eh > MISMATCH_HULL:
                bad = abs(T_used - hyb.T) <= 1
                verdict = "불일치" if bad else "온도가 달라 비교 보류"
            else:
                verdict, bad = "일치", False
            notes.append(f"하이브리드 {hyb.T:.0f} K E_hull {eh * 1000:.1f} meV/atom ({src} 값) — {verdict}")
            if bad:
                item["mismatch"] = out["mismatch"] = True
        if notes:
            item["note"] = " · ".join(notes)
        out["checks"].append(item)
    return out


def _unc(sources: set[str]) -> float | None:
    """하이브리드 곡선 값의 대략 오차(meV/atom) — 반응식에 쓰인 MP 상의 앵커 방식 검증 MAE."""
    if th.MP_GIBBS in sources:
        return float(th.SISSO_MAE)
    return float(th.ANCHOR_MAE["nk"]) if th.MP_0K in sources else None


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
    x = _recipe_x(ctx)
    caveats = ["열역학 판정만 — 반응 속도·핵생성은 반영하지 않음"]
    hyb = gibbs_pd = gibbs_entries = gibbs_only = None
    phases = {c.reduced_formula for c in (c1, c2)}
    T = ctx.T
    if T and T >= 298.0:  # 상온(298.15 K)은 SISSO 하한 300 K 로 계산 — 차이 무시할 만함
        Th = min(max(T, th.T_MIN), th.T_MAX)
        try:
            hyb = th.build(pd, Th)
        except (KeyError, ValueError) as exc:  # pymatgen 원소 G 표에 없는 원소 등
            caveats.append(f"하이브리드 상태도를 만들 수 없어 0 K 만: {type(exc).__name__}: {exc}")
        if hyb is not None:
            ir = InterfacialReactivity(c1, c2, hyb.pd, norm=True,
                                       use_hull_energy=any(hyb.best(c.formula) is None for c in comps))
            phases |= _rxn_phases(ir)
            curves.append({"label": f"{Th:.0f} K (하이브리드: NASA 실험 ΔGf + MP)", "kinks": _kinks(ir)})
            # 참고용 — v0 방식(MP 전 다형 + SISSO)만 쓴 곡선. 탄산염 과안정 문제를 보여 주려고 값으로만 남긴다
            Tg = int(min(th.T_MAX, max(th.T_MIN, round(T / 100) * 100)))
            gibbs_entries = GibbsComputedStructureEntry.from_entries(entries, temp=Tg)
            gibbs_pd = PhaseDiagram(gibbs_entries)
            gibbs_only = {"label": f"{Tg} K (MP Gibbs 근사만, 참고)", "kinks": _kinks(
                InterfacialReactivity(c1, c2, gibbs_pd, norm=True, use_hull_energy=use_hull))}

    sources_by_phase = {f: hyb.source_of(f) for f in sorted(phases) if hyb is not None and hyb.source_of(f)}
    unc = _unc(set(sources_by_phase.values()))
    values = []
    for cv in [*curves, *([gibbs_only] if gibbs_only else [])]:
        u = unc if hyb is not None and cv is curves[-1] else None
        low = min(cv["kinks"], key=lambda k: k["e"])
        values.append(val(f"최저 반응에너지 · {cv['label']}", low["e"], "meV/atom", kind=ValueKind.ENERGY, unc=u))
        if x is not None:
            values.append(val(f"레시피 혼합비에서 · {cv['label']}", round(_interp(cv["kinks"], x), 1),
                              "meV/atom", kind=ValueKind.ENERGY, unc=u))
    main = curves[-1]
    low = min(main["kinks"], key=lambda k: k["e"])
    main_pd = hyb.pd if hyb is not None else pd
    products = _products_at(main_pd, c1, c2, x) if x is not None else None
    if low["e"] < -1.0:
        summary = f"{main['label']}: {low['rxn']} — 반응에너지 {low['e']:.0f} meV/atom, 열역학적으로 유리(ΔG<0)"
        if products:
            summary += (f" · 레시피 혼합비(x={x:.2f})의 평형 생성물 {' + '.join(products)}"
                        f" ({_interp(main['kinks'], x):.0f} meV/atom)")
    else:
        summary = f"{main['label']}: 두 성분 사이에 열역학적으로 유리한 반응이 없음"

    check = _cross_check(ctx, hyb, gibbs_pd, gibbs_entries)
    warnings = []
    status = Status.OK
    if check["mismatch"]:
        status = Status.WARNING
        for c in (c for c in check["checks"] if c.get("mismatch")):
            warnings.append(
                f"교차검증 불일치: NASA 실험 열화학은 {c['formula']} 가 {c['T_used']:.0f} K 에서 "
                f"{'분해된다' if c['decomposed'] else '안정하다'}고 하지만 하이브리드 상태도는 "
                f"E_hull {c['hybrid_e_hull_meV']} meV/atom ({c['hybrid_source']} 값) → 이 온도의 L0 반응 판정은 신뢰 낮음")

    data: dict[str, Any] = {"curves": curves, "recipe_x": x, "labels": [comps[0].formula, comps[1].formula],
                            "crosscheck": check}
    if products is not None:
        data["products_at_recipe"] = products
    if hyb is not None:
        by_src: dict[str, list[str]] = {}
        for f, s in sources_by_phase.items():
            sp = hyb.nasa_species.get(f, f) if s == th.NASA else f
            by_src.setdefault(s, []).append(f if sp == f else f"{f}({sp})")
        caveats += [
            "T 곡선은 에너지 출처 혼합(계획서 3.5절): NASA 실험 ΔGf(T)가 있는 조성은 그 값, 나머지는 MP. "
            "MP 전용 상은 NASA 앵커 상으로부터의 반응에너지로 보정 — 생성물이 모두 고체면 0 K DFT 반응에너지"
            "(노이만–코프, 'MP-0K'), 기체가 끼면 SISSO G(T)('MP-Gibbs')",
            "상별 출처 — " + "; ".join(f"{s}: {', '.join(fs)}" for s, fs in sorted(by_src.items())),
            f"알려진 오차: 앵커 보정 검증 MAE 0 K {th.ANCHOR_MAE['nk']} · SISSO {th.ANCHOR_MAE['sisso']} meV/atom "
            f"(NASA 삼원 산화물 11종), SISSO 자체 약 {th.SISSO_MAE} meV/atom — 그보다 작은 차이는 판정 보류",
        ]
        liquid = [f"{f}({hyb.nasa_species[f]})" for f in (c.reduced_formula for c in (c1, c2))
                  if hyb.nasa_species.get(f, "").endswith("(L)") and hyb.source_of(f) == th.NASA]
        if liquid:
            caveats.append(f"NASA 기준 {hyb.T:.0f} K 에서 액체인 반응물: {', '.join(liquid)} — 고상이 아니라 용융염 반응일 수 있음")
        if abs(T - hyb.T) > 2:
            caveats.append(f"요청 {T:.0f} K 는 MP Gibbs 적용 범위({th.T_MIN:.0f}–{th.T_MAX:.0f} K) 밖이라 {hyb.T:.0f} K 에서 계산")
        if hyb.nasa_out_of_range:
            caveats.append("NASA 에 있으나 이 온도가 유효 구간 밖이라 MP 로 채운 상: " + ", ".join(hyb.nasa_out_of_range[:8]))
        data["sources_by_phase"] = sources_by_phase
        data["hybrid"] = {
            "T": hyb.T, "method": hyb.method, "nasa_data": nasa.nasa9_source(),
            "nasa_species": {f: hyb.nasa_species[f] for f in sources_by_phase if sources_by_phase[f] == th.NASA},
            "e_hull_meV": {c.formula: round(eh * 1000, 1) for c in comps if (eh := hyb.e_hull(c.formula)) is not None},
            "anchoring": {f: hyb.anchoring[f] for f in sources_by_phase if f in hyb.anchoring},
            "coordinate_check_meV": hyb.compare,  # 앵커 상의 NASA − MP Gibbs (원소 기준 맞춘 뒤)
            "element_shift_meV": hyb.element_shift,  # 원소 G: pymatgen 표 − NASA
        }
        if gibbs_only:
            data["gibbs_only"] = gibbs_only
    elif T is None or T < 298.0:
        caveats.append("온도 조건이 없거나 상온 미만이라 0 K DFT 만 — 기체가 나오는 반응(탄산염 분해 등)은 판정이 틀릴 수 있음")

    return Outcome(
        status=status, fidelity=Fidelity.L0, engine="pymatgen", engine_version=_pmg(),
        conditions_basis=f"0 K DFT + {curves[-1]['label']}" if len(curves) > 1 else "0 K DFT",
        energy_reference=(f"0 K 곡선: {mp.ENERGY_REFERENCE} · T 곡선(주 판정): {th.ENERGY_REFERENCE}"
                          if hyb is not None else mp.ENERGY_REFERENCE),
        values=values, sources=[MP_SRC, ("nasa-cea-thermo", nasa.nasa9_source())],
        caveats=caveats, warnings=warnings, summary=summary, data=data,
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
