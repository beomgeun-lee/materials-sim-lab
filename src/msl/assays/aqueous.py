"""A5 수계 안정성·부식 — Pourbaix 도표 (L0, 25 °C, 계획서 3.2절).

대상: 금속(1~3종) + O·H 로만 된 고체 성분이 물·수용액과 함께 있을 때.
판정: 레시피의 pH 와 전위에서 그 성분이 Pourbaix hull 위에 있는지(ΔG_pbx), 아니면 어떤 화학종이 안정한지 보고
면역 / 부동태 / 부식으로 나눈다 (Pourbaix 도감의 관례, 이온 농도 10⁻⁶ mol/kg).

- pH: 레시피 조건 → A4 계산값 → 물·용액 성분만으로 PHREEQC 계산 → 7(중성 가정) 순.
- 전위: 레시피 Eh → 분위기(공기·O₂ = O₂/H₂O 선, 불활성 = H₂/H₂O 선) → 지정이 없으면 두 극단을 모두 본다.
"""

from __future__ import annotations

import math
import re
from importlib.metadata import version

from pymatgen.core import Composition

from msl.assays.base import Context, Outcome, not_applicable, val
from msl.connectors import mp_ion_ref_data
from msl.engines import mp, nasa
from msl.engines import pourbaix as pb
from msl.recipe import METALS_EXCLUDED
from msl.resolve import Resolved
from msl.schema.recipe import State
from msl.schema.refs import Namespace
from msl.schema.result import Fidelity, Status

ENGINE = "pymatgen"
MAX_METALS = 3  # 다원소 Pourbaix 는 원소 수에 따라 조합이 폭증한다
TOL = 1e-3  # eV/atom — 이 이하면 hull 위(안정)로 본다
LIMITS = [[0.0, 14.0], [-2.0, 2.0]]  # 차트 범위: pH, V vs SHE
AIR_PO2 = 0.21
AIR = {"air", "공기"}
OXYGEN = {"o2", "oxygen", "산소"}
INERT = {"n2", "ar", "argon", "nitrogen", "질소", "아르곤", "inert", "불활성", "탈기", "h2", "수소"}
VERDICT_KO = {"immune": "면역", "stable": "안정", "passive": "부동태", "transform": "고체 변환",
              "partial": "일부 용해", "corrosion": "부식", "dissolve": "용해"}

_SUB = str.maketrans("0123456789.", "₀₁₂₃₄₅₆₇₈₉.")
_SUP = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")


def pretty(name: str) -> str:
    """pymatgen Pourbaix 이름 → 읽기 쉬운 표기: 'Fe[+2]' → 'Fe²⁺', 'Mg(HO)2(s)' → 'Mg(OH)₂(s)'."""
    name = name.replace("(HO)", "(OH)")
    subs = lambda body: re.sub(r"(?<=[A-Za-z)])(\d+(?:\.\d+)?)", lambda m: m.group(1).translate(_SUB), body)
    m = re.match(r"^(.*)\[([+-])(\d+(?:\.\d+)?)\]$", name)
    if m:
        body, sign, n = m.groups()
        n = str(int(float(n))) if float(n).is_integer() else n
        return subs(body) + (("" if n == "1" else n) + sign).translate(_SUP)
    body, suffix = re.match(r"^(.*?)(\((?:s|aq)\))?$", name).groups()
    return subs(body) + (suffix or "")


def target(c: Resolved) -> list[str] | None:
    """Pourbaix 대상이면 금속 원소 목록. 금속 + O·H 로만 된 고체 성분."""
    if c.formula is None or c.props or c.is_water or c.state not in (None, State.SOLID):
        return None
    els = {str(e) for e in Composition(c.formula).elements}
    metals = sorted(e for e in els if e not in METALS_EXCLUDED)
    if not metals or not (els - set(metals)) <= {"O", "H"}:
        return None
    return metals


def _pH(ctx: Context, targets: list[Resolved]) -> tuple[float, str]:
    if ctx.recipe.conditions.pH is not None:
        return ctx.recipe.conditions.pH, "레시피 조건"
    if "pH" in ctx.shared:
        return ctx.shared["pH"], "A4 수용액 평형 계산값"
    rest = [c for c in ctx.comps if c not in targets]
    if any(c.is_water or c.state is State.AQUEOUS for c in rest):
        from msl.assays.equilibrium import a4

        sub = Context(ctx.recipe, rest, ctx.registry, {})
        try:
            if a4(sub).status is Status.OK and "pH" in sub.shared:
                return sub.shared["pH"], "물·용액 성분만의 PHREEQC 계산값"
        except Exception:
            pass
    return 7.0, "지정·계산 없음 → 중성 가정"


def _scenarios(ctx: Context, pH: float) -> list[tuple[str, float]]:
    cond = ctx.recipe.conditions
    if cond.Eh is not None:
        return [("지정 전위", cond.Eh)]
    h2, o2 = pb.water_lines(pH)
    air = ("공기 포화 (O₂/H₂O, pO₂ 0.21)", o2 + pb.NERNST / 4 * math.log10(AIR_PO2))
    reduced = ("탈기 (H₂/H₂O)", h2)
    atm = (cond.atmosphere or "").strip().lower()
    if atm in AIR:
        return [air]
    if atm in OXYGEN:
        return [("산소 포화 (O₂/H₂O, pO₂ 1)", o2)]
    if atm in INERT:
        return [reduced]
    return [air, reduced]


def _own_entry(ents: list, comp: Composition):
    same = [e for e in ents if e.phase_type == "Solid"
            and e.entry.composition.reduced_composition.almost_equals(comp.reduced_composition)]
    return min(same, key=lambda e: e.entry.energy_per_atom) if same else None


def _verdict(metal: bool, own: bool, parts: list) -> tuple[str, str]:
    ions = [pretty(p.name) for p in parts if p.phase_type == "Ion"]
    solids = [pretty(p.name) for p in parts if p.phase_type == "Solid"]
    if own:
        return ("immune", "면역 — 금속 그대로 안정") if metal else ("stable", "안정 — 그대로 유지")
    if not ions:
        names = " + ".join(solids)
        return ("passive", f"부동태 — 표면에 {names} 생성") if metal else ("transform", f"{names}(으)로 바뀜")
    if solids:
        return "partial", f"일부 용해 — {', '.join(ions)} 와 {' + '.join(solids)}"
    return ("corrosion" if metal else "dissolve"), f"{'부식' if metal else '용해'} — {', '.join(ions)}(으)로 녹음"


def a5(ctx: Context) -> Outcome:
    if not any(c.is_water or c.state is State.AQUEOUS for c in ctx.comps):
        return not_applicable("물이나 수용액 성분이 없음", ENGINE)
    targets = [(c, m) for c in ctx.comps if (m := target(c))]
    if not targets:
        return not_applicable("금속이나 금속 산화물·수산화물 고체 성분이 없음", ENGINE)
    too_many = [c.name for c, m in targets if len(m) > MAX_METALS]
    targets = [(c, m) for c, m in targets if len(m) <= MAX_METALS]
    if not targets:
        return not_applicable(f"금속 원소가 {MAX_METALS}종을 넘는 성분만 있음: {', '.join(too_many)}", ENGINE)

    pH, pH_how = _pH(ctx, [c for c, _ in targets])
    scenarios = _scenarios(ctx, pH)
    values = [val("pH", round(pH, 2)), val("금속 이온 농도 (판정 기준)", "10⁻⁶", "mol/kg")]
    rows, charts, lines, warnings = [], [], [], []
    mp_phases: set[str] = set()
    for c, metals in targets:
        comp = Composition(c.formula)
        name = c.formula if c.ref.namespace is Namespace.ELEMENT else c.name  # 원소는 기호로 (PubChem 영문명 대신)
        ratio = {m: comp[m] for m in metals}
        total = sum(ratio.values())
        ents, pd, source = pb.diagram(metals, {m: v / total for m, v in ratio.items()})
        own = _own_entry(ents, comp)
        metal = all(str(e) in metals for e in comp.elements)
        points, verdicts = [], []
        for label, E in scenarios:
            stable = pd.get_stable_entry(pH, E)
            parts = getattr(stable, "entry_list", [stable])
            dg = float(pd.get_decomposition_energy(own, pH, E)) if own is not None else None
            code, text = _verdict(metal, dg is not None and dg <= TOL, parts)
            mp_phases |= {p.name for p in parts if p.phase_type == "Solid" and not source.get(p.name, "").startswith(pb.EXPERIMENTAL)}
            stable_txt = " + ".join(pretty(p.name) for p in parts)
            rows.append([name, label, round(pH, 2), round(E, 3), stable_txt, text, None if dg is None else round(dg, 3)])
            values.append(val(f"{name} · {label}", text))
            if dg is not None:
                values.append(val(f"ΔG_pbx {name} · {label}", round(dg, 4), "eV/atom"))
            points.append({"label": label, "pH": round(pH, 3), "E": round(E, 3), "verdict": code, "text": text})
            verdicts.append(f"{label.split(' (')[0]} {VERDICT_KO[code]}({stable_txt})")
        doms = pb.domains(pd, LIMITS)
        for d in doms:
            d["label"] = " + ".join(pretty(p["name"]) for p in d["parts"])
            d["source"] = " · ".join(sorted({source.get(p["name"], "?") for p in d["parts"]}))
        charts.append({"target": name, "formula": c.formula, "elements": metals, "conc": pb.DEFAULT_CONC,
                       "limits": LIMITS, "domains": doms, "points": points})
        lines.append(f"{name}: " + " · ".join(verdicts))

    if too_many:
        warnings.append(f"금속 원소가 {MAX_METALS}종을 넘어 빠진 성분: {', '.join(too_many)}")
    if mp_phases:
        warnings.append(f"실험값이 없어 MP 계산(실험 앵커 보정)에 기댄 상: {', '.join(pretty(n) for n in sorted(mp_phases))} "
                        "— 실험 대비 수십~수백 meV/atom 틀릴 수 있다")
    if ctx.T is not None and abs(ctx.T - pb.T25) > 10:
        warnings.append(f"Pourbaix 는 25 °C 기준 — 입력 온도 {ctx.T - 273.15:.0f} °C 는 반영하지 않음")
    if any("Cl" in {str(e) for e in Composition(x.formula).elements} for x in ctx.comps if x.formula):
        warnings.append("염화물 이온(Cl⁻)은 착이온·공식(孔蝕)으로 부식을 키우지만 이 도표(금속–O–H 계)에는 들어 있지 않음")

    caveats = [
        "평형 열역학 — '부동태'는 산화물·수산화물이 안정하다는 뜻일 뿐 막이 치밀하고 보호성인지는 모름 (예: 철의 녹은 보호성이 약함)",
        "공기 포화 전위는 O₂/H₂O 평형 전위(이론 상한). 실제 자연수의 산화환원 전위는 보통 이보다 0.3~0.6 V 낮음",
        "착화제(Cl⁻·NH₃·유기산 등)와 속도(과전압)는 반영하지 않음",
        f"고체: 실험 ΔGf(NASA·NBS) 우선, 없으면 MP 계산을 실험 앵커에 이어 붙임. MP 상은 실험 관측(ICSD) 근거가 있는 것만 (D15). "
        f"이온: 실험 ΔGf ({len(mp_ion_ref_data.records())}종)",
    ]
    return Outcome(
        status=Status.WARNING if warnings else Status.OK, fidelity=Fidelity.L0, engine=ENGINE,
        engine_version=version("pymatgen"), energy_reference=pb.ENERGY_REFERENCE,
        conditions_basis=f"25 °C, 1 atm, 금속 이온 10⁻⁶ mol/kg, pH {pH:.2f} ({pH_how})",
        values=values, caveats=caveats, warnings=warnings,
        sources=[("materials-project", mp.DB_VERSION), ("mp-ion-ref-data", mp_ion_ref_data.VERSION),
                 ("nasa-cea-thermo", nasa.nasa9_source())],
        summary=" / ".join(lines) + f" — pH {pH:.2f} ({pH_how})",
        data={"pourbaix": charts, "verdict": "; ".join(lines) if len(lines) > 1 else lines[0].split(": ", 1)[1],
              "table": {"columns": ["성분", "조건", "pH", "E (V vs SHE)", "안정한 화학종", "판정", "ΔG_pbx (eV/atom)"], "rows": rows}},
    )
