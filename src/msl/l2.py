"""L2 안정성 확인 — 치환 구조 생성 → uMLIP 이완 → 자기일관 hull → 두 모델 앙상블 (계획서 4단계, D19).

1. 경쟁 상: 대상 원소계의 MP 상(조성마다 바닥 엔트리, MP hull 거리 ≤ 0.1 eV/atom)을 같은 uMLIP 로 다시 이완해 hull 을 만든다.
   현재 MP(v2026) GGA+U 에너지는 LASPH 를 켠 설정이라 MPtrj 로 학습한 uMLIP 과 전이금속 산화물에서 0.05~0.19 eV/atom
   어긋난다 → uMLIP 값을 MP hull 에 바로 대지 않는다. 모든 에너지에는 MP2020 보정을 똑같이 붙인다.
2. 후보 구조 (v2): 대상과 음이온:양이온 비가 같은 같은 원소계 모체의 양이온 자리 배치 — 정전기(Ewald) 에너지가 낮은 순 +
   무작위 몇 개, 그리고 다른 화학계의 같은 조성 패턴(익명 화학식) 안정 구조를 대상 원소로 장식한 원형.
   빠른 모델(ORB)로 모두 이완해 거르고 상위 몇 개만 두 모델로 계산한다. DB 에 있는 조성은 그 구조를 그대로 이완한다.
3. 두 모델(MACE-MPA-0, ORB v3)로 따로 계산해 hull 거리의 평균과 차이를 낸다 — 차이가 크면 불확실.

한계: 배치는 Ewald 순위 상위와 무작위 표본이라 전역 최저를 보장하지 않는다. 원형은 조성 패턴이 정확히 같은 것만 찾는다.
"""

from __future__ import annotations

import json
import time
from fractions import Fraction
from itertools import product
from typing import Any, Callable

import numpy as np
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition, Structure

from msl.engines import mp, umlip
from msl.env import CACHE_DIR

CACHE = CACHE_DIR / "umlip"
REF_MAX_EHULL = 0.05  # eV/atom — 경쟁 상으로 다시 계산할 MP 상의 범위 (uMLIP 이 순서를 바꿀 여지)
REF_MAX_SITES = 40
L3_WINDOW = 0.05  # eV/atom — L2 hull 거리가 이 안이면 DFT(L3)로 확정할 가치가 있다
L3_SPREAD = 0.03  # eV/atom — 두 모델 차이가 이보다 크면 L2 를 믿기 어렵다


def l3_rule(ehull: float, spread: float, phonon: dict | None, known: bool) -> dict[str, Any]:
    """L2 → L3(DFT) 승격 규칙 (계획서 4단계). 새 조성이고 hull 경계 근처(|ehull| ≤ 0.05)면 승격 후보,
    두 모델이 크게 어긋나거나 포논이 불안정하면 그 이유를 붙인다. 명확히 불안정(> 0.1)하면 승격하지 않는다."""
    reasons, verdict = [], "불필요"
    if known:
        return {"verdict": "불필요", "reasons": ["DB 에 DFT 계산이 이미 있음"]}
    if ehull <= L3_WINDOW:
        verdict = "권장"
        reasons.append(f"hull 거리 {ehull:+.3f} eV/atom — 안정·준안정 경계 안 (uMLIP 오차보다 작은 차이는 DFT 로 확정)")
    elif ehull <= 2 * L3_WINDOW:
        verdict = "선택"
        reasons.append(f"hull 거리 {ehull:+.3f} — 준안정 바깥이지만 uMLIP 오차(수십 meV) 안")
    else:
        reasons.append(f"hull 거리 {ehull:+.3f} — 명확히 불안정, DFT 계산 가치 낮음")
    if spread > L3_SPREAD:
        reasons.append(f"두 모델 차이 {spread:.3f} > {L3_SPREAD} — L2 결과 불확실")
        verdict = "권장" if verdict != "불필요" else "선택"
    if phonon and not phonon.get("error") and not phonon.get("dynamically_stable", True):
        reasons.append(f"포논 허수 모드 (최소 {phonon['min_freq_THz']:+.2f} THz) — 구조가 더 낮은 대칭으로 뒤틀릴 수 있음, DFT 로 재이완 필요")
    return {"verdict": verdict, "reasons": reasons}


def _cached_relax(structure: Structure, key: str, model: str) -> dict[str, Any]:
    """이완 결과 캐시 (모델·키별). {energy, structure, steps, converged, seconds}."""
    path = CACHE / model / f"{key}.json"
    if path.exists():
        d = json.loads(path.read_text(encoding="utf-8"))
        d["structure"] = Structure.from_dict(d["structure"])
        return d
    r = umlip.relax(structure, model, fmax=0.05, steps=400)
    d = {"energy": r.energy, "structure": r.structure.as_dict(), "steps": r.steps, "converged": r.converged, "seconds": r.seconds}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(d), encoding="utf-8")
    d["structure"] = r.structure
    return d


def reference_entries(elements: list[str]) -> list:
    """경쟁 상 — 조성마다 MP 바닥 엔트리. MP hull 위(안정) 상은 크기와 관계없이 항상, 준안정(≤ 0.05)은 원자 ≤ 40 만.

    안정 상을 크기로 빼면 hull 이 실제보다 높아져 새 조성이 과하게 안정해 보인다
    (K–Al–O 의 KAlO2 64원자, Li–Al–O 의 LiAl5O8·Li5AlO4 가 빠져 0.1~0.2 eV/atom 과대평가, D24)."""
    ents = mp.entries_in_chemsys(elements)
    pd = PhaseDiagram(ents)
    best: dict[str, Any] = {}
    for e in ents:
        f = e.composition.reduced_formula
        if f not in best or e.energy_per_atom < best[f].energy_per_atom:
            best[f] = e
    out = []
    for f, e in best.items():
        is_element = len(e.composition.elements) == 1
        eh = pd.get_e_above_hull(e)
        if is_element or eh <= 1e-6 or (eh <= REF_MAX_EHULL and len(e.structure) <= REF_MAX_SITES):
            out.append(e)
    return out


def uhull(elements: list[str], model: str, log: Callable[[str], None] = lambda s: None,
          exclude: Composition | None = None) -> tuple[PhaseDiagram, list[str]]:
    """자기일관 hull — 경쟁 상을 모두 같은 uMLIP 로 이완 (MP2020 보정). exclude: 이 조성은 경쟁 상에서 뺀다 (모르는 척 평가)."""
    refs = reference_entries(elements)
    if exclude is not None:
        refs = [e for e in refs if not e.composition.reduced_composition.almost_equals(exclude.reduced_composition)]
    entries, notes = [], []
    t0 = time.time()
    for i, e in enumerate(refs):
        d = _cached_relax(e.structure, str(e.entry_id), model)
        ce = umlip.mp_entry(d["structure"], d["energy"], entry_id=str(e.entry_id))
        if ce is None:
            notes.append(f"{e.composition.reduced_formula}: MP2020 보정 불가 — 제외")
            continue
        entries.append(ce)
        if (i + 1) % 10 == 0:
            log(f"    {model} 경쟁 상 {i + 1}/{len(refs)} ({time.time() - t0:.0f}초)")
    return PhaseDiagram(entries), notes


def _anion(comp: Composition):
    return max(comp.elements, key=lambda e: e.X if e.X == e.X else 0.0)


def parents(target: Composition, limit: int = 3, exclude_same: bool = False) -> list[tuple[Any, Structure]]:
    """대상과 음이온이 같고 음이온:양이온 비가 같은 알려진 구조 (MP hull 가까운 순) — (엔트리, 기본 셀)."""
    an = _anion(target)
    ratio = Fraction(target[an]).limit_denominator(100) / Fraction(target.num_atoms - target[an]).limit_denominator(100)
    ents = mp.entries_in_chemsys(sorted(str(e) for e in target.elements))
    pd = PhaseDiagram(ents)
    cands = []
    for e in ents:
        c = e.composition
        if an not in c.elements or len(c.elements) < 2 or len(e.structure) > 48:
            continue
        if exclude_same and c.reduced_composition.almost_equals(target.reduced_composition):
            continue  # 재발견 시험 — 정답 구조는 모체로 쓰지 않는다
        if str(an) == "O" and e.data.get("oxide_type", "oxide") != "oxide":
            continue  # 과산화물·초산화물(O–O 결합) 구조에 양이온을 치환하면 화학적으로 다른 물질
        r = Fraction(c[an]).limit_denominator(100) / Fraction(c.num_atoms - c[an]).limit_denominator(100)
        if r == ratio:
            cands.append((pd.get_e_above_hull(e), e))
    by_comp: dict[str, list] = {}
    for eh, e in cands:
        by_comp.setdefault(e.composition.reduced_formula, []).append((eh, e))
    ranked = sorted(by_comp.values(), key=lambda v: min(x[0] for x in v))
    out = []
    for group in ranked[:limit]:
        low = min(x[0] for x in group)
        near = [e for eh, e in group if eh <= low + 0.02]  # 에너지가 비슷한 다형 중 기본 셀이 가장 작은 것 (예: 층상 LiNiO2 4원자)
        prims = sorted(((e, e.structure.get_primitive_structure()) for e in near), key=lambda p: len(p[1]))
        out.append(prims[0])
    return out


def _scalings(struct: Structure, max_atoms: int):
    """대각 초격자 배수 — 작은 것부터, 고르게 늘린 것 먼저."""
    for a, b, c in sorted(product(range(1, 6), repeat=3), key=lambda t: (t[0] * t[1] * t[2], max(t) - min(t))):
        if len(struct) * a * b * c <= max_atoms:
            yield (a, b, c)


def orderings(parent: Structure, target: Composition, n: int, max_atoms: int, seed: int = 0) -> list[Structure]:
    """모체의 양이온 자리에 대상 양이온을 무작위 배치한 구조 n 개 (양이온 수가 맞는 가장 작은 초격자)."""
    an = _anion(target)
    cations = [e for e in target.elements if e != an]
    base = Composition(Composition({e: target[e] for e in cations}).get_integer_formula_and_factor()[0])  # 정수 양이온 비
    unit_total = int(round(base.num_atoms))
    rng = np.random.default_rng(seed)
    for scale in _scalings(parent, max_atoms):
        s = parent.copy()
        s.make_supercell(list(scale))
        cat_sites = [i for i, site in enumerate(s) if site.specie != an]
        if not cat_sites or len(cat_sites) % unit_total:
            continue
        m = len(cat_sites) // unit_total
        species = [el for el in cations for _ in range(int(round(base[el] * m)))]
        if len(species) != len(cat_sites):
            continue
        out = []
        for _ in range(n):
            t = s.copy()
            for i, el in zip(cat_sites, rng.permutation(species)):
                t.replace(i, el)
            out.append(t)
        return out
    return []


def charges(target: Composition) -> dict[str, float]:
    """Ewald 순위용 형식 전하. 산화수가 하나뿐인 원소(Li⁺·Mg²⁺ 등)는 고정, 음이온은 가장 흔한 음의 산화수,
    남은 전하는 나머지 원소가 고르게 나눈다 (예: Li6MnNi3O10 → Mn·Ni 평균 +3.5).
    pymatgen 의 자동 추정은 ICSD 빈도로 Mn⁷⁺ 같은 비현실적 배정을 고르기도 해서 쓰지 않는다."""
    an = _anion(target)
    q: dict[str, float] = {str(an): float(min(s for s in an.common_oxidation_states if s < 0)) if any(
        s < 0 for s in an.common_oxidation_states) else -2.0}
    variable = []
    for e in target.elements:
        if e == an:
            continue
        pos = sorted({s for s in (*e.common_oxidation_states, *e.icsd_oxidation_states) if s > 0})
        if len(pos) == 1:  # Ni 는 '흔한' 산화수가 +2 하나뿐이지만 ICSD 에는 +1~+4 → 가변으로 본다
            q[str(e)] = float(pos[0])
        else:
            variable.append(e)
    rest = -sum(target[e] * q[str(e)] for e in target.elements if str(e) in q)
    if variable:
        n_var = sum(target[e] for e in variable)
        for e in variable:
            q[str(e)] = rest / n_var
    return q


def ewald_orderings(parent: Structure, target: Composition, n: int, max_atoms: int) -> list[Structure]:
    """정전기(Ewald) 에너지가 낮은 양이온 배치 n 개 — 대칭상 같은 배치(같은 에너지)는 하나로."""
    from pymatgen.core import Species
    from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation

    an = _anion(target)
    q = charges(target)
    cations = [e for e in target.elements if e != an]
    base = Composition(Composition({e: target[e] for e in cations}).get_integer_formula_and_factor()[0])
    unit_total = int(round(base.num_atoms))
    for scale in _scalings(parent, max_atoms):
        s = parent.copy()
        s.make_supercell(list(scale))
        cat_sites = [i for i, site in enumerate(s) if site.specie != an]
        if not cat_sites or len(cat_sites) % unit_total:
            continue
        occ = {Species(str(e), q[str(e)]): base[e] / unit_total for e in cations}
        for i, site in enumerate(s):
            s.replace(i, Species(str(an), q[str(an)]) if site.specie == an else occ)
        try:
            ranked = OrderDisorderedStructureTransformation(algo=2).apply_transformation(s, return_ranked_list=max(4 * n, 20))
        except Exception:
            return []
        out, seen = [], set()
        for r in ranked:
            e = round(float(r["energy"]), 3)
            if e in seen:
                continue
            seen.add(e)
            st = r["structure"].copy()
            st.remove_oxidation_states()
            out.append(st)
            if len(out) >= n:
                break
        return out
    return []


def prototypes(target: Composition, limit: int = 4, max_sites: int = 40, exclude_chemsys: str | None = None) -> list[tuple[str, Structure]]:
    """다른 화학계의 같은 조성 패턴(익명 화학식) 안정 구조를 대상 원소로 장식한 구조.

    음이온은 음이온에, 나머지는 원자 수가 같은 원소끼리 대응시키되 원자 수가 같은 원소가 여럿이면
    평균 이온 반지름 차이가 가장 작은 배정을 고른다. MP 조회 결과는 캐시한다.
    """
    import itertools
    import pickle

    anon = target.anonymized_formula
    an = _anion(target)
    path = mp.CACHE / f"proto_{anon}_{an}.pkl"
    if path.exists():
        docs = pickle.loads(path.read_bytes())
    else:
        with mp._rester() as mpr:
            got = mpr.materials.summary.search(formula=anon, elements=[str(an)], energy_above_hull=(0, 0.02),
                                               num_sites=(1, max_sites), fields=["material_id", "formula_pretty", "chemsys",
                                                                                "energy_above_hull", "structure"])
        docs = [(str(d.material_id), d.formula_pretty, d.chemsys, float(d.energy_above_hull), d.structure) for d in got]
        mp.CACHE.mkdir(parents=True, exist_ok=True)
        path.write_bytes(pickle.dumps(docs))
    tgt_counts = {e: target.get_el_amt_dict()[str(e)] for e in target.elements}
    rad = lambda e: e.average_ionic_radius or e.atomic_radius or 1.0
    scored = []
    for mid, f, chemsys, eh, st in docs:
        if exclude_chemsys and chemsys == exclude_chemsys:
            continue
        pc = st.composition.reduced_composition
        p_an = _anion(pc)
        if str(p_an) != str(an):
            continue
        pcounts = pc.get_el_amt_dict()
        best = None
        others_t = [e for e in target.elements if e != an]
        others_p = [e for e in pc.elements if e != p_an]
        for perm in itertools.permutations(others_t):
            if any(abs(pcounts[str(pe)] - tgt_counts[te]) > 1e-6 for pe, te in zip(others_p, perm)):
                continue
            cost = sum(abs(float(rad(pe)) - float(rad(te))) for pe, te in zip(others_p, perm))
            if best is None or cost < best[0]:
                best = (cost, dict(zip(map(str, others_p), map(str, perm))))
        if best is None:
            continue
        scored.append((best[0], eh, mid, f, st, best[1]))
    scored.sort(key=lambda x: (x[0], x[1]))
    out = []
    for cost, eh, mid, f, st, mapping in scored[:limit]:
        s = st.copy()
        s.replace_species({**mapping, str(an): str(an)})
        out.append((f"원형 {f}({mid})", s))
    return out


def _key(label: str, s: Structure) -> str:
    """구조 내용으로 만든 캐시 키 (실행마다 같은 값)."""
    import hashlib

    blob = json.dumps([label, [round(x, 4) for x in s.lattice.parameters],
                       [[str(site.specie), *np.round(site.frac_coords, 4).tolist()] for site in s]])
    return "cand-" + hashlib.sha256(blob.encode()).hexdigest()[:20]


def candidates(target: Composition, n_ewald: int = 4, n_random: int = 2, max_atoms: int = 40, use_prototypes: bool = True,
               exclude_same: bool = False) -> list[tuple[str, Structure]]:
    """새 조성의 후보 구조: 같은 원소계 모체의 Ewald 배치 + 무작위 배치, 다른 화학계 구조 원형 장식.

    exclude_same: 대상 조성 자신(과 같은 원소계 원형)은 쓰지 않는다 — 재발견 시험용.
    """
    out: list[tuple[str, Structure]] = []
    for p, prim in parents(target, exclude_same=exclude_same):
        tag = f"{p.composition.reduced_formula}({p.entry_id})"
        out += [(f"{tag} Ewald {k + 1}", s) for k, s in enumerate(ewald_orderings(prim, target, n_ewald, max_atoms))]
        out += [(f"{tag} 무작위 {k + 1}", s) for k, s in enumerate(orderings(prim, target, n_random, max_atoms))]
    if use_prototypes:
        excl = "-".join(sorted(str(e) for e in target.elements)) if exclude_same else None
        try:
            out += prototypes(target, max_sites=max_atoms, exclude_chemsys=excl)
        except Exception:
            pass  # MP 조회 실패 — 같은 원소계 후보만으로
    return out


def evaluate(formula: str, models: tuple[str, ...] = ("mace-mpa-0", "orb-v3"), n_orderings: int = 4, max_atoms: int = 40,
             screen_model: str = "orb-v3", final_top: int = 3, phonon_check: bool = True,
             log: Callable[[str], None] = print) -> dict[str, Any]:
    """조성 하나의 L2 hull 거리 (두 모델). DB 에 있는 조성은 그 구조를, 없는 조성은 후보 구조를 쓴다.

    후보가 여럿이면 빠른 모델(screen_model)로 모두 이완해 상위 final_top 개만 나머지 모델로 계산한다.
    """
    target = Composition(formula)
    els = sorted(str(e) for e in target.elements)
    t0 = time.time()
    ents = mp.entries_in_chemsys(els)
    same = [e for e in ents if e.composition.reduced_composition.almost_equals(target.reduced_composition)]
    if same:
        known = min(same, key=lambda e: e.energy_per_atom)
        structs = [("DB 구조 " + str(known.entry_id), known.structure)]
        mp_ehull = PhaseDiagram(ents).get_e_above_hull(known)
    else:
        known, mp_ehull = None, None
        structs = candidates(target, n_ewald=n_orderings, max_atoms=max_atoms)
        if not structs:
            raise ValueError(f"{target.reduced_formula}: 음이온:양이온 비가 같은 모체·구조 원형을 찾지 못함")
    log(f"{target.reduced_formula}: 후보 구조 {len(structs)}개 · 거르기 {screen_model} → 상위 {final_top}개를 {', '.join(models)}")
    screen: dict[str, Any] = {}
    if len(structs) > final_top:
        pd, _ = uhull(els, screen_model, log=log)
        scored = []
        for label, s in structs:
            d = _cached_relax(s, _key(label, s), screen_model)
            ce = umlip.mp_entry(d["structure"], d["energy"], entry_id=label)
            scored.append((float(pd.get_e_above_hull(ce, allow_negative=True)), label, s))
        scored.sort(key=lambda x: x[0])
        screen = {"model": screen_model, "ranking": [{"label": lb, "ehull": round(eh, 4)} for eh, lb, _ in scored]}
        structs = [(lb, s) for _, lb, s in scored[:final_top]]
    per_model: dict[str, Any] = {}
    for model in models:
        pd, notes = uhull(els, model, log=log)
        best = None
        for label, s in structs:
            d = _cached_relax(s, _key(label, s), model)
            ce = umlip.mp_entry(d["structure"], d["energy"], entry_id=label)
            eh = float(pd.get_e_above_hull(ce, allow_negative=True))
            if best is None or eh < best["ehull"]:
                best = {"ehull": eh, "label": label, "structure": d["structure"], "converged": d["converged"], "steps": d["steps"]}
        dec = pd.get_decomposition(best["structure"].composition)
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        try:
            sg = SpacegroupAnalyzer(best["structure"], symprec=0.1).get_space_group_symbol()
        except Exception:
            sg = "?"
        per_model[model] = {"ehull": round(best["ehull"], 4), "best": best["label"], "space_group": sg, "converged": best["converged"],
                            "decomposition": {e.composition.reduced_formula: round(float(x), 3) for e, x in sorted(dec.items(), key=lambda p: -p[1])},
                            "n_references": len(pd.all_entries), "notes": notes, "engine": umlip.engine_version(model)}
        log(f"  {model}: hull 거리 {best['ehull']:+.3f} eV/atom · 최저 구조 {best['label']} ({sg})")
    vals = [m["ehull"] for m in per_model.values()]
    phon = None
    if phonon_check and "mace-mpa-0" in per_model and float(np.mean(vals)) <= L3_WINDOW * 2:
        try:  # 포논은 float64 인 MACE 로 (ORB float32 는 NaCl 에서 −0.4 THz 잡음)
            from msl import props

            best_s = next(Structure.from_dict(json.loads((CACHE / "mace-mpa-0" / f"{_key(lb, s)}.json").read_text())["structure"])
                          for lb, s in structs if lb == per_model["mace-mpa-0"]["best"]) if not known else known.structure
            phon = props.phonon(best_s, "mace-mpa-0", relax=False)
            log(f"  포논 (MACE): 최소 {phon['min_freq_THz']:+.2f} THz → {'안정' if phon['dynamically_stable'] else '불안정'}")
        except Exception as exc:
            phon = {"error": str(exc)}
    return {"formula": target.reduced_formula, "known": None if known is None else str(known.entry_id),
            "phonon": phon, "l3": l3_rule(float(np.mean(vals)), float(max(vals) - min(vals)), phon, known is not None),
            "mp_ehull": None if mp_ehull is None else round(float(mp_ehull), 4), "n_structures": len(structs) if not screen else len(screen["ranking"]),
            "screen": screen, "models": per_model, "ehull_mean": round(float(np.mean(vals)), 4),
            "ehull_spread": round(float(max(vals) - min(vals)), 4), "seconds": round(time.time() - t0, 1), "method": "v2 Ewald+원형"}


RESULTS = CACHE / "results"


def saved(formula: str) -> dict[str, Any] | None:
    """이미 계산한 L2 결과 (조성별)."""
    path = RESULTS / f"{Composition(formula).reduced_formula}.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def save(res: dict[str, Any]) -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"{res['formula']}.json").write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")


def rediscover(formula: str, models: tuple[str, ...] = ("orb-v3", "mace-mpa-0"), method: str = "v2", final_top: int = 3,
               log: Callable[[str], None] = print) -> dict[str, Any]:
    """재발견 시험 — DB 에 있는 조성을 모른다고 치고 후보 구조를 만들어, 가장 낮은 후보와 DB 구조의 에너지 차이를 잰다.

    method: v1 = 같은 원소계 모체마다 무작위 배치 6개 (개선 전), v2 = Ewald 배치 + 무작위 + 다른 화학계 원형 (개선 후).
    차이(gap)가 0 에 가까우면 DB 구조(또는 그만큼 좋은 구조)를 다시 찾아낸 것이다. 첫 모델로 거르고 상위만 다음 모델로.
    """
    target = Composition(formula)
    els = sorted(str(e) for e in target.elements)
    ents = mp.entries_in_chemsys(els)
    known = min((e for e in ents if e.composition.reduced_composition.almost_equals(target.reduced_composition)),
                key=lambda e: e.energy_per_atom)
    if method == "v1":
        cands = [(f"{p.composition.reduced_formula} 무작위 {k + 1}", s)
                 for p, prim in parents(target, exclude_same=True) for k, s in enumerate(orderings(prim, target, 6, 40))]
    else:
        cands = candidates(target, exclude_same=True)
    t0 = time.time()
    out: dict[str, Any] = {"formula": target.reduced_formula, "known": str(known.entry_id), "method": method, "n_candidates": len(cands), "models": {}}
    pool = cands
    for i, model in enumerate(models):
        pd, _ = uhull(els, model)
        db = _cached_relax(known.structure, str(known.entry_id), model)
        e_db = float(pd.get_e_above_hull(umlip.mp_entry(db["structure"], db["energy"]), allow_negative=True))
        scored = []
        for label, s in pool:
            d = _cached_relax(s, _key(label, s), model)
            scored.append((float(pd.get_e_above_hull(umlip.mp_entry(d["structure"], d["energy"]), allow_negative=True)), label, s))
        scored.sort(key=lambda x: x[0])
        best = scored[0] if scored else (float("nan"), "-", None)
        out["models"][model] = {"ehull_db": round(e_db, 4), "ehull_best": round(best[0], 4), "gap": round(best[0] - e_db, 4), "best": best[1]}
        log(f"  {target.reduced_formula} {method} {model}: DB {e_db:+.3f} · 최선 후보 {best[0]:+.3f} (차이 {best[0] - e_db:+.3f}) — {best[1]}")
        if i == 0:
            pool = [(lb, s) for _, lb, s in scored[:final_top]]  # 다음 모델은 상위만
    out["seconds"] = round(time.time() - t0, 1)
    return out
