"""L2 안정성 확인 — 치환 구조 생성 → uMLIP 이완 → 자기일관 hull → 두 모델 앙상블 (계획서 4단계, D19).

1. 경쟁 상: 대상 원소계의 MP 상(조성마다 바닥 엔트리, MP hull 거리 ≤ 0.1 eV/atom)을 같은 uMLIP 로 다시 이완해 hull 을 만든다.
   현재 MP(v2026) GGA+U 에너지는 LASPH 를 켠 설정이라 MPtrj 로 학습한 uMLIP 과 전이금속 산화물에서 0.05~0.19 eV/atom
   어긋난다 → uMLIP 값을 MP hull 에 바로 대지 않는다. 모든 에너지에는 MP2020 보정을 똑같이 붙인다.
2. 후보 구조: 대상과 음이온:양이온 비가 같은 알려진 구조(모체)의 양이온 자리에 대상 양이온을 무작위로 배치한다
   (초격자 ≤ max_atoms, 배치 n 개). 각 배치를 이완해 가장 낮은 것을 고른다. DB 에 있는 조성은 그 구조를 그대로 이완한다.
3. 두 모델(MACE-MPA-0, ORB v3)로 따로 계산해 hull 거리의 평균과 차이를 낸다 — 차이가 크면 불확실.

한계: 모체는 같은 원소계 안의 구조만 쓴다(다른 화학의 프로토타입 탐색 없음). 배치는 무작위 표본이라 전역 최저는 보장하지 않는다.
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
    """경쟁 상 — 조성마다 MP 바닥 엔트리 (MP hull ≤ 0.05, 원자 ≤ 40). 원소 바닥 상태는 항상 넣는다."""
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
        if is_element or (pd.get_e_above_hull(e) <= REF_MAX_EHULL and len(e.structure) <= REF_MAX_SITES):
            out.append(e)
    return out


def uhull(elements: list[str], model: str, log: Callable[[str], None] = lambda s: None) -> tuple[PhaseDiagram, list[str]]:
    """자기일관 hull — 경쟁 상을 모두 같은 uMLIP 로 이완 (MP2020 보정)."""
    refs = reference_entries(elements)
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


def parents(target: Composition, limit: int = 3) -> list[tuple[Any, Structure]]:
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


def _key(label: str, s: Structure) -> str:
    """구조 내용으로 만든 캐시 키 (실행마다 같은 값)."""
    import hashlib

    blob = json.dumps([label, [round(x, 4) for x in s.lattice.parameters],
                       [[str(site.specie), *np.round(site.frac_coords, 4).tolist()] for site in s]])
    return "cand-" + hashlib.sha256(blob.encode()).hexdigest()[:20]


def evaluate(formula: str, models: tuple[str, ...] = ("mace-mpa-0", "orb-v3"), n_orderings: int = 6, max_atoms: int = 40,
             log: Callable[[str], None] = print) -> dict[str, Any]:
    """조성 하나의 L2 hull 거리 (두 모델). DB 에 있는 조성은 그 구조를, 없는 조성은 치환 구조를 쓴다."""
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
        known, mp_ehull, structs = None, None, []
        for p, prim in parents(target):
            for k, s in enumerate(orderings(prim, target, n_orderings, max_atoms)):
                structs.append((f"{p.composition.reduced_formula}({p.entry_id}) 치환 {k + 1}", s))
        if not structs:
            raise ValueError(f"{target.reduced_formula}: 음이온:양이온 비가 같은 모체 구조를 찾지 못함 (같은 원소계 안)")
    log(f"{target.reduced_formula}: 후보 구조 {len(structs)}개 · 모델 {', '.join(models)}")
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
    return {"formula": target.reduced_formula, "known": None if known is None else str(known.entry_id),
            "mp_ehull": None if mp_ehull is None else round(float(mp_ehull), 4), "n_structures": len(structs),
            "models": per_model, "ehull_mean": round(float(np.mean(vals)), 4), "ehull_spread": round(float(max(vals) - min(vals)), 4),
            "seconds": round(time.time() - t0, 1)}
