"""Pourbaix 도표 — MP 고체 에너지 + 실험 이온 에너지 (Persson 2012 체계, A5).

고체: MP GGA/GGA+U 엔트리에 MP2020 + MaterialsProjectAqueousCompatibility 보정 (mp-api 가 적용).
이온: MPContribs ion_ref_data 스냅샷 (msl.connectors.mp_ion_ref_data). mp-api 가 MPContribs 클라이언트로
받는 부분만 로컬 스냅샷으로 바꾸고, 나머지 조립은 mp-api get_pourbaix_entries 를 그대로 쓴다.
"""

from __future__ import annotations

import json
import pickle
import warnings

from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixEntry

from msl.connectors import mp_ion_ref_data
from msl.engines import mp

DEFAULT_CONC = 1e-6  # mol/kg — 부식 판정 관례 (Pourbaix 도감: 10⁻⁶ M 이하로 녹으면 '부식 없음')
NERNST = 0.0591576  # V/pH, 25 °C (RT ln10 / F)
O2_LINE = 1.229  # V vs SHE, O₂/H₂O 표준 전위
T25 = 298.15
KJ = 96.485  # kJ/mol → eV
ENERGY_REFERENCE = ("혼합 (25 °C): 실험 고체 ΔGf(NASA·이온 기준 고체) 우선 + MP(MP2020·수계 보정) 반응에너지 앵커, "
                    "이온은 실험 ΔGf (MPContribs ion_ref_data)")
EXP_NASA, EXP_ION_REF, EXP_ION = "NASA", "실험 (이온 데이터 기준 고체)", "실험 이온"
MP_ANCHOR, MP_RAW, MP_ION = "MP-앵커", "MP", "실험 이온 + MP 고체 보정"
ELEMENT = "원소 (기준 0)"
EXPERIMENTAL = (EXP_NASA, EXP_ION_REF, EXP_ION, ELEMENT)  # 출처 이름 앞부분


def water_lines(pH: float) -> tuple[float, float]:
    """(H₂/H₂O 선, O₂/H₂O 선) 전위 V vs SHE — 1 atm, 25 °C."""
    return -NERNST * pH, O2_LINE - NERNST * pH


def entries(elements: list[str]) -> list[PourbaixEntry]:
    """화학계(H·O 제외 원소)의 Pourbaix 엔트리. 같은 계는 캐시에서 읽는다."""
    chemsys = "-".join(sorted(set(elements)))
    mp.CACHE.mkdir(parents=True, exist_ok=True)
    path = mp.CACHE / f"pourbaix_{chemsys}_{mp.THERMO_TYPE}_ion{mp_ion_ref_data.VERSION}.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    ions = mp_ion_ref_data.records()
    with mp._rester() as mpr, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mpr.get_ion_reference_data = lambda: ions  # MPContribs 클라이언트 대신 로컬 스냅샷
        out = mpr.get_pourbaix_entries(sorted(set(elements)))
    path.write_bytes(pickle.dumps(out))
    return out


def observed(elements: list[str]) -> set[int]:
    """실험 관측 근거(ICSD, MP theoretical=false)가 있는 고체 material_id (정수). 같은 계는 캐시에서 읽는다."""
    chemsys = "-".join(sorted(set(elements)))
    path = mp.CACHE / f"pourbaix_observed_{chemsys}_{mp.THERMO_TYPE}.json"
    if path.exists():
        return set(json.loads(path.read_text(encoding="utf-8")))
    ids = sorted({str(e.entry_id).split("-GGA")[0].split("-r2SCAN")[0] for e in entries(elements) if e.phase_type == "Solid"})
    docs = mp.summaries(ids, ["theoretical"])
    out = sorted(i for k, d in docs.items() if d.get("theoretical") is False and (i := mp.mp_int(k)) is not None)
    path.write_text(json.dumps(out), encoding="utf-8")
    return set(out)


def _exp_solids(elements: list[str]) -> dict[str, tuple[float, str]]:
    """실험 고체 생성 자유에너지 (25 °C): 약분 화학식 → (eV/atom, 출처). NASA 298.15 K > 이온 기준 고체."""
    from pymatgen.core import Composition

    from msl.engines import thermo_hybrid as th

    els = set(elements) | {"O", "H"}
    out: dict[str, tuple[float, str]] = {}
    for r in mp_ion_ref_data.records():
        d = r["data"]
        c = Composition(d["RefSolid"])
        if len(c.elements) > 1 and {str(e) for e in c.elements} <= els:
            assert d["ΔGᶠRefSolid"]["unit"] == "kJ/mol", d
            out[c.reduced_formula] = (d["ΔGᶠRefSolid"]["value"] / KJ / c.num_atoms, EXP_ION_REF)
    g, _ = th.nasa_gf(sorted(els), T25)
    for f, v in g.items():
        if v.condensed and f != "H2O":
            out[Composition(f).reduced_formula] = (v.dgf_per_atom, f"{EXP_NASA} {v.species}")
    return out


def hybrid_entries(elements: list[str], observed_only: bool = True) -> tuple[list[PourbaixEntry], dict[str, str]]:
    """실험값 우선 Pourbaix 엔트리와 엔트리별 출처 (A2 v1 과 같은 앵커 방식, D15).

    - 고체: NASA 298.15 K ΔGf > 이온 데이터의 기준 고체 ΔGf (둘 다 실험). 같은 조성의 MP 다형은 뺀다.
    - MP 에만 있는 고체 X: 실험값이 있는 조성(앵커)으로의 분해 반응에너지만 MP 에서 가져온다.
        G(X) = G_MP(X) + Σ_p c_p [G_exp(p) − G_MP(p)]   (c_p: 앵커 상태도에서 X 의 분해 원자 분율)
    - 이온: 실험 ΔGf 그대로. 기준 고체도 실험값이라 Persson 2012 의 고체 보정 항이 0 이 된다.
    - observed_only: MP 에만 있는 고체는 실험 관측 구조(ICSD)가 있는 것만 둔다. Pourbaix 도감처럼 실재하는 상으로
      도표를 만든다 (예: GGA 가 과안정화하는 가상의 Cu₂O₃ 가 산소 선 아래에 나타나는 것을 막는다).
    """
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
    from pymatgen.analysis.pourbaix_diagram import MU_H2O, IonEntry
    from pymatgen.core import Composition
    from pymatgen.core.ion import Ion
    from pymatgen.entries.computed_entries import ComputedEntry

    raw = entries(elements)
    exp = _exp_solids(elements)
    seen = observed(elements) if observed_only else None
    solids = [e for e in raw if e.phase_type == "Solid"]
    best: dict[str, PourbaixEntry] = {}
    for e in solids:
        f = e.entry.composition.reduced_formula
        if f not in best or e.entry.energy_per_atom < best[f].entry.energy_per_atom:
            best[f] = e

    base = sorted(set(elements) | {"O", "H"})
    anchors = [PDEntry(Composition(el), 0.0, name=el) for el in base] + [PDEntry(Composition("H2O"), MU_H2O, name="H2O")]
    anchors += [PDEntry(Composition(f), best[f].entry.energy_per_atom * Composition(f).num_atoms, name=f)
                for f in exp if f in best]
    anchor_pd = PhaseDiagram(anchors)

    out: list[PourbaixEntry] = []
    source: dict[str, str] = {}

    def add(comp: Composition, epa: float, eid: str, src: str) -> None:
        ce = ComputedEntry(comp, epa * comp.num_atoms, entry_id=eid)
        pe = PourbaixEntry(ce, entry_id=eid)
        out.append(pe)
        source[pe.name] = src

    for e in solids:
        f = e.entry.composition.reduced_formula
        if len(e.entry.composition.elements) == 1:  # 원소: 생성 에너지 0 이 정의 (바닥 다형만)
            if e is best[f]:
                add(e.entry.composition, 0.0, e.entry_id, ELEMENT)
            continue
        if f in exp:
            if e is best[f]:
                add(e.entry.composition, exp[f][0], e.entry_id, exp[f][1])
            continue
        if seen is not None and mp.mp_int(e.entry_id) not in seen:
            continue  # 실험 관측 근거가 없는 이론 상
        decomp = anchor_pd.get_decomposition(e.entry.composition)
        shift = sum(c * (exp[p.name][0] - p.energy_per_atom) for p, c in decomp.items() if p.name in exp)
        add(e.entry.composition, e.entry.energy_per_atom + shift, e.entry_id, MP_ANCHOR if shift else MP_RAW)
    for f, (epa, src) in exp.items():
        if f not in best:
            add(Composition(f), epa, f"exp:{f}", src)

    allowed = set(elements) | {"O", "H"}
    for n, r in enumerate(mp_ion_ref_data.records()):
        ion = Ion.from_formula(r["formula"])
        if r["data"]["MajElements"] in elements and {str(el) for el in ion.elements} <= allowed:
            pe = PourbaixEntry(IonEntry(ion, r["data"]["ΔGᶠ"]["value"] / KJ), f"ion-{n}")
            out.append(pe)
            source[pe.name] = EXP_ION
    return out, source


def diagram(elements: list[str], comp: dict[str, float] | None = None, conc: float = DEFAULT_CONC,
            method: str = "hybrid") -> tuple[list[PourbaixEntry], PourbaixDiagram, dict[str, str]]:
    """Pourbaix 도표와 엔트리별 에너지 출처. 여러 원소면 comp(원소 → 몰 비)로 조성을 고정한다.

    method: hybrid(실험값 우선, 기본) / mp(MP 원본 — 비교용)
    """
    if method == "hybrid":
        ents, source = hybrid_entries(elements)
    else:
        ents = entries(elements)
        source = {e.name: MP_ION if e.phase_type == "Ion" else MP_RAW for e in ents}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pd = PourbaixDiagram(ents, comp_dict=comp if len(elements) > 1 else None,
                             conc_dict={el: conc for el in elements}, filter_solids=True)
    return ents, pd, source


def domains(pd: PourbaixDiagram, limits: list[list[float]] | None = None, digits: int = 3) -> list[dict]:
    """안정 영역 다각형 [{name, parts: [{name, solid}], vertices: [[pH, V], ...]}] (limits: [[pH0, pH1], [V0, V1]])."""
    _, verts = PourbaixDiagram.get_pourbaix_domains(pd.stable_entries, limits=limits)
    out = []
    for entry, vs in verts.items():
        parts = getattr(entry, "entry_list", [entry])
        out.append({"name": entry.name, "parts": [{"name": e.name, "solid": e.phase_type == "Solid"} for e in parts],
                    "vertices": [[round(float(x), digits), round(float(y), digits)] for x, y in vs]})
    return out
