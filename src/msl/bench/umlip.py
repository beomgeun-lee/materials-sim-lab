"""uMLIP 처리량·정확도 실측 — 계획서 4단계 첫 항목 (노트북 CPU, D5).

대표 구조 20개 (금속·산화물·탄산염·규산염·할로겐화물·페로브스카이트, 1~60 원자)를 MP 캐시에서 골라
모델마다 불러오기·한 번 계산·이완(흐트러진 구조에서) 시간과 MP 원시 DFT 에너지와의 차이를 잰다.
에너지 차이는 GGA / GGA+U(LASPH 켬) / GGA+U(LASPH 끔)로 나눠 본다 — MP v2026 은 GGA+U 를 LASPH 켠 설정으로
다시 계산해, MPtrj(LASPH 끔)로 학습한 uMLIP 과 전이금속 산화물에서 어긋난다 (D19).
"""

from __future__ import annotations

import datetime as dt
import os
import statistics
import time
from typing import Any, Callable

import numpy as np

REPRESENTATIVE = [  # (MP 캐시 화학계, 화학식)
    ("Cu-Ni", "Cu"), ("Al-O", "Al"), ("Al-O-Si", "Si"), ("Al-Li-O", "Li"), ("Au-O-Si", "Au"),
    ("Mg-O-Ti", "MgO"), ("Al-O", "Al2O3"), ("Al-O-Si", "SiO2"), ("Mg-O-Ti", "TiO2"), ("Cl-Fe-Na", "NaCl"),
    ("Co-Li-O", "LiCoO2"), ("Li-Mn-O", "LiMn2O4"), ("Fe-O-Zn", "Fe2O3"), ("Fe-O-Zn", "ZnFe2O4"), ("Li-Ni-O", "NiO"),
    ("C-Ca-O", "CaCO3"), ("C-Ca-O-Si", "CaSiO3"), ("Mg-O-Si", "Mg2SiO4"), ("Ba-C-O-Ti", "BaTiO3"), ("Al-K-O", "KAlO2"),
]


def structures() -> list[dict[str, Any]]:
    """대표 구조 — 각 화학식의 MP GGA_GGA+U 바닥 엔트리."""
    from pymatgen.core import Composition

    from msl.engines import mp

    out = []
    for chemsys, formula in REPRESENTATIVE:
        ents = mp.entries_in_chemsys(chemsys.split("-"))
        target = Composition(formula).reduced_composition  # 문자열로 비교하면 ZnFe2O4 ≠ Zn(FeO2)2
        e = min((x for x in ents if x.composition.reduced_composition.almost_equals(target)), key=lambda x: x.energy_per_atom)
        kind = e.parameters.get("run_type", "?") + (" LASPH" if e.data.get("aspherical") else "")
        out.append({"formula": formula, "entry": e, "kind": kind})
    return out


def _perturb(structure, seed: int = 0, sigma: float = 0.03, strain: float = 0.03):
    s = structure.copy()
    s.apply_strain(strain / 3)  # 부피 약 3% 늘림
    rng = np.random.default_rng(seed)
    for i in range(len(s)):
        s.translate_sites([i], rng.normal(0, sigma, 3), frac_coords=False)
    return s


def bench(models: list[str], log: Callable[[str], None] = print) -> dict[str, Any]:
    import torch

    from msl.engines import umlip

    items = structures()
    res: dict[str, Any] = {"created": dt.datetime.now().isoformat(timespec="seconds"), "device": "cpu",
                           "threads": torch.get_num_threads(), "cpu_count": os.cpu_count(), "torch": torch.__version__,
                           "structures": [{"formula": x["formula"], "material_id": str(x["entry"].entry_id), "sites": len(x["entry"].structure),
                                           "kind": x["kind"]} for x in items], "models": {}}
    for name in models:
        t0 = time.perf_counter()
        umlip.calculator(name)
        load = time.perf_counter() - t0
        umlip.single_point(items[0]["entry"].structure, name)  # 첫 호출 준비 시간은 빼고 잰다
        rows = []
        for x in items:
            e = x["entry"]
            n = len(e.structure)
            en, sp = umlip.single_point(e.structure, name)
            r = umlip.relax(_perturb(e.structure), name, fmax=0.05, steps=400)
            rows.append({"formula": x["formula"], "sites": n, "kind": x["kind"], "sp_s": sp, "sp_ms_per_atom": sp * 1000 / n,
                         "relax_s": r.seconds, "steps": r.steps, "converged": r.converged,
                         "de_raw": en / n - e.uncorrected_energy_per_atom,  # MP 구조에서 한 번 계산한 에너지 − MP 원시 DFT (eV/atom)
                         "de_relaxed": r.energy / n - e.uncorrected_energy_per_atom})
            log(f"  {name:10} {x['formula']:8} {n:>3}원자 · 한 번 {sp * 1000:6.0f} ms · 이완 {r.seconds:5.1f} s {r.steps:>3}단계"
                f"{'' if r.converged else ' (미수렴)'} · ΔE {rows[-1]['de_raw']:+.3f} [{x['kind']}]")
        mae = lambda sel: statistics.mean(abs(r["de_raw"]) for r in rows if sel(r)) if any(sel(r) for r in rows) else None
        res["models"][name] = {
            "engine": umlip.engine_version(name), "load_s": load, "rows": rows,
            "summary": {
                "sp_ms_per_atom_median": statistics.median(r["sp_ms_per_atom"] for r in rows),
                "relax_s_median": statistics.median(r["relax_s"] for r in rows),
                "relax_s_total": sum(r["relax_s"] for r in rows),
                "steps_median": statistics.median(r["steps"] for r in rows),
                "converged": sum(r["converged"] for r in rows), "n": len(rows),
                "relax_per_hour": 3600 / statistics.mean(r["relax_s"] for r in rows),
                "mae_gga": mae(lambda r: r["kind"].startswith("GGA ") or r["kind"] == "GGA"),
                "mae_ggau_lasph": mae(lambda r: r["kind"] == "GGA+U LASPH"),
                "mae_ggau_nolasph": mae(lambda r: r["kind"] == "GGA+U"),
            },
        }
    return res
