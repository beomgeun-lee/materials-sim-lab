"""L2 물성 검증 — uMLIP 탄성률을 MP DFT 탄성 데이터와, 포논 안정성을 알려진 안정/불안정 물질과 비교 (4단계)."""

from __future__ import annotations

import datetime as dt
import pickle
import statistics
from typing import Any, Callable

from msl.engines import mp

STABLE = ["mp-134", "mp-30", "mp-149", "mp-1265", "mp-22862", "mp-1143", "mp-2657"]  # Al Cu Si MgO NaCl Al2O3 TiO2(루틸)
UNSTABLE = ["mp-2998"]  # 입방정 BaTiO3 (Pm-3m) — 0 K 에서 강유전 뒤틀림으로 허수 포논


def _mp_data(ids: list[str]) -> dict[str, dict[str, Any]]:
    path = mp.CACHE / "props_bench.pkl"
    if path.exists():
        cached = pickle.loads(path.read_bytes())
        if set(ids) <= set(cached):
            return cached
    out: dict[str, dict[str, Any]] = {}
    with mp._rester() as mpr:
        for d in mpr.materials.summary.search(material_ids=ids, fields=["material_id", "formula_pretty", "structure", "symmetry"]):
            out[str(d.material_id)] = {"formula": d.formula_pretty, "structure": d.structure, "sg": d.symmetry.symbol, "K": None, "G": None}
        for d in mpr.materials.elasticity.search(material_ids=ids, fields=["material_id", "bulk_modulus", "shear_modulus"]):
            k, g = d.bulk_modulus, d.shear_modulus
            out[str(d.material_id)].update(K=getattr(k, "vrh", None) if k else None, G=getattr(g, "vrh", None) if g else None)
    mp.CACHE.mkdir(parents=True, exist_ok=True)
    path.write_bytes(pickle.dumps(out))
    return out


def bench(models: list[str], log: Callable[[str], None] = print) -> dict[str, Any]:
    from msl import props

    data = _mp_data(STABLE + UNSTABLE)
    res: dict[str, Any] = {"created": dt.datetime.now().isoformat(timespec="seconds"), "models": {}}
    for m in models:
        rows = []
        for mid in STABLE + UNSTABLE:
            d = data[mid]
            row: dict[str, Any] = {"material_id": mid, "formula": d["formula"], "sg": d["sg"], "expect_stable": mid in STABLE,
                                   "K_mp": d["K"], "G_mp": d["G"]}
            if mid in STABLE:
                e = props.elastic(d["structure"], m)
                row.update(K=e["K_vrh"], G=e["G_vrh"], t_elastic=e["seconds"])
            p = props.phonon(d["structure"], m, relax=mid in STABLE)  # 입방정 BaTiO3 는 대칭을 유지한 채 (이완하면 셀만 바뀜)
            row.update(min_freq=p["min_freq_THz"], stable=p["dynamically_stable"], t_phonon=p["seconds"])
            rows.append(row)
            log(f"  {m:10} {d['formula']:6} {d['sg']:8} K {row.get('K', '—')} (MP {row['K_mp'] and round(row['K_mp'])}) · "
                f"G {row.get('G', '—')} (MP {row['G_mp'] and round(row['G_mp'])}) · 포논 최소 {row['min_freq']:+.2f} THz → "
                f"{'안정' if row['stable'] else '불안정'} (기대 {'안정' if row['expect_stable'] else '불안정'})")
        el = [r for r in rows if r.get("K") is not None and r["K_mp"] and r["K_mp"] > 0]
        res["models"][m] = {"rows": rows, "summary": {
            "K_mae": statistics.mean(abs(r["K"] - r["K_mp"]) for r in el),
            "G_mae": statistics.mean(abs(r["G"] - r["G_mp"]) for r in el if r["G_mp"] and r["G_mp"] > 0),  # MP Al G=−14 같은 오류값 제외
            "K_mape": statistics.mean(abs(r["K"] - r["K_mp"]) / r["K_mp"] for r in el),
            "phonon_correct": sum(r["stable"] == r["expect_stable"] for r in rows), "n": len(rows)}}
    return res
