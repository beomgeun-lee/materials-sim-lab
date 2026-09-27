"""L2 안정성 판정 F1 — 계획서 4단계 완료 기준 ("신규 조성 과제에서 MP 참조값 대비 F1 0.85 이상, 자체 세트", D24).

신규 조성 과제: MP 에 있는 3원계 이상 조성을 '모른다고 치고' 판정한다.
- 후보 구조는 대상 조성 자신과 같은 원소계 원형을 빼고 만든다 (candidates(exclude_same=True)).
- 경쟁 상 hull 에서도 대상 조성을 뺀다 → 값은 '분해 에너지' (나머지 상 대비, 음수면 hull 아래).
- 정답도 같은 방식: MP 엔트리로 대상 조성을 뺀 hull 대비 분해 에너지.

판정 기준 (결과를 보기 전에 정함):
- 주 기준: 정답·예측 모두 분해 에너지 ≤ 0 이면 '안정' (Matbench Discovery 관례). 이 F1 을 0.85 와 비교한다.
- 보조: 제품 규칙 — 예측 ≤ L3_WINDOW(0.05) 면 'DFT 로 확인할 후보' (정답은 그대로 ≤ 0). 재현율 위주.
- 구조 탐색 실패와 에너지 오차를 가르려고 'MP 구조를 준 경우'(에너지만)의 F1 도 같이 잰다.

세트: 정해진 화학계 목록에서 hull 위 조성(양성)과 hull 위 0~0.3 eV/atom 조성(음성, 0~0.05 / 0.05~0.3 절반씩)을
MP 구조 원자 수가 작은 순으로 고른다. 후보 구조를 만들 수 없는 조성은 건너뛰고 개수를 기록한다.
"""

from __future__ import annotations

import datetime as dt
import json
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Composition

CHEMSYS = ["Li-Mn-O", "Li-Co-O", "Li-Ni-O", "Fe-Li-O", "Mg-O-Si", "Fe-O-Zn", "Ba-O-Ti", "Al-Li-O", "Mg-O-Ti", "Ca-O-Si",
           "Al-K-O", "Al-Mg-O", "Na-O-Si", "Ca-O-Ti", "Li-O-Ti", "Al-Na-O", "Li-O-Si", "Ca-Fe-O"]
MAX_SITES = 24  # MP 바닥 구조 원자 수 상한 — 후보 구조 크기(40)와 계산 시간을 맞춘다
NEG_MAX = 0.3


def _decomp(pd_entries: list, target: Composition, entry) -> float:
    """대상 조성을 뺀 hull 대비 분해 에너지 (eV/atom)."""
    others = [e for e in pd_entries if not e.composition.reduced_composition.almost_equals(target.reduced_composition)]
    return float(PhaseDiagram(others).get_e_above_hull(entry, allow_negative=True))


def select(n_pos: int = 15, n_neg: int = 15, chemsys: list[str] = CHEMSYS) -> list[dict[str, Any]]:
    from msl.engines import mp

    pos, neg_near, neg_far = [], [], []
    for cs in chemsys:
        ents = mp.entries_in_chemsys(cs.split("-"))
        pd = PhaseDiagram(ents)
        best: dict[str, Any] = {}
        for e in ents:
            f = e.composition.reduced_formula
            if len(e.composition.elements) >= 3 and (f not in best or e.energy_per_atom < best[f].energy_per_atom):
                best[f] = e
        for f, e in best.items():
            if len(e.structure) > MAX_SITES:
                continue
            eh = float(pd.get_e_above_hull(e))
            row = {"formula": f, "chemsys": cs, "mp_id": str(e.entry_id), "mp_ehull": round(eh, 4), "sites": len(e.structure)}
            if eh <= 1e-6:
                pos.append(row)
            elif eh <= 0.05:
                neg_near.append(row)
            elif eh <= NEG_MAX:
                neg_far.append(row)

    def spread(rows: list[dict], k: int) -> list[dict]:
        """화학계가 고르게 섞이도록 — 화학계별로 원자 수 작은 순, 돌아가며 하나씩."""
        by: dict[str, list] = {}
        for r in sorted(rows, key=lambda r: (r["sites"], r["formula"])):
            by.setdefault(r["chemsys"], []).append(r)
        out = []
        while len(out) < k and any(by.values()):
            for cs in list(by):
                if by[cs] and len(out) < k:
                    out.append(by[cs].pop(0))
        return out

    return ([dict(r, label=True) for r in spread(pos, n_pos)] + [dict(r, label=False) for r in spread(neg_near, n_neg // 2)]
            + [dict(r, label=False) for r in spread(neg_far, n_neg - n_neg // 2)])


def run_one(item: dict[str, Any], models: tuple[str, ...] = ("orb-v3", "mace-mpa-0"), final_top: int = 3,
            log: Callable[[str], None] = print) -> dict[str, Any]:
    from msl import l2
    from msl.engines import mp, umlip

    target = Composition(item["formula"])
    els = sorted(str(e) for e in target.elements)
    ents = mp.entries_in_chemsys(els)
    known = next(e for e in ents if str(e.entry_id) == item["mp_id"])
    t0 = time.time()
    truth = _decomp(ents, target, known)
    cands = l2.candidates(target, exclude_same=True)
    out = dict(item, truth=round(truth, 4), n_candidates=len(cands), models={})
    if not cands:
        out["skipped"] = "후보 구조를 만들 수 없음"
        return out
    pool = cands
    for i, model in enumerate(models):
        pd, _ = l2.uhull(els, model, exclude=target)
        d = l2._cached_relax(known.structure, str(known.entry_id), model)
        e_given = float(pd.get_e_above_hull(umlip.mp_entry(d["structure"], d["energy"]), allow_negative=True))
        scored = []
        for label, s in pool:
            d = l2._cached_relax(s, l2._key(label, s), model)
            ce = umlip.mp_entry(d["structure"], d["energy"])
            if ce is not None:
                scored.append((float(pd.get_e_above_hull(ce, allow_negative=True)), label, s))
        scored.sort(key=lambda x: x[0])
        out["models"][model] = {"pred": round(scored[0][0], 4), "given": round(e_given, 4), "best": scored[0][1]}
        if i == 0:
            pool = [(lb, s) for _, lb, s in scored[:final_top]]
    out["pred"] = round(float(np.mean([m["pred"] for m in out["models"].values()])), 4)
    out["given"] = round(float(np.mean([m["given"] for m in out["models"].values()])), 4)
    out["seconds"] = round(time.time() - t0, 1)
    log(f"  {item['formula']:12} 정답 {truth:+.3f} ({'안정' if item['label'] else '불안정'}) · 예측 {out['pred']:+.3f} · "
        f"MP 구조 {out['given']:+.3f} · 후보 {len(cands)} · {out['seconds']:.0f}초")
    return out


def _f1(pairs: list[tuple[bool, bool]]) -> dict[str, Any]:
    tp = sum(t and p for t, p in pairs)
    fp = sum((not t) and p for t, p in pairs)
    fn = sum(t and (not p) for t, p in pairs)
    tn = sum((not t) and (not p) for t, p in pairs)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    return {"f1": round(2 * prec * rec / (prec + rec), 3) if prec + rec else 0.0, "precision": round(prec, 3), "recall": round(rec, 3),
            "accuracy": round((tp + tn) / len(pairs), 3) if pairs else None, "tp": tp, "fp": fp, "fn": fn, "tn": tn}


def score(rows: list[dict[str, Any]], window: float = 0.05) -> dict[str, Any]:
    done = [r for r in rows if "pred" in r]
    truth = [r["truth"] <= 0 for r in done]
    return {
        "n": len(done), "skipped": sum("skipped" in r for r in rows),
        "main": _f1([(t, r["pred"] <= 0) for t, r in zip(truth, done)]),
        "product_rule": _f1([(t, r["pred"] <= window) for t, r in zip(truth, done)]),
        "given_structure": _f1([(t, r["given"] <= 0) for t, r in zip(truth, done)]),
        "mae_pred": round(float(np.mean([abs(r["pred"] - r["truth"]) for r in done])), 4) if done else None,
        "mae_given": round(float(np.mean([abs(r["given"] - r["truth"]) for r in done])), 4) if done else None,
        "search_gap_median": round(float(np.median([r["pred"] - r["given"] for r in done])), 4) if done else None,
    }


def bench(save: Path, n_pos: int = 15, n_neg: int = 15, log: Callable[[str], None] = print) -> dict[str, Any]:
    """세트를 고르고 하나씩 계산해 save 에 이어 쓴다 (중단 후 다시 실행하면 끝난 조성은 건너뜀)."""
    items = select(n_pos, n_neg)
    prev = json.loads(save.read_text(encoding="utf-8")) if save.exists() else {}
    rows = {r["formula"]: r for r in prev.get("rows", [])}
    log(f"세트 {len(items)}개 (안정 {sum(i['label'] for i in items)} · 불안정 {sum(not i['label'] for i in items)}) · "
        f"이미 끝남 {sum('error' not in r for r in rows.values())}")
    for k, item in enumerate(items, 1):
        if "error" not in rows.get(item["formula"], {"error": 1}):
            continue  # 끝난 조성 (실패한 조성은 다시 계산)
        log(f"[{k}/{len(items)}] {item['formula']} ({item['chemsys']})")
        try:
            rows[item["formula"]] = run_one(item, log=log)
        except Exception as exc:
            rows[item["formula"]] = dict(item, error=str(exc)[:300])
            log(f"  실패: {exc}")
        res = {"created": dt.datetime.now().isoformat(timespec="seconds"), "criteria": "분해 에너지 ≤ 0 (정답·예측), 대상 조성 제외 hull",
               "rows": list(rows.values()), "score": score(list(rows.values()))}
        save.write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    return res
