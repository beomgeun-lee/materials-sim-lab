"""제안(BayBE) 되돌아보기 시험 — MP 에 이미 있는 조성을 '아직 안 잰 후보'로 두고 안정 조성을 얼마나 빨리 찾는지 (D25).

후보: L2 F1 세트와 같은 화학계의 3원계 조성(MP 바닥 구조, hull 거리 ≤ 0.3). 정답 = MP hull 거리.
사전값은 넣지 않는다 (L1 은 이 MP 데이터로 학습해 정답이 샌다) — 조성 특성만으로 배우는지 본다.
매 회 batch 개를 '재고'(정답을 알려 주고) 다음을 고른다. 무작위 선택과 비교: 회차별로 찾은 hull 위 조성 비율.
"""

from __future__ import annotations

import random
from typing import Any, Callable

import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram

from msl.bench.l2_f1 import CHEMSYS


def pool(max_ehull: float = 0.3) -> list[dict[str, Any]]:
    from msl.engines import mp

    out = {}
    for cs in CHEMSYS:
        ents = mp.entries_in_chemsys(cs.split("-"))
        pd = PhaseDiagram(ents)
        for e in ents:
            f = e.composition.reduced_formula
            if len(e.composition.elements) < 3:
                continue
            eh = float(pd.get_e_above_hull(e))
            if eh <= max_ehull and (f not in out or eh < out[f]):
                out[f] = eh
    return [{"formula": f, "truth": round(v, 4)} for f, v in sorted(out.items())]


def run(rounds: int = 8, batch: int = 5, seeds: int = 3, log: Callable[[str], None] = print) -> dict[str, Any]:
    from msl.suggest import CampaignSpec, suggest

    items = pool()
    truth = {x["formula"]: x["truth"] for x in items}
    stable = {f for f, v in truth.items() if v <= 1e-6}
    cands = [{"formula": f, "source": "MP", "prior": None} for f in truth]
    log(f"후보 {len(items)}개 · hull 위(안정) {len(stable)}개 · {rounds}회 × {batch}개 · 시드 {seeds}")
    curves: dict[str, list[list[float]]] = {"bo": [], "random": []}
    for seed in range(seeds):
        spec = CampaignSpec(id="camp-bench", name="bench", formulas=["H2O"], batch=batch, seed=seed)
        meas: list[dict[str, Any]] = []
        found = []
        for r in range(rounds):
            res = suggest(spec, candidates=cands, meas=meas)
            meas += [{"formula": s["formula"], "value": truth[s["formula"]], "source": "MP"} for s in res["suggestions"]]
            found.append(len({m["formula"] for m in meas} & stable) / len(stable))
        curves["bo"].append(found)
        rng = random.Random(seed)
        order = rng.sample(sorted(truth), len(truth))
        curves["random"].append([len(set(order[: batch * (r + 1)]) & stable) / len(stable) for r in range(rounds)])
        log(f"  시드 {seed}: BayBE {found[-1]:.0%} · 무작위 {curves['random'][-1][-1]:.0%} (잰 {batch * rounds}개 중 안정 조성 찾은 비율)")
    mean = {k: [round(float(x), 3) for x in np.mean(v, axis=0)] for k, v in curves.items()}
    return {"n_candidates": len(items), "n_stable": len(stable), "rounds": rounds, "batch": batch, "seeds": seeds,
            "found_fraction": mean, "expected_random_final": round(batch * rounds / len(items), 3)}
