"""산화수별 전이금속 에너지 보정 — 화학계 하나씩 빼고 맞춰(LOCO) 3원계 분해 에너지 오차가 줄어드는지."""
import warnings, statistics, json
warnings.filterwarnings("ignore")
import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram
from msl import l2
from msl.engines import mp, umlip

MODELS = ("orb-v3", "mace-mpa-0")
TM = ["Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu"]
Q = list(range(1, 8))

# 화학계별 경쟁 상 (캐시가 모두 있는 산화물 계)
data = {}
for p in sorted(mp.CACHE.glob("entries_*_GGA_GGA+U.pkl")):
    cs = p.name.removeprefix("entries_").split("_GGA")[0]
    els = cs.split("-")
    if "O" not in els or len(els) < 3:
        continue
    refs = l2.reference_entries(els)
    if any(not (l2.CACHE / m / f"{e.entry_id}.json").exists() for e in refs for m in MODELS):
        continue
    u = {m: {} for m in MODELS}
    for m in MODELS:
        for e in refs:
            d = l2._cached_relax(e.structure, str(e.entry_id), m)
            u[m][str(e.entry_id)] = umlip.mp_entry(d["structure"], d["energy"], entry_id=str(e.entry_id))
    data[cs] = (refs, u)
print("화학계", len(data))

def feats(comp):
    """(비 TM 원소 원자 수, TM 원소 × 산화수 모자 함수) — TM 의 평균 형식 산화수는 l2.charges 로."""
    q = l2.charges(comp) if len(comp.elements) > 1 else {str(comp.elements[0]): 0.0}
    f = {}
    for el, n in comp.get_el_amt_dict().items():
        if el in TM and q.get(el, 0) > 0:
            x = min(max(q[el], Q[0]), Q[-1])
            k = int(np.floor(x)); w = x - k
            f[(el, k)] = f.get((el, k), 0) + n * (1 - w)
            if w > 1e-9:
                f[(el, k + 1)] = f.get((el, k + 1), 0) + n * w
        else:
            f[("mu", el)] = f.get(("mu", el), 0) + n
    return f

def fit(css, m, lam=0.05):
    rows, y = [], []
    seen = set()
    for cs in css:
        refs, u = data[cs]
        for e in refs:
            if str(e.entry_id) in seen or u[m][str(e.entry_id)] is None:
                continue
            seen.add(str(e.entry_id))
            rows.append(feats(e.composition)); y.append(e.energy - u[m][str(e.entry_id)].energy)  # MP − uMLIP (총에너지)
    keys = sorted({k for r in rows for k in r}, key=str)
    X = np.array([[r.get(k, 0.0) for k in keys] for r in rows]); y = np.array(y)
    reg = np.array([0.0 if k[0] == "mu" else lam for k in keys])  # TM 산화수 항만 규제
    beta = np.linalg.solve(X.T @ X + np.diag(reg * len(y)), X.T @ y)
    return dict(zip(keys, beta))

def corrected(entry, beta):
    from pymatgen.entries.computed_entries import ComputedEntry
    f = feats(entry.composition)
    return ComputedEntry(entry.composition, entry.energy + sum(beta.get(k, 0.0) * v for k, v in f.items()), entry_id=entry.entry_id)

def dec(entries, target_comp, target):
    others = [x for x in entries if not x.composition.reduced_composition.almost_equals(target_comp)]
    return PhaseDiagram(others).get_e_above_hull(target, allow_negative=True)


def fit_local(refs, u, m, exclude, lam):
    rows, y = [], []
    for e in refs:
        if e.composition.reduced_composition.almost_equals(exclude) or u[m][str(e.entry_id)] is None:
            continue
        rows.append(feats(e.composition)); y.append(e.energy - u[m][str(e.entry_id)].energy)
    keys = sorted({k for r in rows for k in r}, key=str)
    X = np.array([[r.get(k, 0.0) for k in keys] for r in rows]); y = np.array(y)
    reg = np.array([0.0 if k[0] == "mu" else lam for k in keys])
    beta = np.linalg.lstsq(X.T @ X + np.diag(reg * len(y)) + 1e-9 * np.eye(len(keys)), X.T @ y, rcond=None)[0]
    return dict(zip(keys, beta))

for lam in (0.0, 0.01, 0.05, 0.2):
    res = {"raw": [], "corr": []}; per = []
    for cs, (refs, u) in data.items():
        for e in refs:
            c = e.composition.reduced_composition
            if len(c.elements) < 3:
                continue
            dmp = dec(refs, c, e)
            draw, dcor = [], []
            for m in MODELS:
                b = fit_local(refs, u, m, c, lam)
                ue = [x for x in u[m].values() if x is not None]
                t = u[m][str(e.entry_id)]
                draw.append(dec(ue, c, t))
                dcor.append(dec([corrected(x, b) for x in ue], c, corrected(t, b)))
            res["raw"].append(np.mean(draw) - dmp); res["corr"].append(np.mean(dcor) - dmp)
            per.append((cs, c.reduced_formula, round(dmp, 3), round(float(np.mean(draw)), 3), round(float(np.mean(dcor)), 3)))
    print(f"λ={lam}: MAE {np.mean(np.abs(res['raw'])):.4f} → {np.mean(np.abs(res['corr'])):.4f} · 최대 {np.max(np.abs(res['raw'])):.3f} → {np.max(np.abs(res['corr'])):.3f} · 부호 뒤집힘 {sum((a<=0)!=(b<=0) for a,b in [(p[2],p[3]) for p in per])} → {sum((a<=0)!=(b<=0) for a,b in [(p[2],p[4]) for p in per])}")
    for cs in ("Co-Li-O", "Fe-Li-O", "Fe-O-Zn", "Li-Mn-Ni-O", "Li-Mn-O", "Li-Ni-O"):
        s2 = [(abs(p[3]-p[2]), abs(p[4]-p[2])) for p in per if p[0] == cs]
        print(f"    {cs:12} {np.mean([a for a,_ in s2]):.4f} → {np.mean([b for _,b in s2]):.4f}")
