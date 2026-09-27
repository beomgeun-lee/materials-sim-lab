"""시험 카탈로그 구현 — kb/assays.yaml 의 코드 → 실행 함수.

v0 에서 실제로 계산하는 시험: S0, A1, A2, A3, A4, A7, A8, A10.
나머지는 무엇이 준비되면 실행되는지 이유와 함께 '대기(pending)'를 돌려준다.
"""

from __future__ import annotations

from msl.assays import alloy, blend, equilibrium, mp_assays, regulation, safety, supply
from msl.assays.base import AssayFn, Context, Outcome, pending


def _a5(ctx: Context) -> Outcome:
    return pending("3단계에서 붙음 — pymatgen Pourbaix + MP 이온 에너지", "pymatgen")


ASSAYS: dict[str, AssayFn] = {
    "S0": safety.run,
    "A1": mp_assays.a1,
    "A2": mp_assays.a2,
    "A3": equilibrium.a3,
    "A4": equilibrium.a4,
    "A5": _a5,
    "A6": alloy.a6,
    "A7": mp_assays.a7,
    "A8": blend.a8,
    "A9": supply.a9,
    "A10": regulation.a10,
}
