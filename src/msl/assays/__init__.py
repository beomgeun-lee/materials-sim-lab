"""시험 카탈로그 구현 — kb/assays.yaml 의 코드 → 실행 함수.

S0 과 A1~A10 모두 실제로 계산한다 (A5 는 3단계에서 붙음).
준비되지 않은 시험은 무엇이 있으면 실행되는지 이유와 함께 '대기(pending)'를 돌려준다.
"""

from __future__ import annotations

from msl.assays import alloy, aqueous, blend, equilibrium, mp_assays, regulation, safety, supply
from msl.assays.base import AssayFn


ASSAYS: dict[str, AssayFn] = {
    "S0": safety.run,
    "A1": mp_assays.a1,
    "A2": mp_assays.a2,
    "A3": equilibrium.a3,
    "A4": equilibrium.a4,
    "A5": aqueous.a5,
    "A6": alloy.a6,
    "A7": mp_assays.a7,
    "A8": blend.a8,
    "A9": supply.a9,
    "A10": regulation.a10,
}
