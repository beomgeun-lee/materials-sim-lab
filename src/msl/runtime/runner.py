"""레시피 실행기 — 해석 → 라우팅 → 시험 → AssayResult (프로비넌스·캐시, 계획서 5장)."""

from __future__ import annotations

import datetime as dt
import json
import time
import traceback
from dataclasses import dataclass
from functools import cache
from importlib.metadata import version
from pathlib import Path
from typing import Any

from msl import __version__
from msl.assays import ASSAYS
from msl.assays.base import Context, Outcome, pending
from msl.env import CACHE_DIR, KB_DIR, load_dotenv
from msl.recipe import route
from msl.recipe.balance import Balance, balance
from msl.registry.load import Registry, load_registry
from msl.resolve import Resolved, resolve
from msl.schema.provenance import RunProvenance, canonical_hash
from msl.schema.recipe import Recipe
from msl.schema.result import AssayResult, SourceCitation, Status

RESULT_CACHE = CACHE_DIR / "results"
ENGINE_PACKAGES = ["pymatgen", "mp-api", "cantera", "phreeqpython"]


@dataclass
class Report:
    recipe: Recipe
    components: list[Resolved]
    routed: list[tuple[str, str]]
    results: list[AssayResult]
    elapsed: float
    balance: Balance | None = None

    def to_json(self, registry: Registry) -> dict[str, Any]:
        return {
            "recipe": self.recipe.model_dump(mode="json"),
            "components": [
                {"ref": str(c.ref), "name": c.name, "formula": c.formula, "cid": c.cid, "groups": c.groups,
                 "via": c.via, "error": c.error, "amount": str(c.amount) if c.amount else None,
                 "state": c.state.value if c.state else None}
                for c in self.components
            ],
            "routed": [{"code": code, "reason": why, "name": registry.assays[code].name if code in registry.assays else code}
                       for code, why in self.routed],
            "results": [r.model_dump(mode="json") for r in self.results],
            "balance": self.balance.to_json() if self.balance else None,
            "elapsed": round(self.elapsed, 2),
        }


def _versions() -> dict[str, str]:
    return {p: version(p) for p in ENGINE_PACKAGES}


@cache
def _code_hash() -> str:
    """시험·엔진·규칙 코드가 바뀌면 캐시가 무효가 되도록 소스와 kb 규칙 파일을 해시한다."""
    pkg = Path(__file__).resolve().parents[1]
    files = sorted([*pkg.glob("assays/*.py"), *pkg.glob("engines/*.py"), *pkg.glob("resolve/*.py"),
                    *pkg.glob("recipe/*.py"), *pkg.glob("ml/*.py"), *KB_DIR.glob("*.yaml")])
    return canonical_hash({str(f.relative_to(f.parents[1])): f.read_text(encoding="utf-8") for f in files})  # 폴더/파일 (같은 이름 __init__.py 구분)


def _data_fingerprint() -> str:
    """적재 데이터(테이블·소스·버전·행 수)가 바뀌면 캐시가 무효가 되도록 한다. 매 실행 계산 (가볍다)."""
    from msl.db import status

    from msl.ml.train import MODEL_DIR

    models = sorted((p.name, p.read_text(encoding="utf-8")) for p in MODEL_DIR.glob("*.json")) if MODEL_DIR.exists() else []
    return canonical_hash([[(s.table, s.source, s.version, s.rows) for s in status()], models])  # 학습 모델(L1)이 바뀌어도 무효


def _cite(reg: Registry, sid: str, ver: str | None) -> SourceCitation:
    src = reg.sources.get(sid)
    return SourceCitation(source=sid, version=ver or (src.version if src else None),
                          license=src.license if src else "unknown")


def _to_result(code: str, recipe: Recipe, out: Outcome, reg: Registry, h: str, started: dt.datetime,
               params: dict[str, Any]) -> AssayResult:
    return AssayResult(
        assay=code, recipe_id=recipe.id, status=out.status, fidelity=out.fidelity, engine=out.engine,
        engine_version=out.engine_version, model_weights=out.model_weights, energy_reference=out.energy_reference,
        conditions_basis=out.conditions_basis, values=out.values,
        sources=[_cite(reg, s, v) for s, v in out.sources], caveats=out.caveats, warnings=out.warnings,
        data={**out.data, "summary": out.summary},
        provenance=RunProvenance(input_hash=h, msl_version=__version__, engine_versions=_versions(),
                                 parameters=params, started_at=started, finished_at=dt.datetime.now(dt.UTC)),
    )


def run_recipe(recipe: Recipe, use_cache: bool = True, registry: Registry | None = None) -> Report:
    load_dotenv()
    reg = registry or load_registry()
    t0 = time.monotonic()
    comps = [resolve(c) for c in recipe.components]
    routed = route(recipe, comps)
    bal = balance(comps)
    ctx = Context(recipe=recipe, comps=comps, registry=reg, shared={"balance": bal})
    results: list[AssayResult] = []
    data_fp = _data_fingerprint()
    for code, why in routed:
        params = {"route_reason": why}
        h = canonical_hash({"recipe": recipe.model_dump(mode="json"), "assay": code, "versions": _versions(),
                            "msl": __version__, "code": _code_hash(), "data": data_fp, "components": [(str(c.ref), c.formula, c.groups) for c in comps]})
        cached = RESULT_CACHE / f"{h}.json"
        if use_cache and cached.exists():
            results.append(AssayResult.model_validate_json(cached.read_text(encoding="utf-8")))
            continue
        started = dt.datetime.now(dt.UTC)
        fn = ASSAYS.get(code)
        try:
            out = fn(ctx) if fn else pending(f"{code}: 구현 없음", "—")
        except Exception as exc:  # 실패도 결과로 남긴다 — 리포트에서 원인을 볼 수 있게
            out = pending(f"실행 실패: {type(exc).__name__}: {exc}", "—")
            out.status = Status.FAILED
            out.warnings = [traceback.format_exc(limit=3)]
        res = _to_result(code, recipe, out, reg, h, started, params)
        if res.status in (Status.OK, Status.WARNING, Status.NOT_APPLICABLE):
            RESULT_CACHE.mkdir(parents=True, exist_ok=True)
            cached.write_text(res.model_dump_json(), encoding="utf-8")
        results.append(res)
    return Report(recipe=recipe, components=comps, routed=routed, results=results,
                  elapsed=time.monotonic() - t0, balance=bal)


def save_report(report: Report, registry: Registry, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report.to_json(registry), ensure_ascii=False, indent=1), encoding="utf-8")
