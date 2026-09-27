"""MCP 서버 — 다른 AI 도구(Claude Code·Claude Desktop 등)가 물질 검색·레시피 실행·리포트 조회를 부를 수 있게 한다 (계획서 5단계, D23).

stdio 로만 돈다 (`msl mcp`). 웹 작업대와 같은 함수를 부르므로 결과·검증 규칙이 같다.
파일 접근은 예제·내 레시피 폴더와 data/reports 로 한정한다 — 임의 경로를 읽거나 쓰지 않는다.
L2(uMLIP) 계산은 몇 분~몇십 분 걸려 도구로 열지 않고, 저장된 결과 조회만 연다.
"""

from __future__ import annotations

import datetime as dt
import json
import re
from functools import cache
from pathlib import Path
from typing import Any

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations

from msl.env import DATA_DIR
from msl.web import workbench as wb

REPORTS = DATA_DIR / "reports"
REPORT_ID = re.compile(r"^[0-9]{14}-[a-z0-9][a-z0-9\-]*$")
READ = ToolAnnotations(readOnlyHint=True, openWorldHint=False)
COMPUTE = ToolAnnotations(readOnlyHint=True, idempotentHint=True, openWorldHint=True)  # MP·PubChem 조회가 섞일 수 있음

INSTRUCTIONS = """materials-sim-lab — 물질 조합(레시피)을 가상 시험 11종(S0 안전 게이트 + A1~A10)으로 계산한다.
흐름: search_substance 로 성분 확인 → check_recipe 로 레시피 YAML 검증 → run_recipe 로 실행(리포트 id 반환) → get_report.
레시피 형식은 get_recipe 로 예제를 한 번 읽어 보면 된다. 결과마다 충실도(L0 DB 조회, L1 예측, L2 uMLIP, T 평형)와
출처·라이선스가 붙는다. S0 가 위험을 찾으면 다른 시험보다 먼저 알린다 — 사용자에게 반드시 전달할 것.
새 조성의 물성은 predict_properties(L1, 80% 구간), 목표로 조성 찾기는 recommend,
평가 비용이 큰 단계(L2·DFT·실험)에서 다음에 잴 조성은 suggest_next (측정값은 add_measurement)."""

mcp = MCPServer("materials-sim-lab", instructions=INSTRUCTIONS)


@cache
def _registry():
    from msl.registry.load import load_registry

    return load_registry()


def _recipe_files() -> list[tuple[Path, str]]:
    return wb.list_files(wb.EXAMPLE_RECIPES, "example") + wb.list_files(wb.USER_RECIPES, "user")


def _parse(recipe_yaml: str):
    data, errors = wb.parse_yaml(recipe_yaml)
    if data is None:
        return None, errors
    return wb.validate_recipe(data, _registry())


# ── 조회 ──────────────────────────────────────────────────────────────────


@mcp.tool(annotations=READ)
def list_assays() -> list[dict[str, Any]]:
    """가상 시험 목록 — 코드(S0, A1~A10), 이름, 묻는 질문, 적용 대상, 충실도."""
    return [{"code": a.code, "name": a.name, "question": a.question, "applies_to": a.applies_to,
             "fidelity": [f.value for f in a.fidelity]} for a in _registry().assays.values()]


@mcp.tool(annotations=COMPUTE)
def search_substance(text: str) -> dict[str, Any]:
    """성분 하나를 찾는다. 이름(국문 관용명 포함: 염산·가성소다)·CAS·화학식·광물명·원소 기호를 받는다.
    돌려주는 ref 를 레시피 components[].ref 에 그대로 쓰면 된다."""
    from msl.web.app import ResolveRequest, resolve_text

    return resolve_text(ResolveRequest(text=text))


@mcp.tool(annotations=READ)
def list_recipes() -> list[dict[str, str]]:
    """예제(example)와 사용자가 저장한(user) 레시피 파일 목록."""
    import yaml

    out = []
    for path, source in _recipe_files():
        try:
            d = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except Exception:
            continue
        out.append({"file": path.name, "source": source, "id": str(d.get("id")), "name": str(d.get("name"))})
    return out


@mcp.tool(annotations=READ)
def get_recipe(file: str) -> str:
    """레시피 파일의 YAML (list_recipes 의 file 이름)."""
    path = wb.find_file(_recipe_files(), file)
    if path is None:
        raise ToolError(f"레시피 파일이 없음: {file}")
    return path.read_text(encoding="utf-8")


@mcp.tool(annotations=READ)
def check_recipe(recipe_yaml: str) -> dict[str, Any]:
    """레시피 YAML 을 실행하지 않고 검증한다. ok 면 정리된 YAML 과 자동으로 고를 시험 목록을 준다."""
    recipe, errors = _parse(recipe_yaml)
    if recipe is None:
        return {"ok": False, "errors": errors}
    from msl.recipe import route
    from msl.resolve import resolve

    comps = [resolve(c) for c in recipe.components]
    return {"ok": True, "errors": [], "yaml": wb.dump_yaml(recipe),
            "components": [{"ref": str(c.ref), "name": c.name, "formula": c.formula, "error": c.error} for c in comps],
            "assays": [{"code": code, "reason": why} for code, why in route(recipe, comps)]}


# ── 실행 ──────────────────────────────────────────────────────────────────


def _save_report(body: dict[str, Any]) -> str:
    REPORTS.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d%H%M%S")
    rid = f"{stamp}-{re.sub(r'[^a-z0-9-]', '-', body['recipe']['id'].lower())}"
    (REPORTS / f"{rid}.json").write_text(json.dumps(body, ensure_ascii=False), encoding="utf-8")
    return rid


def _render(body: dict[str, Any], fmt: str) -> str:
    from msl.report import to_markdown

    if fmt == "json":
        slim = {k: v for k, v in body.items() if k not in ("licenses", "source_names", "recipe_yaml")}
        return json.dumps(slim, ensure_ascii=False, indent=1)
    return to_markdown(body)


@mcp.tool(annotations=COMPUTE)
def run_recipe(recipe_yaml: str | None = None, file: str | None = None, fmt: str = "markdown") -> str:
    """레시피를 실행한다 — 해석 → S0 안전 게이트 → 시험들. recipe_yaml(직접 쓴 YAML) 또는 file(list_recipes 이름) 중 하나.
    결과는 data/reports 에 저장되고 첫 줄에 report_id 가 나온다. fmt: markdown(기본) | json.
    보통 수 초, 처음 보는 화학계는 MP 조회로 수십 초 걸릴 수 있다."""
    from msl.report import payload
    from msl.runtime.runner import run_recipe as _run

    if (recipe_yaml is None) == (file is None):
        raise ToolError("recipe_yaml 과 file 중 하나만 주세요")
    if file is not None:
        recipe_yaml = get_recipe(file)
    recipe, errors = _parse(recipe_yaml)
    if recipe is None:
        raise ToolError("레시피 오류:\n" + "\n".join(errors))
    body = payload(_run(recipe, registry=_registry()), _registry(), recipe_yaml=wb.dump_yaml(recipe))
    rid = _save_report(body)
    return f"report_id: {rid}\n\n" + _render(body, fmt)


@mcp.tool(annotations=READ)
def list_reports(limit: int = 20) -> list[dict[str, Any]]:
    """MCP 로 실행해 저장한 리포트 목록 (최근 순)."""
    out = []
    for p in sorted(REPORTS.glob("*.json"), reverse=True)[: max(1, min(limit, 200))]:
        try:
            b = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        out.append({"report_id": p.stem, "recipe": b["recipe"]["id"], "name": b["recipe"]["name"],
                    "assays": [f"{x['assay']}:{x['status']}" for x in b.get("results", [])]})
    return out


@mcp.tool(annotations=READ)
def get_report(report_id: str, fmt: str = "markdown") -> str:
    """저장된 리포트. fmt: markdown(기본) | json."""
    if not REPORT_ID.match(report_id) or not (REPORTS / f"{report_id}.json").exists():
        raise ToolError(f"리포트가 없음: {report_id}")
    return _render(json.loads((REPORTS / f"{report_id}.json").read_text(encoding="utf-8")), fmt)


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, openWorldHint=False))
def save_recipe(recipe_yaml: str, overwrite: bool = False) -> dict[str, str]:
    """레시피를 recipes/{id}.yaml 로 저장한다 (웹 작업대의 '내 레시피'에 나온다). 예제와 같은 id 는 거부,
    같은 id 파일이 있으면 overwrite=True 일 때만 덮어쓴다."""
    import yaml

    recipe, errors = _parse(recipe_yaml)
    if recipe is None:
        raise ToolError("레시피 오류:\n" + "\n".join(errors))
    examples = {str((yaml.safe_load(p.read_text(encoding="utf-8")) or {}).get("id")) for p in wb.EXAMPLE_RECIPES.glob("*.yaml")}
    if recipe.id in examples:
        raise ToolError(f"예제와 같은 id: {recipe.id} — 다른 id 로 저장해 주세요")
    path, err = wb.save_text(wb.USER_RECIPES, recipe.id, wb.dump_yaml(recipe), overwrite)
    if err:
        raise ToolError(err)
    return {"file": path.name, "path": wb.shown(path)}


# ── 예측·추천 ─────────────────────────────────────────────────────────────


@mcp.tool(annotations=COMPUTE)
def predict_properties(formulas: list[str]) -> list[dict[str, Any]]:
    """L1 조성 대리모델 — 화학식(또는 이름·CAS)만으로 형성에너지·hull 거리·안정성·밴드갭·밀도 예측 (80% 구간).
    DB(MP)에 있는 조성은 DFT 값도 함께 준다. 한 번에 50개까지. 새 조성의 형성에너지 MAE 약 0.07 eV/atom."""
    from msl.ml.predict import predict
    from msl.web.app import _formula_of

    if len(formulas) > 50:
        raise ToolError("한 번에 50개까지")
    out = []
    for text in formulas:
        try:
            out.append({"input": text, **predict(_formula_of(text))})
        except Exception as exc:
            out.append({"input": text, "error": str(exc)})
    return out


@mcp.tool(annotations=COMPUTE)
def recommend(goal_yaml: str, top: int = 15) -> dict[str, Any]:
    """목표 기반 조성 추천. goal_yaml 예:
    id: goal-my-search / name: 이름 / required: [Li, Mn, O] / optional: [Ni, Co] / max_elements: 4 / max_atoms: 10 /
    constraints: [{prop: ehull, max: 0.05}] / objective: {prop: ehull, direction: min}.
    DB 조성은 DFT, 새 조성은 L1 예측으로 평가한다. 상위 top 개와 전체 실행 기록 파일 이름을 준다."""
    import yaml

    from msl.recommend import Goal, save_run
    from msl.recommend import recommend as _recommend

    from pydantic import ValidationError

    try:
        goal = Goal.model_validate(yaml.safe_load(goal_yaml))
        res = _recommend(goal)
    except ValidationError as exc:
        raise ToolError("목표 오류:\n" + "\n".join(wb.errors_of(exc))) from exc
    except ValueError as exc:
        raise ToolError(str(exc)) from exc
    file = save_run(res).name
    return {"file": file, "summary": {k: v for k, v in res.items() if k not in ("results", "scatter")},
            "results": res["results"][: max(1, min(top, 100))]}


@mcp.tool(annotations=READ)
def l2_result(formula: str) -> dict[str, Any] | None:
    """이미 계산해 둔 L2(uMLIP 두 모델) 안정성 결과 — hull 거리, 모델 차이, 포논 안정성, L3(DFT) 승격 판단.
    없으면 null. 새 계산은 오래 걸려 웹 /recommend·/predict 화면이나 `msl l2` 로 돌린다."""
    from msl import l2

    return l2.saved(formula)


CAMPAIGNS = [wb.ROOT / "examples" / "campaigns", wb.ROOT / "campaigns"]


def _campaign(file: str):
    from msl.suggest import load_campaign

    for folder in CAMPAIGNS:
        path = folder / file
        if re.fullmatch(r"[A-Za-z0-9._\-]+\.yaml", file) and path.exists():
            return load_campaign(path)
    raise ToolError(f"캠페인 파일이 없음: {file} (examples/campaigns · campaigns)")


@mcp.tool(annotations=READ)
def list_campaigns() -> list[dict[str, Any]]:
    """베이지안 최적화 캠페인 목록 — 파일 이름, 목표, 측정값 개수."""
    from msl.suggest import load_campaign, measurements

    out = []
    for folder in CAMPAIGNS:
        for p in sorted(folder.glob("*.yaml")):
            try:
                c = load_campaign(p)
            except Exception:
                continue
            out.append({"file": p.name, "id": c.id, "name": c.name, "target": c.target.text(), "n_measured": len(measurements(c))})
    return out


@mcp.tool(annotations=COMPUTE)
def suggest_next(campaign_file: str, batch: int | None = None) -> dict[str, Any]:
    """다음에 잴(L2·DFT·실험) 조성 제안 — BayBE 베이지안 최적화. 측정값이 없으면 고르게, 쌓이면 대리모델 평균±표준편차와 함께.
    notes 의 경고(측정 출처 섞임 등)는 사용자에게 전할 것. 선택 설치 bo 가 필요하다."""
    try:
        from msl.suggest import suggest
    except ImportError as exc:
        raise ToolError("BayBE 가 없음 — uv sync --extra bo") from exc
    try:
        return suggest(_campaign(campaign_file), batch)
    except ValueError as exc:
        raise ToolError(str(exc)) from exc


@mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False, openWorldHint=False))
def add_measurement(campaign_file: str, formula: str, value: float, source: str, note: str = "") -> dict[str, Any]:
    """캠페인에 측정값 하나를 더한다. source: L2 · DFT · 실험 · 기타. 사용자가 준 값만 넣을 것 — 추정값을 넣지 않는다."""
    from msl.suggest import add_measurement as _add

    if source not in ("L2", "DFT", "실험", "기타"):
        raise ToolError("source 는 L2 · DFT · 실험 · 기타 중 하나")
    try:
        return _add(_campaign(campaign_file), formula, value, source, note)  # type: ignore[arg-type]
    except ValueError as exc:
        raise ToolError(str(exc)) from exc


def main() -> None:
    mcp.run("stdio")
