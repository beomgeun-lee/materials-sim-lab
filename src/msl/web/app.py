"""로컬 웹 화면 — 레시피를 고르고 [실행]하면 실제 엔진으로 계산한 리포트를 보여 준다."""

from __future__ import annotations

import html
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel

from msl.env import ROOT
from msl.registry.load import load_registry
from msl.report import payload
from msl.runtime.runner import run_recipe
from msl.schema.recipe import RecipeError, load_recipe

RECIPES = ROOT / "examples" / "recipes"
SPACES = ROOT / "examples" / "spaces"
STATIC = Path(__file__).parent / "static"

app = FastAPI(title="가상 조합 실험실", docs_url=None, redoc_url=None)
registry = load_registry()


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (STATIC / "index.html").read_text(encoding="utf-8")


@app.get("/api/recipes")
def recipes() -> list[dict]:
    out = []
    for path in sorted(RECIPES.glob("*.yaml")):
        try:
            r = load_recipe(path, known_assays=set(registry.assays))
        except RecipeError:
            continue
        out.append({
            "file": path.name, "id": r.id, "name": r.name, "mode": r.mode.value,
            "yaml": path.read_text(encoding="utf-8"),
            "components": [{"ref": str(c.ref), "amount": str(c.amount) if c.amount else None,
                            "state": c.state.value if c.state else None} for c in r.components],
            "conditions": {k: str(getattr(r.conditions, k)) for k in type(r.conditions).model_fields
                           if getattr(r.conditions, k) is not None},
            "assays": r.assays,
        })
    return out


class RunRequest(BaseModel):
    file: str
    no_cache: bool = False


@app.post("/api/run")
def run(req: RunRequest) -> dict:
    path = (RECIPES / req.file).resolve()
    if path.parent != RECIPES.resolve() or not path.exists():
        raise HTTPException(404, "레시피 파일이 없음")
    recipe = load_recipe(path, known_assays=set(registry.assays))
    report = run_recipe(recipe, use_cache=not req.no_cache, registry=registry)
    return payload(report, registry)


# ── 조합 공간 (3단계) ─────────────────────────────────────────────────────


@app.get("/spaces", response_class=HTMLResponse)
def spaces() -> str:
    from msl.recipe.space import expand, load_space
    from msl.report.space import describe_space, page

    items = []
    for path in sorted(SPACES.glob("*.yaml")):
        try:
            sp = load_space(path)
            n = len(expand(sp))
        except RecipeError as exc:
            items.append(f'<li><b>{html.escape(path.name)}</b> — 정의 오류: {html.escape(str(exc))}</li>')
            continue
        items.append(f'<li><a href="/space/{html.escape(path.stem)}"><b>{html.escape(sp.name)}</b></a> '
                     f'<span class="sub">· {html.escape(describe_space(sp))} · 변형 {n}개</span>'
                     f'<br><span class="sub">{html.escape(sp.description or "")}</span></li>')
    body = ('<p class="eyebrow"><a href="/">← 가상 조합 실험실</a></p><h1>조합 공간</h1>'
            '<p class="sub">레시피 틀 하나에서 혼합비·양·조합을 바꾼 레시피를 여러 개 만들어 한 번에 시험합니다. 처음 실행은 계산 시간이 걸리고, 이후는 캐시에서 읽습니다.</p>'
            f'<ul class="spaces">{"".join(items)}</ul>')
    return page("조합 공간", body)


@app.get("/space/{name}", response_class=HTMLResponse)
def space(name: str, no_cache: bool = False) -> str:
    from msl.recipe.space import load_space, run_space
    from msl.report.space import to_html

    path = (SPACES / f"{name}.yaml").resolve()
    if path.parent != SPACES.resolve() or not path.exists():
        raise HTTPException(404, "공간 파일이 없음")
    try:
        res = run_space(load_space(path), registry, use_cache=not no_cache)
    except RecipeError as exc:
        raise HTTPException(400, str(exc)) from exc
    return to_html(res, nav='<a href="/spaces">← 조합 공간 목록</a>')
