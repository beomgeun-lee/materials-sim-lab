"""로컬 웹 화면 — 레시피를 고르고 [실행]하면 실제 엔진으로 계산한 리포트를 보여 준다."""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse
from pydantic import BaseModel

from msl.env import ROOT
from msl.registry.load import load_registry
from msl.runtime.runner import run_recipe
from msl.schema.recipe import RecipeError, load_recipe

RECIPES = ROOT / "examples" / "recipes"
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
    body = report.to_json(registry)
    body["licenses"] = {k: {"name": v.name, "partition": v.partition.value} for k, v in registry.licenses.items()}
    body["source_names"] = {k: v.name for k, v in registry.sources.items()}
    body["assay_names"] = {k: v.name for k, v in registry.assays.items()}
    return body
