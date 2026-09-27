"""로컬 웹 화면 — 작업대에서 레시피를 만들거나 고쳐 [실행]하면 실제 엔진으로 계산한 리포트를 보여 준다."""

from __future__ import annotations

import html
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, PlainTextResponse
from pydantic import BaseModel

from msl.registry.load import load_registry
from msl.report import payload
from msl.runtime.runner import run_recipe
from msl.schema.recipe import RecipeError, load_recipe
from msl.web import workbench as wb

STATIC = Path(__file__).parent / "static"
GASES = {"O2", "N2", "H2", "CH4", "CO2", "CO", "Ar", "He", "Ne", "NH3", "Cl2", "F2", "C2H6", "C3H8", "C2H2", "C2H4", "SO2", "NO2", "NO", "H2S"}

app = FastAPI(title="가상 조합 실험실", docs_url=None, redoc_url=None)
registry = load_registry()


def _recipe_files() -> list[tuple[Path, str]]:
    return wb.list_files(wb.EXAMPLE_RECIPES, "example") + wb.list_files(wb.USER_RECIPES, "user")


def _space_files() -> list[tuple[Path, str]]:
    return wb.list_files(wb.EXAMPLE_SPACES, "example") + wb.list_files(wb.USER_SPACES, "user")


def _example_ids(kind: str) -> set[str]:
    """예제 파일의 id (파일 이름과 id 가 다르다: cu-ni-alloy.yaml ↔ rcp-cu-ni-alloy)."""
    import yaml

    folder = wb.EXAMPLE_RECIPES if kind == "recipe" else wb.EXAMPLE_SPACES
    out = set()
    for p in folder.glob("*.yaml"):
        try:
            out.add(str(yaml.safe_load(p.read_text(encoding="utf-8")).get("id")))
        except Exception:
            continue
    return out


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (STATIC / "index.html").read_text(encoding="utf-8")


@app.get("/api/recipes")
def recipes() -> list[dict]:
    """예제·내 레시피 목록 — 작업대에 바로 올릴 수 있는 폼 형태로."""
    out = []
    for path, source in _recipe_files():
        try:
            r = load_recipe(path, known_assays=set(registry.assays))
        except RecipeError:
            continue
        out.append({"file": path.name, "source": source, "id": r.id, "name": r.name,
                    "yaml": path.read_text(encoding="utf-8"), "form": wb.to_form(r)})
    return out


class RunRequest(BaseModel):
    file: str
    no_cache: bool = False


@app.post("/api/run")
def run(req: RunRequest) -> dict:
    """파일 이름으로 실행 (예제·내 레시피)."""
    path = wb.find_file(_recipe_files(), req.file)
    if path is None:
        raise HTTPException(404, "레시피 파일이 없음")
    recipe = load_recipe(path, known_assays=set(registry.assays))
    report = run_recipe(recipe, use_cache=not req.no_cache, registry=registry)
    return payload(report, registry, recipe_yaml=path.read_text(encoding="utf-8"))


# ── 작업대: 성분 찾기 · 레시피 확인 · 실행 · 저장 ─────────────────────────


class ResolveRequest(BaseModel):
    text: str


@app.post("/api/resolve")
def resolve_text(req: ResolveRequest) -> dict[str, Any]:
    """자유 입력(이름·CAS·화학식·광물명·원소)을 참조로 바꿔 해석한다."""
    from msl.resolve import resolve
    from msl.schema.recipe import Component

    ref = wb.guess_ref(req.text)
    try:
        comp = Component.model_validate({"ref": ref})
    except Exception as exc:  # 참조 형식 오류 (예: 존재하지 않는 원소 기호)
        return {"ref": ref, "error": str(exc).splitlines()[-1]}
    r = resolve(comp)
    state = "solid"
    if r.formula:
        from pymatgen.core import Composition

        try:
            c = Composition(r.formula)
            state = ("gas" if any(c.almost_equals(Composition(g)) for g in GASES)
                     else "liquid" if c.almost_equals(Composition("H2O")) else "solid")
        except Exception:
            pass
    if comp.ref.namespace.value == "name":  # 국문 관용명 표에 기본 상태가 있으면 (예: 염산 → 수용액)
        from msl.resolve import korean_aliases

        state = (korean_aliases().get(comp.ref.key.replace(" ", "")) or {}).get("state", state)
    return {"ref": str(r.ref), "name": r.name, "formula": r.formula, "via": r.via, "cid": r.cid, "cas": r.cas,
            "groups": r.groups, "error": r.error, "state": state}


class RecipeRequest(BaseModel):
    recipe: dict[str, Any]
    no_cache: bool = False
    overwrite: bool = False


def _validated(form: dict[str, Any]):
    recipe, errors = wb.validate_recipe(form, registry)
    if recipe is None:
        raise HTTPException(422, {"errors": errors})
    return recipe


@app.post("/api/recipe/check")
def recipe_check(req: RecipeRequest) -> dict[str, Any]:
    """검증 + 성분 해석 + 자동으로 고를 시험 + 질량수지 미리보기 (계산은 하지 않음)."""
    from msl.recipe import route
    from msl.recipe.balance import balance
    from msl.resolve import resolve

    recipe, errors = wb.validate_recipe(req.recipe, registry)
    if recipe is None:
        return {"ok": False, "errors": errors}
    comps = [resolve(c) for c in recipe.components]
    return {
        "ok": True, "errors": [], "yaml": wb.dump_yaml(recipe), "form": wb.to_form(recipe),
        "components": [{"ref": str(c.ref), "name": c.name, "formula": c.formula, "via": c.via, "error": c.error,
                        "groups": c.groups} for c in comps],
        "route": [{"code": code, "reason": why, "name": registry.assays[code].name if code in registry.assays else code}
                  for code, why in route(recipe, comps)],
        "balance": balance(comps).to_json(),
    }


class YamlRequest(BaseModel):
    yaml: str
    no_cache: bool = False
    overwrite: bool = False


@app.post("/api/recipe/parse")
def recipe_parse(req: YamlRequest) -> dict[str, Any]:
    """YAML 편집 → 폼. 검증 오류는 줄 단위 메시지로."""
    data, errors = wb.parse_yaml(req.yaml)
    if data is None:
        return {"ok": False, "errors": errors}
    recipe, errors = wb.validate_recipe(data, registry)
    if recipe is None:
        return {"ok": False, "errors": errors}
    return {"ok": True, "errors": [], "form": wb.to_form(recipe), "yaml": wb.dump_yaml(recipe)}


@app.post("/api/recipe/run")
def recipe_run(req: RecipeRequest) -> dict[str, Any]:
    recipe = _validated(req.recipe)
    report = run_recipe(recipe, use_cache=not req.no_cache, registry=registry)
    return payload(report, registry, recipe_yaml=wb.dump_yaml(recipe))


@app.post("/api/recipe/save")
def recipe_save(req: RecipeRequest) -> dict[str, Any]:
    recipe = _validated(req.recipe)
    if recipe.id in _example_ids("recipe"):
        raise HTTPException(409, f"예제와 같은 id: {recipe.id} — 다른 id 로 저장해 주세요 (예제는 고치지 않음)")
    path, err = wb.save_text(wb.USER_RECIPES, recipe.id, wb.dump_yaml(recipe), req.overwrite)
    if err:
        raise HTTPException(409, err)
    return {"file": path.name, "source": "user", "path": wb.shown(path)}


# ── 조합 공간 ─────────────────────────────────────────────────────────────


def _load_space_text(text: str):
    from msl.recipe.space import Space

    data, errors = wb.parse_yaml(text)
    if data is None:
        return None, errors
    try:
        return Space.model_validate(data), []
    except Exception as exc:
        from pydantic import ValidationError

        return None, wb.errors_of(exc) if isinstance(exc, ValidationError) else [str(exc)]


@app.get("/spaces", response_class=HTMLResponse)
def spaces() -> str:
    from msl.recipe.space import expand, load_space
    from msl.report.space import describe_space, page

    groups: dict[str, list[str]] = {"example": [], "user": []}
    for path, source in _space_files():
        try:
            sp = load_space(path)
            n = len(expand(sp))
        except RecipeError as exc:
            groups[source].append(f'<li><b>{html.escape(path.name)}</b> — 정의 오류: {html.escape(str(exc))}</li>')
            continue
        groups[source].append(
            f'<li><a href="/space/{html.escape(path.stem)}"><b>{html.escape(sp.name)}</b></a> '
            f'<span class="sub">· {html.escape(describe_space(sp))} · 변형 {n}개</span> '
            f'<button type="button" class="linkish" data-edit="{html.escape(path.stem)}">편집기로</button>'
            f'<br><span class="sub">{html.escape(sp.description or "")}</span></li>')
    user = "".join(groups["user"]) or '<li class="sub">아직 없음 — 아래 편집기에서 만들어 [저장]하면 여기에 나옵니다.</li>'
    body = ('<p class="eyebrow"><a href="/">← 가상 조합 실험실</a></p><h1>조합 공간</h1>'
            '<p class="sub">레시피 틀 하나에서 혼합비·양·조합을 바꾼 레시피를 여러 개 만들어 한 번에 시험합니다. '
            '처음 실행은 계산 시간이 걸리고, 이후는 캐시에서 읽습니다.</p>'
            f'<h2 class="sec">예제</h2><ul class="spaces">{"".join(groups["example"])}</ul>'
            f'<h2 class="sec">내 조합 공간</h2><ul class="spaces">{user}</ul>'
            + (STATIC / "space-editor.html").read_text(encoding="utf-8"))
    return page("조합 공간", body)


@app.get("/space/{name}", response_class=HTMLResponse)
def space(name: str, no_cache: bool = False) -> str:
    from msl.recipe.space import load_space, run_space
    from msl.report.space import to_html

    path = wb.find_file(_space_files(), name)
    if path is None:
        raise HTTPException(404, "공간 파일이 없음")
    try:
        res = run_space(load_space(path), registry, use_cache=not no_cache)
    except RecipeError as exc:
        raise HTTPException(400, str(exc)) from exc
    return to_html(res, nav='<a href="/spaces">← 조합 공간 목록</a>')


@app.get("/api/space/{name}/yaml", response_class=PlainTextResponse)
def space_yaml(name: str) -> str:
    path = wb.find_file(_space_files(), name)
    if path is None:
        raise HTTPException(404, "공간 파일이 없음")
    return path.read_text(encoding="utf-8")


@app.post("/api/space/check")
def space_check(req: YamlRequest) -> dict[str, Any]:
    from msl.recipe.space import expand
    from msl.report.space import describe_space

    sp, errors = _load_space_text(req.yaml)
    if sp is None:
        return {"ok": False, "errors": errors}
    try:
        variants = expand(sp)
    except (RecipeError, ValueError) as exc:
        return {"ok": False, "errors": [str(exc)]}
    return {"ok": True, "errors": [], "n": len(variants), "describe": describe_space(sp), "name": sp.name,
            "variants": [{"id": v.recipe.id, "label": v.label} for v in variants[:60]]}


@app.post("/api/space/run")
def space_run(req: YamlRequest) -> dict[str, Any]:
    from msl.recipe.space import run_space
    from msl.report.space import to_html

    sp, errors = _load_space_text(req.yaml)
    if sp is None:
        raise HTTPException(422, {"errors": errors})
    try:
        res = run_space(sp, registry, use_cache=not req.no_cache)
    except (RecipeError, ValueError) as exc:
        raise HTTPException(422, {"errors": [str(exc)]}) from exc
    return {"html": to_html(res), "rows": len(res.rows), "elapsed": round(res.elapsed, 1)}


@app.post("/api/space/save")
def space_save(req: YamlRequest) -> dict[str, Any]:
    sp, errors = _load_space_text(req.yaml)
    if sp is None:
        raise HTTPException(422, {"errors": errors})
    if sp.id in _example_ids("space"):
        raise HTTPException(409, f"예제와 같은 id: {sp.id} — 다른 id 로 저장해 주세요 (예제는 고치지 않음)")
    path, err = wb.save_text(wb.USER_SPACES, sp.id, req.yaml.rstrip() + "\n", req.overwrite)
    if err:
        raise HTTPException(409, err)
    return {"file": path.name, "path": wb.shown(path)}


# ── 물성 예측 (L1 조성 대리모델) ──────────────────────────────────────────


@app.get("/predict", response_class=HTMLResponse)
def predict_page() -> str:
    from msl.report.space import page

    return page("물성 예측", (STATIC / "predict.html").read_text(encoding="utf-8"))


@app.get("/api/predict/info")
def predict_info() -> dict[str, Any]:
    from msl.ml.predict import ModelMissing, info

    try:
        return info()
    except ModelMissing as exc:
        return {"error": str(exc)}


class PredictRequest(BaseModel):
    texts: list[str]


def _formula_of(text: str) -> str:
    """화학식이면 그대로, 아니면 해석기로 화학식을 찾는다 (이름·CAS·광물명)."""
    from pymatgen.core import Composition

    t = text.strip()
    if wb.FORMULA_RE.match(t):
        Composition(t)
        return t
    from msl.resolve import resolve
    from msl.schema.recipe import Component

    r = resolve(Component.model_validate({"ref": wb.guess_ref(t)}))
    if not r.formula:
        raise ValueError(f"화학식을 찾지 못함 — {r.error or t}")
    return r.formula


@app.post("/api/predict")
def predict_api(req: PredictRequest) -> dict[str, Any]:
    from msl.ml.predict import ModelMissing, predict

    if len(req.texts) > 50:
        raise HTTPException(422, "한 번에 50개까지")
    results = []
    for text in req.texts:
        try:
            p = predict(_formula_of(text))
            p["input"] = text
            results.append(p)
        except ModelMissing as exc:
            raise HTTPException(503, str(exc)) from exc
        except Exception as exc:
            results.append({"input": text, "error": str(exc)})
    return {"results": results}
