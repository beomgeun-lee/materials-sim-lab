"""웹 작업대 — 화면에서 만든 레시피·조합 공간을 검증·변환·저장한다.

화면은 레시피를 '폼' 형태(수량은 "25 °C" 같은 문자열)로 주고받는다. 검증은 CLI 와 같은 스키마(Recipe·Space)로 한다.
사용자 파일은 예제(examples/)와 섞지 않고 recipes/ · spaces/ 에 둔다.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml
from pydantic import ValidationError

from msl.env import ROOT
from msl.registry.load import Registry
from msl.schema.recipe import Recipe
from msl.schema.refs import ELEMENTS, Namespace

EXAMPLE_RECIPES = ROOT / "examples" / "recipes"
EXAMPLE_SPACES = ROOT / "examples" / "spaces"
USER_RECIPES = ROOT / "recipes"
USER_SPACES = ROOT / "spaces"

CAS_RE = re.compile(r"^\d{2,7}-\d{2}-\d$")
FORMULA_RE = re.compile(r"^(?:[A-Z][a-z]?\d*(?:\.\d+)?|\((?:[A-Z][a-z]?\d*)+\)\d*)+(?:[·.]\d*(?:[A-Z][a-z]?\d*)+)*$")
HANGUL = re.compile(r"[가-힣]")
NAMESPACES = {n.value for n in Namespace}
QUANTITY_KEYS = ("T", "P", "time")


def guess_ref(text: str) -> str:
    """자유 입력 → 참조 문자열. 'cas:…' 처럼 네임스페이스가 있으면 그대로 쓴다."""
    t = text.strip()
    ns = t.split(":", 1)[0].lower() if ":" in t else ""
    if ns in NAMESPACES:
        return f"{ns}:{t.split(':', 1)[1].strip()}"
    if CAS_RE.match(t):
        return f"cas:{t}"
    if re.fullmatch(r"KE-\d+", t, re.I):
        return f"ke:{t.upper()}"
    if t in ELEMENTS:
        return f"element:{t}"
    if HANGUL.search(t):
        return f"name:{t}"
    if FORMULA_RE.match(t):
        return f"formula:{t}"
    from msl.resolve import _mineral_key, mineral_db

    if _mineral_key(t) in mineral_db():
        return f"mineral:{t.lower()}"
    return f"name:{t}"


def _clean(v: Any) -> Any:
    """빈 문자열·None·빈 목록을 뺀다 (폼에서 비워 둔 칸)."""
    if isinstance(v, dict):
        out = {k: _clean(x) for k, x in v.items()}
        return {k: x for k, x in out.items() if x not in (None, "", [], {})}
    if isinstance(v, list):
        return [_clean(x) for x in v]
    if isinstance(v, str):
        return v.strip()
    return v


FIELD_KO = {"id": "id", "name": "이름", "components": "성분", "ref": "참조", "amount": "양", "state": "상태",
            "conditions": "조건", "T": "온도", "P": "압력", "pH": "pH", "Eh": "전위 Eh", "sweep": "스윕", "assays": "시험",
            "mode": "모드", "generator": "생성기", "collect": "수집 열", "base": "레시피 틀"}


def _loc(loc: tuple) -> str:
    parts = []
    for x in loc:
        if isinstance(x, int):
            parts.append(f"{x + 1}번째")
        elif x in FIELD_KO:
            parts.append(FIELD_KO[x])
        elif not str(x).startswith(("function-", "tagged-union")):
            parts.append(str(x))
    return " · ".join(parts) or "(레시피)"


def _msg(e: dict) -> str:
    t, ctx = e["type"], e.get("ctx") or {}
    if t == "missing":
        return "꼭 적어야 하는 항목이 없음" + (" (성분은 하나 이상)" if e["loc"][-1:] == ("components",) else "")
    if t == "string_pattern_mismatch":
        pat = ctx.get("pattern", "")
        if pat.startswith("^rcp-"):
            return "형식이 맞지 않음 — rcp- 로 시작하고 영소문자·숫자·하이픈만 (예: rcp-my-test)"
        if pat.startswith("^spc-"):
            return "형식이 맞지 않음 — spc- 로 시작하고 영소문자·숫자·하이픈만 (예: spc-my-space)"
        return f"형식이 맞지 않음 ({pat})"
    if t in ("too_short", "too_long"):
        return f"개수가 맞지 않음 (최소 {ctx.get('min_length', '?')}개)" if t == "too_short" else f"너무 많음 (최대 {ctx.get('max_length', '?')}개)"
    if t in ("greater_than_equal", "greater_than", "less_than_equal", "less_than"):
        words = {"ge": "최소", "gt": "초과", "le": "최대", "lt": "미만"}
        return "값이 허용 범위를 벗어남 (" + ", ".join(f"{words.get(k, k)} {v:g}" if isinstance(v, (int, float)) else f"{k} {v}"
                                          for k, v in ctx.items()) + ")"
    if t in ("float_parsing", "int_parsing", "float_type", "int_type"):
        return "숫자가 아님"
    if t in ("enum", "literal_error"):
        return f"허용하는 값이 아님 — {ctx.get('expected', '')}"
    if t == "extra_forbidden":
        return "알 수 없는 항목 (철자 확인)"
    return e["msg"].removeprefix("Value error, ")


def errors_of(exc: ValidationError) -> list[str]:
    return [f"{_loc(e['loc'])}: {_msg(e)}" for e in exc.errors()]


def validate_recipe(form: dict[str, Any], registry: Registry) -> tuple[Recipe | None, list[str]]:
    data = _clean(form)
    try:
        recipe = Recipe.model_validate(data)
    except ValidationError as exc:
        return None, errors_of(exc)
    if recipe.assays != "auto":
        unknown = [c for c in recipe.assays if c not in registry.assays]
        if unknown:
            return None, [f"assays: 시험 카탈로그에 없는 코드 {unknown}"]
    return recipe, []


def to_form(recipe: Recipe) -> dict[str, Any]:
    """Recipe → 화면 폼 (수량은 문자열)."""
    cond = recipe.conditions
    form: dict[str, Any] = {
        "id": recipe.id, "name": recipe.name, "mode": recipe.mode.value,
        "components": [{"ref": str(c.ref), "amount": str(c.amount) if c.amount else "",
                        "state": c.state.value if c.state else "", "label": c.label or ""} for c in recipe.components],
        "conditions": {k: (str(v) if k in QUANTITY_KEYS else v) for k in type(cond).model_fields
                       if (v := getattr(cond, k)) is not None},
        "assays": recipe.assays,
        "notes": recipe.notes or "",
    }
    if recipe.sweep:
        sw = recipe.sweep
        form["sweep"] = {"param": sw.param, "start": str(sw.start), "stop": str(sw.stop), "steps": sw.steps}
    return form


def dump_yaml(recipe: Recipe) -> str:
    """사람이 읽기 좋은 순서의 레시피 YAML."""
    f = to_form(recipe)
    data: dict[str, Any] = {"id": f["id"], "name": f["name"]}
    if f["mode"] != "mixture":
        data["mode"] = f["mode"]
    data["components"] = [_clean(c) for c in f["components"]]
    if f["conditions"]:
        data["conditions"] = f["conditions"]
    if "sweep" in f:
        data["sweep"] = f["sweep"]
    data["assays"] = f["assays"]
    if f["notes"]:
        data["notes"] = f["notes"]
    return yaml.safe_dump(data, allow_unicode=True, sort_keys=False, width=120)


def parse_yaml(text: str) -> tuple[dict[str, Any] | None, list[str]]:
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        return None, [f"YAML 형식 오류: {exc}"]
    if not isinstance(data, dict):
        return None, ["YAML 최상위가 매핑(키: 값)이 아님"]
    return data, []


def list_files(folder: Path, source: str) -> list[tuple[Path, str]]:
    return [(p, source) for p in sorted(folder.glob("*.yaml"))] if folder.exists() else []


def find_file(folder_pairs: list[tuple[Path, str]], name: str) -> Path | None:
    """파일 이름(확장자 포함 또는 제외)으로 찾는다. 경로 이탈 방지."""
    stem = Path(name).stem
    if not re.fullmatch(r"[A-Za-z0-9_.\-]+", stem):
        return None
    return next((p for p, _ in folder_pairs if p.stem == stem), None)


def shown(path: Path) -> str:
    """화면에 보여 줄 경로 — 프로젝트 안이면 상대 경로."""
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def save_text(folder: Path, stem: str, text: str, overwrite: bool) -> tuple[Path | None, str | None]:
    """사용자 폴더에 저장. stem 은 레시피·공간 id (스키마가 이미 영소문자·숫자·하이픈으로 제한)."""
    if not re.fullmatch(r"(rcp|spc)-[a-z0-9][a-z0-9\-]*", stem):
        return None, f"파일 이름으로 쓸 수 없는 id: {stem}"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{stem}.yaml"
    if path.exists() and not overwrite:
        return None, f"이미 있는 파일: {shown(path)} — 덮어쓰려면 확인 필요"
    path.write_text(text, encoding="utf-8")
    return path, None
