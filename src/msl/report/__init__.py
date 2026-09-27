"""통합 리포트 — 화면(JSON)·HTML 파일·Markdown 파일 (계획서 6장 2단계)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from msl.registry.load import Registry
from msl.runtime.runner import Report

STATIC = Path(__file__).resolve().parents[1] / "web" / "static" / "index.html"
STATUS_KO = {"ok": "완료", "warning": "주의", "not-applicable": "해당 없음", "pending": "대기", "failed": "실패"}


def payload(report: Report, registry: Registry, recipe_yaml: str | None = None) -> dict[str, Any]:
    """웹 화면이 그리는 JSON — 결과에 라이선스·소스 이름·시험 이름을 붙인다."""
    body = report.to_json(registry)
    body["licenses"] = {k: {"name": v.name, "partition": v.partition.value} for k, v in registry.licenses.items()}
    body["source_names"] = {k: v.name for k, v in registry.sources.items()}
    body["assay_names"] = {k: v.name for k, v in registry.assays.items()}
    if recipe_yaml is not None:
        body["recipe_yaml"] = recipe_yaml
    return body


def to_html(body: dict[str, Any]) -> str:
    """웹 화면 템플릿에 결과를 넣은 단일 HTML (서버 없이 열린다)."""
    html = STATIC.read_text(encoding="utf-8")
    data = json.dumps(body, ensure_ascii=False).replace("</", "<\\/")
    return html.replace("<script>\nconst $ =", f"<script>window.__REPORT__ = {data};</script>\n<script>\nconst $ =", 1)


def _fmt(v: Any) -> str:
    if isinstance(v, float):
        return f"{v:,.0f}" if v.is_integer() else f"{v:.4g}"
    return str(v)


def to_markdown(body: dict[str, Any]) -> str:
    r = body["recipe"]
    lines = [f"# {r['name']}", "", f"`{r['id']}` · 서버 처리 {body['elapsed']}초", "", "## 성분", "",
             "| 참조 | 이름 | 화학식 | 양 | 해석 |", "|---|---|---|---|---|"]
    for c in body["components"]:
        lines.append(f"| `{c['ref']}` | {c['name']} | {c['formula'] or '—'} | {c['amount'] or '—'} | {c['error'] or c['via']} |")
    bal = body.get("balance") or {}
    if bal.get("rows"):
        lines += ["", "## 질량수지", "", "| 화학식 | mol | g | wt% |", "|---|---|---|---|"]
        for row in bal["rows"]:
            wt = f"{row['mass_fraction'] * 100:.4g}" if row.get("mass_fraction") is not None else "—"
            lines.append(f"| {row['formula']} | {_fmt(row['moles'])} | {_fmt(row['mass_g']) if row.get('mass_g') is not None else '—'} | {wt} |")
    lines += ["", "## 시험 결과", ""]
    for x in body["results"]:
        name = body["assay_names"].get(x["assay"], x["assay"])
        lines += [f"### {x['assay']} {name} — {STATUS_KO.get(x['status'], x['status'])} · 충실도 {x['fidelity']}", "",
                  (x["data"] or {}).get("summary", ""), ""]
        for w in x["warnings"]:
            if x["status"] != "failed":
                lines.append(f"> ⚠ {w}")
        if x["values"]:
            lines += ["", "| 항목 | 값 |", "|---|---|"]
            lines += [f"| {v['name']} | {_fmt(v['value'])}{' ' + v['unit'] if v.get('unit') else ''} |" for v in x["values"]]
        if x["status"] in ("ok", "warning"):
            src = ", ".join(f"{body['source_names'].get(s['source'], s['source'])} ({s['license']})" for s in x["sources"])
            lines += ["", f"- 엔진: {x['engine']} {x['engine_version']}" + (f" · 에너지 기준: {x['energy_reference']}" if x.get("energy_reference") else ""),
                      f"- 조건: {x['conditions_basis']}", f"- 출처: {src or '—'}"]
            lines += [f"- 한계: {c}" for c in x["caveats"]]
        lines.append("")
    return "\n".join(lines)
