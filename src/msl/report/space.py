"""조합 공간 리포트 — CSV · Markdown · HTML (웹 화면과 같은 스타일의 단일 파일)."""

from __future__ import annotations

import csv
import html
import io
import json
import math
import re
from typing import Any

from msl.recipe.space import AmountGen, RatioGen, SpaceResult, SubsetGen
from msl.report import STATIC

VERDICT_CLASS = {"부적합": "bad", "주의": "warn", "규칙 해당 없음": "none", "해당 없음": "none"}


def _cell(v: Any) -> str:
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:,.0f}" if v.is_integer() and abs(v) >= 100 else f"{v:.4g}"
    return str(v)


def describe(res: SpaceResult) -> str:
    return describe_space(res.space)


def describe_space(space) -> str:
    g = space.generator
    if isinstance(g, RatioGen):
        return f"혼합비 격자 — {', '.join(c.short for c in g.components)} · 간격 {g.step:g}" + (f" · 총량 {g.total}" if g.total else f" ({g.unit})")
    if isinstance(g, AmountGen):
        return f"양 스윕 — {g.component.short} {g.start.value:g}→{g.stop.value:g} {g.start.unit} · {g.steps}점"
    assert isinstance(g, SubsetGen)
    req = f" + 항상 {', '.join(c.short for c in g.require)}" if g.require else ""
    return f"부분집합 열거 — {len(g.pool)}개 중 {', '.join(map(str, g.sizes))}개씩{req}"


def to_csv(res: SpaceResult) -> str:
    buf = io.StringIO()
    keys = list(res.rows[0]["params"]) if res.rows else []
    w = csv.writer(buf)
    w.writerow(["변형", *keys, *res.columns, "recipe_id"])
    for r in res.rows:
        params = [" + ".join(v) if isinstance(v, list) else v for v in r["params"].values()]
        w.writerow([r["label"], *params, *(r["values"].get(c) for c in res.columns), r["recipe_id"]])
    return buf.getvalue()


def to_markdown(res: SpaceResult) -> str:
    s = res.space
    lines = [f"# {s.name}", "", f"`{s.id}` · {describe(res)} · 변형 {len(res.rows)}개 · {res.elapsed:.1f}초", ""]
    if s.description:
        lines += [s.description, ""]
    lines += ["| 변형 | " + " | ".join(res.columns) + " |", "|---|" + "---|" * len(res.columns)]
    for r in res.rows:
        lines.append(f"| {r['label']} | " + " | ".join(_cell(r["values"].get(c)).replace("|", "\\|") for c in res.columns) + " |")
    return "\n".join(lines) + "\n"


FONT_LINKS = '<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin><link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Hahmlet:wght@500;700&family=IBM+Plex+Sans+KR:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">'  # 웹 화면과 같은 글꼴 (없으면 하위 화면 제목이 명조 대체 글꼴로 보였다)


def _style() -> str:
    m = re.search(r"<style>(.*?)</style>", STATIC.read_text(encoding="utf-8"), re.S)
    return m.group(1) if m else ""


def _nice(raw: float) -> float:
    p = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    n = raw / p
    return (1 if n < 1.5 else 2 if n < 3 else 5 if n < 7 else 10) * p


SERIES = ["var(--accent)", "var(--f-T)", "var(--ok)", "var(--f-L2)"]


def _line_chart(xname: str, xs: list[float], title: str, series: list[tuple[str, list[Any]]]) -> str:
    lines = [(name, [(x, y) for x, y in zip(xs, ys) if isinstance(y, (int, float))]) for name, ys in series]
    lines = [(n, pts) for n, pts in lines if len(pts) >= 2]
    if not lines:
        return ""
    W, H, L, R, T, B = 620, 240, 58, 18, 22, 40
    pw, ph = W - L - R, H - T - B
    allp = [p for _, pts in lines for p in pts]
    x0, x1 = min(p[0] for p in allp), max(p[0] for p in allp)
    y0, y1 = min(p[1] for p in allp), max(p[1] for p in allp)
    pad = (y1 - y0) * 0.08 or 1
    y0, y1 = y0 - pad, y1 + pad
    X = lambda v: L + (v - x0) / ((x1 - x0) or 1) * pw
    Y = lambda v: T + (y1 - v) / (y1 - y0) * ph
    g = ""
    step = _nice((y1 - y0) / 4)
    v = math.ceil(y0 / step) * step
    while v <= y1 + 1e-12:
        g += f'<line class="grid" x1="{L}" x2="{L + pw}" y1="{Y(v):.1f}" y2="{Y(v):.1f}"/><text x="{L - 6}" y="{Y(v) + 4:.1f}" text-anchor="end">{_cell(round(v, 6))}</text>'
        v += step
    xstep = _nice((x1 - x0) / 5 or 1)
    v = math.ceil(x0 / xstep) * xstep
    while v <= x1 + 1e-12:
        g += f'<text x="{X(v):.1f}" y="{H - B + 16}" text-anchor="middle">{_cell(round(v, 6))}</text>'
        v += xstep
    legend = ""
    for i, (name, pts) in enumerate(lines):
        col = SERIES[i % len(SERIES)]
        g += f'<polyline fill="none" stroke="{col}" stroke-width="2" points="{" ".join(f"{X(x):.1f},{Y(y):.1f}" for x, y in pts)}"/>'
        g += "".join(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="3" fill="{col}"><title>{html.escape(name)} · {html.escape(xname)} {x:g} → {y:g}</title></circle>' for x, y in pts)
        legend += f'<span><i style="background:{col}"></i>{html.escape(name)}</span>'
    g += f'<text x="{L}" y="{T - 7}">{html.escape(title)}</text><text x="{L + pw}" y="{H - 5}" text-anchor="end">{html.escape(xname)} →</text>'
    return (f'<div class="chart"><svg viewBox="0 0 {W} {H}" role="img" aria-label="{html.escape(title)}">{g}</svg>'
            + (f'<div class="legend-line">{legend}</div>' if len(lines) > 1 else "") + "</div>")


def _charts(res: SpaceResult) -> str:
    ax = res.axis()
    if not ax:
        return ""
    groups: dict[str, list[tuple[str, list[Any]]]] = {}
    for c in res.space.collect:
        groups.setdefault(c.chart or c.title, []).append((c.title, [r["values"].get(c.title) for r in res.rows]))
    charts = [_line_chart(ax[0], ax[1], title, series) for title, series in groups.items()]
    charts = [c for c in charts if c]
    return '<section class="charts">' + "".join(charts) + "</section>" if charts else ""


def _matrix(m: dict[str, Any]) -> str:
    head = "".join(f"<th>{html.escape(n)}</th>" for n in m["names"])
    rows = ""
    for name, cells in zip(m["names"], m["cells"]):
        tds = ""
        for v in cells:
            cls = VERDICT_CLASS.get(str(v), "")
            tds += f'<td class="mx {cls}">{html.escape(_cell(v)) if v is not None else ""}</td>'
        rows += f"<tr><th>{html.escape(name)}</th>{tds}</tr>"
    note = ""
    if any(v == "규칙 해당 없음" for row in m["cells"] for v in row):
        note = ('<p class="sub">‘규칙 해당 없음’은 S0 자체 규칙 v0 에 없는 조합이라는 뜻이며 안전하다는 뜻이 아니다 '
                '(예: 락스 + 에탄올은 클로로폼 등을 만들 수 있으나 v0 규칙에 없음). NOAA 쌍별 판정표 확보 후 보강.</p>')
    return (f'<div class="chart-title">{html.escape(m["title"])} — 쌍별 행렬</div>'
            f'<div class="tbl"><table class="t mxt"><thead><tr><th></th>{head}</tr></thead><tbody>{rows}</tbody></table></div>{note}')


def to_html(res: SpaceResult, nav: str = "") -> str:
    s = res.space
    body = [f'<p class="eyebrow">{nav + " · " if nav else ""}조합 공간 · {html.escape(s.id)}</p><h1>{html.escape(s.name)}</h1>',
            f'<p class="sub">{html.escape(describe(res))} · 변형 {len(res.rows)}개 · {res.elapsed:.1f}초</p>']
    if s.description:
        body.append(f"<p>{html.escape(s.description)}</p>")
    body.append(_charts(res))
    m = res.matrix()
    if m:
        body.append(_matrix(m))
    head = "".join(f"<th>{html.escape(c)}</th>" for c in ["변형", *res.columns])
    trs = ""
    for r in res.rows:
        tds = "".join(f'<td class="{VERDICT_CLASS.get(str(r["values"].get(c)), "")}">{html.escape(_cell(r["values"].get(c)))}</td>' for c in res.columns)
        trs += f'<tr><td>{html.escape(r["label"])}</td>{tds}</tr>'
    body.append(f'<div class="tbl"><table class="t"><thead><tr>{head}</tr></thead><tbody>{trs}</tbody></table></div>')
    data = json.dumps(res.to_json(), ensure_ascii=False).replace("</", "<\\/")
    return page(s.name, "".join(body), f'<script type="application/json" id="space-data">{data}</script>')


EXTRA_CSS = """
.page { max-width: 1040px; margin: 0 auto; padding: 24px 16px 48px; display: grid; gap: 14px; }
.page h1 { font-family: var(--font-display); font-size: 26px; font-weight: 700; }
.page .sub { color: var(--muted); font-size: 13px; }
.page a { color: var(--accent); }
.charts { display: grid; grid-template-columns: repeat(auto-fill, minmax(min(100%, 460px), 1fr)); gap: 12px; }
.charts .chart { background: var(--surface); border: 1px solid var(--line); border-radius: 8px; padding: 8px; }
table.t td.bad, table.t td.mx.bad { background: var(--bad-soft); color: var(--bad); font-weight: 600; }
table.t td.warn, table.t td.mx.warn { background: var(--warn-soft); color: var(--warn); font-weight: 600; }
table.mxt td.mx { text-align: center; min-width: 64px; }
table.t td.none { color: var(--muted); }
table.t { background: var(--surface); }
ul.spaces { display: grid; gap: 10px; list-style: none; padding: 0; }
ul.spaces li { background: var(--surface); border: 1px solid var(--line); border-radius: 8px; padding: 10px 12px; }
"""


def page(title: str, body: str, tail: str = "") -> str:
    """웹 화면과 같은 스타일의 단일 HTML 페이지."""
    return (f'<!doctype html><html lang="ko"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<title>{html.escape(title)}</title><link rel="icon" href="data:,">{FONT_LINKS}<style>{_style()}{EXTRA_CSS}</style></head>'
            f'<body><main class="page">{body}</main>{tail}</body></html>')
