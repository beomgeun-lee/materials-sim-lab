"""msl 명령줄 도구."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from msl.registry.load import KB_DIR, check_registry, load_registry, summarize

app = typer.Typer(help="materials-sim-lab — 조합 시뮬레이션 프레임워크", no_args_is_help=True)
registry_app = typer.Typer(help="지식 레지스트리(kb/*.yaml)", no_args_is_help=True)
recipe_app = typer.Typer(help="조합 레시피", no_args_is_help=True)
app.add_typer(registry_app, name="registry")
app.add_typer(recipe_app, name="recipe")

KbOption = Annotated[Path, typer.Option("--kb", help="kb 디렉토리")]


@registry_app.command("check")
def registry_check(kb: KbOption = KB_DIR) -> None:
    """레지스트리 형식·교차 참조·라이선스 정책을 검사한다."""
    reg = load_registry(kb)
    typer.echo(
        f"라이선스 {len(reg.licenses)} · 소스 {len(reg.sources)} · "
        f"엔진 {len(reg.engines)} · 시험 {len(reg.assays)}"
    )
    for title, counter in summarize(reg).items():
        if counter:
            items = ", ".join(f"{k} {v}" for k, v in sorted(counter.items()))
            typer.echo(f"  {title}: {items}")
    issues = check_registry(reg)
    if issues:
        typer.secho(f"\n문제 {len(issues)}건", fg=typer.colors.RED, err=True)
        for issue in issues:
            typer.echo(f"  - {issue}", err=True)
        raise typer.Exit(1)
    typer.secho("\n통과", fg=typer.colors.GREEN)


STATUS_MARK = {"ok": "●", "warning": "▲", "failed": "✗", "not-applicable": "–", "pending": "…"}


@app.command("run")
def run(
    path: Path,
    no_cache: Annotated[bool, typer.Option("--no-cache", help="캐시를 쓰지 않고 다시 계산")] = False,
    json_out: Annotated[Path | None, typer.Option("--json", help="리포트를 JSON 으로 저장")] = None,
    report_out: Annotated[Path | None, typer.Option("--report", help="리포트 파일 (.html 또는 .md)")] = None,
) -> None:
    """레시피 하나를 실행한다: 해석 → 안전 게이트 → 라우팅된 시험들."""
    from msl.runtime.runner import run_recipe, save_report
    from msl.schema.recipe import RecipeError, load_recipe

    reg = load_registry()
    try:
        recipe = load_recipe(path, known_assays=set(reg.assays))
    except RecipeError as exc:
        typer.secho(f"레시피 오류\n{exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    report = run_recipe(recipe, use_cache=not no_cache, registry=reg)
    typer.secho(f"{recipe.name}  [{recipe.id}]", bold=True)
    for c in report.components:
        info = c.error or f"{c.formula or '-'} · {c.via}" + (f" · 반응성 그룹 {', '.join(c.groups)}" if c.groups else "")
        typer.echo(f"  {str(c.ref):26} {c.name[:30]:30} {info}")
    typer.echo("")
    for r in report.results:
        name = reg.assays[r.assay].name if r.assay in reg.assays else r.assay
        typer.echo(f"{STATUS_MARK.get(r.status.value, '?')} {r.assay:4} {name:12} [{r.fidelity.value}] {r.data.get('summary', '')}")
        for v in r.values:
            unit = f" {v.unit}" if v.unit else ""
            shown = int(v.value) if isinstance(v.value, float) and v.value.is_integer() else v.value
            typer.echo(f"        {v.name}: {shown}{unit}")
        for w in r.warnings:
            typer.secho(f"        ! {w.splitlines()[0]}", fg=typer.colors.YELLOW)
    typer.echo(f"\n{report.elapsed:.1f}초")
    if json_out:
        save_report(report, reg, json_out)
        typer.echo(f"저장: {json_out}")
    if report_out:
        from msl.report import payload, to_html, to_markdown

        body = payload(report, reg, recipe_yaml=Path(path).read_text(encoding="utf-8"))
        text = to_markdown(body) if report_out.suffix.lower() == ".md" else to_html(body)
        report_out.parent.mkdir(parents=True, exist_ok=True)
        report_out.write_text(text, encoding="utf-8")
        typer.echo(f"리포트: {report_out}")


@app.command("serve")
def serve(port: Annotated[int, typer.Option(help="포트")] = 8000) -> None:
    """웹 화면을 띄운다 (http://127.0.0.1:포트)."""
    import uvicorn

    typer.echo(f"가상 조합 실험실 — http://127.0.0.1:{port}  (종료: Ctrl+C)")
    uvicorn.run("msl.web.app:app", host="127.0.0.1", port=port, log_level="warning")


db_app = typer.Typer(help="데이터 적재 (1단계 데이터 코어)", no_args_is_help=True)
app.add_typer(db_app, name="db")


def _pick(source: str) -> dict:
    from msl.connectors import available

    mods = available()
    if source == "all":
        return mods
    if source not in mods:
        typer.secho(f"커넥터 없음: {source} (있는 것: {', '.join(mods)})", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    return {source: mods[source]}


@db_app.command("fetch")
def db_fetch(source: str) -> None:
    """원본 스냅샷을 받는다 (source id 또는 all)."""
    for sid, mod in _pick(source).items():
        typer.echo(f"↓ {sid} … ", nl=False)
        path = mod.fetch()
        typer.secho(f"{path}", fg=typer.colors.GREEN)


@db_app.command("load")
def db_load(source: str) -> None:
    """스냅샷을 정규화해 적재한다 (source id 또는 all)."""
    for sid, mod in _pick(source).items():
        typer.echo(f"⇢ {sid} … ", nl=False)
        result = mod.load()
        typer.secho(", ".join(f"{t} {parts}" for t, parts in result.items()), fg=typer.colors.GREEN)


@db_app.command("status")
def db_status() -> None:
    """적재 현황: 테이블 · 소스 · 파티션 · 행 수 · 버전."""
    from msl.connectors import available
    from msl.db import status

    rows = status()
    typer.echo(f"{'테이블':18} {'소스':24} {'파티션':10} {'행':>9}  버전")
    for r in rows:
        typer.echo(f"{r.table:18} {r.source:24} {r.partition:10} {r.rows:>9,}  {r.version[:40]}")
    loaded = {r.source for r in rows}
    waiting = [s for s in available() if s not in loaded]
    typer.echo(f"\n적재 {len(loaded)}개 소스 · {sum(r.rows for r in rows):,}행" + (f" · 커넥터만 있음: {', '.join(waiting)}" if waiting else ""))


@app.command("validate")
def validate(
    path: Annotated[Path, typer.Option("--set", help="검증 세트")] = KB_DIR / "validation" / "assays_v0.yaml",
    no_cache: Annotated[bool, typer.Option("--no-cache", help="캐시 없이 다시 계산")] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="통과 사례의 검사 내용도 출력")] = False,
) -> None:
    """검증 세트를 돌려 채점한다 (계획서 8절: S0 재현율 100%, 나머지 90% 이상)."""
    from msl.validation import run_set, score

    results = run_set(path, registry=load_registry(), use_cache=not no_cache)
    for r in results:
        mark = typer.style("✓", fg=typer.colors.GREEN) if r.passed else typer.style("✗", fg=typer.colors.RED)
        typer.echo(f"{mark} {r.assay:4} {r.id:30} {r.summary[:70]}")
        if verbose or not r.passed:
            for d in r.detail:
                typer.echo(f"        {d}")
    s = score(results)
    typer.echo("\n" + " · ".join(f"{a} {p}/{t}" for a, (p, t) in s["by_assay"].items()))
    color = typer.colors.GREEN if s["passed"] else typer.colors.RED
    typer.secho(f"S0 위험 재현율 {s['s0_recall']:.0%} (기준 100%) · 나머지 통과율 {s['other_rate']:.0%} (기준 90%)", fg=color)
    if not s["passed"]:
        raise typer.Exit(1)


space_app = typer.Typer(help="조합 공간 생성기 — 레시피 틀에서 여러 레시피를 만들어 한 표로 (3단계)", no_args_is_help=True)
app.add_typer(space_app, name="space")


def _load_space(path: Path):
    from msl.recipe.space import load_space
    from msl.schema.recipe import RecipeError

    try:
        return load_space(path)
    except RecipeError as exc:
        typer.secho(f"공간 정의 오류\n{exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc


@space_app.command("expand")
def space_expand(path: Path) -> None:
    """공간이 만들 레시피 목록만 보여 준다 (실행하지 않음)."""
    from msl.recipe.space import expand
    from msl.report.space import describe
    from msl.schema.recipe import RecipeError

    space = _load_space(path)
    try:
        variants = expand(space)
    except RecipeError as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    typer.secho(f"{space.name}  [{space.id}] — 변형 {len(variants)}개", bold=True)
    for v in variants:
        comps = ", ".join(f"{c.ref}{' ' + str(c.amount.value) + ' ' + c.amount.unit if c.amount else ''}" for c in v.recipe.components)
        typer.echo(f"  {v.recipe.id:34} {v.label:28} {comps}")


@space_app.command("run")
def space_run(
    path: Path,
    no_cache: Annotated[bool, typer.Option("--no-cache", help="캐시 없이 다시 계산")] = False,
    csv_out: Annotated[Path | None, typer.Option("--csv", help="결과 표를 CSV 로 저장")] = None,
    report_out: Annotated[Path | None, typer.Option("--report", help="리포트 파일 (.html 또는 .md)")] = None,
    json_out: Annotated[Path | None, typer.Option("--json", help="결과를 JSON 으로 저장")] = None,
) -> None:
    """공간의 모든 레시피를 실행하고 collect 값을 한 표로 모은다."""
    import json

    from msl.recipe.space import run_space
    from msl.report.space import describe, to_csv, to_html, to_markdown
    from msl.schema.recipe import RecipeError

    space = _load_space(path)
    progress = lambda i, n, v: typer.echo(f"\r  {i + 1}/{n} {v.label[:50]:50}", nl=False)
    try:
        res = run_space(space, load_registry(), use_cache=not no_cache, progress=progress)
    except RecipeError as exc:
        typer.secho(f"\n{exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    typer.echo("\r" + " " * 60 + "\r", nl=False)
    typer.secho(f"{space.name}  [{space.id}]", bold=True)
    typer.echo(f"{describe(res)} · 변형 {len(res.rows)}개 · {res.elapsed:.1f}초\n")
    width = max(len(r["label"]) for r in res.rows) if res.rows else 10
    typer.echo(f"{'변형':{width}}  " + "  ".join(f"{c[:22]:>22}" for c in res.columns))
    for r in res.rows:
        cells = []
        for c in res.columns:
            v = r["values"].get(c)
            v = "—" if v is None else (f"{v:.4g}" if isinstance(v, float) else str(v))
            cells.append(f"{v[:22]:>22}")
        typer.echo(f"{r['label']:{width}}  " + "  ".join(cells))
    for out, text in ((csv_out, to_csv), (json_out, None), (report_out, None)):
        if out is None:
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        if out is csv_out:
            out.write_text(to_csv(res), encoding="utf-8-sig")
        elif out is json_out:
            out.write_text(json.dumps(res.to_json(), ensure_ascii=False, indent=1), encoding="utf-8")
        else:
            out.write_text(to_markdown(res) if out.suffix.lower() == ".md" else to_html(res), encoding="utf-8")
        typer.echo(f"저장: {out}")


bench_app = typer.Typer(help="엔진 기준 재현 (3단계 완료 기준)", no_args_is_help=True)
app.add_typer(bench_app, name="bench")


@bench_app.command("phreeqc")
def bench_phreeqc(
    dist: Annotated[Path, typer.Option("--dist", help="USGS PHREEQC 배포본 폴더 (phreeqc-3.8.6-17100)")],
    binary: Annotated[Path, typer.Option("--binary", help="같은 배포본으로 빌드한 공식 phreeqc 실행 파일")],
    only: Annotated[str | None, typer.Option("--only", help="예제 이름 (쉼표로)")] = None,
    save: Annotated[Path | None, typer.Option("--save", help="공식 실행 결과(기준값)를 JSON 으로 저장")] = None,
) -> None:
    """PHREEQC 공식 예제를 공식 실행 파일과 A4 경로(IPhreeqc)로 돌려 비교한다."""
    import json

    from msl.bench.phreeqc import bench

    res = bench(dist, binary, only.split(",") if only else None,
                progress=lambda n: typer.echo(f"\r  {n:8}", nl=False))
    typer.echo("\r" + " " * 20 + "\r", nl=False)
    ok = n_cmp = 0
    for name, r in res["examples"].items():
        if "error" in r:
            typer.secho(f"✗ {name:6} {r['error']}", fg=typer.colors.RED)
            continue
        c, fc = r["compare"], r["files"]
        ok += r["reproduced"]
        n_cmp += r["comparable"]
        mark = (typer.style("✓", fg=typer.colors.GREEN) if r["reproduced"] else
                typer.style("✗", fg=typer.colors.RED) if r["comparable"] else typer.style("–", fg=typer.colors.YELLOW))
        eff = r.get("db_version_effect")
        detail = (f"용액 {c['solutions'][0]:>3}/{c['solutions'][1]:<3} ΔpH {c['max_dpH']:.3f} Δpe {c['max_dpe']:.3f} "
                  f"ΔI {c['max_rel_dI']:.0e} ΔSI {c['max_dSI']:.2f}" if c["solutions"][0] else "용액 기술 출력 없음")
        if fc["files"][0]:
            detail += f" · 파일 {','.join(fc['files'][0])} 값 {fc['cells']}개 최대 상대차 {fc['max_rel']:.0e}"
        elif not r["comparable"]:
            detail += " · 선택 출력 파일도 없음 (USER_GRAPH 전용) → 비교 불가"
        typer.echo(f"{mark} {name:6} {detail}  [{r['db']}]"
                   + (f"  (3.8.6 DB 와 차이: ΔpH {eff['max_dpH']:.3f}, ΔSI {eff['max_dSI']:.2f})" if eff and eff["compared"] else ""))
        if r["error_ours"] and not r["reproduced"]:
            typer.echo(f"        IPhreeqc 오류: {r['error_ours'].strip().splitlines()[-1][:120]}")
    typer.secho(f"\n재현 {ok}/{n_cmp} (비교 가능한 예제) · 비교 불가 {len(res['examples']) - n_cmp}",
                fg=typer.colors.GREEN if ok == n_cmp else typer.colors.YELLOW)
    if save:
        save.write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
        typer.echo(f"저장: {save}")


@app.command("resolve-check")
def resolve_check(
    path: Annotated[Path, typer.Option("--list", help="시험 목록")] = KB_DIR / "validation" / "resolver_v1.yaml",
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="통과 항목도 출력")] = False,
) -> None:
    """해석기 시험 목록으로 성공률을 잰다 (계획서 1단계 완료 기준: 95% 이상)."""
    import yaml
    from pymatgen.core import Composition

    from msl.resolve import resolve
    from msl.schema.recipe import Component

    items = yaml.safe_load(path.read_text(encoding="utf-8"))
    ok, by_ns, fails = 0, {}, []
    for it in items:
        r = resolve(Component(ref=it["ref"]))
        ns = it["ref"].split(":", 1)[0]
        good = r.formula is not None and r.error is None
        if good and it.get("expect_formula"):
            try:
                good = Composition(r.formula).reduced_composition.almost_equals(
                    Composition(it["expect_formula"]).reduced_composition)
            except Exception:
                good = False
        ok += good
        passed, total = by_ns.get(ns, (0, 0))
        by_ns[ns] = (passed + good, total + 1)
        if not good:
            fails.append(f"  ✗ {it['ref']:32} → {r.formula or '-'} (기대 {it.get('expect_formula', '화학식')}) {r.error or r.via}")
        elif verbose:
            typer.echo(f"  ✓ {it['ref']:32} → {r.formula} · {r.via}")
    for line in fails:
        typer.echo(line)
    rate = ok / len(items)
    typer.echo("\n" + " · ".join(f"{ns} {p}/{t}" for ns, (p, t) in by_ns.items()))
    color = typer.colors.GREEN if rate >= 0.95 else typer.colors.RED
    typer.secho(f"해석 성공 {ok}/{len(items)} = {rate:.1%} (기준 95%)", fg=color)
    if rate < 0.95:
        raise typer.Exit(1)


@recipe_app.command("check")
def recipe_check(paths: list[Path], kb: KbOption = KB_DIR) -> None:
    """레시피 YAML 을 스키마와 시험 카탈로그로 검증한다."""
    from msl.schema.recipe import RecipeError, load_recipe

    known = set(load_registry(kb).assays)
    failed = 0
    for path in paths:
        try:
            recipe = load_recipe(path, known_assays=known)
        except RecipeError as exc:
            failed += 1
            typer.secho(f"✗ {path}", fg=typer.colors.RED, err=True)
            for line in str(exc).splitlines():
                typer.echo(f"    {line}", err=True)
            continue
        refs = ", ".join(str(c.ref) for c in recipe.components)
        typer.secho(f"✓ {path}  [{recipe.id}] {refs}", fg=typer.colors.GREEN)
    if failed:
        raise typer.Exit(1)
