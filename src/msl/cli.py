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


@app.command("serve")
def serve(port: Annotated[int, typer.Option(help="포트")] = 8000) -> None:
    """웹 화면을 띄운다 (http://127.0.0.1:포트)."""
    import uvicorn

    typer.echo(f"가상 조합 실험실 — http://127.0.0.1:{port}  (종료: Ctrl+C)")
    uvicorn.run("msl.web.app:app", host="127.0.0.1", port=port, log_level="warning")


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
