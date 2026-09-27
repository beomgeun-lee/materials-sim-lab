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


@bench_app.command("umlip")
def bench_umlip(
    models: Annotated[str, typer.Option("--models", help="쉼표로 (mace-mpa-0, orb-v3)")] = "mace-mpa-0,orb-v3",
    save: Annotated[Path | None, typer.Option("--save", help="결과 JSON")] = None,
) -> None:
    """uMLIP 처리량·정확도 실측 — 대표 구조 20개, 노트북 CPU (계획서 4단계)."""
    import json

    from msl.bench.umlip import bench

    res = bench(models.split(","), log=typer.echo)
    typer.echo(f"\nCPU 스레드 {res['threads']} · torch {res['torch']}")
    typer.echo(f"{'모델':12} {'불러오기':>8} {'한 번(ms/원자)':>14} {'이완 중앙(s)':>12} {'단계 중앙':>9} {'수렴':>6} {'이완/시간':>9} "
               f"{'ΔE GGA':>8} {'GGA+U LASPH':>12} {'GGA+U 끔':>9}")
    for name, m in res["models"].items():
        s = m["summary"]
        f = lambda v: "—" if v is None else f"{v:.3f}"
        typer.echo(f"{name:12} {m['load_s']:7.1f}s {s['sp_ms_per_atom_median']:14.1f} {s['relax_s_median']:12.2f} {s['steps_median']:9.0f} "
                   f"{s['converged']:>3}/{s['n']:<2} {s['relax_per_hour']:9.0f} {f(s['mae_gga']):>8} {f(s['mae_ggau_lasph']):>12} {f(s['mae_ggau_nolasph']):>9}")
    if save:
        save.write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
        typer.echo(f"저장: {save}")


@bench_app.command("l2")
def bench_l2(
    formulas: Annotated[str, typer.Option("--formulas", help="쉼표로")] = "LiMn2O4,LiCoO2,Li2MnO3,LiFeO2,Mg2SiO4,ZnFe2O4,BaTiO3",
    methods: Annotated[str, typer.Option("--methods", help="v1(무작위) · v2(Ewald+원형)")] = "v1,v2",
    save: Annotated[Path | None, typer.Option("--save", help="결과 JSON")] = None,
) -> None:
    """L2 재발견 시험 — DB 물질을 모른다고 치고 구조를 생성해 DB 구조와의 에너지 차이를 잰다 (후보 생성 방식 비교)."""
    import json
    import statistics

    from msl import l2

    rows = []
    for f in formulas.split(","):
        for m in methods.split(","):
            try:
                rows.append(l2.rediscover(f, method=m, log=typer.echo))
            except Exception as exc:
                typer.secho(f"  {f} {m}: 실패 — {exc}", fg=typer.colors.RED)
    typer.echo(f"\n{'조성':10} {'방식':4} {'후보':>4} {'ORB 차이':>9} {'MACE 차이':>10}  최선 후보 (ORB)")
    for r in rows:
        o, m_ = r["models"].get("orb-v3", {}), r["models"].get("mace-mpa-0", {})
        typer.echo(f"{r['formula']:10} {r['method']:4} {r['n_candidates']:>4} {o.get('gap', float('nan')):+9.3f} {m_.get('gap', float('nan')):+10.3f}  {o.get('best', '')}")
    for m in methods.split(","):
        mine = [r for r in rows if r["method"] == m and "orb-v3" in r["models"]]
        gaps = [r["models"]["orb-v3"]["gap"] for r in mine if r["n_candidates"]]
        none = sum(1 for r in mine if not r["n_candidates"])
        if mine:
            typer.echo(f"{m}: ORB 차이 중앙값 {statistics.median(gaps) if gaps else float('nan'):+.3f} · "
                       f"0.01 이내 재발견 {sum(g <= 0.01 for g in gaps)}/{len(mine)}" + (f" · 후보를 못 만듦 {none}" if none else ""))
    if save:
        save.write_text(json.dumps(rows, ensure_ascii=False, indent=1), encoding="utf-8")
        typer.echo(f"저장: {save}")


@bench_app.command("props")
def bench_props(
    models: Annotated[str, typer.Option("--models")] = "mace-mpa-0,orb-v3",
    save: Annotated[Path | None, typer.Option("--save")] = None,
) -> None:
    """uMLIP 물성 검증 — 탄성률(MP DFT 대비)과 포논 안정성(알려진 안정·불안정 물질)."""
    import json

    from msl.bench.props import bench

    res = bench(models.split(","), log=typer.echo)
    for m, r in res["models"].items():
        s = r["summary"]
        typer.echo(f"{m}: 체적탄성률 MAE {s['K_mae']:.1f} GPa ({s['K_mape']:.0%}) · 전단 MAE {s['G_mae']:.1f} GPa · 포논 판정 {s['phonon_correct']}/{s['n']}")
    if save:
        save.write_text(json.dumps(res, ensure_ascii=False, indent=1, default=str), encoding="utf-8")
        typer.echo(f"저장: {save}")


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


ml_app = typer.Typer(help="L1 조성 대리모델 (4단계 선행)", no_args_is_help=True)
app.add_typer(ml_app, name="ml")


@ml_app.command("train")
def ml_train(
    stability_n: Annotated[int, typer.Option("--stability-n", help="안정성 점검에 쓸 시험 조성 수")] = 1500,
) -> None:
    """MP 바닥 다형으로 조성 → 형성에너지·밴드갭·금속 여부·부피 모델을 학습하고 평가한다."""
    from msl.ml.train import train

    train(log=typer.echo, stability_n=stability_n)


@app.command("predict")
def predict_cmd(
    formulas: Annotated[list[str], typer.Argument(help="화학식 (예: Li1.2Ni0.6Mn0.2O2 LiFePO4)")],
    json_out: Annotated[Path | None, typer.Option("--json", help="결과를 JSON 으로 저장")] = None,
) -> None:
    """조성만으로 물성을 예측한다 (L1 대리모델 — 형성에너지·안정성·밴드갭·밀도)."""
    import json

    from msl.ml.predict import ModelMissing, info, predict

    try:
        meta = info()
    except ModelMissing as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    r, g = meta["metrics"]["random"], meta["metrics"]["chemsys"]
    typer.echo(f"L1 조성 모델 {meta['version']} · MP {meta['data']['version']} {meta['data']['n']:,}개 조성으로 학습")
    typer.echo(f"  형성에너지 MAE {r['ef_mae']:.3f} (새 화학계 {g['ef_mae']:.3f}) eV/atom · "
               f"안정성 판정 정확도 {r['stability_accuracy']:.0%} (hull 거리 MAE {r['ehull_mae']:.3f})\n")
    out = []
    for f in formulas:
        try:
            p = predict(f)
        except Exception as exc:
            typer.secho(f"✗ {f}: {exc}", fg=typer.colors.RED)
            continue
        out.append(p)
        lo, hi = p["ef_interval"]
        typer.secho(f"{p['formula']}" + ("" if p["in_domain"] else "  ⚠ 학습 범위 밖 — 신뢰 낮음"), bold=True)
        typer.echo(f"  형성에너지 {p['ef']:+.3f} eV/atom (80% 구간 {lo:+.3f} ~ {hi:+.3f})")
        if p["ehull"] is not None:
            dec = " + ".join(p["decomposition"]) if p["decomposition"] else ""
            typer.echo(f"  hull 거리 {p['ehull']:+.3f} eV/atom → {p['stability']}" + (f" [분해: {dec}]" if p["ehull"] > 0 and dec else ""))
        typer.echo(f"  밴드갭 {p['gap']:.2f} eV (금속일 확률 {p['p_metal']:.0%}) · 밀도 {p['density']} g/cm³")
        if p["known"]:
            k = p["known"]
            typer.secho(f"  DB 에 있음 {k['material_id']}: 형성에너지 {k['ef']:+.3f}, hull {k['ehull']:.3f}, 밴드갭 {k['gap']:.2f} — DFT 값을 쓰는 것이 낫다", fg=typer.colors.CYAN)
        typer.echo("  비슷한 알려진 물질: " + ", ".join(f"{n['formula']} ({n['ef']:+.2f})" for n in p["nearest"]))
    if json_out:
        json_out.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
        typer.echo(f"저장: {json_out}")


@app.command("recommend")
def recommend_cmd(
    path: Path,
    csv_out: Annotated[Path | None, typer.Option("--csv", help="상위 후보를 CSV 로 저장")] = None,
    show: Annotated[int, typer.Option("--show", help="터미널에 보여 줄 후보 수")] = 20,
) -> None:
    """목표(원소·조건·정렬)를 주면 조성 후보를 만들어 평가하고 순위를 매긴다 (DB 는 DFT, 새 조성은 L1 예측)."""
    import yaml
    from pydantic import ValidationError

    from msl.ml.predict import ModelMissing
    from msl.recommend import Goal, recommend, save_run, to_csv

    try:
        goal = Goal.model_validate(yaml.safe_load(path.read_text(encoding="utf-8")))
    except ValidationError as exc:
        from msl.web.workbench import errors_of

        typer.secho("목표 정의 오류\n" + "\n".join(errors_of(exc)), fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    try:
        res = recommend(goal, progress=typer.echo)
    except (ModelMissing, ValueError) as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    saved = save_run(res)
    c = res["counts"]
    typer.secho(f"\n{goal.name}", bold=True)
    typer.echo(f"조건: {', '.join(res['constraints_text']) or '없음'} · 정렬: {res['objective_text']}")
    typer.echo(f"후보 {c['candidates']:,} → 전하 균형 {c['charge_balanced']:,} → 평가 {c['evaluated']:,} (DB {c['known']} · 새 조성 {c['new']}) · "
               f"충족 {c['충족']} · 가능성 있음 {c['가능성 있음']} · 불충족 {c['불충족']}\n")
    typer.echo(f"{'':3} {'화학식':14} {'출처':4} {'판정':8} {'hull':>8} {'밴드갭':>7} {'밀도':>6}  비고")
    for i, r in enumerate(res["results"][:show], 1):
        v = r["values"]
        note = "; ".join(r["flags"] + ([f"불확실: {', '.join(r['uncertain'])}"] if r["uncertain"] else []))
        typer.echo(f"{i:>3} {r['formula']:14} {r['source']:4} {r['status']:8} {v['ehull']:+8.3f} {v['gap']:7.2f} "
                   f"{(v['density'] or 0):6.2f}  {r['material_id'] or ''} {note}")
    typer.echo(f"\n저장: {saved}")
    if csv_out:
        csv_out.write_text(to_csv(res), encoding="utf-8-sig")
        typer.echo(f"CSV: {csv_out}")


@app.command("l2")
def l2_cmd(
    formula: Annotated[str, typer.Argument(help="화학식 (예: Li6MnNi3O10, LiCoO2)")],
    models: Annotated[str, typer.Option("--models", help="쉼표로 (mace-mpa-0, orb-v3)")] = "mace-mpa-0,orb-v3",
    orderings: Annotated[int, typer.Option("--orderings", help="모체마다 Ewald 순위 배치 수 (무작위 2개는 별도)")] = 4,
    max_atoms: Annotated[int, typer.Option("--max-atoms", help="후보 구조 원자 수 상한")] = 40,
    json_out: Annotated[Path | None, typer.Option("--json", help="결과 JSON")] = None,
) -> None:
    """L2 안정성 확인 — 치환 구조 → uMLIP 이완 → 자기일관 hull (두 모델). 처음 계는 경쟁 상 이완에 몇 분."""
    import json

    from msl import l2
    from msl.engines.umlip import UmlipUnavailable

    try:
        res = l2.evaluate(formula, models=tuple(models.split(",")), n_orderings=orderings, max_atoms=max_atoms, log=typer.echo)
    except (UmlipUnavailable, ValueError) as exc:
        typer.secho(str(exc), fg=typer.colors.RED, err=True)
        raise typer.Exit(1) from exc
    typer.secho(f"\n{res['formula']} — L2 hull 거리 {res['ehull_mean']:+.3f} eV/atom (두 모델 차이 {res['ehull_spread']:.3f})", bold=True)
    if res["known"]:
        typer.echo(f"  DB 에 있음 {res['known']} · MP(현재 DB) hull 거리 {res['mp_ehull']:.3f}")
    for m, r in res["models"].items():
        dec = " + ".join(f"{k} {v:.0%}" for k, v in r["decomposition"].items())
        typer.echo(f"  {m:10} {r['ehull']:+.3f} · 최저 구조 {r['best']} ({r['space_group']}) · 경쟁 상 {r['n_references']}개 · 분해 {dec}")
    ph = res.get("phonon")
    if ph and not ph.get("error"):
        typer.echo(f"  포논 (MACE): 최소 {ph['min_freq_THz']:+.2f} THz → {'동역학적 안정' if ph['dynamically_stable'] else '허수 모드 — 불안정'}")
    if res.get("l3"):
        typer.secho(f"  L3(DFT) 승격: {res['l3']['verdict']}", bold=True)
        for r in res["l3"]["reasons"]:
            typer.echo(f"    - {r}")
    typer.echo(f"  후보 구조 {res['n_structures']}개 · {res['seconds']:.0f}초")
    if json_out:
        json_out.write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
        typer.echo(f"저장: {json_out}")


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
