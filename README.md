# materials-sim-lab — 소재 조합 시뮬레이션 프레임워크 (가칭)

원소·광물·화학물질을 **조합**하고, 그 조합을 **여러 각도로 가상 시험**하는 프레임워크.
공공데이터로 통합 DB를 만들고, 검증된 오픈 계산 엔진을 연결한다.

<details>
<summary><b>English summary</b></summary>

**materials-sim-lab** is a framework for combining elements, minerals and chemicals into a *recipe*
and running that recipe through a battery of **virtual tests**. It builds an integrated database from
public data sources and wraps proven open-source engines instead of re-implementing them.
Documentation and UI are in Korean; code identifiers and CLI are in English.

**Virtual tests (11)** — every result carries its fidelity level, sources and licenses:

| Code | Test | Engine / data |
|---|---|---|
| S0 | Mixing safety gate (runs first) | Reactivity-group rules (CAMEO-style) |
| A1 | Composition stability, polymorphs | pymatgen + Materials Project |
| A2 | Solid-state / interface reactions | NASA thermo + MP hybrid |
| A3 | Gas equilibrium, adiabatic flame | Cantera (NASA) |
| A4 | Aqueous equilibrium (pH, saturation) | PHREEQC, cross-checked with Reaktoro |
| A5 | Aqueous corrosion (Pourbaix) | Experimental ΔGf + MP |
| A6 | Alloy phase equilibria | pycalphad + open TDBs (research use) |
| A7 | Property estimates | MP, L1 model |
| A8 | Blends / composites | Mixing rules, bounds |
| A9 | Supply risk and price | USGS MCS, World Bank |
| A10 | Korean chemical regulations | data.go.kr services |

**Funnel beyond the database**

- **L1** — composition-only surrogate (112 features, gradient boosting, conformal 80% intervals):
  formation energy MAE 0.07 eV/atom on unseen chemical systems.
- **L2** — universal ML interatomic potentials (MACE-MPA-0, ORB v3): candidate structures from Ewald-ranked
  substitutions and cross-chemistry prototypes, relaxed and scored on a self-consistent convex hull.
  Stability F1 **0.90** on a held-out set of 30 ternary oxides (`msl bench l2-f1`).
- **Goal-based recommendation** — enumerate charge-balanced compositions and rank them (DFT where known, L1 otherwise).
- **Next-experiment suggestion** — Bayesian optimization with BayBE over a candidate pool (`msl suggest`).

**Interfaces** — CLI (`msl`), local web workbench (`msl serve`), and an MCP server (`msl mcp`, stdio, 15 tools)
so AI assistants can search substances, validate and run recipes, and read reports.

**Quick start**

```bash
uv sync                                        # Python 3.12
cp .env.example .env                           # add your Materials Project API key
uv run msl run examples/recipes/mgo-alumina.yaml
uv run msl serve                               # http://127.0.0.1:8000
uv sync --extra l2 --extra bo                  # optional: uMLIP (PyTorch) and BayBE
```

Status: stages 0–4 of the plan are complete; stage 5 (interfaces) is in progress.
Raw data (`data/`) and API keys are not part of the repository; each source's license is recorded in
`kb/sources.yaml`. Code is Apache-2.0 — see [LICENSE](LICENSE) and [NOTICE](NOTICE).

</details>

---

## 현황

| | |
|---|---|
| 단계 | 0~4단계 완료 (4단계 완료 기준: L2 안정성 판정 F1 0.90, D24) · **5단계(인터페이스) 진행** — MCP 서버 완료, 사용성 테스트 남음 · [3단계 보고](docs/07-3단계-평형-트랙.md) |
| 검증 | 검증 세트 40개 사례 모두 통과 (S0 위험 재현율 100%) — `msl validate` · PHREEQC 공식 예제 28/28 재현 — `msl bench phreeqc` |
| 실제 계산 시험 | **11종 모두** — S0 안전 · A1 안정성(광물 다형) · A2 고상반응(고온 하이브리드) · A3 기체 · A4 수용액 · A5 수계 부식(Pourbaix) · A6 합금 · A7 물성 · A8 블렌드 · A9 공급 · A10 국내 규제 |
| 조합 공간 | 혼합비 격자 · 양 스윕 · 부분집합 열거 → 표·곡선·쌍별 행렬 — `msl space`, 웹 `/spaces` |
| 웹 작업대 | 성분 검색(국문 관용명 포함)·시약 선반·양/조건/시험 편집·YAML 편집 → 실행·저장, 조합 공간 편집기 — [사용 안내](docs/08-사용-안내.md) |
| 물성 예측 (L1) | 화학식만으로 형성에너지(MAE 0.07 eV/atom, 80% 구간)·안정성(판정 73%)·밴드갭·밀도 — `msl predict`, 웹 `/predict`, MP 에 없는 성분은 A1·A7 이 자동 사용 |
| 목표 기반 추천 | 원소·조건(hull 거리·밴드갭·밀도 등)·정렬을 주면 조성 후보를 전수 평가해 순위 (DB 는 DFT, 새 조성은 L1) — `msl recommend`, 웹 `/recommend` |
| 다음에 잴 조성 | 잰 값(L2·DFT·실험)으로 다음 조성 제안 — BayBE, 되돌아보기에서 무작위의 약 2.7배 — `msl suggest` |
| MCP 서버 | 다른 AI 도구에서 성분 검색·레시피 검증/실행·리포트 조회·예측·추천 (도구 15개, stdio) — `msl mcp`, [사용 안내 9절](docs/08-사용-안내.md) |
| L2 안정성 확인 | uMLIP(MACE-MPA-0·ORB v3) — Ewald 배치·구조 원형으로 후보 구조 생성 → 이완 → 자기일관 hull, 두 모델 평균·차이. 재발견 시험 7/7 — `msl l2`, 웹 [L2 확인] |
| 적재 데이터 | 15개 소스 · 277,479행 — MP 163k · COD 21k · IMA 6,239종(MP 다형 매칭) · USGS · 가격 · NASA · PHREEQC · CAMEO · 합금 TDB 2,774 · 수용액 이온 362 |
| 해석기 | 원소·화학식·IMA 광물명·CAS·영문/국문 이름·국문 관용명·KE 번호 — 99.5% |
| 테스트 | 404개 통과 |

### 실행

```bash
uv sync                                               # 환경 설치 (Python 3.12)
uv run msl db load all                                # 데이터 적재 (원본은 msl db fetch 로 받음)
uv run msl db status                                  # 적재 현황
uv sync --extra bo && uv run msl suggest next examples/campaigns/li-mn-o-cathode.yaml   # 다음에 잴 조성 (BayBE)
uv run msl mcp                                        # MCP 서버 (stdio) — claude mcp add msl -- uv --directory <폴더> run msl mcp
uv run msl serve                                      # 웹 작업대 → http://127.0.0.1:8000 (레시피 만들기·고치기·실행·저장)
uv run msl run examples/recipes/mgo-alumina.yaml      # 터미널에서 레시피 하나 실행
uv run msl run examples/recipes/cu-ni-alloy.yaml --report out.html   # 리포트 파일 (.html/.md)
uv run msl validate                                   # 검증 세트 채점
uv run msl ml train                                   # L1 조성 모델 학습 (처음 한 번, 약 15분)
uv run msl predict Li1.2Ni0.6Mn0.2O2 LiFePO4          # 조성 → 물성 예측 (웹: /predict)
uv run msl recommend examples/goals/li-mn-o-cathode.yaml   # 목표 기반 추천 (웹: /recommend)
uv sync --extra l2 && uv run msl l2 Li6MnNi3O10         # L2 안정성 확인 (uMLIP, 처음 원소계는 몇 분)
uv run msl space run examples/spaces/cu-ni-composition.yaml --report out.html   # 조합 공간 (웹: /spaces)
uv run msl bench phreeqc --dist <배포본> --binary <phreeqc>   # PHREEQC 공식 예제 재현
uv run msl resolve-check                              # 해석기 시험 목록 (기준 95%)
uv run pytest                                         # 테스트
```

`.env`에 `MP_API_KEY`(A1·A2·A7)와 `DATA_GO_KR_SERVICE_KEY`(A10, 국문명 해석)가 필요합니다. 없으면 '대기'로 표시됩니다.
IMA PDF 적재에는 시스템의 `pdftotext`(poppler)가 필요합니다. USGS MCS 원본은 브라우저로 받아 둡니다([D12](docs/03-결정-기록.md)).

### 데모 레시피 (`examples/recipes/`)

| 레시피 | 실제로 계산되는 것 |
|---|---|
| 마그네시아 + 알루미나 | MgAl₂O₄(스피넬) 생성 반응에너지 −54 meV/atom(0 K) — 실험값(약 −52)과 일치 |
| 석회석 + 석영 | 하이브리드 고온 열역학(NASA + MP 앵커): 1,273 K 에서 −127 meV/atom, CaSiO₃ + CO₂ |
| 표백제 + 염산 | S0 **부적합** — 염소(Cl₂) 발생. A10: 염화수소 사고대비물질·인체등유해성물질, GHS '위험'. PHREEQC 기본 DB 는 차아염소산을 다루지 못한다고 명시 |
| 방해석 + 물 | pH 8.21, Ca 0.54 mmol/kgw (대기 CO₂ 평형) |
| 메탄 + 공기 | 단열 평형 온도 2,225 K |
| 에폭시 + 유리섬유 | 밀도 1.61 g/cm³, 영률 경계 4.9~24.1 GPa |
| Li–Co–O 계 | 안정 화합물 11개 + L2 재확인 후보 33개(표에는 가까운 8개), MP 물성 |
| Cu–Ni 합금 | 고상선 1,193.6 °C · 액상선 1,239.6 °C · Scheil 응고 종료 (A6, 연구용 TDB) |
| 철 + 물 (공기) | A5 Pourbaix 도표 — 공기 포화 전위에서 Fe₂O₃ 부동태, pH 5.61(대기 CO₂) |
| 알루미늄 + NaOH 수용액 | S0 **부적합**(수소 발생) · A5 pH 12.88 에서 Al(OH)₄⁻ 로 부식 |

---

## 핵심 아이디어

기존 도구는 도메인별로 끊겨 있다.

| 도메인 | 도구 |
|---|---|
| 결정 | Materials Project |
| 수용액 | PHREEQC |
| 합금 | CALPHAD |
| 혼합위험 | CAMEO |

**계산 엔진은 새로 만들지 않는다.** 대신 그 사이를 잇는 계층을 만든다.
- 공통 조합 레시피
- 도메인 라우터
- 안전 게이트
- 출처 추적

조합 하나를 넣으면 적용할 수 있는 시험이 자동으로 골라져 실행된다. 아래는 `msl run` 실제 출력을 줄인 것이다.

```
마그네시아 + 알루미나 → 스피넬  [rcp-mgo-alumina]
● S0   혼합 안전 게이트  [L0] 반응성 그룹 정보가 있는 성분이 2개 미만
● A1   조성 안정성     [L0] 반응물 안정성 — MgO 0 meV/atom, Al2O3 0 meV/atom
● A2   고상·계면 반응   [L0] 1700 K (Gibbs 근사): 0.5 MgO + 0.5 Al2O3 -> 0.5 MgAl2O4 — 반응에너지 -28 meV/atom
        최저 반응에너지 · 0 K (DFT): -54.4 meV/atom
… A9   공급·경제성     [L0] 1단계에서 USGS MCS 2026·World Bank 가격 자료를 적재한 뒤 실행
```

시험 종류(S0, A1~A10)와 충실도 계층(L0~L3, 평형 트랙)은 [계획서 3장](docs/02-계획서.md)에 있다.

---

## 문서

| | |
|---|---|
| [00 조사 · 공공데이터 소스](docs/00-조사/01-공공데이터-소스.md) | 약 70개 소스, 라이선스, MVP 적재 Top 10, 국내 공공데이터 |
| [00 조사 · 시뮬레이션 엔진](docs/00-조사/02-시뮬레이션-엔진.md) | 충실도 계층, uMLIP 현황(Matbench Discovery), 조합 유형별 엔진 매핑 |
| [00 조사 · 플랫폼 벤치마킹](docs/00-조사/03-플랫폼-벤치마킹.md) | 비교 매트릭스, 빈틈 분석, 차용 Top 10, 포지셔닝 |
| [01 벤치마킹 종합](docs/01-벤치마킹-종합.md) | 조사 결론을 우리 시스템 관점에서 정리 (만들 것·래핑할 것·참고만 할 것) |
| [02 계획서](docs/02-계획서.md) | 목표·설계·DB 구축·아키텍처·로드맵·검증·리스크·결정 사항 |
| [03 결정 기록](docs/03-결정-기록.md) | 확정된 결정 D1~D16과 파생 제약 |
| [04 기관 문의 초안](docs/04-기관-문의-초안.md) | NOAA·KOSHA·국가소재연구데이터센터·KRISS 문의문 |
| [05 1단계 적재 보고](docs/05-1단계-데이터-코어.md) | 소스별 적재, 완료 기준 점검, 확인된 사실 |
| [06 2단계 결과 보고](docs/06-2단계-조합-시험.md) | 검증 세트, A2 하이브리드, A6 합금, 광물 다형, 리포트 |
| [07 3단계 결과 보고](docs/07-3단계-평형-트랙.md) | A5 Pourbaix, 조합 공간 생성기, PHREEQC 공식 예제 재현 |
| [08 사용 안내](docs/08-사용-안내.md) | 만들어진 것 목록, 웹·터미널 사용법, 레시피·조합 공간 작성법 |
| [화면 목업](prototype/virtual-lab-mockup.html) | 교육용 첫 화면 목업 (실제 화면은 `msl serve`) |

---

## 다음

1. GPU(MPS) 환경에서 uMLIP 처리량 재측정 — 4단계 마지막 항목
2. SGTE 부록 라이선스 확인 → A6 배포 가능 여부
3. 고온·고압 수계 DB 라이선스 정리 (SUPCRT98/07 확인 또는 SUPCRTBL 연구 전용) — D22
