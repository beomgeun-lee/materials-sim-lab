---
name: materials-sim-lab
description: 물질 조합(레시피)의 안전성·수용액 평형·부식·안정성·물성·공급·국내 규제를 가상으로 계산할 때. materials-sim-lab MCP 서버(msl) 도구로 성분 검색, 레시피 검증·실행, 리포트 조회, 조성 물성 예측, 목표 기반 조성 추천을 한다.
---

# materials-sim-lab 사용법

MCP 서버 `msl`(stdio)이 연결돼 있어야 한다. 등록:

```bash
claude mcp add msl -- uv --directory <materials-sim-lab 폴더> run msl mcp
```

## 기본 흐름

1. **성분 확인** — `search_substance("염산")` → `ref`(예: `cas:7647-01-0`)·화학식·해석 경로.
   국문 관용명·CAS·화학식·광물명·원소 기호를 받는다. `error` 가 있으면 사용자에게 다른 이름을 물어본다.
2. **레시피 작성** — 형식이 헷갈리면 `list_recipes` → `get_recipe("calcite-water.yaml")` 로 예제를 본다.
3. **검증** — `check_recipe(yaml)` → 오류는 한국어 줄 단위. ok 면 자동으로 고를 시험(`assays`)과 이유가 나온다.
4. **실행** — `run_recipe(recipe_yaml=...)` 또는 `run_recipe(file=...)`. 첫 줄의 `report_id` 를 기억해 두고
   나중에 `get_report(report_id)` 로 다시 연다. 구조화된 값이 필요하면 `fmt="json"`.
5. **저장**(사용자가 원할 때만) — `save_recipe(yaml)`. 같은 id 가 있으면 사용자 확인 뒤 `overwrite=True`.

## 레시피 뼈대

```yaml
id: rcp-my-test            # rcp- 로 시작, 영소문자·숫자·하이픈
name: 구리 + 묽은 염산
components:
  - {ref: element:Cu, amount: 5 g, state: solid}          # solid · liquid · gas · aqueous
  - {ref: cas:7647-01-0, amount: 0.1 mol, state: aqueous}
  - {ref: cas:7732-18-5, amount: 1 L, state: liquid}
conditions: {T: 25 °C, atmosphere: air}                   # 선택: P, pH, Eh(V vs SHE), time
assays: auto                                              # 또는 [S0, A4, A5]
```

양 단위: g·kg·mg·mol·mmol·L·mL·mol/L·mol/kg·wt%·at%·mol%·vol%. 원소계 탐색은 `mode: system`(양 없이 원소만).

## 결과 읽는 법 — 사용자에게 전할 때 지킬 것

- **S0(혼합 안전 게이트)가 위험을 찾으면 가장 먼저, 분명하게 알린다.** 계산 결과보다 우선한다.
- 결과마다 **충실도**가 있다: L0 = DB 조회(DFT·실험값), L1 = 조성 모델 예측(구간 포함), L2 = uMLIP, T = 열역학 평형 계산.
  L1 예측은 "예측"이라고 말하고 80% 구간을 함께 전한다.
- `caveats`(주의)를 빼지 말고 전한다 — 예: phreeqc.dat 는 대략 0–100 °C 검증, PHREEQC·Reaktoro pH 차이 경고.
- `not-applicable` 은 실패가 아니라 "이 조합에 맞지 않는 시험"이다. 이유를 짧게 전한다.
- 출처·라이선스가 결과에 붙는다. restricted(연구용) 데이터가 쓰인 값은 상업 용도로 쓰면 안 된다고 알린다.
- 평형 계산은 반응 속도·점화 여부를 모른다. "일어날 수 있다"와 "일어난다"를 구분해서 말한다.

## 그 밖의 도구

| 도구 | 쓰임 |
|---|---|
| `list_assays` | 시험 11종의 질문·적용 대상 |
| `predict_properties([화학식…])` | L1 — 형성에너지·hull 거리·안정성·밴드갭·밀도 (DB 에 있으면 DFT 값도). 50개까지 |
| `recommend(goal_yaml)` | 원소·조건·정렬을 주면 조성 후보 순위. 예: `required: [Li, Mn, O]`, `constraints: [{prop: ehull, max: 0.05}]` |
| `l2_result(화학식)` | 이미 계산한 uMLIP 안정성·포논·DFT 승격 판단. 새 계산은 사용자가 웹이나 `msl l2` 로 돌린다 |
| `list_reports` | 저장된 리포트 목록 |
