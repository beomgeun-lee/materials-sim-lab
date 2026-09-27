# 04. 공개 CALPHAD 열역학 DB(TDB) — 라이선스·커버리지

- 작성일: 2026-09-27 (모든 사실은 이 날짜 기준 웹·파일 확인)
- 범위: A6 합금 상평형(pycalphad)에 쓸 수 있는 공개 TDB. 계획서 10절 #6 "공개 TDB가 커버하는 합금계"에 답한다. 수계·기체 DB는 01 문서 6.5·6.7절.
- 방법: 배포 페이지와 TDB 파일 머리 주석에서 라이선스 문구를 직접 인용했다. 파일은 받아서 pycalphad 0.11.2로 읽고 평형을 계산해 봤다. Zenodo·GitHub 부가 후보(4.3절)는 저장소 API 메타데이터로 확인했고, 파일 내용까지 모두 열어 보지는 않았다. 라이선스를 찾지 못한 것은 "미확인"으로 적었다.
- 결정 규칙: D2(연구 단계는 비상업·미확인도 적재하되 `restricted`). `kb/sources.yaml` 에는 확인되지 않은 라이선스를 `unknown` 으로 적는다. 법률 자문이 아니다.

---

## 0. 요약

### 0.1 핵심 결론

1. **Cu–Ni 는 된다.** 평가된 Cu–Ni(an Mey 1992, Calphad 16 255)가 두 공개 DB에 있다. 하나는 SGTE 2원계 모음(Hallstedt 2025)의 `CuNi-92Mey-LB.tdb`, 다른 하나는 COST 507 경합금 DB다. 데모 레시피(Cu 70·Ni 30 at%)의 평형 고상선은 **1,193.6 °C**, 액상선은 **1,239.6 °C** 다. C71500 용융 범위 1,170–1,240 °C 와 비교하면 고상선 +23.6 K, 액상선 −0.4 K 로 둘 다 ±30 K 안이다(5절).
2. **가장 넓은 소스는 SGTE 2원계 모음이다.** Calphad 89 (2025) 102833 논문 부록에 평가된 2원계 TDB 2,773개가 있다(1,241계). pycalphad 로 2,754개가 읽힌다. Al–Cu, Fe–Cr, Al–Zn, Cu–Zn, Fe–Ni, Cu–Ni 모두 들어 있고, 계마다 Landolt-Börnstein 수록본(-LB)이 있다.
3. **적재한 TDB는 모두 라이선스 미확인(`unknown`) → `restricted` 다.**
   - SGTE 모음은 논문이 CC BY 4.0 이지만, 부록 zip 과 TDB 파일에는 문구가 없다. CC BY 가 부록 데이터에도 적용된다는 명시를 확인하지 못했다(본문 페이지 403).
   - COST 507 은 "It may be used freely but on your own risk" 뿐이고, 상업 이용·재배포는 명시되지 않았다.
   - 그래서 지금 A6 결과는 모두 연구용이다. 결과 caveats 에 '연구용 TDB' 가 붙는다.
4. **재배포 가능(open) 라이선스가 명시된 후보는 둘이지만 바로는 못 쓴다.**
   - MatCalc 공개 DB(ODbL 1.0)는 평가 DB이고 Cu–Ni·Fe–Ni·Fe–Cr·Al–Cu 를 포함한다. 그러나 MatCalc 전용 명령이 있어 pycalphad 가 그대로 읽지 못한다(변환기 필요).
   - LibreCalphad `mf-steel.tdb`(MIT, 파일에 직접 명시)는 pycalphad 로 읽힌다. 그러나 평형 결과가 물리적으로 틀린다: Cu–Ni 에서 FCC_L12 가 2,100 K 에서도 남는다.
   - 둘 다 적재를 보류했다(`status: planned`).
5. **쓰면 안 되는 것도 분명하다.**
   - NIMS CPDDB 는 재배포가 금지된다. pycalphad 동봉 테스트 DB 중 9개가 그 사본이다.
   - NIST Solder 는 data.nist.gov 에 SRD 139 로 등록돼 있다.
   - OpenCalphad GitHub 예제는 상용 SSOL2 에서 발췌한 것이다.
   - Equilipy 의 AlCuMgSi 파일은 "must not be distributed" 다.

### 0.2 비교표

범례: 재배포 ○ 가능 / △ 추정·미확인 / × 불가. 판정은 이번 단계의 적재 여부다.

| 소스 | 계(원소) | 출처 | 라이선스 근거(원문) | 재배포 | 품질 | 판정 · sources.yaml id |
|---|---|---|---|---|---|---|
| SGTE 2원계 모음 (Hallstedt 2025) | 2원계 1,241계 (2,773파일) | Calphad 89 (2025) 102833 부록 mmc3.zip [T1] | 논문 CC BY 4.0 (출판사 메타데이터 `openaccessUserLicense`, Crossref `vor`) [T2]. 파일·zip 에는 문구 없음 | △ | 평가 (SGTE 검수, -LB = Landolt-Börnstein 수록) | **적재** · `sgte-binary-collection` (unknown → restricted) |
| COST 507 Round II (1999) | 29원소 경합금, 평가 2원계·3원계 | EUR 18499, OpenCalphad 배포 [T3] | 파일: "It may be used freely but on your own risk". 보고서: "Reproduction is authorised provided the source is acknowledged" | △ | 평가 (동결판) | **적재** · `cost507` (unknown → restricted) |
| MatCalc 공개 DB (mc_fe·mc_ni·mc_al) | 25·22·12원소 | TU Wien [T4] | 파일 머리: "made available under the Open Database License … Database Contents License" | ○ (동일조건) | 평가 (공학용 다원계) | 보류 — pycalphad 변환 필요 · `matcalc-open-tdb` |
| LibreCalphad `mf-steel.tdb` | 약 90원소 (2원계는 일부) | M. Frichtl, GitHub [T5] | 파일 머리에 MIT 허락 문구 전문 | ○ | 문헌 파라미터 모음, 진행 중 | 보류 — 평형 결과 이상 · `librecalphad` |
| NIST Solder (SRD 139) | Ag·Bi·Cu·Pb·Sb·Sn | U.R. Kattner, 2017-08-04 [T6] | 페이지·파일 문구 없음. data.nist.gov 에 `nrd:SRD` 로 등록 (doi:10.18434/T4759N) | × (SRD Act) | 평가+외삽 | 보류 · `nist-solder-tdb` (NIST-SRD) |
| SGTE Unary 5.0 | 순원소 87 | SGTE [T7] | "can only be used for extracting data for assessment work or to tabulate or plot data for the pure elements" | × (용도 제한) | 평가 | 적재 안 함 (2원계 파일에 순원소 값 포함) · `sgte-unary` |
| NIMS CPDDB | 700+ 계 | NIMS [T8] | "No reproduction, republication or distribution to third parties of any content is permitted without written permission of NIMS." | × | 평가 | 제외 (내부 참조) · `nims-cpddb` |
| TDBDB (Brown대) | 색인 약 811건 | van de Walle 외 [T9] | 링크만 제공. "you agree not to copy over our database's content" | — | 색인 | 적재 대상 아님 · `tdbdb` |
| pycalphad 동봉 테스트 DB | 파일마다 다름 (TDB 38개) | pycalphad 0.11.2 `tests/databases` [T10] | 저장소 MIT 는 "Software" 대상. 파일마다 별도 표기 (4.1절) | 파일별 | 대부분 테스트용 | 테스트에만 사용 (COST507 수정본은 계산용) |
| kawin `CuNi.tdb` | Cu–Ni (FCC·LIQUID) | kawin 예제 [T11] | 저장소 MIT (데이터 파일 명시 없음) | △ | an Mey 1992 사본 | 후보 (open 대체안) |

---

## 1. 적재한 TDB

### 1.1 SGTE 2원계 데이터셋 모음 (Hallstedt 2025) — `sgte-binary-collection`

- **출처**: B. Hallstedt, "The SGTE collection of binary datasets", *Calphad* 89 (2025) 102833, doi:10.1016/j.calphad.2025.102833.
  - 부록 `mmc3.zip`(7.6 MB, 2025-05 온라인)에 TDB 2,773개가 있다 [T1]. SHA-256 은 `6a0d53de…747d` 다.
  - 같은 파일을 pycalphad 팀이 phasediagrams.org 에도 게시한다(2,833파일·1,258계, 라이선스 표기 없음) [T12].
- **라이선스**
  - 논문은 CC BY 4.0 이다. Elsevier API `openaccessUserLicense` 와 Crossref `license[vor]` 가 모두 `creativecommons.org/licenses/by/4.0/` 이다 [T2].
  - 부록 zip 과 파일 2,773개에는 저작권·라이선스 문구가 한 줄도 없다(전수 검색).
  - 출판사 공개 라이선스 페이지에도 부록 데이터에 관한 문장이 없다. 논문 본문(데이터 제공 문구)은 403 으로 읽지 못했다.
  - → **`unknown`**. 부록에 CC BY 가 적용됨이 확인되면 `CC-BY-4.0` 으로 올린다. 저자(SGTE)에게 문의할 대상이다.
- **품질**
  - 파일마다 원 평가 논문, 오식 수정 메모, "Checked against paper and author's tdb" 확인 기록이 머리 주석에 있다.
  - `-LB` 는 Landolt-Börnstein 수록본, `-3g` 는 3세대 순원소 기술이다.
- **pycalphad 호환**: 2,754개가 읽힌다. 못 읽는 19개는 `C_S`·`STATUS`·`COMP_SETS`·`NEVER_DIS_PAR` 같은 비표준 TYPE_DEF, 규칙상 정의 누락, 문법 오류 때문이다. 이 파일들은 `parse_ok=false` 로 남기고 A6 은 건너뛴다.
- **대표 계와 대표 평가(preferred)**

| 계 | 파일 수 | preferred |
|---|---|---|
| Cu–Ni | 4 | CuNi-92Mey-LB |
| Al–Cu | 5 | AlCu-98Sau-LB |
| Fe–Cr | 7 | CrFe-87And-LB |
| Al–Zn | 5 | AlZn-93Mey-LB |
| Cu–Zn | 6 | CuZn-93Kow-LB |
| Fe–Ni | 12 | FeNi-03Dup-LB |
| Al–Mg · Al–Si · Cu–Sn · Ag–Cu · Cr–Ni · Ni–Ti · Fe–Mn | 6–13 | 모두 -LB 있음 |

### 1.2 COST 507 경합금 DB — `cost507`

- **출처**: COST 507 최종 DB, Round II(1999-01, 동결판). I. Ansara, A.T. Dinsdale, M.H. Rand 편, EUR 18499 [T3].
  - OpenCalphad 가 `COST507.zip`(2017-04-21) 과 보고서 PDF 를 배포한다.
  - 29원소다. 평가 2원계는 Al–Zn·Cu–Ni·Cu–Si·Cu–Zn·Mg–Ni·Al–Cu·Al–Fe·Cu–Fe 등이다. 3원계는 Al–Cu–Mg·Al–Mg–Si·Al–Mg–Zn·Al–Si–Zn·Cu–Mg–Zn 등 25개다.
  - **Fe–Ni 는 없다.**
- **라이선스**
  - 파일 `DATABASE_INFO`: "It may be used freely but on your own risk, no one takes responsibilty that the calculated results are correct, even for the assessed systems." 같은 곳에 "Commercial versions of this database can be obtained from various organisations" 도 있다.
  - 보고서 PDF: "© European Communities, 1998 / Reproduction is authorised provided the source is acknowledged."
  - 상업 이용·재배포가 명시되지 않았으므로 → **`unknown`**.
- **pycalphad 호환**
  - OpenCalphad 원본은 pycalphad 가 못 읽는다. 문법 오류가 있기 때문이다: Y 원자량 `8.89059+01`, `REF: 0`(공백), 중복 `FUNCTION R`, `TYPE_DEF R IF(…)`, `!.0E-4`.
  - pycalphad 팀이 고친 사본이 `tests/databases/COST507.tdb` 로 동봉돼 있다. 원본 zip 과 이 사본을 함께 스냅샷하고, 계산에는 사본을 쓴다.
- **주의**
  - Hallstedt 메모에 따르면 Cu–Ni 의 `BM(FCC_A1,CU,NI)` 에 오류 가능성이 있다. 자기 항이라 고온 결과에는 영향이 없다: 고상선·액상선은 SGTE 파일과 0.1 K 안에서 같다.
  - A6 순위 규칙상 2원계는 SGTE 모음이 먼저 쓰인다. COST 507 은 3원계 경합금에 쓰인다.

## 2. 보류한 후보 (라이선스는 명시, 기술 문제)

### 2.1 MatCalc 공개 DB — `matcalc-open-tdb`

- **파일**: mc_fe 2.062(2024-11-08), mc_ni 2.036(2024-08-20), mc_al 2.037(2025 개정) [T4]. 확산(.ddb)·물리(.pdb) DB 도 있다. mc_cu·mc_mg 공개판은 없다.
- **라이선스**
  - 파일 머리: "This database … is made available under the Open Database License: http://opendatacommons.org/licenses/odbl/1.0/. Any rights in individual contents of the database are licensed under the Database Contents License".
  - 상업 이용·재배포가 가능하다. 공개 이용하는 파생 DB 는 같은 라이선스여야 한다(ODbL 4.4).
  - `licenses.yaml` 에 ODbL id 가 없어 성격이 같은 `CC-BY-SA-4.0` 으로 매핑했다. **ODbL-1.0 id 추가가 필요하다.**
- **커버리지**: 세 파일 모두 Cu–Ni 가 있다. 액상은 Jansson 1987, fcc 는 an Mey 값이다. mc_fe 에는 Fe–Ni·Fe–Cr 도 있다.
- **막힌 곳**: pycalphad 0.11.2 가 그대로 읽지 못한다.
  - 파일이 latin-1 인코딩이다.
  - MatCalc 전용 명령이 있다: `REFERENCE_ELEMENT`, `ATTACH_CONTRIBUTION … ORDER_DISORDER`, `ADD_COMPOSITION_SET`.
  - `HMVA` 파라미터가 있다.
  - 규칙-불규칙 기여를 `TYPE_DEF … DIS_PART` 로 바꾸는 변환기가 필요하다. pycalphad 예제의 `mc_fe_v2.059.pycalphad.tdb` 가 변환 선례다.
  - **open 파티션 A6 의 1순위 후속 과제**다.

### 2.2 LibreCalphad `mf-steel.tdb` — `librecalphad`

- **라이선스**: 파일 머리에 "This file is part of LibreCalphad. Copyright (c) 2024 Matthew Frichtl / Permission is hereby granted, free of charge…"(MIT 전문)가 있다 [T5]. 커밋 `a7dace2`(2024-12-23)로 고정해 받았다.
- **pycalphad 로 계산해 본 결과**
  - `FCC2_L10` 상은 모델 생성이 실패한다(규칙-불규칙 부격자 비 불일치).
  - 이 상을 빼고 계산하면, Cu–Ni 70/30 에서 FCC_L12(규칙상)가 2,100 K 에서도 남는다. 순 Ni 융점은 1,728 K 다. Fe–Ni 90/10 과 Fe–Cr 80/20 도 1,900 K 에서 고체다.
  - 같은 DB 에서 `LIQUID`+`FCC_A1` 만 넣으면 Cu–Ni 가 정상이다. 따라서 규칙상 모델 조합 문제로 추정한다.
- **파일 대조로 확인한 Cu–Ni 차이 (an Mey 원문 대비)**
  - fcc 자기 상호작용이 순 Ni 값으로 들어가 있다: `TC +633`, `BMAGN +0.52`. 원문은 −935.5/−594.9, −0.7316/−0.3174 다.
  - 액상 `L0` 의 T 계수가 1.29093 이다. 원문은 1.29893 이다.
- **판단**: 적재하면 A6 이 강(Fe–Cr–Ni)에 이 DB 를 골라 틀린 값을 낸다. 그래서 커넥터를 만들었다가 뺐다. pycalphad 호환이 해결되면 다시 검토한다.

## 3. 쓰지 않는 소스

- **NIST Solder (SRD 139)** [T6]
  - Ag·Bi·Cu·Pb·Sb·Sn 이다(Ni 없음). `NIST-solder.tdb` 는 2017-08-04 갱신판이다.
  - 페이지와 파일에 저작권 문구는 없다. 그러나 data.nist.gov 레코드가 `@type: nrd:SRD`, 제목 "… Solder Systems - SRD 139" 다.
  - NIST 저작권 페이지는 SRD 에 "the Secretary of Commerce to secure copyright on behalf of the United States in Standard Reference Data" 를 적용한다. 비-SRD 데이터는 17 U.S.C. §105 로 미국 내 저작권이 없다고 적는다.
  - 1차 조사(01 문서 6.6)의 '퍼블릭 도메인 추정'을 바꿔 **NIST-SRD** 로 둔다.
- **SGTE Unary 5.0**: 다운로드 페이지에 용도 제한 문구가 있다(0.2 표). 순원소 값은 적재한 2원계 파일에 이미 들어 있다.
- **NIMS CPDDB**
  - 로그인이 필요하고 대량 다운로드가 금지된다. 재배포 금지 원문은 0.2 표에 있다.
  - TDBDB 로 찾은 Cu–Ni 파일(cuni_mey·cuni_mey2·cuni_sri·cuni_sha)도 여기 있다.
- **TDBDB**
  - 현재 주소는 https://avdwgroup.engin.brown.edu/ 다. 약 811건이다(NIMS 412, Calphad 부록 321, NIST 72).
  - "This site provides links to the data files and acticles [sic] rather than data/article itself (for copyright and fairness reasons)." 최신 수록 논문이 2019년이다.
  - Calphad 저널 부록 TDB 는 Crossref 에 Elsevier TDM 약관만 있다.

## 4. 기타 후보 조사 결과

### 4.1 pycalphad 동봉 테스트 DB (0.11.2, `tests/databases`)

- **저장소 LICENSE**: MIT 다. "this software and associated documentation files" 만 대상이고, DB 파일 README 는 없다.
- **파일별 표기**
  - **NIMS 사본 9개**: `Copyright (C) NIMS 20xx`. 예: cumg, alzn_mey, pbsn, alfe, cuo, crtiv_ghosh, nbre_liu, Al-Mg_Zhong, al2o3_nd2o3_zro2.
  - **테스트 전용 표기**: "FOR TESTING PURPOSES ONLY -- NOT FOR RESEARCH". alcocrni, femn, diffusion, FeNi_deep_branching 이다.
  - **Hallstedt 데이터셋**: AuSn-13Don, CoV-20Wan 이다. 1.1절 모음과 같은 계열이다.
  - **MatCalc 파생**: `mc_fecocrnbti.tdb` 에 "Based on Matcalc steel database 2.060; used under ODBL" 이 있다.
  - **논문 기반**: Al-Fe_sundman2009(Acta Mater. 57 2896), femns(Dilner 2015, Calphad 48 95), alfeo(Lindwall 2015) 등이다.
- **판단**: 적재하지 않는다. `tests/test_alloy.py` 는 COST507(수정본)로 A6 로직을 네트워크 없이 검사한다.

### 4.2 Cu–Ni 가 들어 있는 다른 파일

- **kawin** `examples/mobility_fitting/databases/CuNi.tdb` [T11]
  - an Mey 1992 사본이고 FCC·LIQUID 만 있다.
  - 저장소는 MIT 지만 데이터 파일을 명시하지 않는다. open 대체안으로 둔다.
- **LibreCalphad**: 2.2절.
- **MatCalc**: 2.1절.

### 4.3 Zenodo·GitHub 부가 후보 (API 메타데이터 기준)

| 후보 | 계 | 라이선스 표기 | 비고 |
|---|---|---|---|
| Zenodo 15266690 (Shi·Wei·LLorca 2025) | Ni–Co | CC BY 4.0 | 제1원리+CALPHAD 연구 모델 |
| Zenodo 22161809 (Kashyrina 2025) | Al–Ni–Sn | CC BY-SA 4.0 | 문헌 2원계 조합 |
| Zenodo 20848011 · 20181685 | Mo–Si–Y–Hf · Co–Se | CC BY 4.0 | 2026 평가, 목표 계 아님 |
| Zenodo 7767663 `feni_cac.tdb` | Fe–Ni | CC BY 4.0 표기 | 파일 머리 "Copyright (C) NIMS 2012" → 충돌, 쓰지 않음 |
| Zenodo 11198868 · 10568692 | Cr–Cu–Fe–Ni · Cr–Fe–Mo–Nb–Ni | CC BY 4.0 표기 | 파일 접근 제한 |
| ESPEI 문서 `Cr-Ni_mcmc.tdb` | Cr–Ni | MIT | 데모 피팅 |
| OpenCalphad GitHub `FENI.TDB`·`crfe.TDB` | Fe–Ni·Cr–Fe | LICENSE 파일 없음 (GPL-v3 PDF만) | "From database: SSOL2"(상용) → 제외 |
| Equilipy `AlCuMgSi_SK.tdb` | Al–Cu–Mg–Si | 저장소 BSD-3 | "private collection … must not be distributed without permission" → 제외 |

- Materials Commons 공개 데이터셋 2,345개 메타데이터에서 TDB 를 찾지 못했다. Figshare 에서도 결과가 없었다.

## 5. A6 검증 (pycalphad 0.11.2 + scheil 0.3.0, 노트북)

| 레시피 | TDB | 평형 고상선 | 평형 액상선 | 문헌 | 차이 |
|---|---|---|---|---|---|
| **Cu 70 · Ni 30 at%** (데모) | CuNi-92Mey-LB | 1,466.8 K (1,193.6 °C) | 1,512.8 K (1,239.6 °C) | C71500 1,170–1,240 °C | +23.6 K / −0.4 K |
| Cu 70 · Ni 30 wt% (= Ni 31.7 at%) | 〃 | 1,199.5 °C | 1,246.7 °C | 〃 (C71500 은 wt% 기준, Fe·Mn 소량 포함) | +29.5 K / +6.7 K |
| Cu 70 · Zn 30 wt% | CuZn-93Kow-LB | 918 °C | 949 °C | C26000 915–955 °C | +3 K / −6 K |
| Pb 38.1 · Sn 61.9 wt% | PbSn-95Oht-LB | 183.0 °C | 183.4 °C | 공정 183 °C | 0 K |
| Al 87.4 · Si 12.6 wt% | AlSi-97Feu-LB | 577 °C | 579 °C | 공정 577 °C | 0 K |
| Ag 72 · Cu 28 wt% | AgCu-04Wit-LB | 780 °C | 780 °C | 공정 779–780 °C | ≤1 K |

- **계산 시간**
  - Cu–Ni 데모: 31점 스윕 + 고상선·액상선 세분 + 상온 1점 + Scheil 을 합쳐 약 1.0 s. 레시피 전체 실행은 3.1 s.
  - 3원계(COST507 Al–Cu–Mg): 약 5 s.
  - 적재(`msl db load sgte-binary-collection`): 파일 2,773개를 pycalphad 로 읽어 색인하는 데 약 52 s.
- **데모 레시피의 다른 결과**
  - 상온(298 K) 평형은 FCC_A1 두 조성 세트(혼화 갭)다. 실제 재료는 확산이 느려 단상으로 남는다.
  - Scheil 응고 종료 온도는 1,359.8 K 다. 고상 확산이 없다고 본 극단값이라 순 Cu 융점 근처까지 내려간다.

## 6. A6 선택 규칙 (코드: `src/msl/assays/alloy.py`)

1. `tdb_systems` 에서 레시피 원소를 모두 포함하고 pycalphad 로 읽히는 파일만 후보로 삼는다.
2. 레시피 원소의 **모든 쌍이 `binaries`**(G·L 파라미터에 함께 나오는 쌍)에 있어야 한다. 없는 쌍을 이상용액으로 계산하면 그럴듯한 오답이 나오기 때문이다. 조건을 못 채우면 `not_applicable("공개 TDB 에 <계> 없음 …")` 을 낸다.
3. 순위: 평가 DB > 계가 딱 맞음(2원계 파일) > preferred(-LB·최신) > open 파티션 > 원소 수 적음.
4. TDB 라이선스가 restricted 면 결과는 내되 caveats 첫 줄에 '연구용 TDB' 를 적는다.

## 7. 미확인·후속

1. **SGTE 2원계 모음 부록의 라이선스**: CC BY 4.0 이 부록 데이터에도 적용되는가. 논문 본문의 데이터 제공 문구를 확인하거나 저자에게 문의한다. 확인되면 A6 대부분이 open 이 된다.
2. **COST 507 재배포·상업 이용 조건**: SGTE 또는 OpenCalphad(B. Sundman)에 문의한다.
3. **`licenses.yaml` 에 ODbL-1.0 추가**: 현재 MatCalc 는 CC-BY-SA-4.0 으로 매핑해 두었다.
4. **MatCalc → pycalphad 변환기**: 재배포 가능한 평가 DB 로 open 파티션 A6 을 만드는 길이다.
5. **LibreCalphad FCC_L12 문제**: 저자에게 제보할지 판단한다.
6. **NIST Solder**: SRD 139 재배포 허가 여부를 확인한다.

## 출처 (조회일 2026-09-27)

- [T1] Hallstedt 2025 부록 zip: https://ars.els-cdn.com/content/image/1-s2.0-S0364591625000367-mmc3.zip · 논문 https://doi.org/10.1016/j.calphad.2025.102833
- [T2] Crossref https://api.crossref.org/works/10.1016/j.calphad.2025.102833 · Elsevier API https://api.elsevier.com/content/article/pii/S0364591625000367 · 출판사 라이선스 안내 https://www.elsevier.com/about/policies-and-standards/open-access-licenses
- [T3] OpenCalphad https://www.opencalphad.com/databases.php · https://www.opencalphad.com/databases/COST507.zip · 보고서 https://www.opencalphad.com/databases/CGNA18499ENC_001.pdf
- [T4] MatCalc https://www.matcalc.at/index.php/databases/open-databases · ODbL https://opendatacommons.org/licenses/odbl/1-0/
- [T5] LibreCalphad https://github.com/mfrichtl/librecalphad (커밋 a7dace2)
- [T6] NIST Solder https://www.metallurgy.nist.gov/phase/solder/solder.html · https://doi.org/10.18434/T4759N · NIST 저작권 https://www.nist.gov/open/copyright-fair-use-and-licensing-statements-srd-data-software-and-technical-series-publications
- [T7] SGTE https://www.sgte.net/en/free-pure-elements-database
- [T8] NIMS CPDDB https://cpddb.nims.go.jp/en/
- [T9] TDBDB https://doi.org/10.1016/j.calphad.2018.04.003 · https://avdwgroup.engin.brown.edu/
- [T10] pycalphad https://github.com/pycalphad/pycalphad (LICENSE.txt, pycalphad/tests/databases)
- [T11] kawin https://github.com/materialsgenomefoundation/kawin/blob/main/examples/mobility_fitting/databases/CuNi.tdb
- [T12] phasediagrams.org https://phasediagrams.org/about
- 4.3절: https://zenodo.org/records/15266690 · https://zenodo.org/records/22161809 · https://zenodo.org/records/20848011 · https://zenodo.org/records/20181685 · https://zenodo.org/records/7767663 · https://github.com/PhasesResearchLab/ESPEI · https://github.com/sundmanbo/opencalphad
