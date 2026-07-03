# pipeline_fixes.md — `bsa_docking_desolv_pipeline_v2.py` 수정 지침

> **실제 적용은 이렇게**: FIX①④ = **`python apply_fixes.py`** 가 자동 적용(유일성 assert + ast.parse 게이트 +
> 멱등 가드). 이 문서는 **rationale + anchored old→new 원문**이며, FIX②③(재실행 필요)의 수동 참조용이다.

> 전제: 재해석 자체는 **`reinterpret.py` 하나로 끝**(docking 재실행 불필요, 기존 CSV만 사용).
> 아래 패치는 **미래 full-run을 정직·재현 가능하게** 만드는 영구 수정. 각 패치는 anchor 문자열로
> 위치를 잡고, `old -> new`를 **정확히 1회** 치환. 치환 전 `grep -n`으로 유일성 확인, 후 `python -c "import ast; ast.parse(open('bsa_docking_desolv_pipeline_v2.py').read())"`로 파싱 검증.
> 라인번호는 참고값(현재 파일 1260줄 기준). anchor 문자열이 진짜 기준.

---

## 검증으로 확정된 사실 [Fact — CSV 실측]

| 항목 | 값 |
|---|---|
| A×none bias (exp−pred) | **+0.578** (pred 일괄 과소예측) |
| A×none absR² (raw, offset −1.823) | −0.675 |
| A×none absR² (offset 자유화, slope=1) | **−0.226**  ← ① 기여 제거해도 여전히 음수 |
| A×none calR² (=pearson², 완전 선형보정) | **+0.219** (pearson 0.468) ← 진짜 성능 상한 |
| desolv 셀 calR² | 0.041~0.147 (전부 none 미만) |
| docking_poses.csv `box` 열 | **전부 0** (L813 버그로 유실) → ②③ 오프라인 불가 |

**결론**: ①④②는 pred를 **선형으로만** 이동 → pearson·calR² **불변**. absR²만 왜곡. 정직한 지표 calR²는
**A×none이 최고(0.219), 모든 desolv·B1·B2가 그 미만**. → 버그 수정은 absR² '절대크기'를 정상화할 뿐,
**"desolvation 무익 · A>B1>B2 · absR²<0" 결론은 버그와 무관하게 유지.**

---

## FIX ① — offset 물리 검증 (코드 변경 아님, 확인 항목)

- L67 `LOG_M_BSA = math.log10(66.463)  # 1.823` 은 **몰→질량 기준 변환**(logK_pw[L/kg] = logK_A − log10 M_BSA).
- **검증 필요**: `compounds_input_table1.csv`의 `albumin_exp`가 **질량기준 log K_BSA/w (L/kg)** 인가?
  출처 Endo & Goss 2011 / UFZ 정의 확인. 질량기준이면 −1.823 **정답**(bias +0.58은 진짜 docking 과소예측).
  몰기준이면 offset 제거 대상.
- **패치 없음** — 결론이 어느 쪽이든 calR²·순위 불변이므로 재해석엔 영향 0. 문서화만:

```
# anchor (L1010 근처)
old:
                dG_corr = dG_contact - gsolv
                logK_pred = -dG_corr / (2.302585 * RT) - LOG_M_BSA
new:
                # [FIX①] LOG_M_BSA(-1.823)는 몰->질량기준 변환. albumin_exp가 질량기준
                #        log K_BSA/w(L/kg)일 때만 정합(Endo&Goss2011/UFZ 정의 확인).
                #        offset은 순수 가산상수 -> pearson/calR² 불변, absR²/bias만 이동.
                # [FIX④] Vina ΔG는 이미 수용액 결합에 대해 경험적 학습됨 -> xTB Gsolv 1:1 가산은
                #        수상 desolvation을 이중계상(reference state 불일치). GSOLV_SCALE로 감도만 봄.
                dG_corr = dG_contact - GSOLV_SCALE * gsolv
                logK_pred = -dG_corr / (2.302585 * RT) - LOG_M_BSA
```

## FIX ④ — Gsolv 스케일 knob (기본값 1.0 = 기존과 동일, backward-compatible)

L67 블록 끝(`LOG_M_BSA = ...` 다음 줄)에 상수 추가:

```
# anchor (L67)
old:
LOG_M_BSA = math.log10(M_BSA)  # = 1.823
new:
LOG_M_BSA = math.log10(M_BSA)  # = 1.823
# [FIX④] xTB Gsolv 가산 스케일. 1.0=기존(이중계상 의심). desolvation 감도분석용 knob.
#        검증결과 |Gsolv|mean~4.4 kcal/mol을 Vina ΔG에 1:1 가산 -> pred -3.6 shift -> absR² 폭락.
#        calR²(스케일불변)로 보면 desolvation은 어차피 none보다 나쁨.
GSOLV_SCALE = 1.0
```
> 위 FIX① 치환에서 `dG_contact - gsolv` → `dG_contact - GSOLV_SCALE * gsolv` 이미 반영됨.
> 스모크(L621,625)도 동일하게 바꾸려면 그쪽 `best_dG - gsolv`/`best_dG + gsolv`에도 적용(선택).

## FIX ② — B1/B2 pool 오염 (오프라인 부분수정 + 근원은 재실행 필요)

- **오프라인 가능**: B2(전포즈 산술평균)는 화합물별 포즈수(13~66) 편차로 계상 불공정 → calR² 0.154로 최악.
  B1(Boltzmann)은 저에너지 포즈 지배라 A와 사실상 동급(calR² 0.208). **권고: B2 비중 낮추고 A/B1만 신뢰.**
- **근원(재실행 필요)**: 뜬 박스(Box7·8) 포즈 필터는 **box 식별자 필요**하나 `docking_poses.csv`엔 유실.
  L813이 원인:

```
# anchor (L813) — box 유실 근원. pose_idx < len 이라 항상 0.
old:
                    'box': pose_idx // len(unique_coords),  # Re-index after dedup
                    'pose': pose_idx % len(unique_coords),
new:
                    # [FIX②] 원 box 식별자 보존해야 뜬 박스 필터 가능. dock_compound가
                    #        (dG, coords, box_id) 3-tuple 반환하도록 바꾼 뒤 여기서 box_id 저장.
                    #        (현재 dock_results는 (dG,coords) 2-tuple -> box 복원 불가 = 재실행 필요)
                    'box': pose_box_ids[pose_idx] if 'pose_box_ids' in dir() else 0,
                    'pose': pose_idx,
```
> **주의**: 이 패치는 `dock_compound`(L392~)가 box_id를 함께 반환하도록 고쳐야 실효.
> ②를 진짜로 하려면 docking 재실행 불가피 — **우선순위 낮음**(calR² 결론 불변).

## FIX ③ — dedup RMSD (재실행 필요, 보류 권장)

```
# anchor (L382)
old:
                diff = coords[i] - coords[j]
                rmsd = np.sqrt(np.mean(diff**2))
new:
                # [FIX③] naive 좌표차 RMSD는 대칭·정렬 무시. RDKit GetBestRMS로 교체해야 정확하나
                #        coords가 docking_poses.csv에 없음 -> docking 재실행 시에만 적용 가능. (보류)
                diff = coords[i] - coords[j]
                rmsd = np.sqrt(np.mean(diff**2))
```
> A(최저1)·calR² 결론엔 무영향. B1/B2 정밀도용. **보류.**

---

## 적용 순서 (Claude Code)

1. `conda activate rapids-cuml` (없으면 안내·멈춤).
2. **`python reinterpret.py`** 먼저 실행 → 재해석 테이블·figure 산출(수정 前에도 동작). 결과 확인.
3. 위 FIX①④ 치환(2곳: L67 상수, L1010 블록) → `ast.parse` 검증 → git commit.
4. ②③는 **재실행 결정 후에만**. 안 할 거면 문서 주석만(위 patch의 comment 부분) 반영.
5. `hand_off.md` 갱신 상태 확인.
