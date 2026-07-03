# 실행 프롬프트 — 버그 ②③ 실제 수정 + docking 재실행 + 결과 확인 (도구 무관·서버 직접 실행)

너는 서버 `/home1/s9383/autodock_vina_research` 에서 **bash를 직접 실행**하는 agentic 어시스턴트다.
STEP을 **순서대로** 실행. 각 STEP: [명령] 그대로 → [기대] 대조 → 불일치면 [실패시]대로 **멈추고 보고**.
명령 임의 변경 금지. 경로는 작업폴더 기준. 규약: [Fact](코드·데이터·본문) vs [추론] 분리. 논문 인용 시 title·authors·journal·요약 + read여부.

## 이번 세션 목적 (사용자 확정)
"버그 실제 수정 → 재실행 → 결과 확인. **기존 논리대로**, 새 결정 추가 금지. 결과는 그대로 받는다(안 좋아도 됨)."
- **②③만 실제 수정 + docking 재실행 + factorial 재계산 + 결과 확인.**
- **④·①은 그대로 유지.** ④ GSOLV_SCALE=1.0(기존과 동일), ① offset −1.823.

## 확정 사실 [Fact — 코드·논문 본문 확인]
- **① offset −1.823 = not-a-bug [Fact 본문 확인 ✓]**: Endo & Goss, *Chem. Res. Toxicol.* 2011, 24, 2293–2301 (+SI).
  Table 1 = log K_BSA/w **[L/kg], 37°C**(우리 83개 albumin_exp 출처, 값 일치). Appendix eq12 K_a=K_BSA/w·MW_BSA(~67 kg/mol),
  본문 **"log K_a = log K_BSA/w + 1.83"** → 코드 `logKpw=logKA−1.823` 정합(log10 66.463=1.8226≈1.83). **손대지 않음.**
- **④ Gsolv 1:1 가산**: [추론] Vina 스코어가 수용액 결합 학습이라 xTB Gsolv 1:1 가산은 이중계상 의심. 단 desolv 셀이
  애초 none보다 나쁨 → 스케일 고쳐도 무익. 처리방식 변경은 새 결정이므로 **기존 논리대로 손대지 않음(1.0 유지).**
- **② box 유실 [Fact]**: `dock_compound`가 (dG,coords)만 반환 → 저장 시 `box: pose_idx//len` = **항상 0**. 뜬 박스 식별 불가.
- **③ RMSD 공식 오류 [Fact]**: `np.sqrt(np.mean(diff**2))`는 3N으로 나눠 RMSD **√3배 과소평가** → 2.0Å 컷오프서 과잉 병합.

## 지표 정책 [사용자 확정 — 반드시 준수]
- **primary = absolute R²(원 설계)**. 재산출·정직 보고. **음수여도 그대로.** calR²로 결론을 **대체하지 마라.**
- **calR²·Pearson = 진단 병기**(absR² 저하가 스케일 왜곡인지 형태 불일치인지 분해용).
- **벤치마크 비교 시 [Fact 본문]**: 논문 R²=0.76(logKow 회귀)/0.78–0.79(PP-LFER)는 **fitted 회귀/MLR 값 = 우리 calR²와 동급**.
  → 벤치마크는 **우리 calR²(≈0.22) vs 논문 R²(0.76–0.79)** 로 비교. **우리 absR²를 논문 R²와 직접 비교하지 마라**(무보정 물리예측 vs fitted, 지표 종류 다름).

## 불변 규칙
- 치환은 **`apply_fixes_23.py`로만**(assert+ast.parse+멱등). 손 sed 금지. 재실행 전 **백업 필수**. 폴링 금지(squeue·tail 1회). pip 임의설치 금지(먼저 conda activate). Gsolv는 `gsolv_fixed.py`로만. push는 사용자 승인 후.

---

## STEP 0 — 환경·파일
```bash
cd /home1/s9383/autodock_vina_research
source ~/.bashrc 2>/dev/null; conda activate rapids-cuml && echo "ENV_OK: $(python -c 'import sys;print(sys.version.split()[0])')"
ls -1 apply_fixes_23.py box_filter_check.py reinterpret.py bsa_docking_desolv_pipeline_v2.py \
      docking_poses.csv gsolv.csv predictions_all_cells.csv results_factorial.csv compounds_input_table1.csv
python -c "import pandas,numpy,scipy,sklearn; print('DEPS_OK')"
```
[기대] `ENV_OK: 3.x`, 9개 파일 존재, `DEPS_OK`.
[실패시] env명 다르면 멈추고 질의. 스크립트 없으면 "번들 업로드 필요" 보고(임의 생성 금지).

## STEP 1 — 현재(수정 前) 스냅샷 보존
```bash
python reinterpret.py 2>&1 | sed -n '1,30p'
cp results_factorial.csv results_factorial.PRE.csv
cp predictions_all_cells.csv predictions_all_cells.PRE.csv
cp docking_poses.csv docking_poses.PRE.csv
```
[기대] A×none absR²=−0.675 / calR²=0.219 등 기존 값 출력, `.PRE.csv` 3개 생성.
[용도] 재실행 후 PRE/POST 비교, 특히 A×none 불변 확인.

## STEP 2 — 코드 수정 (②③) · 안전 스크립트로만
```bash
cp bsa_docking_desolv_pipeline_v2.py bsa_docking_desolv_pipeline_v2.py.bak23
python apply_fixes_23.py
python -c "import ast; ast.parse(open('bsa_docking_desolv_pipeline_v2.py').read()); print('ast.parse OK')"
grep -n "results.append((float(dG), coords, box_idx))" bsa_docking_desolv_pipeline_v2.py
grep -n "np.sum(diff\*\*2, axis=1)" bsa_docking_desolv_pipeline_v2.py
grep -n "'box': box_id" bsa_docking_desolv_pipeline_v2.py
grep -c "pose_idx // len" bsa_docking_desolv_pipeline_v2.py   # 0 이어야
```
[기대] `[OK] FIX2/3 적용 완료 + ast.parse 통과.`, `ast.parse OK`, box_idx/np.sum/'box':box_id 각 1줄, 마지막 `0`.
[실패시] `[STOP] anchor count=...` → `cp *.bak23 원본` 복원 후 멈추고 보고. `[SKIP]` → 이미 적용됨(정상).

## STEP 3 — 재실행 준비: 기존 docking 산출물 제거 (⚠ 필수)
[Fact] `run_full_docking`은 **체크포인트+append**(L755, L821). 안 지우면 "완료"로 스킵하거나 중복 append됨.
```bash
rm -f docking_poses.csv _checkpoint_docking.json     # PRE는 STEP1에서 이미 백업됨
ls -la docking_poses.csv 2>&1 | head -1              # No such file 이어야
ls -la _checkpoint_docking.json 2>&1 | head -1
```
[기대] 둘 다 "No such file". `_checkpoint_receptor.json`·`gsolv.csv`·gsolv 체크포인트는 **남겨둠**(재사용).

## STEP 4 — docking 재실행 (~10분, 83 cpd × 9 box)
```bash
nohup python bsa_docking_desolv_pipeline_v2.py --docking-only > rerun_docking.log 2>&1 &
echo "PID=$!"
```
[대기·1회만 확인] (폴링 금지)
```bash
sleep 600; tail -n 15 rerun_docking.log
python -c "import pandas as pd; d=pd.read_csv('docking_poses.csv'); print('cpd=%d rows=%d box_uniq=%s'%(d.compound.nunique(), len(d), sorted(d.box.unique())))"
```
[기대] `cpd=83`, `box_uniq=[0,1,2,3,4,5,6,7,8]`(전부 0 아님 = ② 반영). rows는 ③(dedup 교정)로 기존 2867과 다를 수 있음.
[실패시] cpd<83 또는 box_uniq=[0] → 로그 확인·보고. 실행 중이면 sleep 1회 더만.

## STEP 5 — factorial 재계산 (gsolv 재사용) + primary(absR²) 재산출
```bash
python bsa_docking_desolv_pipeline_v2.py --analyze > rerun_analyze.log 2>&1
tail -n 5 rerun_analyze.log
python reinterpret.py
```
[기대]
- **primary = absR²** 표(전 셀 음수, A×none≈−0.675) + 진단 calR²(A×none≈0.219).
- **A×none absR²·calR² 불변**(Vina seed 42 고정 → 최저 dG 불변 → A 불변)이 재실행 정상 기준.
- B1/B2만 ③(dedup 교정)로 변할 수 있음.
[실패시] A×none이 크게 달라지면 `.PRE.csv`와 대조·보고 후 멈춤.

## STEP 6 — 결과 확인: PRE/POST 비교 + 뜬 박스 필터(②의 목적)
```bash
echo "=== PRE vs POST (primary=absR², 진단 calR²) ==="
python - <<'PY'
import pandas as pd
pre=pd.read_csv("results_factorial.PRE.csv"); post=pd.read_csv("results_factorial.csv")
k=['pose_method','desolv_method']
m=pre.merge(post,on=k,suffixes=('_pre','_post'))
print(m[k+['absR2_pre','absR2_post','calR2_pre','calR2_post']].to_string(index=False,
      float_format=lambda x:f"{x:+.3f}"))
PY
echo "=== box별 dG + 뜬박스 필터 ==="
python box_filter_check.py
```
[판정]
- `box_filter_check.py` (1)에서 **평균 dG가 0 근처인 box**가 뜬 박스. 그 번호를 상단 `EXCLUDE_BOXES=[...]`에 넣고 재실행:
```bash
# 예: (1)에서 box 7,8이 평균 dG≈0 이면
sed -i 's/^EXCLUDE_BOXES: list = \[\]/EXCLUDE_BOXES: list = [7, 8]/' box_filter_check.py
python box_filter_check.py | sed -n '/(2)/,$p'
```
[해석 규칙 — 결과 안 바꿈]
- primary=absR²로 보고. calR²는 진단 병기. **결과가 여전히 약해도(absR²<0) 그대로 보고.**
- 벤치마크: 우리 **calR²(≈0.22) vs 논문 R²(0.76–0.79, fitted)** → docking이 fitted-ceiling에 크게 못 미침(정합).
- ②③ 수정 목표 = A×none 불변·B1/B2만 변동 확인. 필터 후에도 A≥B1≥B2면 "A 최고 유지"[Fact], 역전하면 그 수치로 보고.
- **뜬 박스 번호는 데이터로 판별(추측 금지).**

## STEP 7 — 커밋 + hand_off 갱신
```bash
git add apply_fixes_23.py box_filter_check.py reinterpret.py bsa_docking_desolv_pipeline_v2.py \
        docking_poses.csv gsolv.csv predictions_all_cells.csv results_factorial.csv \
        results_reinterpreted.csv scatter_A_none.png calR2_heatmap.png absR2_heatmap.png hand_off.md
git status --short
git commit -m "FIX2/3: real box id + corrected heavy-atom RMSD; re-dock + recompute; primary=absR2 (calR2 diagnostic); FIX1 confirmed from Endo&Goss2011 body"
```
[STEP7 hand_off.md 갱신] §4b(새 absR²/calR²), §4c②③(수정완료), §5, §7에 1~2줄 + commit hash. A×none 불변 여부·B1/B2 변화·뜬 박스 필터 결과·EXCLUDE_BOXES 값 기록. **push는 사용자 승인 후.** [보고 후 STOP]

---

## 부록 A — 재실행 정상성 self-check
- A×none absR²≈−0.675, calR²≈0.219 **불변**(seed 고정). 벗어나면 이상.
- `docking_poses.csv` box_uniq=[0..8](전부 0 아님). pose 수는 ③로 변동 가능.
- gsolv.csv·`_checkpoint_receptor.json` 재사용(건드리지 않음).

## 부록 B — 무엇을 바꿨나 (기존 논리 유지)
- ② dock_compound가 box_idx 반환 → dedup 통과 → 저장 시 실제 box. (유실 복구, 새 결정 없음)
- ③ RMSD를 per-atom 공식(sum(axis=1))으로 교정. heavy-atom RMSD<2.0Å 의도 그대로. (새 파라미터 없음)
- ④·① 미변경. 뜬 박스 번호는 데이터로 판별(하드코딩 없음).

## 부록 C — 한 줄
**apply_fixes_23.py로 ②③ 수정 → 기존 docking_poses.csv·docking 체크포인트 제거 → --docking-only 재실행(box 복원·dedup 교정)
→ --analyze 재계산 → reinterpret로 primary=absR²(음수 그대로) 보고·A×none 불변 확인, calR²는 진단·논문 R²(fitted)와 calR² 비교
→ box_filter_check로 뜬 박스 데이터 판별 → 커밋. ④·①은 기존 논리대로 불변. FIX① offset은 논문 본문 확인으로 not-a-bug 확정.**
