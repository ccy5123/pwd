# hand_off — BSA docking + xTB desolvation factorial (재해석 완료·상세 인계)

> **이전 `hand_off_session_bugfix.md` 대체판.** 차이: 버그 ①④②의 영향을 제공 CSV로 **실측 검증**해
> 결론과 올바른 fix를 **확정**함. 핵심 산출물 = **`reinterpret.py`**(docking 재실행 없이 재해석 완결).
> 표기 규약: **[Fact]** = 코드·데이터·논문 본문 직접 확인 · **[추론]** = Claude 해석. 논문은 read 여부 명시.
> 스타일: 한국어 + 영어 도메인 용어, terse, table 선호. 착수 전 모호하면 focused 1문항(객관식).
>
> **각 수치는 아래 [재현] 블록으로 서버에서 그대로 재산출 가능**(도구 무관, bash 직접 실행 전제).

작업폴더: `/home1/s9383/autodock_vina_research`  ·  env: `conda activate rapids-cuml`  ·  Rocky Linux/GPFS/SLURM, no sudo.

---

## 0. 불변 규칙
- **docking 재실행 금지**(`docking_poses.csv` 보존). 폴링 금지(`squeue`·`tail` 1회). pip 임의설치 금지(먼저 conda activate).
- 파일 치환은 **`apply_fixes.py`로만**(유일성 assert + ast.parse 게이트). 손 sed 금지.
- Gsolv는 `gsolv_fixed.py`로만. git은 이전 세션 관례(push는 사용자 승인 후).
- 규약: [Fact] vs [추론] 분리. 논문 인용 시 title·authors·journal·요약 + read 여부.

---

## 1. 연구 목표 [Fact]
논문 docking 프로토콜을 **중성 유기물 · BSA**에 이식하여, docking ΔG에 **xTB 수상 desolvation** 을 가산하는 보정이
albumin–water 분배계수(**log K_BSA/w**) 예측을 개선하는지, **포즈 3방식 × desolvation 5수준 = 15셀 factorial**로 측정.
**primary 지표 = absolute R²(원 설계) 그대로 유지·보고**(음수여도 그대로). calR²·Pearson은 **진단 병기**
(absR² 저하가 스케일 왜곡 때문인지 형태 불일치 때문인지 구분용). calR²로 결론을 **대체하지 않음**.
프레임: PFAS→중성유기물 / HSA→BSA(1AO6→4F5S) / solvation 미고려(none) ↔ 고려(xTB) 비교.

## 2. 논문 (docking 프로토콜 출처) [read: 본문·SI 실제로 읽음 ✓]
- **Title**: Deriving Membrane–Water and Protein–Water Partition Coefficients from In Vitro Experiments for Per- and Polyfluoroalkyl Substances (PFAS)
- **Authors**: Ruiwen Chen, Derek Muensterman, Jennifer Field, Carla Ng
- **Journal**: *Environ. Sci. Technol.* 2025, 59, 82–91 (+ SI)
- **핵심 프로토콜 [Fact]**: HSA 1AO6, AutoDock Vina, 9-grid box(각 26³Å, spacing 1.0Å), blind, exhaustiveness=32,
  20 modes, 최저에너지 채택. ΔG = −RT·ln K_A.
- **[Fact] 중요**: 논문은 9-center **생성법을 미기재**(SI Table S4에 좌표만). 진단상 그 center 중 5개가 1AO6에서도
  표면 밖(반경 10Å 내 Cα 0개) → "결합부위 표적"이 아니라 "전체 encompass grid". 그래서 박스는 **BSA(4F5S) chainA
  bounding box 3×3×1 = 9개**로 논문 방식을 재현(HSA→BSA 좌표이식은 폐기). **이 결정은 확정.**

## 3. 확정 설정 [Fact]
- 단백질: BSA **4F5S** chain A, HETATM/water 제거, PDB2PQR/PROPKA @pH7.4.
- 박스: **4F5S chainA bounding box 3×3×1 = 9개**, 각 26³Å, spacing=1.0Å. Box7(16Å)·Box8(8.2Å)이 단백질서 뜸
  → A(최저1)엔 무해, B1/B2엔 영향(§4c ②).
- ligand: RDKit(SMILES→AddHs→ETKDGv3 seed42→MMFF)→Meeko. **중성 유지(탈양성자화 X).**
- 포즈 3: **A**(전체 최저1) / **B1**(전포즈 Boltzmann 가중, w∝exp(−ΔG/RT)) / **B2**(전포즈 산술평균). B1·B2 RMSD<2.0Å dedup.
- desolv 5: **none** + {GFN1,GFN2}×{GBSA,ALPB}(water).
- 보정식: `dG_corr = dG_contact − gsolv`. 변환: `logKpw = −dG/(2.302585·RT) − 1.823`.
  상수 [Fact]: RT=0.616 kcal/mol(T=310.15K), LOG_M_BSA=log10(66.463)=1.823(몰→질량 kg/mol 기준 변환).

---

## 4. 결과 + 버그 검증 [Fact — CSV 실측 완료]

### 4a. 정합성 [Fact]
Job 718126 COMPLETED. `docking_poses.csv`=83 cpd, 2867 poses(13~66/cpd, mean 34.5).
`predictions_all_cells.csv`=1245행(3×5×83, **albumin_exp 병합됨**). `gsolv.csv`=332행(83×4). `results_factorial.csv`=15셀.

[재현]
```bash
python - <<'PY'
import pandas as pd
d=pd.read_csv("docking_poses.csv"); g=d.groupby("compound").size()
print("cpd=%d poses=%d /cpd: min=%d max=%d mean=%.1f  box uniq=%s"%(
      d.compound.nunique(), len(d), g.min(), g.max(), g.mean(), sorted(d.box.unique())))
print("predictions rows:", len(pd.read_csv("predictions_all_cells.csv")))
PY
```

### 4b. 15셀 재해석 [Fact — `python reinterpret.py` 실측]

**primary — absR² (원 설계 지표)**
| pose | none | GFN1/gbsa | GFN1/alpb | GFN2/gbsa | GFN2/alpb |
|---|---|---|---|---|---|
| **A** | **−0.675** | −24.8 | −24.3 | −25.5 | −30.3 |
| B1 | −1.553 | −29.4 | −28.8 | −30.1 | −35.1 |
| B2 | −5.244 | −41.9 | −41.1 | −42.4 | −47.6 |

전 셀 음수. 최고 = A×none −0.675 (offset 자유화해도 −0.226, slope≠1). desolv 열은 ④ 스케일로 폭락.

**진단 — calR² (=pearson², 선형보정 후; absR² 저하 원인 분해용)**
| pose | none | GFN1/gbsa | GFN1/alpb | GFN2/gbsa | GFN2/alpb |
|---|---|---|---|---|---|
| **A** | **0.219** | 0.125 | 0.147 | 0.083 | 0.128 |
| B1 | 0.208 | 0.116 | 0.138 | 0.077 | 0.121 |
| B2 | 0.154 | 0.069 | 0.088 | 0.041 | 0.083 |

- **최고 셀 = A×none**: calR²=0.219, pearson=0.468, absR²(raw)=−0.675, absR²(offset자유)=−0.226,
  **bias(exp−pred)=+0.578**(pred가 exp보다 일괄 낮음), n=83.

[재현] `python reinterpret.py` → 위 3개 표 + 검정 2줄 출력, `results_reinterpreted.csv` 저장.

### 4c. 버그 4개 — 검증 결과 [Fact] / 해석 [추론]

**① 변환 offset (−1.823)**
- [Fact] 위치: `bsa_docking_desolv_pipeline_v2.py` L1010(변환) · L67(상수 정의).
- [Fact] offset은 **순수 가산상수** → pred에 상수를 더해도 `pearson`·`calR²` 불변(실측 확인:
  pred±2.0에도 pearson=0.4676 고정). bias(exp−pred)=+0.578.
- [Fact — 본문 확인 ✓ (Claude가 albumin.pdf 본문 + albumin_sup.pdf SI 직접 열람)] 확정:
  Endo & Goss 2011 (Chem. Res. Toxicol. 24, 2293–2301) Table 1의 albumin_exp = **log K_BSA/w [L/kg], 37°C**
  (우리 83개 값 출처, 일치). Appendix eq12: K_a = K_BSA/w·MW_BSA(≈67 kg/mol). 본문 명시 **"log K_a = log K_BSA/w + 1.83"**.
  ∴ 코드 `logKpw = logKA − 1.823`는 정합(산술 검산: log10 66.463 = 1.8226 ≈ 1.83). **offset은 not-a-bug, 유지.**
  bias(exp−pred)=+0.578은 docking의 진짜 과소예측. (어느 경우든 calR²·순위 불변은 그대로.)

**② B1/B2 pool 오염**
- [Fact] L957 `dock_pivot = groupby('compound')['dG'].apply(list)` = 전 포즈 pool. L990~994에서 B1/B2 산출.
- [Fact 정정 — 중요] `docking_poses.csv` 의 `box` 열이 **전부 0**(L813 `pose_idx // len(unique_coords)`가 항상 0,
  dedup 시 box 태그 폐기). coords도 미저장. → **뜬 박스 필터·좌표 RMSD 오프라인 불가**, docking 재실행 필요.
- [Fact] B2(산술평균)는 화합물별 포즈수(13~66) 편차로 calR² 0.154(최악). B1(Boltzmann)은 A와 사실상 동급(0.208).
- [추론] pose 방식은 calR²를 거의 안 바꿈. absR²의 A>B1>B2 격차는 대부분 같은 스케일/bias 아티팩트. **B2 신뢰 낮음.**

**③ dedup RMSD**
- [Fact] L382 `rmsd = sqrt(mean(diff**2))` = naive 좌표차(대칭·정렬 무시). coords 미저장 → 재실행 필요.
- [추론] A·calR² 결론 무영향. B1/B2 정밀도용. **보류.**

**④ Gsolv 스케일**
- [Fact] L1009 `dG_corr = dG_contact − gsolv`, gsolv 음수(|mean|≈4.4 kcal/mol). → dG_corr = ΔG + |gsolv| →
  logKpw 음의 방향 이동. desolv 셀 pred가 none 대비 **−3.6** shift → absR² −24~−48.
- [추론] Vina 스코어는 이미 수용액 결합에 경험적 학습됨 → xTB Gsolv 1:1 가산은 **수상 desolvation 이중계상**
  (reference-state 불일치). 부호 방향은 맞으나 **스케일이 물리적으로 부정합.** 단 desolv 셀 calR²가 애초 none 미만
  이므로 **스케일을 고쳐도 desolvation은 무익.**

[재현] (④·② 스케일/pool 효과)
```bash
python - <<'PY'
import pandas as pd, numpy as np
from scipy.stats import pearsonr; from sklearn.metrics import r2_score
p=pd.read_csv("predictions_all_cells.csv"); g=pd.read_csv("gsolv.csv")
for dm in ['none','GFN1/gbsa','GFN2/alpb']:
    c=p[(p.pose_method=='A')&(p.desolv_method==dm)]; e,q=c.albumin_exp.values,c.logK_pred.values
    r=pearsonr(e,q)[0]; gm=0 if dm=='none' else g[g.method==dm.upper()].gsolv.abs().mean()
    print(f"A×{dm:<10} absR2={r2_score(e,q):+8.2f} calR2={r*r:.3f} pearson={r:+.3f} bias={np.mean(e-q):+.2f} |gsolv|={gm:.2f}")
PY
```

### 4d. 재해석 결론 [Fact 기반 — 확정]  (primary = absR², calR²는 진단)
- **primary(absR²) 결과**: 전 15셀 absR² < 0. 최고 = A×none −0.675(offset 자유화해도 −0.226, slope≠1). desolv 셀은
  −24~−48로 폭락(④ 스케일 아티팩트). **원 설계 지표로 "이 접근은 absR²>0 미달"이 1차 결론.**
- **진단(calR²/pearson, 선형변환 불변)**: absR² 저하가 "스케일 왜곡"인지 "형태 불일치"인지 분해. A×none calR²=0.219
  (pearson 0.47)이 최고이고 모든 desolvation·B1·B2가 그 미만 → desolv의 absR² 폭락은 스케일 탓이지만 **calR²로도 none
  미만**이라, 스케일을 고쳐도 desolvation은 무익. B1/B2도 A 못 넘음. (calR²는 결론을 대체하지 않고 원인을 설명.)
- **결론(primary 기준, 진단으로 뒷받침)**: 1) desolvation 보정 무익, 2) pose A>B1>B2, 3) docking+xTB는 이 계열에서
  absR²>0 미달(약함). 버그 수정은 absR² '절대크기'만 정상화할 뿐 이 결론을 안 바꿈.

**정직한 예측 상한 벤치마크 [Fact — 본문 확인 ✓ (Claude가 albumin.pdf/SI 열람)]**
Endo & Goss 2011: log Kow 회귀(eq3) **R²=0.76, SD=0.43**(n=76); PP-LFER(Table 2, eq4/5) **R²=0.78–0.79, SD=0.41–0.42**(n=82).
**[중요 — 지표 정합]** 논문 R²는 **fitted 회귀/MLR 값 = 우리 calR²(선형보정 허용)와 동급**. 따라서 벤치마크는
**우리 calR²(≈0.22) vs 논문 R²(0.76–0.79)** 로 비교 → docking이 fitted-ceiling에 크게 못 미침. (우리 absR²는 무보정
물리예측이라 더 엄격한 별도 지표이므로 논문 R²와 직접 비교하지 말 것.) 논문 자체도 PP-LFER SD 0.41–0.42로 "solvent
dissolution 모델로 albumin 결합을 정확히 못 잡는다 + 결합부위 size limitation"을 결론 → 우리 "약함" 판정과 정합.

**지표 정책 [사용자 확정]**: primary = **absolute R²**(원 설계) 그대로 정직 보고(음수여도 그대로), calR²·Pearson 병기(진단).
결과를 좋게 만들려 하지 않음 — ②③ 수정 후 A×none 불변·B1/B2만 변동 확인이 목표.

---

## 5. 작업 상태 + 남은 작업 [확정] — 상세 순서는 `NEW_SESSION_PROMPT.md` STEP 0~7

**완료(이전 세션)**
- FIX① 확정 [Fact — 본문 확인 ✓ (Claude가 albumin.pdf/SI 열람)]: albumin_exp = log K_BSA/w [L/kg], 37°C (Endo & Goss 2011 Table 1, 83개 값 일치). 본문 "log K_a = log K_BSA/w + 1.83" + eq12(K_a=K_BSA/w·MW_BSA, ~67 kg/mol) → offset −1.823 논문 일치. **not-a-bug, 유지.**
- FIX①④: `apply_fixes.py` 적용(GSOLV_SCALE=1.0 = 기존과 동일, caveat 주석). **수치 불변**, ast.parse OK.

**이번 세션 = FIX②③ 실제 수정 + 재실행 + 결과 확인** (사용자 지시: "버그 수정후 결과 확인, 기존 논리대로")
1. `apply_fixes_23.py` → ② box 식별자 보존(dock_compound→dedup→save) + ③ RMSD 공식 교정(sum(axis=1)). (STEP2)
2. 기존 `docking_poses.csv` + `_checkpoint_docking.json` 백업·제거 (append/체크포인트 함정 회피). (STEP3)
3. `--docking-only` **재실행**(~10분) → box 0..8 복원 + dedup 교정. (STEP4)
4. `--analyze` 재계산(gsolv.csv 재사용) → `reinterpret.py`로 새 15셀. **A×none 불변(≈0.219) 확인**이 정상 기준. (STEP5)
5. `box_filter_check.py` → box별 dG로 뜬 박스 **데이터 판별**(추측 금지) → EXCLUDE_BOXES 필터 B1/B2 비교. (STEP6)
6. git commit → hand_off 갱신. push는 사용자 승인 후. (STEP7)

**④·①은 기존 논리대로 불변**(GSOLV_SCALE=1.0, offset −1.823). 재실행 시 자동 재계산만.

---

## 6. 파일 인벤토리 (번들 = 서버 작업폴더에 있어야 함)
| 파일 | 역할 |
|---|---|
| `bsa_docking_desolv_pipeline_v2.py` | 파이프라인(1260줄). 패치 대상(FIX①④). |
| **`reinterpret.py`** | **재해석(신규, 검증완료). docking 없이 CSV로 재해석·figure 산출.** |
| **`apply_fixes.py`** | FIX①④ 적용기(GSOLV_SCALE knob + caveat, no-op). 이미 적용됨. |
| **`apply_fixes_23.py`** | **FIX②③ 실제 적용기(box 보존 + RMSD 공식 교정). full-block 치환 + assert + ast.parse.** |
| **`box_filter_check.py`** | **② 목적 확인: 재실행 후 box별 dG로 뜬 박스 판별 + B1/B2 필터 비교(EXCLUDE_BOXES).** |
| `pipeline_fixes.md` | 패치 rationale + anchored old→new(②③ 포함). |
| `gsolv_fixed.py` · `run_gsolv.py` | Gsolv 계산(정상, 그대로). |
| `docking_poses.csv` | docking 산출(83 cpd, 2867 poses, **box 전부 0**). **보존·재사용.** |
| `gsolv.csv` · `predictions_all_cells.csv` · `results_factorial.csv` | 현 결과(exp 병합). |
| `compounds_input_table1.csv` | 입력(83 유효, albumin_exp). |
| 산출: `results_reinterpreted.csv` · `scatter_A_none.png` · `calR2_heatmap.png` · `absR2_heatmap.png` | reinterpret.py 출력. |
| `hand_off.md`(본 파일) · `NEW_SESSION_PROMPT.md` | 인계·실행 프롬프트. |

## 7. 버그 이력 (해결됨 — 재발방지) [Fact]
CSV append 유실 / Meeko ATOM coords 빈배열 / v.poses() 'M' 반환 / dG=energies[:,0] / 박스 99개·coverage 오판 /
Kabsch R 방향(→9-center 이식 폐기, 3×3×1 확정) / box save 유실(L813) → **이번 세션 `apply_fixes_23.py`로 수정+재실행** / RMSD 공식 √3 과소평가(L382, 과잉병합) → **이번 세션 교정**.

---

### 착수 요약 (한 줄)
**이번 세션 = ②③ 실제 수정: `apply_fixes_23.py`로 box 식별자 보존 + RMSD 공식 교정 → 기존 docking 산출물·체크포인트
제거 → `--docking-only` 재실행 → `--analyze` 재계산 → `reinterpret.py`(A×none 불변 확인)·`box_filter_check.py`(뜬 박스
데이터 판별)로 결과 확인 → 커밋. ④·①은 기존 논리대로 불변. [추론] A×none은 Vina seed 고정으로 불변 예상,
B1/B2만 dedup 교정으로 변동 — 결과는 재실행으로 확정.**
