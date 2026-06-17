# Background
Tissue/blood partition coefficients (logK) drive bioaccumulation and PBPK models.
The partition-toolkit screens cheap GFN1/GFN2-xTB single-point energies in 25 ALPB
+ 12–14 GBSA implicit solvents (78 single points per molecule) and uses each
solvent's `logK_solvent/water` (from the body-temperature transfer free energy) as a
*surrogate* for a biological phase. The surrogates are validated against 440 UFZ-LSER
compounds across five phases: phospholipid membrane, storage lipid, albumin, and
muscle protein (chicken / fish). The toolkit README reports a per-phase R² table; this
proposal re-states that result as a frozen, autonomously judged, re-verifiable claim.

# Goal
Establish, per phase, whether the **pre-registered** xTB surrogate is a *useful direct
predictor* of the measured partition coefficient — operationalised as coefficient of
determination **R² ≥ 0.70** over the UFZ test set.

# Method
For each phase, pair the FROZEN pre-selected surrogate `logK_*` column (the column the
toolkit's prior 78-single-point screening identified, README table) with the measured
UFZ `*_exp` column, join on CAS, drop missing rows, and compute

    R² = 1 − SS_res / SS_tot

with the surrogate used **directly** as the prediction (the toolkit's reported metric;
it can be negative when the uncalibrated surrogate is worse than the mean). The frozen
surrogate per phase:

| Phase | Measured column | Frozen surrogate |
|-------|-----------------|------------------|
| phospholipid membrane–water | `membrane_exp` | `logK_GFN2_gbsa_acetonitrile` |
| storage lipid–water | `storage_exp` | `logK_GFN1_alpb_hexane` |
| albumin–water | `albumin_exp` | `logK_GFN1_alpb_woctanol` |
| muscle protein–water (chicken) | `muscle_chicken_exp` | `logK_GFN1_alpb_octanol` |
| muscle protein–water (fish) | `muscle_fish_exp` | `logK_GFN1_alpb_octanol` |

Each R² is judged against the frozen 0.70 threshold autonomously by the sci-adk
DecisionEngine (numeric threshold rule, no LLM). Each hypothesis is **confirmatory**
and **empirical**; the Evidence is `data_source="measured"` (the dependent variable is
real measured UFZ logK), so the evidence-validity adequacy gate admits the verdict.

**Limitation (pre-registered honestly):** the surrogate per phase was selected on this
same dataset, so these are *resubstitution* R² — optimistic upper bounds, not
cross-validated. Pearson r² (the variance a linear calibration would explain) is
recorded alongside each verdict for transparency.

# Expected Output
A per-phase SUPPORTED / REFUTED verdict. From prior screening we expect the lipid
phases (membrane, storage) to be **SUPPORTED** and the protein phases (albumin, muscle
chicken/fish) to be **REFUTED** — i.e. the xTB single-conformer implicit-solvent
surrogate is a useful direct predictor of lipid partitioning but not of protein
partitioning. The whole run must re-derive from its record (`sci-adk verify` exit 0).

---

## 연구 배경 (Korean summary)
조직/혈액 분배계수(logK)는 생체축적·PBPK 모델의 핵심 입력이다. partition-toolkit은
GFN1/GFN2-xTB 단일점 에너지(분자당 78 SP, 25 ALPB + 12–14 GBSA 용매)로 각 용매의
`logK_solvent/water`를 생체 상(phase)의 *대리(surrogate)* 로 사용하고, UFZ-LSER 440
화합물 5개 상에 대해 검증한다. 본 제안은 README의 R² 표를 **고정된(frozen)·자동
판정되는·재검증 가능한** 주장으로 다시 진술한다.

## 연구 목표
상(phase)별로, 사전등록된 xTB surrogate가 측정된 분배계수의 *유용한 직접 예측자*인지를
**R² ≥ 0.70**(결정계수)로 판정한다.

## 연구 방법
상별로 사전선택·고정된 surrogate `logK_*` 컬럼과 측정값 `*_exp`를 CAS로 결합하고 결측을
제거한 뒤 `R² = 1 − SS_res/SS_tot`(surrogate를 직접 예측값으로 사용)를 계산해 고정된 0.70
임계값에 대해 sci-adk DecisionEngine이 자동 판정한다. 각 가설은 confirmatory·empirical이며
Evidence는 `data_source="measured"`(종속변수가 실제 측정 UFZ logK). **한계:** surrogate가
동일 데이터에서 선택되었으므로 resubstitution R²(낙관적 상한)이며 교차검증이 아니다.

## 기대 산출물
상별 SUPPORTED/REFUTED 판정. 예상: 지질상(membrane, storage)=SUPPORTED, 단백질상
(albumin, muscle)=REFUTED. 전체 run은 기록으로부터 재도출 가능해야 한다(`sci-adk verify`
exit 0).
