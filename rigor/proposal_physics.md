# Background
Tissue/blood partition coefficients (log K_phase/water) are core inputs to PBPK and
bioaccumulation models. The established empirical route — PP-LFER (Endo/Brown/Goss
2013, this dataset's source) — is **not an a priori prediction**: it requires (i)
*measured* Abraham solute descriptors (E, S, A, B, V, L) and (ii) phase system
coefficients **fitted to experimental partition data**. It therefore cannot predict for
a solute or a phase outside its fitted calibration.

A first-principles alternative computes log K **purely from quantum-chemical
thermodynamics, with zero parameters fitted to the partition data**: represent each
biological phase as an *explicit molecular model* (storage lipid, DMPC membrane,
albumin / muscle protein as an amino-acid pseudo-solvent at the tissue's composition),
obtain each solute's σ-profile from a quantum calculation, and use COSMO-RS to compute
the solute's infinite-dilution activity coefficient ln γ∞ in water and in each phase.
The partition coefficient then follows from the free-energy difference,

    log K_phase/water = (ln γ∞_water − ln γ∞_phase) / ln 10 + c_unit ,

where c_unit is a fixed molar-volume unit conversion, **not a fitted parameter**. This
is the method already drafted in the repo (`*/cosmo_kpw_openrs_*_v4.py`,
openCOSMO-RS / openCOSMORS24a, body temperature 310.15 K).

This proposal is explicitly **NOT** a descriptor-based / machine-learning surrogate
model. The target value is *derived from physical law* (statistical-thermodynamic
activity coefficients over quantum σ-profiles), not regressed from features. (The prior
`partition-model-cv` run did the latter — using the physically computed logK values as
ML features — which is the opposite of this goal and is superseded here.)

# Goal
Determine how well **parameter-free, first-principles COSMO-RS thermodynamics** — phase
= explicit molecular pseudo-solvent, σ-profile = quantum calculation, **no fit to
partition data** — reproduces the measured tissue–water partition coefficients across
the five biological phases, and identify *where the physical model systematically
fails* (which missing physics it reveals).

# Method
Phase models (explicit molecules, not descriptors):
- storage lipid → a triacylglycerol model (e.g. triolein) as pure pseudo-solvent;
- phospholipid membrane → DMPC (1,2-dimyristoyl-PC, zwitterionic) as pure pseudo-solvent;
- albumin → BSA amino-acid composition as a 20-residue pseudo-solvent (capped Ace-X-Nme);
- muscle protein → actin+myosin amino-acid composition pseudo-solvent (capped residues).

Per solute and per phase residue/molecule:
1. 3D structure (RDKit ETKDG) → geometry optimization → COSMO/CPCM σ-profile from a
   quantum calculation (the repo's cached route: xTB-optimized geometry → CPCM
   BP86/def2-TZVPD → `.orcacosmo`); parameterization openCOSMORS24a.
2. Compute ln γ∞ of the solute (a) in pure water and (b) in the phase pseudo-solvent
   (amino acids mixed at the tissue mole fractions) at 310.15 K.
3. log K_pw = (ln γ∞_water − ln γ∞_phase)/ln 10 + c_unit. **No parameter is fitted to
   the experimental K.**

Validation against the measured UFZ-LSER sets (membrane 207, storage 247, albumin 83,
muscle chicken/fish ≈46, plus the Endo-46 muscle check). Report **RMSE, mean signed
error (bias), R², and the predicted-vs-measured slope** per phase. Decision thresholds
are **field-standard and pre-committed, NOT tuned to the observed numbers**:
- *useful for screening*: RMSE ≤ 1.0 log units;
- *quantitative / near chemical accuracy*: RMSE ≤ 0.5 log units;
- *unbiased*: |mean signed error| ≤ 0.3 log units (no systematic offset).

The decisive distinction from PP-LFER and from any ML model: this prediction uses **no
measured descriptors and no parameters fit to partition data** — it is fully a priori.

# Expected Output
Per phase, a parameter-free physical-prediction accuracy (RMSE, bias) with a
SUPPORTED/REFUTED verdict against the pre-set bars. Honest prior expectation from the
literature (recorded as prior work, not as a tuned target): a *homogeneous* pseudo-
solvent COSMO-RS is known to be limited for membranes — Endo 2011 reports RMSE ≈ 1.0
for the homogeneous-PC "Approach 2" because a bilayer is anisotropic and needs the
depth-resolved COSMOmic treatment; the amino-acid pseudo-solvent for proteins is
likewise an approximation. The real deliverable is therefore twofold: (1) the a priori
accuracy physics alone reaches with zero fitting, and (2) *the systematic error pattern*
(e.g. PAHs under-predicted in membrane; H-bond-donor compounds in protein) that
pinpoints the missing physics — the mechanistic insight an ML surrogate cannot give.

NOTE (anti-HARKing): prior COSMO-RS results already exist in the repo, so this is an
EXPLORATORY re-evaluation; the bars above are field-standard, set independently of those
numbers, and a true confirmatory test would pre-register on an untouched solute set.

---

## 연구 배경 (Korean summary)
조직/혈액 분배계수는 PBPK·생체축적 모델의 핵심 입력이다. 기존 PP-LFER는 *측정* Abraham
descriptor와 *분배데이터에 피팅된* 상 계수를 요구하므로 진정한 a priori 예측이 아니다.
1st-principles 대안은 각 생체 상을 **명시적 분자 모델**(저장지질, DMPC 막, 알부민·근육은
아미노산 pseudo-solvent)로 표현하고, 용질의 σ-profile을 **양자계산**으로 얻어 COSMO-RS로
물·상에서의 무한희석 활동계수 lnγ∞를 계산, `log K=(lnγ_w−lnγ_p)/ln10`으로 **분배데이터에
피팅하는 파라미터 0개**로 logK를 *물리적으로 도출*한다 (repo의 `cosmo_kpw_openrs_*_v4.py`).
**이것은 descriptor 기반 ML surrogate가 아니다** — 직전 `partition-model-cv`(물리 결과를
ML 특징으로 회귀)는 목표의 정반대였고 여기서 폐기된다.

## 연구 목표
무피팅 1st-principles COSMO-RS 열역학(상=명시적 분자, σ-profile=양자계산, 분배데이터 무피팅)이
5개 생체상의 측정 분배계수를 얼마나 재현하는지, 그리고 **물리모델이 어디서 체계적으로 실패하는지**
(빠진 물리)를 규명한다.

## 연구 방법
상 모델(명시적 분자): 저장지질=트리올레인, 막=DMPC, 알부민=BSA 아미노산 조성, 근육=actin+myosin
조성(capped Ace-X-Nme). 용질·상: 3D→최적화→COSMO σ-profile(xTB 기하+CPCM BP86/def2-TZVPD)→
openCOSMORS24a. 310.15 K에서 lnγ∞(물)·lnγ∞(상) 계산 → `log K=(lnγ_w−lnγ_p)/ln10+c_unit`
(c_unit=몰부피비, 피팅 아님). 검증: 측정 UFZ셋에 RMSE·bias·R²·기울기. **임계값은 분야표준 사전고정**
(유용 RMSE≤1.0; 정량 RMSE≤0.5; 무편향 |bias|≤0.3). PP-LFER·ML과의 본질 차이 = **측정 descriptor·
분배데이터 피팅 0**.

## 기대 산출물
상별 무피팅 물리예측 정확도 + 판정. 문헌 기반 정직한 예상: 균질 pseudo-solvent COSMO-RS는 막에서
부정확(~RMSE 1.0; 깊이분해 COSMOmic 필요, Endo 2011), 단백질 AA-pseudo-solvent도 근사. 진짜
산출물은 (1) 무피팅 물리가 도달하는 정확도와 (2) *체계적 오차 패턴*(막의 PAH 과소예측, 단백질의
H결합 화합물)으로 **빠진 물리를 짚는 것** — ML surrogate가 줄 수 없는 통찰.
