# Establish, per phase, whether the pre-registered xTB surrogate is a useful direct predictor of the measured partition coefficient (coefficient of determination R^2 >= 0.70).

> Draft compiled by sci-adk from Spec `partition-xtb-ufz` (v1). Belief state is revisable as Evidence accrues.

## Goal
Establish, per phase, whether the pre-registered xTB surrogate is a useful direct predictor of the measured partition coefficient (coefficient of determination R^2 >= 0.70).

## Background
Tissue/blood partition coefficients drive bioaccumulation and PBPK models. The partition-toolkit screens cheap GFN1/GFN2-xTB single-point energies in 25 ALPB + 12-14 GBSA implicit solvents (78 SP/molecule) and uses each solvent's logK_solvent/water as a surrogate for a biological phase, validated against 440 UFZ-LSER compounds across five phases.

## Method
For each phase, pair the FROZEN pre-selected surrogate logK column with the measured UFZ logK, drop missing rows, and compute R^2 = 1 - SS_res/SS_tot using the surrogate directly as the prediction (the toolkit's reported metric). Judge each R^2 against the frozen 0.70 threshold autonomously. NOTE (limitation): the surrogate was selected on this same set, so these are resubstitution R^2 (optimistic upper bounds), not cross-validated.

Planned approaches:
- GFN1/GFN2-xTB single-point energies with ALPB and GBSA implicit solvation
- solvent-as-phase surrogacy: logK_solvent/water from the body-temperature transfer free energy
- coefficient-of-determination validation against measured UFZ-LSER logK

## Hypotheses and findings

### The pre-registered single-conformer xTB surrogate 'logK_GFN2_gbsa_acetonitrile' predicts phospholipid membrane-water logK with coefficient of determination R^2 >= 0.70 over the UFZ test set
- Hypothesis id: `hyp-membrane` (confirmatory)
- Decision rule (threshold): coefficient of determination R^2 (surrogate logK as a direct, uncalibrated predictor of the measured phase logK) >= 0.70 over the UFZ test set => support (useful predictor for this phase); R^2 < 0.70 => refute
- **Status: supported** — confidence 0.102 (credence)
- Basis: threshold rule: statistic 'point'=0.808065 >= 0.7 is met (combine='latest', margin=0.108065)
- Evidence validity: referent=empirical; data_source(s)=measured

### The pre-registered single-conformer xTB surrogate 'logK_GFN1_alpb_hexane' predicts storage lipid-water logK with coefficient of determination R^2 >= 0.70 over the UFZ test set
- Hypothesis id: `hyp-storage` (confirmatory)
- Decision rule (threshold): coefficient of determination R^2 (surrogate logK as a direct, uncalibrated predictor of the measured phase logK) >= 0.70 over the UFZ test set => support (useful predictor for this phase); R^2 < 0.70 => refute
- **Status: supported** — confidence 0.128 (credence)
- Basis: threshold rule: statistic 'point'=0.836644 >= 0.7 is met (combine='latest', margin=0.136644)
- Evidence validity: referent=empirical; data_source(s)=measured

### The pre-registered single-conformer xTB surrogate 'logK_GFN1_alpb_woctanol' predicts albumin-water logK with coefficient of determination R^2 >= 0.70 over the UFZ test set
- Hypothesis id: `hyp-albumin` (confirmatory)
- Decision rule (threshold): coefficient of determination R^2 (surrogate logK as a direct, uncalibrated predictor of the measured phase logK) >= 0.70 over the UFZ test set => support (useful predictor for this phase); R^2 < 0.70 => refute
- **Status: refuted** — confidence 0.325 (credence)
- Basis: threshold rule: statistic 'point'=0.307376 >= 0.7 is not met (combine='latest', margin=0.392624)
- Evidence validity: referent=empirical; data_source(s)=measured

### The pre-registered single-conformer xTB surrogate 'logK_GFN1_alpb_octanol' predicts muscle protein-water (chicken) logK with coefficient of determination R^2 >= 0.70 over the UFZ test set
- Hypothesis id: `hyp-muscle_chicken` (confirmatory)
- Decision rule (threshold): coefficient of determination R^2 (surrogate logK as a direct, uncalibrated predictor of the measured phase logK) >= 0.70 over the UFZ test set => support (useful predictor for this phase); R^2 < 0.70 => refute
- **Status: refuted** — confidence 0.517 (credence)
- Basis: threshold rule: statistic 'point'=-0.0274823 >= 0.7 is not met (combine='latest', margin=0.727482)
- Evidence validity: referent=empirical; data_source(s)=measured

### The pre-registered single-conformer xTB surrogate 'logK_GFN1_alpb_octanol' predicts muscle protein-water (fish) logK with coefficient of determination R^2 >= 0.70 over the UFZ test set
- Hypothesis id: `hyp-muscle_fish` (confirmatory)
- Decision rule (threshold): coefficient of determination R^2 (surrogate logK as a direct, uncalibrated predictor of the measured phase logK) >= 0.70 over the UFZ test set => support (useful predictor for this phase); R^2 < 0.70 => refute
- **Status: refuted** — confidence 0.507 (credence)
- Basis: threshold rule: statistic 'point'=-0.00787159 >= 0.7 is not met (combine='latest', margin=0.707872)
- Evidence validity: referent=empirical; data_source(s)=measured

## Evidence
- `evi-partition-membrane-20260617-044218-f80e5f72` (experiment_run): point=0.808065, finding={"phase": "membrane", "statistic": "coef_det_R2", "coef_det_R2": 0.808065, "pearson_r2": 0.829772, "rmse_log_units": 0.8
- `evi-partition-storage-20260617-044218-a7ce4dcc` (experiment_run): point=0.836644, finding={"phase": "storage", "statistic": "coef_det_R2", "coef_det_R2": 0.836644, "pearson_r2": 0.859761, "rmse_log_units": 0.75
- `evi-partition-albumin-20260617-044218-f9996827` (experiment_run): point=0.307376, finding={"phase": "albumin", "statistic": "coef_det_R2", "coef_det_R2": 0.307376, "pearson_r2": 0.501277, "rmse_log_units": 0.71
- `evi-partition-muscle_chicken-20260617-044218-1264e6e9` (experiment_run): point=-0.0274823, finding={"phase": "muscle_chicken", "statistic": "coef_det_R2", "coef_det_R2": -0.027482, "pearson_r2": 0.62001, "rmse_log_units
- `evi-partition-muscle_fish-20260617-044218-2eefd2d2` (experiment_run): point=-0.00787159, finding={"phase": "muscle_fish", "statistic": "coef_det_R2", "coef_det_R2": -0.007872, "pearson_r2": 0.564251, "rmse_log_units":
