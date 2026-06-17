# Test whether a regularized multi-surrogate model, evaluated by honest repeated cross-validation, (a) makes each phase usefully predictable and (b) does so robustly -- in particular whether it rescues the protein phases the single-surrogate baseline declared unpredictable.

> Draft compiled by sci-adk from Spec `partition-model-cv` (v1). Belief state is revisable as Evidence accrues.

## Goal
Test whether a regularized multi-surrogate model, evaluated by honest repeated cross-validation, (a) makes each phase usefully predictable and (b) does so robustly -- in particular whether it rescues the protein phases the single-surrogate baseline declared unpredictable.

## Background
The partition-toolkit's single-surrogate xTB screening predicts lipid phases well (R^2 ~ 0.8) but reports protein phases (albumin, muscle) as unpredictable (R^2 ~ 0). That used ONE uncalibrated solvent surrogate per phase, scored by resubstitution. 72 xTB logK surrogates per molecule are already computed for all 440 UFZ compounds.

## Method
Features = all 72 logK_* columns (median-imputed, standardized in-pipeline). Model = RidgeCV(alphas=logspace(-2.0, 4.0, 25)). Evaluate by 20x repeated 5-fold CV; statistic = out-of-fold coefficient of determination R^2 per seed. Judge mean R^2 (>=0.60, useful) and min R^2 (>=0.00, robust) against frozen thresholds. Baseline = the toolkit's uncalibrated single-surrogate R^2, recorded per phase.

Planned approaches:
- multi-surrogate ridge regression over the 72 xTB logK_solvent/water columns
- honest model selection: alpha chosen per training fold (RidgeCV), median imputation + standardization inside the CV pipeline (no leakage)
- repeated k-fold cross-validation; out-of-fold coefficient of determination
- resampling-stability bar (min over repeats) to flag data-limited phases

## Hypotheses and findings

### A cross-validated multi-surrogate ridge model predicts phospholipid membrane-water logK with mean out-of-fold R^2 >= 0.60
- Hypothesis id: `hyp-membrane-useful` (confirmatory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.60 => support (useful predictor on average); < 0.60 => refute
- **Status: supported** — confidence 0.264 (credence)
- Basis: threshold rule: statistic 'point'=0.905935 >= 0.6 is met (combine='latest', margin=0.305935)
- Evidence validity: referent=empirical; data_source(s)=measured

### That phospholipid membrane-water model is robust to resampling: min out-of-fold R^2 over 20 CV repeats >= 0.00
- Hypothesis id: `hyp-membrane-robust` (confirmatory)
- Decision rule (threshold): min out-of-fold R^2 across the 20 CV resamples >= 0.00 => support (never worse than the no-skill mean on any resample -- robust); < 0.00 => refute (unstable / data-limited)
- **Status: supported** — confidence 0.577 (credence)
- Basis: threshold rule: statistic 'point'=0.860854 >= 0 is met (combine='latest', margin=0.860854)
- Evidence validity: referent=empirical; data_source(s)=measured

### A cross-validated multi-surrogate ridge model predicts storage lipid-water logK with mean out-of-fold R^2 >= 0.60
- Hypothesis id: `hyp-storage-useful` (confirmatory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.60 => support (useful predictor on average); < 0.60 => refute
- **Status: supported** — confidence 0.258 (credence)
- Basis: threshold rule: statistic 'point'=0.89816 >= 0.6 is met (combine='latest', margin=0.29816)
- Evidence validity: referent=empirical; data_source(s)=measured

### That storage lipid-water model is robust to resampling: min out-of-fold R^2 over 20 CV repeats >= 0.00
- Hypothesis id: `hyp-storage-robust` (confirmatory)
- Decision rule (threshold): min out-of-fold R^2 across the 20 CV resamples >= 0.00 => support (never worse than the no-skill mean on any resample -- robust); < 0.00 => refute (unstable / data-limited)
- **Status: supported** — confidence 0.504 (credence)
- Basis: threshold rule: statistic 'point'=0.701233 >= 0 is met (combine='latest', margin=0.701233)
- Evidence validity: referent=empirical; data_source(s)=measured

### A cross-validated multi-surrogate ridge model predicts albumin-water logK with mean out-of-fold R^2 >= 0.60
- Hypothesis id: `hyp-albumin-useful` (confirmatory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.60 => support (useful predictor on average); < 0.60 => refute
- **Status: supported** — confidence 0.0429 (credence)
- Basis: threshold rule: statistic 'point'=0.643803 >= 0.6 is met (combine='latest', margin=0.0438025)
- Evidence validity: referent=empirical; data_source(s)=measured

### That albumin-water model is robust to resampling: min out-of-fold R^2 over 20 CV repeats >= 0.00
- Hypothesis id: `hyp-albumin-robust` (confirmatory)
- Decision rule (threshold): min out-of-fold R^2 across the 20 CV resamples >= 0.00 => support (never worse than the no-skill mean on any resample -- robust); < 0.00 => refute (unstable / data-limited)
- **Status: supported** — confidence 0.44 (credence)
- Basis: threshold rule: statistic 'point'=0.579336 >= 0 is met (combine='latest', margin=0.579336)
- Evidence validity: referent=empirical; data_source(s)=measured

### A cross-validated multi-surrogate ridge model predicts muscle protein-water (chicken) logK with mean out-of-fold R^2 >= 0.60
- Hypothesis id: `hyp-muscle_chicken-useful` (confirmatory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.60 => support (useful predictor on average); < 0.60 => refute
- **Status: supported** — confidence 0.0448 (credence)
- Basis: threshold rule: statistic 'point'=0.645867 >= 0.6 is met (combine='latest', margin=0.0458675)
- Evidence validity: referent=empirical; data_source(s)=measured

### That muscle protein-water (chicken) model is robust to resampling: min out-of-fold R^2 over 20 CV repeats >= 0.00
- Hypothesis id: `hyp-muscle_chicken-robust` (confirmatory)
- Decision rule (threshold): min out-of-fold R^2 across the 20 CV resamples >= 0.00 => support (never worse than the no-skill mean on any resample -- robust); < 0.00 => refute (unstable / data-limited)
- **Status: refuted** — confidence 0.476 (credence)
- Basis: threshold rule: statistic 'point'=-0.646688 >= 0 is not met (combine='latest', margin=0.646688)
- Evidence validity: referent=empirical; data_source(s)=measured

### A cross-validated multi-surrogate ridge model predicts muscle protein-water (fish) logK with mean out-of-fold R^2 >= 0.60
- Hypothesis id: `hyp-muscle_fish-useful` (confirmatory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.60 => support (useful predictor on average); < 0.60 => refute
- **Status: supported** — confidence 0.0467 (credence)
- Basis: threshold rule: statistic 'point'=0.647872 >= 0.6 is met (combine='latest', margin=0.0478717)
- Evidence validity: referent=empirical; data_source(s)=measured

### That muscle protein-water (fish) model is robust to resampling: min out-of-fold R^2 over 20 CV repeats >= 0.00
- Hypothesis id: `hyp-muscle_fish-robust` (confirmatory)
- Decision rule (threshold): min out-of-fold R^2 across the 20 CV resamples >= 0.00 => support (never worse than the no-skill mean on any resample -- robust); < 0.00 => refute (unstable / data-limited)
- **Status: supported** — confidence 0.213 (credence)
- Basis: threshold rule: statistic 'point'=0.239669 >= 0 is met (combine='latest', margin=0.239669)
- Evidence validity: referent=empirical; data_source(s)=measured

## Evidence
- `evi-pmodel-membrane-useful-20260617-052019-66004764` (experiment_run): point=0.905935, finding={"phase": "membrane", "n": 207, "n_features": 72, "cv_r2_mean": 0.905935, "cv_r2_std": 0.013007, "cv_r2_min": 0.860854, 
- `evi-pmodel-membrane-robust-20260617-052019-1a28b7b0` (experiment_run): point=0.860854, finding={"phase": "membrane", "n": 207, "n_features": 72, "cv_r2_mean": 0.905935, "cv_r2_std": 0.013007, "cv_r2_min": 0.860854, 
- `evi-pmodel-storage-useful-20260617-052021-7a4faa16` (experiment_run): point=0.89816, finding={"phase": "storage", "n": 247, "n_features": 72, "cv_r2_mean": 0.89816, "cv_r2_std": 0.065187, "cv_r2_min": 0.701233, "c
- `evi-pmodel-storage-robust-20260617-052021-ac735069` (experiment_run): point=0.701233, finding={"phase": "storage", "n": 247, "n_features": 72, "cv_r2_mean": 0.89816, "cv_r2_std": 0.065187, "cv_r2_min": 0.701233, "c
- `evi-pmodel-albumin-useful-20260617-052021-cf7f0b87` (experiment_run): point=0.643803, finding={"phase": "albumin", "n": 83, "n_features": 72, "cv_r2_mean": 0.643803, "cv_r2_std": 0.029585, "cv_r2_min": 0.579336, "c
- `evi-pmodel-albumin-robust-20260617-052021-6f0f86f4` (experiment_run): point=0.579336, finding={"phase": "albumin", "n": 83, "n_features": 72, "cv_r2_mean": 0.643803, "cv_r2_std": 0.029585, "cv_r2_min": 0.579336, "c
- `evi-pmodel-muscle_chicken-useful-20260617-052022-23067a1c` (experiment_run): point=0.645867, finding={"phase": "muscle_chicken", "n": 46, "n_features": 72, "cv_r2_mean": 0.645867, "cv_r2_std": 0.299649, "cv_r2_min": -0.64
- `evi-pmodel-muscle_chicken-robust-20260617-052022-2975e218` (experiment_run): point=-0.646688, finding={"phase": "muscle_chicken", "n": 46, "n_features": 72, "cv_r2_mean": 0.645867, "cv_r2_std": 0.299649, "cv_r2_min": -0.64
- `evi-pmodel-muscle_fish-useful-20260617-052023-4a6ad8b4` (experiment_run): point=0.647872, finding={"phase": "muscle_fish", "n": 45, "n_features": 72, "cv_r2_mean": 0.647872, "cv_r2_std": 0.120152, "cv_r2_min": 0.239669
- `evi-pmodel-muscle_fish-robust-20260617-052023-ddac8fe4` (experiment_run): point=0.239669, finding={"phase": "muscle_fish", "n": 45, "n_features": 72, "cv_r2_mean": 0.647872, "cv_r2_std": 0.120152, "cv_r2_min": 0.239669
