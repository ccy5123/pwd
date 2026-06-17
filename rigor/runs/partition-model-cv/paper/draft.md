# Test, with honest repeated cross-validation, whether a regularized model over the 72 xTB logK surrogates makes each phase (a) useful (mean CV R^2 >= 0.50) and (b) reliable (95% across-seed interval above the no-skill null).

> Draft compiled by sci-adk from Spec `partition-model-cv` (v1). Belief state is revisable as Evidence accrues.

## Goal
Test, with honest repeated cross-validation, whether a regularized model over the 72 xTB logK surrogates makes each phase (a) useful (mean CV R^2 >= 0.50) and (b) reliable (95% across-seed interval above the no-skill null).

## Background
The partition-toolkit's single-surrogate xTB screening predicts lipid phases well but reports protein phases (albumin, muscle) as unpredictable (R^2 ~ 0), using ONE uncalibrated solvent surrogate scored by resubstitution. PRIOR WORK (recorded as Evidence): multi-descriptor PP-LFER models of these exact phases are established -- Endo/Brown/Goss 2013 (this dataset's source) reaches RMSE 0.32-0.53 log units. So 'predictable by a multi-descriptor model' is NOT novel; the open question is whether a fully in-silico xTB-only surrogate model (no measured Abraham descriptors) can recover the protein phases.

## Method
Features = all 72 logK_* columns (median-imputed, standardized in-pipeline). Model = RidgeCV(alphas=logspace(-2.0, 4.0, 25)), alpha chosen per training fold. Evaluate by 20x repeated 5-fold CV; statistic = out-of-fold coefficient of determination R^2 per seed. Bars are PRINCIPLED and pre-committed (0.50 = majority variance; null=0 = no-skill), NOT tuned to the observed numbers.

Planned approaches:
- multi-surrogate ridge regression over the 72 xTB logK_solvent/water columns
- per-fold alpha selection + in-pipeline imputation/standardization (no leakage)
- repeated k-fold CV; out-of-fold coefficient of determination
- across-seed 95% interval vs a no-skill null as a principled reliability test

## Hypotheses and findings

### A cross-validated multi-surrogate ridge model predicts phospholipid membrane-water logK with mean out-of-fold R^2 >= 0.50
- Hypothesis id: `hyp-membrane-useful` (exploratory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.50 (majority of out-of-sample variance) => support; < 0.50 => refute
- **Status: supported** — confidence 0.334 (credence)
- Basis: threshold rule: statistic 'point'=0.905935 >= 0.5 is met (combine='latest', margin=0.405935)
- Evidence validity: referent=empirical; data_source(s)=measured

### That phospholipid membrane-water model is reliably better than no-skill: the 95% across-seed CV-R^2 interval lies entirely above 0.00
- Hypothesis id: `hyp-membrane-reliable` (exploratory)
- Decision rule (interval): the 95% across-seed interval [p2.5, p97.5] of the CV R^2 lies entirely above the no-skill null 0.00 => support (reliably better than no-skill); interval contains the null => neutral (reliability inconclusive); entirely below => refute
- **Status: supported** — confidence 1 (credence)
- Basis: interval rule: CI=[0.875066, 0.92078] above null_value=0 (support_side='above', combine='latest')
- Evidence validity: referent=empirical; data_source(s)=measured

### A cross-validated multi-surrogate ridge model predicts storage lipid-water logK with mean out-of-fold R^2 >= 0.50
- Hypothesis id: `hyp-storage-useful` (exploratory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.50 (majority of out-of-sample variance) => support; < 0.50 => refute
- **Status: supported** — confidence 0.328 (credence)
- Basis: threshold rule: statistic 'point'=0.89816 >= 0.5 is met (combine='latest', margin=0.39816)
- Evidence validity: referent=empirical; data_source(s)=measured

### That storage lipid-water model is reliably better than no-skill: the 95% across-seed CV-R^2 interval lies entirely above 0.00
- Hypothesis id: `hyp-storage-reliable` (exploratory)
- Decision rule (interval): the 95% across-seed interval [p2.5, p97.5] of the CV R^2 lies entirely above the no-skill null 0.00 => support (reliably better than no-skill); interval contains the null => neutral (reliability inconclusive); entirely below => refute
- **Status: supported** — confidence 0.956 (credence)
- Basis: interval rule: CI=[0.720383, 0.951173] above null_value=0 (support_side='above', combine='latest')
- Evidence validity: referent=empirical; data_source(s)=measured

### A cross-validated multi-surrogate ridge model predicts albumin-water logK with mean out-of-fold R^2 >= 0.50
- Hypothesis id: `hyp-albumin-useful` (exploratory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.50 (majority of out-of-sample variance) => support; < 0.50 => refute
- **Status: supported** — confidence 0.134 (credence)
- Basis: threshold rule: statistic 'point'=0.643803 >= 0.5 is met (combine='latest', margin=0.143803)
- Evidence validity: referent=empirical; data_source(s)=measured

### That albumin-water model is reliably better than no-skill: the 95% across-seed CV-R^2 interval lies entirely above 0.00
- Hypothesis id: `hyp-albumin-reliable` (exploratory)
- Decision rule (interval): the 95% across-seed interval [p2.5, p97.5] of the CV R^2 lies entirely above the no-skill null 0.00 => support (reliably better than no-skill); interval contains the null => neutral (reliability inconclusive); entirely below => refute
- **Status: supported** — confidence 0.997 (credence)
- Basis: interval rule: CI=[0.584901, 0.682786] above null_value=0 (support_side='above', combine='latest')
- Evidence validity: referent=empirical; data_source(s)=measured

### A cross-validated multi-surrogate ridge model predicts muscle protein-water (chicken) logK with mean out-of-fold R^2 >= 0.50
- Hypothesis id: `hyp-muscle_chicken-useful` (exploratory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.50 (majority of out-of-sample variance) => support; < 0.50 => refute
- **Status: supported** — confidence 0.136 (credence)
- Basis: threshold rule: statistic 'point'=0.645867 >= 0.5 is met (combine='latest', margin=0.145867)
- Evidence validity: referent=empirical; data_source(s)=measured

### That muscle protein-water (chicken) model is reliably better than no-skill: the 95% across-seed CV-R^2 interval lies entirely above 0.00
- Hypothesis id: `hyp-muscle_chicken-reliable` (exploratory)
- Decision rule (interval): the 95% across-seed interval [p2.5, p97.5] of the CV R^2 lies entirely above the no-skill null 0.00 => support (reliably better than no-skill); interval contains the null => neutral (reliability inconclusive); entirely below => refute
- **Status: proposed** — confidence 0 (credence)
- Basis: interval rule: CI=[-0.0445381, 0.77377] contains null_value=0 (support_side='above', combine='latest')
- Evidence validity: referent=empirical; data_source(s)=measured

### A cross-validated multi-surrogate ridge model predicts muscle protein-water (fish) logK with mean out-of-fold R^2 >= 0.50
- Hypothesis id: `hyp-muscle_fish-useful` (exploratory)
- Decision rule (threshold): mean out-of-fold R^2 over 20x repeated 5-fold CV of the multi-surrogate ridge model >= 0.50 (majority of out-of-sample variance) => support; < 0.50 => refute
- **Status: supported** — confidence 0.137 (credence)
- Basis: threshold rule: statistic 'point'=0.647872 >= 0.5 is met (combine='latest', margin=0.147872)
- Evidence validity: referent=empirical; data_source(s)=measured

### That muscle protein-water (fish) model is reliably better than no-skill: the 95% across-seed CV-R^2 interval lies entirely above 0.00
- Hypothesis id: `hyp-muscle_fish-reliable` (exploratory)
- Decision rule (interval): the 95% across-seed interval [p2.5, p97.5] of the CV R^2 lies entirely above the no-skill null 0.00 => support (reliably better than no-skill); interval contains the null => neutral (reliability inconclusive); entirely below => refute
- **Status: supported** — confidence 0.584 (credence)
- Basis: interval rule: CI=[0.358402, 0.76698] above null_value=0 (support_side='above', combine='latest')
- Evidence validity: referent=empirical; data_source(s)=measured

## Evidence
- `evi-pmodel-priorwork-20260617-061636-9344de80` (literature): finding={"searched_via": "agent web_search (June 2026)", "key_prior_art": ["Endo, Brown, Goss 2013, Environ. Sci. Technol. 47, 6
- `evi-pmodel-membrane-surrogate-scan-20260617-061636-66deced0` (observation): finding={"phase": "membrane", "what": "full single-surrogate direct-R2 field (all 72)", "purpose": "anti method-shopping: record
- `evi-pmodel-membrane-useful-20260617-061636-2649d5d4` (experiment_run): point=0.905935, finding={"phase": "membrane", "n": 207, "n_features": 72, "cv_r2_mean": 0.905935, "cv_r2_std": 0.013007, "cv_r2_p2.5": 0.875066,
- `evi-pmodel-membrane-reliable-20260617-061636-20cada82` (experiment_run): point=0.905935, ci=[0.875066475, 0.920779875], finding={"phase": "membrane", "n": 207, "n_features": 72, "cv_r2_mean": 0.905935, "cv_r2_std": 0.013007, "cv_r2_p2.5": 0.875066,
- `evi-pmodel-storage-surrogate-scan-20260617-061636-e8ded909` (observation): finding={"phase": "storage", "what": "full single-surrogate direct-R2 field (all 72)", "purpose": "anti method-shopping: record 
- `evi-pmodel-storage-useful-20260617-061636-44b75b31` (experiment_run): point=0.89816, finding={"phase": "storage", "n": 247, "n_features": 72, "cv_r2_mean": 0.89816, "cv_r2_std": 0.065187, "cv_r2_p2.5": 0.720383, "
- `evi-pmodel-storage-reliable-20260617-061636-7a81647f` (experiment_run): point=0.89816, ci=[0.7203826249999999, 0.95117335], finding={"phase": "storage", "n": 247, "n_features": 72, "cv_r2_mean": 0.89816, "cv_r2_std": 0.065187, "cv_r2_p2.5": 0.720383, "
- `evi-pmodel-albumin-surrogate-scan-20260617-061636-9b766c26` (observation): finding={"phase": "albumin", "what": "full single-surrogate direct-R2 field (all 72)", "purpose": "anti method-shopping: record 
- `evi-pmodel-albumin-useful-20260617-061636-2f380d39` (experiment_run): point=0.643803, finding={"phase": "albumin", "n": 83, "n_features": 72, "cv_r2_mean": 0.643803, "cv_r2_std": 0.029585, "cv_r2_p2.5": 0.584901, "
- `evi-pmodel-albumin-reliable-20260617-061636-3c700375` (experiment_run): point=0.643803, ci=[0.584900625, 0.682786125], finding={"phase": "albumin", "n": 83, "n_features": 72, "cv_r2_mean": 0.643803, "cv_r2_std": 0.029585, "cv_r2_p2.5": 0.584901, "
- `evi-pmodel-muscle_chicken-surrogate-scan-20260617-061636-9c4c4625` (observation): finding={"phase": "muscle_chicken", "what": "full single-surrogate direct-R2 field (all 72)", "purpose": "anti method-shopping: 
- `evi-pmodel-muscle_chicken-useful-20260617-061636-0be82e80` (experiment_run): point=0.645867, finding={"phase": "muscle_chicken", "n": 46, "n_features": 72, "cv_r2_mean": 0.645867, "cv_r2_std": 0.299649, "cv_r2_p2.5": -0.0
- `evi-pmodel-muscle_chicken-reliable-20260617-061636-65f33506` (experiment_run): point=0.645867, ci=[-0.044538100000000025, 0.773770125], finding={"phase": "muscle_chicken", "n": 46, "n_features": 72, "cv_r2_mean": 0.645867, "cv_r2_std": 0.299649, "cv_r2_p2.5": -0.0
- `evi-pmodel-muscle_fish-surrogate-scan-20260617-061636-c8df0dbc` (observation): finding={"phase": "muscle_fish", "what": "full single-surrogate direct-R2 field (all 72)", "purpose": "anti method-shopping: rec
- `evi-pmodel-muscle_fish-useful-20260617-061636-c71a042c` (experiment_run): point=0.647872, finding={"phase": "muscle_fish", "n": 45, "n_features": 72, "cv_r2_mean": 0.647872, "cv_r2_std": 0.120152, "cv_r2_p2.5": 0.35840
- `evi-pmodel-muscle_fish-reliable-20260617-061636-100c8d94` (experiment_run): point=0.647872, ci=[0.3584019, 0.766979725], finding={"phase": "muscle_fish", "n": 45, "n_features": 72, "cv_r2_mean": 0.647872, "cv_r2_std": 0.120152, "cv_r2_p2.5": 0.35840
