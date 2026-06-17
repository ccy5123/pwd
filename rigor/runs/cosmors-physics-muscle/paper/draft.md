# Quantify how well parameter-free COSMO-RS reproduces measured muscle logK, and diagnose WHERE the physics systematically fails (the missing physics).

> Draft compiled by sci-adk from Spec `cosmors-physics-muscle` (v1). Belief state is revisable as Evidence accrues.

## Goal
Quantify how well parameter-free COSMO-RS reproduces measured muscle logK, and diagnose WHERE the physics systematically fails (the missing physics).

## Background
Tissue partition coefficients are needed a priori for PBPK/bioaccumulation. PP-LFER needs measured Abraham descriptors + coefficients FITTED to partition data. A first-principles route computes log K purely from quantum sigma-profiles + COSMO-RS activity coefficients, with the protein as an explicit amino-acid pseudo-solvent and NO parameter fit to partition data. This run validates that a-priori physical prediction for muscle protein (the one phase with a committed COSMO-RS result; the others need ORCA sigma-profiles, unavailable here).

## Method
Protein = actin+myosin amino-acid composition pseudo-solvent (capped Ace-X-Nme residues); solute sigma-profiles from quantum calc; openCOSMORS24a; ln_gamma_inf in water and in the protein mix at 310.15 K; log K = (ln_gamma_w - ln_gamma_p)/ln10, NO fit to partition data. Validate vs the Endo measured set (n=46): RMSE, mean signed error (bias), R^2, slope, and per-class residuals (polar O/N vs nonpolar). Bars are field-standard and pre-committed: useful RMSE<=1.0, unbiased |bias|<=0.3, localized-gap>=0.5.

Planned approaches:
- first-principles COSMO-RS: quantum sigma-profile -> activity coefficient -> log K
- protein as amino-acid pseudo-solvent at tissue composition (no partition-data fit)
- validation vs measured logK: RMSE / bias / R^2 / slope + per-class residual diagnosis

## Hypotheses and findings

### Parameter-free first-principles COSMO-RS (amino-acid pseudo-solvent) predicts muscle protein-water logK with RMSE <= 1.0 log units
- Hypothesis id: `hyp-muscle-useful` (exploratory)
- Decision rule (threshold): parameter-free COSMO-RS RMSE vs measured logK <= 1.0 log units => support (useful screening); > 1.0 => refute
- **Status: supported** — confidence 0.0565 (credence)
- Basis: threshold rule: statistic 'point'=0.941892 <= 1 is met (combine='latest', margin=0.0581077)
- Evidence validity: referent=empirical; data_source(s)=measured

### That physical prediction has no systematic offset: |mean signed error| <= 0.3 log units
- Hypothesis id: `hyp-muscle-unbiased` (exploratory)
- Decision rule (threshold): |mean signed error (pred-exp)| <= 0.3 => support (no systematic offset); > 0.3 => refute (biased)
- **Status: refuted** — confidence 0.38 (credence)
- Basis: threshold rule: statistic 'point'=0.777505 <= 0.3 is not met (combine='latest', margin=0.477505)
- Evidence validity: referent=empirical; data_source(s)=measured

### The prediction error is a localized systematic over-prediction of polar/H-bond solutes (mean polar residual exceeds mean nonpolar residual by >= 0.5), diagnosing the homogeneous amino-acid pseudo-solvent as the missing physics (exposed polar residues over-counted vs a folded protein)
- Hypothesis id: `hyp-muscle-missing-physics` (exploratory)
- Decision rule (threshold): mean(polar residual) - mean(nonpolar residual) >= 0.5 => support (over-prediction localized to polar/H-bond solutes -- diagnoses the homogeneous-pseudo-solvent missing physics); < 0.5 => refute
- **Status: supported** — confidence 0.0531 (credence)
- Basis: threshold rule: statistic 'point'=0.554595 >= 0.5 is met (combine='latest', margin=0.0545946)
- Evidence validity: referent=empirical; data_source(s)=measured

## Evidence
- `evi-cosmors-priorwork-20260617-065538-923d9d40` (literature): finding={"searched_via": "agent web_search + repo data/README", "key_prior_art": ["Endo, Bauerfeind, Goss 2012, Environ. Sci. Te
- `evi-cosmors-validation-record-20260617-065538-e49e81d1` (observation): finding={"phase": "muscle_protein", "method": "first-principles openCOSMO-RS (no fit)", "n": 46, "rmse": 0.9419, "bias_mean_sign
- `evi-cosmors-hyp-muscle-useful-20260617-065538-ba0e9ee2` (experiment_run): point=0.941892, finding={"phase": "muscle_protein", "method": "first-principles openCOSMO-RS (no fit)", "n": 46, "rmse": 0.9419, "bias_mean_sign
- `evi-cosmors-hyp-muscle-unbiased-20260617-065538-f74730c3` (experiment_run): point=0.777505, finding={"phase": "muscle_protein", "method": "first-principles openCOSMO-RS (no fit)", "n": 46, "rmse": 0.9419, "bias_mean_sign
- `evi-cosmors-hyp-muscle-missing-physics-20260617-065538-a4aa3a9c` (experiment_run): point=0.554595, finding={"phase": "muscle_protein", "method": "first-principles openCOSMO-RS (no fit)", "n": 46, "rmse": 0.9419, "bias_mean_sign
