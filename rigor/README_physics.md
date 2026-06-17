# First-principles COSMO-RS partition prediction (physics run)

This is the **physics** run (`proposal_physics.md`), as opposed to the ML surrogate run.
The partition coefficient is **computed from physical law** — quantum σ-profiles →
COSMO-RS infinite-dilution activity coefficients → `log K = (lnγ_water − lnγ_phase)/ln10`
— with the protein as an explicit amino-acid pseudo-solvent and **zero parameters fitted
to the partition data**. It does NOT regress descriptors; it validates an *a-priori*
physical prediction and diagnoses where the physics fails.

**Scope (honest):** only the muscle-protein phase has a committed first-principles result
in the repo (`structural_protein/cosmo_work_v4/combined_results_46.json`, n=46). The other
four phases need σ-profile generation (ORCA QM), unavailable in this sandbox.

## Result (sci-adk DecisionEngine; `sci-adk verify` exit 0)

Pre-committed, field-standard bars (not data-tuned):

| Hypothesis | statistic | verdict |
|---|---|---|
| **predictive** (coef-det R² ≥ 0) | **coefficient of determination R² = −0.06** | **REFUTED** |
| **useful** (RMSE ≤ 1.0) | RMSE = **0.94** | **SUPPORTED** (lax bar — see note) |
| **unbiased** (\|bias\| ≤ 0.3) | mean signed error = **+0.78** | **REFUTED** |
| **missing-physics** (polar−nonpolar residual ≥ 0.5) | gap = **+0.55** | **SUPPORTED** |

> **Two different "R²" — do not confuse them.** The **coefficient of determination**
> R² = 1 − SS_res/SS_tot (the uncalibrated prediction vs the 1:1 line) is **−0.06**:
> SS_res (40.8) > SS_tot (38.4), so the a-priori prediction is *worse than predicting
> the mean* as a direct quantitative predictor. The **Pearson r²** (correlation only,
> scale/offset-invariant) is **0.66** — the physics captures the right *trend/ranking*
> but not the *absolute values*. The repo's cosmo script (and my first report) called
> the Pearson r² "R²", which over-stated accuracy; corrected here. Also slope = 0.63
> (range compression). The "useful" RMSE ≤ 1.0 bar is lax: the data spread is std = 0.91,
> so RMSE 0.94 barely exceeds it — consistent with the negative coefficient of
> determination. The honest predictive verdict is the **REFUTED** coef-det.

## What the physics says

Parameter-free COSMO-RS **captures the trend** of muscle-protein partitioning (Pearson
r² = 0.66) but, **uncalibrated, is not a quantitatively accurate direct predictor** —
its coefficient of determination is **negative (−0.06)**, i.e. no better than the mean
baseline. The cause is a **systematic over-prediction localized to polar / H-bond
solutes** (polar residual +1.04 vs nonpolar +0.39; worst: TBP +2.56, dibutyl ether
+1.75, N,N-diethylaniline, diazepam, 1-octanol).

**Diagnosed missing physics:** the amino-acid pseudo-solvent models the protein as a
*homogeneous mixture of fully-exposed side chains*. A real folded protein buries most
polar residues or satisfies them by intramolecular H-bonding, so the model **over-counts
available polar/H-bond interactions** → over-predicts protein affinity for polar solutes,
while non-polar solutes (alkanes, PAHs, partitioning by dispersion) are predicted well.
This matches Endo 2011's finding that a homogeneous pseudo-solvent (~RMSE 1.0) is
inadequate and a structure-resolved treatment is needed. **This mechanistic diagnosis is
the deliverable an ML surrogate cannot give.**

## Reproduce

```bash
pip install "git+https://github.com/ccy5123/sci-adk@6f8376d407a49240153e115a94a3dbe7bc14ed0c"
pip install numpy
python rigor/run_physics.py                         # drives `sci-adk run --capability cosmors-physics-muscle`
sci-adk verify rigor/runs/cosmors-physics-muscle    # all REPRODUCED; exit 0
```

## Limitations / next
- **EXPLORATORY** (COSMO-RS results pre-exist in the repo); a confirmatory test needs an
  untouched solute set, and the field-standard bars were set independently of the numbers.
- **Only muscle**; the other four phases need ORCA σ-profiles (not available here). Running
  `cosmo_kpw_openrs_*_v4.py` on a QM-equipped host extends this to all five phases.
- Next physics step: a **structure/depth-resolved** protein model (buried-residue
  weighting or COSMOmic-style) to test whether correcting the over-counted polar exposure
  removes the polar bias.
