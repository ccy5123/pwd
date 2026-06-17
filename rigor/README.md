# rigor/ — advancing the partition research with sci-adk

This directory uses the [**sci-adk**](https://github.com/ccy5123/sci-adk) rigor engine
(*agents propose; the engine judges by frozen criteria*) to **advance** the
partition-toolkit's research, not just re-record it. It contains two runs:

1. **`partition-model-cv`** — the research advance (headline below).
2. **`partition-xtb-ufz`** — the baseline: a faithful, verifiable restatement of the
   toolkit's existing single-surrogate README numbers (the thing the advance improves on).

Both produce a run directory a third party can re-derive **without an LLM**
(`sci-adk verify`, exit 0) and that is **tamper-evident** (editing any recorded value
flips a claim to `DIVERGED` and changes the record digest).

---

## Headline result — protein phases are predictable after all

**Open problem:** the toolkit README concluded protein phases (albumin, muscle) are
*unpredictable* (R² ≈ 0). That used **one uncalibrated** xTB solvent surrogate per
phase, scored by **resubstitution**. The 72 xTB `logK_*` surrogates are already
computed for all 440 UFZ compounds — so we fit a **regularized multi-surrogate ridge
model** and evaluate it by **honest repeated cross-validation** (20× 5-fold,
out-of-fold R²). sci-adk judges two frozen, pre-registered properties per phase:
**useful** (mean CV R² ≥ 0.60) and **robust** (min CV R² ≥ 0.00 across the 20 repeats).

| Phase | n | Baseline R² (toolkit, single surrogate) | **CV R² mean ± std** (this model) | CV R² min | Δ mean | useful / robust |
|-------|--:|--:|--:|--:|--:|---|
| phospholipid membrane–water | 207 | 0.808 | **0.906 ± 0.013** | 0.861 | +0.10 | ✅ / ✅ |
| storage lipid–water | 247 | 0.837 | **0.898 ± 0.065** | 0.701 | +0.06 | ✅ / ✅ |
| albumin–water | 83 | 0.307 | **0.644 ± 0.030** | 0.579 | **+0.34** | ✅ / ✅ |
| muscle protein–water (fish) | 45 | −0.008 | **0.648 ± 0.120** | 0.240 | **+0.66** | ✅ / ✅ |
| muscle protein–water (chicken) | 46 | −0.027 | 0.646 ± **0.300** | **−0.647** | +0.67 | ✅ / ❌ |

**What sci-adk judged (autonomous, reproduced from the record):**

- **The "proteins are unpredictable" conclusion is overturned.** A multi-surrogate CV
  model **rescues albumin** (0.31 → **0.64 ± 0.03**, robust) and makes **muscle (fish)**
  useful and robust (−0.01 → **0.65**, min 0.24) — both REFUTED phases in the baseline
  are now SUPPORTED. The toolkit's near-zero R² was an artifact of the single
  uncalibrated surrogate, not a real ceiling.
- **Lipid phases hold up out-of-fold** (membrane 0.91, storage 0.90 under CV — not just
  resubstitution), so the improvement is honest.
- **muscle (chicken) is flagged, not hidden.** It is `useful` on average (0.65) but
  **`robust` is REFUTED**: with n = 46 the model is unstable across resamples (min
  −0.65). The engine records this as a first-class REFUTED — the honest "promising but
  data-limited" verdict the single-number README could not express.

**Net:** 9/10 frozen hypotheses SUPPORTED; the one REFUTED (muscle-chicken robustness)
is a real, actionable finding (needs more data or stronger regularization/feature
reduction), not a failure of the run.

---

## Files

- `proposal.md` — four-pane pre-registration of the baseline study (EN + KR).
- `partition_capability.py` — baseline capability (single-surrogate validation → `partition-xtb-ufz`).
- `partition_model_capability.py` — **the research advance**: the multi-surrogate
  cross-validated model capability (→ `partition-model-cv`). Frozen Spec (per phase:
  `useful` mean-CV-R² and `robust` min-CV-R² confirmatory hypotheses) + an `ExperimentFn`
  that actually fits the ridge model and runs repeated CV on the real data, emitting one
  `data_source="measured"` Evidence item per hypothesis.
- `run_rigor.py` / `run_model.py` — zero-LLM compile drivers for the two runs.
- `runs/partition-model-cv/` and `runs/partition-xtb-ufz/` — generated, verifiable
  artifacts (`spec.json` · `evidence/*.json` · `claims/*.json` · `paper/draft.md`).

## Frozen evaluation protocol (model run)

```
features  = all 72 logK_* columns   (median-imputed + standardized inside the CV pipeline)
model     = RidgeCV(alphas = logspace(-2, 4, 25))   # alpha chosen on each train fold only
CV        = 5-fold, repeated over seeds 0..19; statistic = out-of-fold R² = 1 - SS_res/SS_tot
judge     = useful: mean R² ≥ 0.60   |   robust: min R² ≥ 0.00   (thresholds frozen in the Spec)
```

Each Evidence item records the full per-seed R² distribution, the single-surrogate
baseline, and the improvement — so the verdict is fully auditable.

## Reproduce

```bash
pip install "git+https://github.com/ccy5123/sci-adk@6f8376d407a49240153e115a94a3dbe7bc14ed0c"
pip install numpy pandas scikit-learn

python rigor/run_model.py                       # the research advance
sci-adk verify rigor/runs/partition-model-cv    # REPRODUCED per claim; exit 0

python rigor/run_rigor.py                        # the baseline restatement
sci-adk verify rigor/runs/partition-xtb-ufz
```

`sci-adk verify` re-applies each frozen `DecisionRule` to the recorded Evidence and
checks the re-derived SUPPORTED/REFUTED matches what was recorded; it re-runs no model,
calls no LLM, and writes nothing. (Reproducing the *numbers* needs the pinned
scikit-learn; re-*verifying* the recorded run does not.)

## Honesty / limitations (recorded in the Specs)

- The model run is **cross-validated** (out-of-fold), which fixes the baseline's
  resubstitution optimism — but the muscle sets (n ≈ 45) remain small, hence the
  robustness REFUTE for muscle-chicken. More measured data or feature reduction is the
  next step.
- Features are the 72 xTB transfer-`logK` surrogates only; descriptor expansion
  (PP-LFER E/S/A/B/V/L, raw σ-profile moments) is unexplored and a natural next run.
- Why route through sci-adk: it separates **record** (append-only measured Evidence)
  from **belief** (revisable Claims), forbids self-certification (the model proposes a
  direction; the **DecisionEngine** renders the binding verdict against the frozen
  threshold), and enforces an **evidence-validity gate** (an empirical claim needs
  `measured` data). The muscle-chicken instability surfaces as a REFUTED claim instead
  of being averaged away.
