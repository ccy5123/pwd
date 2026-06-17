# rigor/ — sci-adk verification of the partition-coefficient findings

This directory runs the partition-toolkit's central empirical question through the
[**sci-adk**](https://github.com/ccy5123/sci-adk) rigor/verification kernel — a
*referee, not a player*. Instead of a hand-maintained R² table in the top-level
README, the per-phase result is a **frozen pre-registration (Spec)** → **measured
Evidence** → an **autonomously judged Claim**, packaged as a run directory a third
party can re-derive **without an LLM** (`sci-adk verify`, exit 0) and that is
**tamper-evident** (a record digest changes if any recorded value is edited).

## The question

Does the pre-registered cheap xTB single-conformer implicit-solvent surrogate predict
each tissue/blood logK well, operationalised as coefficient of determination
**R² ≥ 0.70** over the 440-compound UFZ-LSER set? (Full pre-registration:
[`proposal.md`](proposal.md).)

## Result (autonomous verdicts, reproduced from the record)

| Phase | Frozen surrogate | R² (coef. det.) | n | Verdict |
|-------|------------------|----------------:|--:|---------|
| phospholipid membrane–water | `logK_GFN2_gbsa_acetonitrile` | **0.808** | 207 | **SUPPORTED** |
| storage lipid–water | `logK_GFN1_alpb_hexane` | **0.837** | 247 | **SUPPORTED** |
| albumin–water | `logK_GFN1_alpb_woctanol` | 0.307 | 83 | **REFUTED** |
| muscle protein–water (chicken) | `logK_GFN1_alpb_octanol` | −0.027 | 46 | **REFUTED** |
| muscle protein–water (fish) | `logK_GFN1_alpb_octanol` | −0.008 | 45 | **REFUTED** |

**Headline claim (judged, not asserted):** the xTB implicit-solvent surrogate is a
useful *direct* predictor of **lipid** partitioning (membrane, storage) but **fails for
protein** phases (albumin, muscle) — the uncalibrated surrogate is worse than the mean
there (negative R²). These R² reproduce the top-level README table exactly.

> Caveat, recorded in the Spec: the surrogate per phase was selected on this same set,
> so these are *resubstitution* R² (optimistic upper bounds), not cross-validated. Each
> Evidence item also records the Pearson r² (calibrated linear signal) for transparency.

## Files

- `proposal.md` — the human-readable four-pane pre-registration.
- `partition_capability.py` — the sci-adk capability: builds the frozen `Spec` (one
  confirmatory, empirical hypothesis per phase with a numeric R² ≥ 0.70 threshold rule)
  and the `ExperimentFn` that computes each R² from the **real** repo data and emits one
  `data_source="measured"` Evidence item per phase.
- `run_rigor.py` — driver that compiles the run (zero-LLM deterministic path).
- `runs/partition-xtb-ufz/` — the generated, verifiable artifacts:
  `spec.json` (frozen) · `evidence/*.json` (append-only, measured) ·
  `claims/*.json` (judged belief) · `paper/draft.md` (auto-rendered).

## Reproduce

```bash
# 1. Install the sci-adk kernel + numeric deps (kernel pinned to the commit this run used)
pip install "git+https://github.com/ccy5123/sci-adk@6f8376d407a49240153e115a94a3dbe7bc14ed0c"
pip install numpy pandas

# 2. Re-compile the run from the real data (regenerates runs/partition-xtb-ufz/)
python rigor/run_rigor.py

# 3. Headless, read-only, LLM-free re-verification (third party, CI-style)
sci-adk verify rigor/runs/partition-xtb-ufz     # prints REPRODUCED per claim; exit 0
```

`sci-adk verify` re-applies each frozen `DecisionRule` to the recorded Evidence and
checks the re-derived SUPPORTED/REFUTED matches what was recorded. It re-runs no
experiment, calls no model, and writes nothing. Editing any recorded R² or claim
status makes it report `DIVERGED` and exit 1 (verified: a claim-status flip and an
evidence-value edit are both caught, and the latter changes the record digest).

## Why route this through sci-adk

The toolkit's prior README reduced each phase to a single R² number with no machinery
to stop "build passed = result true." sci-adk separates **record** (append-only,
measured Evidence) from **belief** (revisable Claims), forbids self-certification
(the experiment proposes a direction; the **DecisionEngine** renders the binding
verdict against the frozen threshold), and enforces an **evidence-validity gate**
(an *empirical* claim cannot be affirmed without `measured` data). The negative-R²
protein phases are recorded as first-class **REFUTED** results, not hidden — which is
exactly the honest scientific reporting the toolkit needs.
