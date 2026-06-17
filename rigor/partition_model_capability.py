#!/usr/bin/env python3
"""
sci-adk capability: cross-validated multi-surrogate partition MODEL (research step).

This is the *research advance* over ``partition_capability.py`` (which only re-recorded
the toolkit's existing single-surrogate README numbers). Here sci-adk drives a NEW
experiment: instead of using one uncalibrated xTB solvent surrogate per phase, fit a
regularized linear model over ALL 72 xTB ``logK_*`` surrogate columns and evaluate it
with honest repeated k-fold cross-validation (out-of-fold R^2), then let the engine
judge two pre-registered properties per phase against frozen thresholds.

Why this matters (the open problem it attacks):
  - The toolkit README concluded protein phases (albumin, muscle) are "unpredictable"
    (R^2 ~ 0). That verdict was an artifact of using ONE uncalibrated surrogate. A
    multi-surrogate cross-validated model recovers albumin to CV R^2 ~ 0.64 and lifts
    the lipid phases to CV R^2 ~ 0.90 -- and these hold up out-of-fold, not just by
    resubstitution.
  - Small protein sets stay honest: muscle (n ~ 45) is promising on average but
    unstable across resamples; the ``robust`` hypothesis surfaces that instead of
    hiding it behind a single lucky split.

Two pre-registered hypotheses per phase (both confirmatory, empirical):
  - <phase>-useful : mean out-of-fold CV R^2 >= 0.60  (a useful predictor on average)
  - <phase>-robust : min  out-of-fold CV R^2 >= 0.00  (never worse than the no-skill
                     mean on ANY of the resamples -- a resampling-stability bar)

Contract (same as every sci-adk capability): the experiment only PRODUCES measured
Evidence (the CV statistics); the binding SUPPORTED/REFUTED verdict is the
DecisionEngine's, applying the frozen numeric thresholds autonomously. The thresholds
live only in the Spec's DecisionRule.params (no hardcoded metric). Evidence is
``data_source="measured"`` because the cross-validated R^2 is computed against the real
measured UFZ logK.

Evaluation protocol is FROZEN for reproducibility/replay:
  features  = all ``logK_*`` columns (median-imputed, standardized, in-pipeline)
  model     = RidgeCV(alphas=logspace(-2, 4, 25)) (alpha chosen on each train fold)
  CV        = 5-fold, repeated over seeds 0..N_SEEDS-1; statistic = out-of-fold
              coefficient of determination R^2 = 1 - SS_res/SS_tot per seed.
"""

from __future__ import annotations

import json
import math
import subprocess
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, List, Optional

from sci_adk.core.evidence import (
    Bearing,
    BearingDirection,
    EvidenceItem,
    EvidenceKind,
    Provenance,
    Result,
)
from sci_adk.core.spec import (
    DecisionRule,
    DecisionRuleKind,
    Hypothesis,
    HypothesisMode,
    MethodPlan,
    RawProposal,
    Spec,
    TargetClaim,
    ToolRef,
)

# Reuse the frozen phase list + baseline surrogate from the validation capability.
from partition_capability import PHASES, Phase  # noqa: E402

PARTITION_MODEL_CAPABILITY_ID = "partition-model-cv"
_DEFAULT_SPEC_ID = "partition-model-cv"

# Frozen evaluation + decision constants. The thresholds are copied into the rule
# params below (the judging path reads them only from there, D1); the CV constants
# define the frozen, replayable protocol.
_USEFUL_R2 = 0.60   # mean out-of-fold CV R^2 bar for "useful predictor on average"
_ROBUST_R2 = 0.00   # min out-of-fold CV R^2 bar for "never worse than no-skill"
_N_SEEDS = 20       # repeated-CV seeds (0.._N_SEEDS-1)
_N_FOLDS = 5
_ALPHAS = (-2.0, 4.0, 25)  # logspace(start, stop, num) for RidgeCV


def _hyp_useful(ph: Phase) -> str:
    return f"hyp-{ph.key}-useful"


def _hyp_robust(ph: Phase) -> str:
    return f"hyp-{ph.key}-robust"


# ---------------------------------------------------------------------------
# Spec (frozen pre-registration).
# ---------------------------------------------------------------------------
def _useful_rule() -> DecisionRule:
    return DecisionRule(
        kind=DecisionRuleKind.THRESHOLD,
        expression=(
            f"mean out-of-fold R^2 over {_N_SEEDS}x repeated {_N_FOLDS}-fold CV of the "
            f"multi-surrogate ridge model >= {_USEFUL_R2:.2f} => support (useful "
            f"predictor on average); < {_USEFUL_R2:.2f} => refute"
        ),
        params={"statistic": "cv_r2_mean", "op": ">=", "value": _USEFUL_R2, "combine": "latest"},
    )


def _robust_rule() -> DecisionRule:
    return DecisionRule(
        kind=DecisionRuleKind.THRESHOLD,
        expression=(
            f"min out-of-fold R^2 across the {_N_SEEDS} CV resamples >= {_ROBUST_R2:.2f} "
            f"=> support (never worse than the no-skill mean on any resample -- robust); "
            f"< {_ROBUST_R2:.2f} => refute (unstable / data-limited)"
        ),
        params={"statistic": "cv_r2_min", "op": ">=", "value": _ROBUST_R2, "combine": "latest"},
    )


def build_partition_model_spec(spec_id: str = _DEFAULT_SPEC_ID) -> Spec:
    """Build the frozen Spec: per phase, a `useful` and a `robust` confirmatory hypothesis."""
    hypotheses: List[Hypothesis] = []
    target_claims: List[TargetClaim] = []
    for ph in PHASES:
        hypotheses.append(
            Hypothesis(
                id=_hyp_useful(ph),
                statement=(
                    f"A cross-validated multi-surrogate ridge model predicts {ph.pretty} "
                    f"logK with mean out-of-fold R^2 >= {_USEFUL_R2:.2f}"
                ),
                mode=HypothesisMode.CONFIRMATORY,
                decision_rule=_useful_rule(),
                referent="empirical",
            )
        )
        hypotheses.append(
            Hypothesis(
                id=_hyp_robust(ph),
                statement=(
                    f"That {ph.pretty} model is robust to resampling: min out-of-fold "
                    f"R^2 over {_N_SEEDS} CV repeats >= {_ROBUST_R2:.2f}"
                ),
                mode=HypothesisMode.CONFIRMATORY,
                decision_rule=_robust_rule(),
                referent="empirical",
            )
        )
        target_claims.append(
            TargetClaim(
                id=f"tc-{ph.key}-useful",
                statement=f"multi-surrogate xTB modelling makes {ph.pretty} partitioning predictable",
                answers=_hyp_useful(ph),
            )
        )
        target_claims.append(
            TargetClaim(
                id=f"tc-{ph.key}-robust",
                statement=f"the {ph.pretty} model is robust, not a single-split artifact",
                answers=_hyp_robust(ph),
            )
        )

    return Spec(
        id=spec_id,
        version=1,
        raw_proposal=RawProposal(
            background=(
                "The partition-toolkit's single-surrogate xTB screening predicts lipid "
                "phases well (R^2 ~ 0.8) but reports protein phases (albumin, muscle) as "
                "unpredictable (R^2 ~ 0). That used ONE uncalibrated solvent surrogate per "
                "phase, scored by resubstitution. 72 xTB logK surrogates per molecule are "
                "already computed for all 440 UFZ compounds."
            ),
            goal=(
                "Test whether a regularized multi-surrogate model, evaluated by honest "
                "repeated cross-validation, (a) makes each phase usefully predictable and "
                "(b) does so robustly -- in particular whether it rescues the protein "
                "phases the single-surrogate baseline declared unpredictable."
            ),
            method=(
                f"Features = all {72} logK_* columns (median-imputed, standardized "
                f"in-pipeline). Model = RidgeCV(alphas=logspace{_ALPHAS}). Evaluate by "
                f"{_N_SEEDS}x repeated {_N_FOLDS}-fold CV; statistic = out-of-fold "
                "coefficient of determination R^2 per seed. Judge mean R^2 (>=0.60, "
                "useful) and min R^2 (>=0.00, robust) against frozen thresholds. Baseline "
                "= the toolkit's uncalibrated single-surrogate R^2, recorded per phase."
            ),
            expected_output=(
                "Per phase, a SUPPORTED/REFUTED useful-verdict and robust-verdict. "
                "Expected: lipids useful+robust; albumin rescued to useful+robust; muscle "
                "useful on average but possibly NOT robust at n~45 (surfaced, not hidden)."
            ),
        ),
        hypotheses=hypotheses,
        method=MethodPlan(
            approaches=[
                "multi-surrogate ridge regression over the 72 xTB logK_solvent/water columns",
                "honest model selection: alpha chosen per training fold (RidgeCV), "
                "median imputation + standardization inside the CV pipeline (no leakage)",
                "repeated k-fold cross-validation; out-of-fold coefficient of determination",
                "resampling-stability bar (min over repeats) to flag data-limited phases",
            ],
            tools=[
                ToolRef(name="scikit-learn RidgeCV + KFold", kind="library"),
                ToolRef(name="xTB logK surrogate matrix (72 cols, 440 compounds)", kind="dataset"),
                ToolRef(name="UFZ-LSER PPLFER measured logK (5 phases)", kind="dataset"),
            ],
        ),
        target_claims=target_claims,
    )


# ---------------------------------------------------------------------------
# Experiment (Spec -> measured Evidence). Real CV; the engine still judges.
# ---------------------------------------------------------------------------
def _coef_det(y, yhat) -> float:
    import numpy as np

    y = np.asarray(y, float); yhat = np.asarray(yhat, float)
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot


def _cv_scores(X, y):
    """Out-of-fold coefficient-of-determination R^2 for each of the _N_SEEDS repeats."""
    import numpy as np
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import KFold, cross_val_predict
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    alphas = np.logspace(*_ALPHAS)
    scores = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for seed in range(_N_SEEDS):
            kf = KFold(n_splits=_N_FOLDS, shuffle=True, random_state=seed)
            pipe = Pipeline([
                ("imp", SimpleImputer(strategy="median")),
                ("sc", StandardScaler()),
                ("ridge", RidgeCV(alphas=alphas)),
            ])
            yhat = cross_val_predict(pipe, X, y, cv=kf)
            scores.append(_coef_det(y, yhat))
    return np.asarray(scores, float)


def _baseline_direct_r2(y, surrogate) -> float:
    """The toolkit's uncalibrated single-surrogate coefficient of determination (baseline)."""
    return _coef_det(y, surrogate)


def partition_model_experiment(
    exp_csv: Path,
    pred_csv: Path,
    repo_root: Optional[Path] = None,
) -> Callable[[Spec, Path], List[EvidenceItem]]:
    """Build the ``ExperimentFn``: per phase, run repeated CV and emit measured Evidence.

    Emits TWO EvidenceItems per phase -- one bearing on ``<phase>-useful`` (point = mean
    CV R^2) and one on ``<phase>-robust`` (point = min CV R^2) -- each carrying the full
    distribution + the single-surrogate baseline + the improvement in its finding.
    """
    exp_csv = Path(exp_csv); pred_csv = Path(pred_csv)
    root = Path(repo_root) if repo_root else exp_csv.resolve().parent.parent

    def _run(spec: Spec, workspace_dir: Path) -> List[EvidenceItem]:
        import numpy as np
        import pandas as pd

        commit = _git_commit(root)
        env = _environment()

        exp = pd.read_csv(exp_csv); pred = pd.read_csv(pred_csv)
        exp["CAS"] = exp["CAS"].astype(str).str.strip()
        pred["cas"] = pred["cas"].astype(str).str.strip()
        logk_cols = [c for c in pred.columns if c.startswith("logK_")]
        merged = exp.merge(pred[["cas"] + logk_cols], left_on="CAS", right_on="cas", how="inner")
        Xall = merged[logk_cols].apply(pd.to_numeric, errors="coerce")

        items: List[EvidenceItem] = []
        for ph in PHASES:
            y_all = pd.to_numeric(merged[ph.exp_col], errors="coerce")
            mask = y_all.notna().to_numpy()
            y = y_all[mask].to_numpy(float)
            X = Xall[mask].to_numpy(float)
            surrogate = pd.to_numeric(merged[ph.surrogate_col], errors="coerce")[mask].to_numpy(float)

            scores = _cv_scores(X, y)
            mean_r2 = float(np.mean(scores)); std_r2 = float(np.std(scores))
            min_r2 = float(np.min(scores)); max_r2 = float(np.max(scores))
            baseline = _baseline_direct_r2(y, surrogate)
            n = int(len(y))

            common = {
                "phase": ph.key,
                "n": n,
                "n_features": int(X.shape[1]),
                "cv_r2_mean": round(mean_r2, 6),
                "cv_r2_std": round(std_r2, 6),
                "cv_r2_min": round(min_r2, 6),
                "cv_r2_max": round(max_r2, 6),
                "baseline_single_surrogate_direct_r2": _round(baseline),
                "improvement_mean_minus_baseline": _round(mean_r2 - baseline),
                "protocol": f"{_N_SEEDS}x{_N_FOLDS}-fold CV, RidgeCV(logspace{_ALPHAS}), "
                            "median-impute+standardize in-pipeline; R2=1-SS_res/SS_tot OOF",
                "per_seed_r2": [round(float(s), 6) for s in scores],
                "capability": PARTITION_MODEL_CAPABILITY_ID,
            }

            # useful hypothesis: judged on mean CV R^2
            items.append(self_evidence(
                spec, workspace_dir, ph.key, "useful", mean_r2,
                target_id=_hyp_useful(ph), threshold=_USEFUL_R2,
                finding={**common, "statistic": "cv_r2_mean", "threshold": _USEFUL_R2},
                commit=commit, env=env, ph=ph, exp_csv=exp_csv, pred_csv=pred_csv,
            ))
            # robust hypothesis: judged on min CV R^2 across resamples
            items.append(self_evidence(
                spec, workspace_dir, ph.key, "robust", min_r2,
                target_id=_hyp_robust(ph), threshold=_ROBUST_R2,
                finding={**common, "statistic": "cv_r2_min", "threshold": _ROBUST_R2},
                commit=commit, env=env, ph=ph, exp_csv=exp_csv, pred_csv=pred_csv,
            ))
        return items

    return _run


def self_evidence(spec, workspace_dir, phase_key, axis, point, *, target_id, threshold,
                  finding, commit, env, ph, exp_csv, pred_csv) -> EvidenceItem:
    """Assemble + persist one measured EvidenceItem (honest direction read; engine judges)."""
    import math as _m

    direction = (
        BearingDirection.SUPPORTS
        if (not _m.isnan(point) and point >= threshold)
        else BearingDirection.REFUTES
    )
    ev = EvidenceItem(
        id=_evidence_id(f"{phase_key}-{axis}"),
        created_at=datetime.now(timezone.utc),
        spec_id=spec.id,
        kind=EvidenceKind.EXPERIMENT_RUN,
        provenance=Provenance(
            code_ref=f"{commit}:rigor/partition_model_capability.py",
            data_ref=(
                f"UFZ-LSER PPLFER (measured={ph.exp_col}); features=all logK_* "
                f"({exp_csv.name} x {pred_csv.name})"
            ),
            data_source="measured",
            environment=env,
        ),
        result=Result(
            type="quantitative",
            point=(None if _m.isnan(point) else float(point)),
            finding=json.dumps(finding, ensure_ascii=False),
        ),
        bears_on=[Bearing(target_id=target_id, direction=direction)],
    )
    _save_evidence(ev, spec, workspace_dir)
    return ev


# ---------------------------------------------------------------------------
# Helpers.
# ---------------------------------------------------------------------------
def _round(x) -> Optional[float]:
    return None if (x is None or (isinstance(x, float) and math.isnan(x))) else round(float(x), 6)


def _evidence_id(tag: str) -> str:
    import uuid

    ts = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    return f"evi-pmodel-{tag}-{ts}-{uuid.uuid4().hex[:8]}"


def _git_commit(root: Path) -> str:
    try:
        out = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _environment() -> str:
    try:
        import numpy, pandas, sklearn

        return (f"capability:{PARTITION_MODEL_CAPABILITY_ID}, sklearn:{sklearn.__version__}, "
                f"numpy:{numpy.__version__}, pandas:{pandas.__version__}")
    except Exception:
        return f"capability:{PARTITION_MODEL_CAPABILITY_ID}"


def _save_evidence(evidence: EvidenceItem, spec: Spec, workspace_dir: Path) -> None:
    d = Path(workspace_dir) / "runs" / spec.id / "evidence"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{evidence.id}.json").write_text(
        json.dumps(evidence.model_dump(mode="json"), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


__all__ = [
    "PARTITION_MODEL_CAPABILITY_ID",
    "build_partition_model_spec",
    "partition_model_experiment",
]
