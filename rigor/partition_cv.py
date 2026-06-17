#!/usr/bin/env python3
"""
Pure partition cross-validation science — NO sci-adk dependency.

This module is deliberately free of any `sci_adk` import so it can run INSIDE the
Docker experiment container (which has only numpy/pandas/scikit-learn, not the sci-adk
kernel). It is the "science" half of the execution seam (mirroring sci-adk's split of
pure `t1_encoding.py` from the wiring in `t1_capability.py`): the container and the
in-process fallback both call `compute_cv_payload`, so the statistics are identical and
only the *container* is a seam — the science is never faked.

The frozen evaluation protocol (kept in sync with partition_model_capability.py):
  features = all logK_* columns (median-imputed + standardized in-pipeline)
  model    = RidgeCV(alphas=logspace(-2, 4, 25)), alpha chosen per training fold
  CV       = N_FOLDS-fold, repeated over seeds 0..N_SEEDS-1; out-of-fold R2 per seed.
"""

from __future__ import annotations

import json
import sys
import warnings
from typing import Dict, List, Tuple

# Frozen CV protocol constants (the single source of truth; partition_model_capability
# imports these so the Spec text and the computation cannot drift).
N_SEEDS = 20
N_FOLDS = 5
ALPHAS = (-2.0, 4.0, 25)  # logspace(start, stop, num)

# The five phases: (key, measured column, frozen single-surrogate column for the baseline).
PHASE_COLUMNS: Tuple[Tuple[str, str, str], ...] = (
    ("membrane", "membrane_exp", "logK_GFN2_gbsa_acetonitrile"),
    ("storage", "storage_exp", "logK_GFN1_alpb_hexane"),
    ("albumin", "albumin_exp", "logK_GFN1_alpb_woctanol"),
    ("muscle_chicken", "muscle_chicken_exp", "logK_GFN1_alpb_octanol"),
    ("muscle_fish", "muscle_fish_exp", "logK_GFN1_alpb_octanol"),
)


def coef_det(pred, exp) -> float:
    """R^2 = 1 - SS_res/SS_tot (can be negative; the toolkit's reported metric)."""
    import numpy as np

    pred = np.asarray(pred, float)
    exp = np.asarray(exp, float)
    ss_res = float(np.sum((exp - pred) ** 2))
    ss_tot = float(np.sum((exp - np.mean(exp)) ** 2))
    return float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot


def cv_scores(X, y) -> Tuple[list, list]:
    """Out-of-fold R^2 and RMSE for each of the N_SEEDS repeated CV runs."""
    import numpy as np
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import KFold, cross_val_predict
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    alphas = np.logspace(*ALPHAS)
    r2s, rmses = [], []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for seed in range(N_SEEDS):
            kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
            pipe = Pipeline([
                ("imp", SimpleImputer(strategy="median")),
                ("sc", StandardScaler()),
                ("ridge", RidgeCV(alphas=alphas)),
            ])
            yhat = cross_val_predict(pipe, X, y, cv=kf)
            r2s.append(coef_det(yhat, y))
            rmses.append(float(np.sqrt(np.mean((np.asarray(y, float) - yhat) ** 2))))
    return r2s, rmses


def compute_cv_payload(exp_csv: str, pred_csv: str) -> Dict[str, dict]:
    """Compute every per-phase statistic the capability needs, as JSON-serializable data.

    Returns ``{phase_key: {n, n_features, r2s, rmses, baseline_direct_r2, field}}`` where
    ``field`` is the FULL single-surrogate direct-R^2 list (all logK_* columns), sorted
    desc — the anti-method-shopping record of every attempt, not just the winner.
    """
    import numpy as np
    import pandas as pd

    exp = pd.read_csv(exp_csv)
    pred = pd.read_csv(pred_csv)
    exp["CAS"] = exp["CAS"].astype(str).str.strip()
    pred["cas"] = pred["cas"].astype(str).str.strip()
    logk_cols = [c for c in pred.columns if c.startswith("logK_")]
    merged = exp.merge(pred[["cas"] + logk_cols], left_on="CAS", right_on="cas", how="inner")
    Xall = merged[logk_cols].apply(pd.to_numeric, errors="coerce")

    out: Dict[str, dict] = {}
    for key, exp_col, surrogate_col in PHASE_COLUMNS:
        y_all = pd.to_numeric(merged[exp_col], errors="coerce")
        mask = y_all.notna().to_numpy()
        y = y_all[mask].to_numpy(float)
        X = Xall[mask].to_numpy(float)
        surrogate = pd.to_numeric(merged[surrogate_col], errors="coerce")[mask].to_numpy(float)

        r2s, rmses = cv_scores(X, y)

        field: List[Tuple[str, float]] = []
        for col in logk_cols:
            s = pd.to_numeric(merged[col], errors="coerce")[mask].to_numpy(float)
            if np.isfinite(s).all() and np.std(s) > 0:
                field.append((col, round(coef_det(s, y), 6)))
        field.sort(key=lambda t: t[1], reverse=True)

        out[key] = {
            "n": int(len(y)),
            "n_features": int(X.shape[1]),
            "r2s": [round(float(v), 6) for v in r2s],
            "rmses": [round(float(v), 6) for v in rmses],
            "baseline_direct_r2": round(coef_det(surrogate, y), 6),
            "field": field,
        }
    return out


# Container entry point: `python partition_cv.py <exp_csv> <pred_csv>` prints the payload
# JSON on the last stdout line (this is what the Docker executor parses back).
if __name__ == "__main__":
    payload = compute_cv_payload(sys.argv[1], sys.argv[2])
    print(json.dumps(payload))
