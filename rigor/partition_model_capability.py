#!/usr/bin/env python3
"""
sci-adk capability: cross-validated multi-surrogate partition model (v2, corrected).

This supersedes the v1 design (which set its thresholds AFTER inspecting the CV
numbers -- a HARKing violation of the frozen-Spec / anti-HARKing principle, and
labelled the result ``confirmatory``). The corrections here:

  1. PRINCIPLED, NON-TUNED bars (not reverse-engineered from the data):
       - useful   (threshold): mean out-of-fold CV R^2 >= 0.50  (explains the
                   MAJORITY of out-of-sample variance -- a standard round bar, not
                   the v1 0.60 that hugged albumin/muscle's observed 0.644/0.646).
       - reliable (interval) : the 95% across-seed interval [p2.5, p97.5] of the CV
                   R^2 lies entirely ABOVE the no-skill null 0 (a textbook
                   significance criterion; null=0, not the v1 "min >= 0" picked to
                   single out muscle-chicken).
  2. HONEST mode = EXPLORATORY. The data has already been examined on this fixed
     dataset, so genuine confirmatory pre-registration is no longer possible here;
     true confirmation needs an untouched external/validation set (stated as a
     limitation in the Spec). The numeric rules are still judged autonomously.
  3. PRIOR WORK is recorded (a LITERATURE EvidenceItem) and the novelty is reframed:
     multi-descriptor modelling of these exact phases is established PP-LFER work
     (Endo/Brown/Goss 2013, the dataset's own source paper, RMSE 0.32-0.53). This
     run does NOT beat that; its honest contribution is that a fully in-silico
     xTB-only surrogate model (no measured Abraham descriptors) recovers the protein
     phases the toolkit's SINGLE-surrogate route reported as unpredictable.

Contract unchanged: the experiment only PRODUCES measured Evidence (the CV
statistics); the DecisionEngine renders the binding verdict against the frozen rules.
Thresholds live only in DecisionRule.params (no hardcoded metric). Evidence is
``data_source="measured"`` (CV R^2 computed against real measured UFZ logK).
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

from partition_capability import PHASES, Phase  # noqa: E402

PARTITION_MODEL_CAPABILITY_ID = "partition-model-cv"
_DEFAULT_SPEC_ID = "partition-model-cv"

# Principled, pre-committed bars (copied into rule params; judging path reads them
# only from there, D1). 0.50 = majority of out-of-sample variance; 0.0 = no-skill null.
_USEFUL_R2 = 0.50
_NULL_R2 = 0.00
_N_SEEDS = 20
_N_FOLDS = 5
_ALPHAS = (-2.0, 4.0, 25)

# Prior art recorded at Spec time (resolves the discovery checkpoint substance).
_PRIOR_WORK = {
    "searched_via": "agent web_search (June 2026)",
    "key_prior_art": [
        "Endo, Brown, Goss 2013, Environ. Sci. Technol. 47, 6630-6639, "
        "doi:10.1021/es401772m -- the DATASET'S OWN SOURCE: PP-LFER models for "
        "membrane/storage lipid, albumin, and structural protein vs water; "
        "composition model RMSE 0.32-0.53 log units (the real benchmark).",
        "Streamlining Linear Free Energy Relationships of Proteins (JCIM 2024, "
        "doi:10.1021/acs.jcim.4c01289) -- linear modelling of protein LFERs.",
        "PP-LFER for storage lipids (Environ. Sci. Technol. 2024, "
        "doi:10.1021/acs.est.4c01994).",
        "Ehlert, Stahn, Spicher, Grimme 2021, JCTC 17, 4250-4261, "
        "doi:10.1021/acs.jctc.1c00471 -- GFN2-xTB(ALPB) logKow MAD ~0.65 (basis of "
        "the surrogate approach).",
    ],
    "novelty_assessment": (
        "Multi-descriptor modelling of these phases is NOT novel (established "
        "PP-LFER work, better RMSE). Honest contribution: a fully in-silico xTB-only "
        "surrogate model (no measured Abraham descriptors) recovers protein-phase "
        "predictability that the toolkit's single uncalibrated surrogate missed."
    ),
}


def _hyp_useful(ph: Phase) -> str:
    return f"hyp-{ph.key}-useful"


def _hyp_reliable(ph: Phase) -> str:
    return f"hyp-{ph.key}-reliable"


def _useful_rule() -> DecisionRule:
    return DecisionRule(
        kind=DecisionRuleKind.THRESHOLD,
        expression=(
            f"mean out-of-fold R^2 over {_N_SEEDS}x repeated {_N_FOLDS}-fold CV of the "
            f"multi-surrogate ridge model >= {_USEFUL_R2:.2f} (majority of out-of-sample "
            f"variance) => support; < {_USEFUL_R2:.2f} => refute"
        ),
        params={"statistic": "cv_r2_mean", "op": ">=", "value": _USEFUL_R2, "combine": "latest"},
    )


def _reliable_rule() -> DecisionRule:
    return DecisionRule(
        kind=DecisionRuleKind.INTERVAL,
        expression=(
            f"the 95% across-seed interval [p2.5, p97.5] of the CV R^2 lies entirely "
            f"above the no-skill null {_NULL_R2:.2f} => support (reliably better than "
            f"no-skill); interval contains the null => neutral (reliability "
            f"inconclusive); entirely below => refute"
        ),
        params={"null_value": _NULL_R2, "support_side": "above", "combine": "latest"},
    )


def build_partition_model_spec(spec_id: str = _DEFAULT_SPEC_ID) -> Spec:
    hypotheses: List[Hypothesis] = []
    target_claims: List[TargetClaim] = []
    for ph in PHASES:
        hypotheses.append(Hypothesis(
            id=_hyp_useful(ph),
            statement=(
                f"A cross-validated multi-surrogate ridge model predicts {ph.pretty} "
                f"logK with mean out-of-fold R^2 >= {_USEFUL_R2:.2f}"
            ),
            mode=HypothesisMode.EXPLORATORY,
            decision_rule=_useful_rule(),
            referent="empirical",
        ))
        hypotheses.append(Hypothesis(
            id=_hyp_reliable(ph),
            statement=(
                f"That {ph.pretty} model is reliably better than no-skill: the 95% "
                f"across-seed CV-R^2 interval lies entirely above {_NULL_R2:.2f}"
            ),
            mode=HypothesisMode.EXPLORATORY,
            decision_rule=_reliable_rule(),
            referent="empirical",
        ))
        target_claims.append(TargetClaim(
            id=f"tc-{ph.key}-useful",
            statement=f"in-silico xTB surrogate modelling makes {ph.pretty} partitioning predictable",
            answers=_hyp_useful(ph),
        ))
        target_claims.append(TargetClaim(
            id=f"tc-{ph.key}-reliable",
            statement=f"the {ph.pretty} model's skill is reliable across resampling",
            answers=_hyp_reliable(ph),
        ))

    return Spec(
        id=spec_id,
        version=1,
        raw_proposal=RawProposal(
            background=(
                "The partition-toolkit's single-surrogate xTB screening predicts lipid "
                "phases well but reports protein phases (albumin, muscle) as "
                "unpredictable (R^2 ~ 0), using ONE uncalibrated solvent surrogate scored "
                "by resubstitution. PRIOR WORK (recorded as Evidence): multi-descriptor "
                "PP-LFER models of these exact phases are established -- Endo/Brown/Goss "
                "2013 (this dataset's source) reaches RMSE 0.32-0.53 log units. So "
                "'predictable by a multi-descriptor model' is NOT novel; the open "
                "question is whether a fully in-silico xTB-only surrogate model (no "
                "measured Abraham descriptors) can recover the protein phases."
            ),
            goal=(
                "Test, with honest repeated cross-validation, whether a regularized "
                "model over the 72 xTB logK surrogates makes each phase (a) useful "
                "(mean CV R^2 >= 0.50) and (b) reliable (95% across-seed interval above "
                "the no-skill null)."
            ),
            method=(
                f"Features = all 72 logK_* columns (median-imputed, standardized "
                f"in-pipeline). Model = RidgeCV(alphas=logspace{_ALPHAS}), alpha chosen "
                f"per training fold. Evaluate by {_N_SEEDS}x repeated {_N_FOLDS}-fold CV; "
                "statistic = out-of-fold coefficient of determination R^2 per seed. "
                "Bars are PRINCIPLED and pre-committed (0.50 = majority variance; null=0 "
                "= no-skill), NOT tuned to the observed numbers."
            ),
            expected_output=(
                "Per phase, a useful-verdict (SUPPORTED/REFUTED) and a reliable-verdict "
                "(SUPPORTED/NEUTRAL/REFUTED). LIMITATION: the dataset has been examined, "
                "so these hypotheses are EXPLORATORY, not confirmatory -- genuine "
                "confirmation needs an untouched external validation set. Performance is "
                "also below the PP-LFER benchmark (RMSE 0.32-0.53)."
            ),
        ),
        hypotheses=hypotheses,
        method=MethodPlan(
            approaches=[
                "multi-surrogate ridge regression over the 72 xTB logK_solvent/water columns",
                "per-fold alpha selection + in-pipeline imputation/standardization (no leakage)",
                "repeated k-fold CV; out-of-fold coefficient of determination",
                "across-seed 95% interval vs a no-skill null as a principled reliability test",
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
# Experiment.
# ---------------------------------------------------------------------------
def _coef_det(y, yhat) -> float:
    import numpy as np

    y = np.asarray(y, float); yhat = np.asarray(yhat, float)
    ss_res = float(np.sum((y - yhat) ** 2)); ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot


def _cv_scores(X, y):
    import numpy as np
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import KFold, cross_val_predict
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    alphas = np.logspace(*_ALPHAS)
    r2s, rmses = [], []
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
            r2s.append(_coef_det(y, yhat))
            rmses.append(float(np.sqrt(np.mean((np.asarray(y, float) - yhat) ** 2))))
    return np.asarray(r2s, float), np.asarray(rmses, float)


def partition_model_experiment(
    exp_csv: Path,
    pred_csv: Path,
    repo_root: Optional[Path] = None,
) -> Callable[[Spec, Path], List[EvidenceItem]]:
    exp_csv = Path(exp_csv); pred_csv = Path(pred_csv)
    root = Path(repo_root) if repo_root else exp_csv.resolve().parent.parent

    def _run(spec: Spec, workspace_dir: Path) -> List[EvidenceItem]:
        import numpy as np
        import pandas as pd

        commit = _git_commit(root); env = _environment()
        items: List[EvidenceItem] = []

        # (0) Prior-work as a LITERATURE record (append-only; resolves discovery substance).
        lit = EvidenceItem(
            id=_evidence_id("priorwork"),
            created_at=datetime.now(timezone.utc),
            spec_id=spec.id,
            kind=EvidenceKind.LITERATURE,
            provenance=Provenance(
                code_ref=f"{commit}:rigor/partition_model_capability.py",
                data_ref="prior-art web_search (June 2026)",
                environment=env,
            ),
            result=Result(type="qualitative", finding=json.dumps(_PRIOR_WORK, ensure_ascii=False)),
            bears_on=[],
        )
        _save_evidence(lit, spec, workspace_dir); items.append(lit)

        exp = pd.read_csv(exp_csv); pred = pd.read_csv(pred_csv)
        exp["CAS"] = exp["CAS"].astype(str).str.strip()
        pred["cas"] = pred["cas"].astype(str).str.strip()
        logk_cols = [c for c in pred.columns if c.startswith("logK_")]
        merged = exp.merge(pred[["cas"] + logk_cols], left_on="CAS", right_on="cas", how="inner")
        Xall = merged[logk_cols].apply(pd.to_numeric, errors="coerce")

        for ph in PHASES:
            y_all = pd.to_numeric(merged[ph.exp_col], errors="coerce")
            mask = y_all.notna().to_numpy()
            y = y_all[mask].to_numpy(float); X = Xall[mask].to_numpy(float)
            surrogate = pd.to_numeric(merged[ph.surrogate_col], errors="coerce")[mask].to_numpy(float)

            # Anti-method-shopping (F3): record the FULL single-surrogate field -- every
            # one of the 72 surrogates' direct R^2 -- as one measured OBSERVATION, so the
            # toolkit's "best surrogate per phase" cannot be cherry-picked: the whole field
            # (winners AND losers) is on the append-only record, not just the winner.
            field = []
            for col in logk_cols:
                s = pd.to_numeric(merged[col], errors="coerce")[mask].to_numpy(float)
                if np.isfinite(s).all() and np.std(s) > 0:
                    field.append((col, round(_coef_det(s, y), 6)))
            field.sort(key=lambda t: t[1], reverse=True)
            scan_ev = EvidenceItem(
                id=_evidence_id(f"{ph.key}-surrogate-scan"),
                created_at=datetime.now(timezone.utc), spec_id=spec.id,
                kind=EvidenceKind.OBSERVATION,
                provenance=_prov(commit, env, ph, exp_csv, pred_csv),
                result=Result(type="qualitative", finding=json.dumps({
                    "phase": ph.key, "what": "full single-surrogate direct-R2 field (all 72)",
                    "purpose": "anti method-shopping: record every attempt, not just the best",
                    "n": int(len(y)), "n_surrogates": len(field),
                    "best": field[0] if field else None,
                    "worst": field[-1] if field else None,
                    "all_surrogates_r2": field,
                }, ensure_ascii=False)),
                bears_on=[],
            )
            _save_evidence(scan_ev, spec, workspace_dir); items.append(scan_ev)

            r2s, rmses = _cv_scores(X, y)
            mean_r2 = float(np.mean(r2s)); std_r2 = float(np.std(r2s))
            p2_5, p97_5 = (float(v) for v in np.percentile(r2s, [2.5, 97.5]))
            baseline = _coef_det(y, surrogate)
            n = int(len(y))

            common = {
                "phase": ph.key, "n": n, "n_features": int(X.shape[1]),
                "cv_r2_mean": round(mean_r2, 6), "cv_r2_std": round(std_r2, 6),
                "cv_r2_p2.5": round(p2_5, 6), "cv_r2_p97.5": round(p97_5, 6),
                "cv_rmse_mean": round(float(np.mean(rmses)), 6),
                "baseline_single_surrogate_direct_r2": _round(baseline),
                "ppLFER_benchmark_rmse": "0.32-0.53 (Endo 2013, measured Abraham descriptors)",
                "per_seed_r2": [round(float(s), 6) for s in r2s],
                "protocol": f"{_N_SEEDS}x{_N_FOLDS}-fold CV, RidgeCV(logspace{_ALPHAS}), "
                            "median-impute+standardize in-pipeline; R2=1-SS_res/SS_tot OOF",
                "capability": PARTITION_MODEL_CAPABILITY_ID,
            }

            # useful (threshold on mean CV R^2)
            d_use = BearingDirection.SUPPORTS if mean_r2 >= _USEFUL_R2 else BearingDirection.REFUTES
            ev_use = EvidenceItem(
                id=_evidence_id(f"{ph.key}-useful"),
                created_at=datetime.now(timezone.utc), spec_id=spec.id,
                kind=EvidenceKind.EXPERIMENT_RUN,
                provenance=_prov(commit, env, ph, exp_csv, pred_csv),
                result=Result(type="quantitative", point=mean_r2,
                              finding=json.dumps({**common, "statistic": "cv_r2_mean",
                                                  "threshold": _USEFUL_R2}, ensure_ascii=False)),
                bears_on=[Bearing(target_id=_hyp_useful(ph), direction=d_use)],
            )
            _save_evidence(ev_use, spec, workspace_dir); items.append(ev_use)

            # reliable (interval [p2.5, p97.5] vs null 0, support_side=above)
            if p2_5 > _NULL_R2:
                d_rel = BearingDirection.SUPPORTS
            elif p97_5 < _NULL_R2:
                d_rel = BearingDirection.REFUTES
            else:
                d_rel = BearingDirection.NEUTRAL
            ev_rel = EvidenceItem(
                id=_evidence_id(f"{ph.key}-reliable"),
                created_at=datetime.now(timezone.utc), spec_id=spec.id,
                kind=EvidenceKind.EXPERIMENT_RUN,
                provenance=_prov(commit, env, ph, exp_csv, pred_csv),
                result=Result(type="quantitative", point=mean_r2, ci=[p2_5, p97_5],
                              finding=json.dumps({**common, "statistic": "cv_r2_95pct_interval",
                                                  "null_value": _NULL_R2}, ensure_ascii=False)),
                bears_on=[Bearing(target_id=_hyp_reliable(ph), direction=d_rel)],
            )
            _save_evidence(ev_rel, spec, workspace_dir); items.append(ev_rel)
        return items

    return _run


# ---------------------------------------------------------------------------
# Helpers.
# ---------------------------------------------------------------------------
def _prov(commit, env, ph, exp_csv, pred_csv) -> Provenance:
    return Provenance(
        code_ref=f"{commit}:rigor/partition_model_capability.py",
        data_ref=(f"UFZ-LSER PPLFER (measured={ph.exp_col}); features=all logK_* "
                  f"({exp_csv.name} x {pred_csv.name})"),
        data_source="measured",
        environment=env,
    )


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


def partition_model_provider(exp_csv: Path, pred_csv: Path, repo_root: Path):
    """Build the adapter ``CapabilityProvider`` for the real `sci-adk run --capability` path.

    Mirrors sci-adk's own ``_t1_provider`` (adapter/registry.py): the kernel is driven
    through its designed seam -- the CLI (composition root) registers this provider, then
    ``sci-adk run --capability partition-model-cv`` resolves it from the registry and runs
    its built-in demo (Spec + options). The kernel never imports this module; the seam
    direction stays adapter -> kernel.
    """
    from sci_adk.adapter.registry import CapabilityProvider

    def experiment_fn(**options):
        return partition_model_experiment(
            options["exp_csv"], options["pred_csv"], options.get("repo_root")
        )

    return CapabilityProvider(
        id=PARTITION_MODEL_CAPABILITY_ID,
        experiment_fn=experiment_fn,
        demo_spec=build_partition_model_spec,
        demo_options=lambda: {"exp_csv": exp_csv, "pred_csv": pred_csv, "repo_root": repo_root},
    )


__all__ = [
    "PARTITION_MODEL_CAPABILITY_ID",
    "build_partition_model_spec",
    "partition_model_experiment",
    "partition_model_provider",
]
