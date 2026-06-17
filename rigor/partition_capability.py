#!/usr/bin/env python3
"""
sci-adk capability: UFZ partition-coefficient surrogate validation.

This is a domain plugin for the sci-adk rigor/verification kernel
(https://github.com/ccy5123/sci-adk). It turns the partition-toolkit's central
empirical question -- "do the cheap xTB implicit-solvent surrogates predict each
tissue-water partition coefficient well?" -- into a *frozen Spec* + *measured
Evidence* + an *autonomously judged Claim*, so the result is reproducible and
tamper-evident rather than a hand-maintained number in a README table.

Design contract (mirrors sci_adk.adapter.t1_capability):
  - Agents propose, the engine judges. This module only PRODUCES Evidence
    (one coefficient-of-determination R^2 per phase -> Result.point). The binding
    SUPPORTED/REFUTED verdict is the DecisionEngine's, applying the frozen
    numeric threshold rule autonomously (no LLM, no self-certification).
  - No hardcoded metric. The usefulness threshold (R^2 >= 0.70) lives only in
    the Spec's DecisionRule.params, never as a constant in the judging path.
  - Evidence validity. Each phase hypothesis is ``empirical`` (a partition
    coefficient is a physical phenomenon), and its Evidence is
    ``data_source="measured"`` because the R^2 is computed against the real
    measured UFZ-LSER experimental logK. The adequacy gate therefore admits the
    binding verdict honestly.

The statistic is the coefficient of determination of the *uncalibrated* surrogate
used directly as the predictor (R^2 = 1 - SS_res/SS_tot), which is the exact
metric the partition-toolkit README reports. The per-phase surrogate is the one
identified by the toolkit's prior 78-single-point screening; it is FROZEN here
(pre-registration), so this run is a confirmatory test of that pre-selected
method, not a fresh best-of-72 scan. The optimism of selecting on the same set is
recorded as a limitation in the Spec and the paper draft (resubstitution R^2, not
cross-validated).
"""

from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Tuple

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

# The capability id the runtime selector resolves to. HOW, not WHAT: it is never a
# frozen Spec field -- it travels only in Evidence provenance.
PARTITION_CAPABILITY_ID = "partition-xtb-ufz"
_DEFAULT_SPEC_ID = "partition-xtb-ufz"

# Frozen usefulness threshold for the coefficient of determination. A surrogate that
# explains >= 70% of the variance in the measured logK is treated as a useful direct
# predictor for that phase. Lives in the rule params (D1); copied here only to build
# those params + the human-readable expression, never read by the judging path.
_R2_THRESHOLD = 0.70


@dataclass(frozen=True)
class Phase:
    """One tissue/blood phase: its measured column + the pre-registered surrogate.

    Attributes:
        key: short phase key (also the hypothesis id stem).
        pretty: human-readable phase name for statements/paper.
        exp_col: experimental logK column in the UFZ input table (the MEASURED
            dependent variable).
        surrogate_col: the FROZEN xTB surrogate logK column used as the direct
            predictor (chosen by the toolkit's prior screening, README table).
    """

    key: str
    pretty: str
    exp_col: str
    surrogate_col: str

    @property
    def hyp_id(self) -> str:
        return f"hyp-{self.key}"

    @property
    def target_id(self) -> str:
        return f"tc-{self.key}"


# The five phases the partition-toolkit reports, each paired with the surrogate the
# README pre-identified. This list is the SINGLE source of truth: build_partition_spec
# derives one hypothesis per row, and the experiment evaluates exactly these rows.
PHASES: Tuple[Phase, ...] = (
    Phase("membrane", "phospholipid membrane-water", "membrane_exp", "logK_GFN2_gbsa_acetonitrile"),
    Phase("storage", "storage lipid-water", "storage_exp", "logK_GFN1_alpb_hexane"),
    Phase("albumin", "albumin-water", "albumin_exp", "logK_GFN1_alpb_woctanol"),
    Phase("muscle_chicken", "muscle protein-water (chicken)", "muscle_chicken_exp", "logK_GFN1_alpb_octanol"),
    Phase("muscle_fish", "muscle protein-water (fish)", "muscle_fish_exp", "logK_GFN1_alpb_octanol"),
)


# ---------------------------------------------------------------------------
# Spec builder (the frozen pre-registration contract).
# ---------------------------------------------------------------------------
def _phase_rule() -> DecisionRule:
    """The shared frozen threshold rule: coefficient of determination R^2 >= 0.70.

    The engine reads the statistic from ``Result.point`` (its threshold handler is
    hardwired to ``point``); ``statistic`` documents which measurement that is, and
    ``op``/``value`` are the machine-usable threshold (no hardcoded metric, D1).
    """
    return DecisionRule(
        kind=DecisionRuleKind.THRESHOLD,
        expression=(
            "coefficient of determination R^2 (surrogate logK as a direct, uncalibrated "
            "predictor of the measured phase logK) >= 0.70 over the UFZ test set => "
            "support (useful predictor for this phase); R^2 < 0.70 => refute"
        ),
        params={"statistic": "coef_det_R2", "op": ">=", "value": _R2_THRESHOLD, "combine": "latest"},
    )


def build_partition_spec(spec_id: str = _DEFAULT_SPEC_ID) -> Spec:
    """Build the frozen partition-validation Spec (one confirmatory hypothesis/phase).

    Each hypothesis is ``confirmatory`` (a pre-registered method + threshold) and
    ``empirical`` (a partition coefficient is a measured physical phenomenon, so the
    adequacy gate requires measured Evidence to affirm it). The surrogate per phase is
    frozen from the toolkit's prior screening; the rule is identical across phases so
    the only thing under test is whether each phase clears the usefulness threshold.
    """
    hypotheses = [
        Hypothesis(
            id=ph.hyp_id,
            statement=(
                f"The pre-registered single-conformer xTB surrogate "
                f"'{ph.surrogate_col}' predicts {ph.pretty} logK with coefficient of "
                f"determination R^2 >= {_R2_THRESHOLD:.2f} over the UFZ test set"
            ),
            mode=HypothesisMode.CONFIRMATORY,
            decision_rule=_phase_rule(),
            referent="empirical",
        )
        for ph in PHASES
    ]
    target_claims = [
        TargetClaim(
            id=ph.target_id,
            statement=(
                f"xTB implicit-solvent surrogacy is {'useful' if ph.key in ('membrane', 'storage') else 'evaluated'} "
                f"for {ph.pretty} partitioning"
            ),
            answers=ph.hyp_id,
        )
        for ph in PHASES
    ]
    return Spec(
        id=spec_id,
        version=1,
        raw_proposal=RawProposal(
            background=(
                "Tissue/blood partition coefficients drive bioaccumulation and PBPK "
                "models. The partition-toolkit screens cheap GFN1/GFN2-xTB single-point "
                "energies in 25 ALPB + 12-14 GBSA implicit solvents (78 SP/molecule) and "
                "uses each solvent's logK_solvent/water as a surrogate for a biological "
                "phase, validated against 440 UFZ-LSER compounds across five phases."
            ),
            goal=(
                "Establish, per phase, whether the pre-registered xTB surrogate is a "
                "useful direct predictor of the measured partition coefficient "
                "(coefficient of determination R^2 >= 0.70)."
            ),
            method=(
                "For each phase, pair the FROZEN pre-selected surrogate logK column with "
                "the measured UFZ logK, drop missing rows, and compute R^2 = 1 - "
                "SS_res/SS_tot using the surrogate directly as the prediction (the "
                "toolkit's reported metric). Judge each R^2 against the frozen 0.70 "
                "threshold autonomously. NOTE (limitation): the surrogate was selected on "
                "this same set, so these are resubstitution R^2 (optimistic upper bounds), "
                "not cross-validated."
            ),
            expected_output=(
                "A per-phase SUPPORTED/REFUTED verdict. Expected from prior screening: "
                "lipid phases (membrane, storage) SUPPORTED; protein phases (albumin, "
                "muscle chicken/fish) REFUTED -- i.e. xTB surrogacy works for lipids but "
                "not proteins."
            ),
        ),
        hypotheses=hypotheses,
        method=MethodPlan(
            approaches=[
                "GFN1/GFN2-xTB single-point energies with ALPB and GBSA implicit solvation",
                "solvent-as-phase surrogacy: logK_solvent/water from the body-temperature "
                "transfer free energy",
                "coefficient-of-determination validation against measured UFZ-LSER logK",
            ],
            tools=[
                ToolRef(name="tblite (GFN1-xTB/GFN2-xTB)", kind="solver"),
                ToolRef(name="RDKit ETKDGv3+MMFF94 (3D conformer)", kind="library"),
                ToolRef(name="UFZ-LSER PPLFER partition dataset (440 compounds)", kind="dataset"),
            ],
        ),
        target_claims=target_claims,
    )


# ---------------------------------------------------------------------------
# Experiment (Spec -> measured Evidence). The honest read; the engine still judges.
# ---------------------------------------------------------------------------
def _coefficient_of_determination(pred, exp) -> float:
    """R^2 = 1 - SS_res/SS_tot with the surrogate used DIRECTLY as the prediction.

    This is the toolkit's reported metric: it can be negative when the uncalibrated
    surrogate predicts worse than the mean of the measured values (the protein phases),
    which is exactly the signal that the surrogate is not a usable direct predictor.
    """
    import numpy as np

    pred = np.asarray(pred, dtype=float)
    exp = np.asarray(exp, dtype=float)
    ss_res = float(np.sum((exp - pred) ** 2))
    ss_tot = float(np.sum((exp - np.mean(exp)) ** 2))
    if ss_tot == 0.0:
        return float("nan")
    return 1.0 - ss_res / ss_tot


def _pearson_r2(pred, exp) -> float:
    """Squared Pearson correlation -- the variance a LINEAR CALIBRATION would explain.

    Recorded alongside the coefficient of determination for transparency: a phase can
    have a moderate Pearson r^2 (some linear signal) yet a negative coefficient of
    determination (the uncalibrated surrogate is still a bad direct predictor).
    """
    import numpy as np

    pred = np.asarray(pred, dtype=float)
    exp = np.asarray(exp, dtype=float)
    if np.std(pred) == 0 or np.std(exp) == 0:
        return float("nan")
    r = float(np.corrcoef(pred, exp)[0, 1])
    return r * r


def _rmse(pred, exp) -> float:
    import numpy as np

    pred = np.asarray(pred, dtype=float)
    exp = np.asarray(exp, dtype=float)
    return float(np.sqrt(np.mean((exp - pred) ** 2)))


def _load_paired(exp_csv: Path, pred_csv: Path, phase: Phase):
    """Return (measured, predicted) numpy arrays for one phase, NaN-dropped, joined on CAS."""
    import pandas as pd

    exp = pd.read_csv(exp_csv)
    pred = pd.read_csv(pred_csv)
    exp["CAS"] = exp["CAS"].astype(str).str.strip()
    pred["cas"] = pred["cas"].astype(str).str.strip()
    if phase.surrogate_col not in pred.columns:
        raise KeyError(f"surrogate column '{phase.surrogate_col}' not in {pred_csv}")
    if phase.exp_col not in exp.columns:
        raise KeyError(f"experimental column '{phase.exp_col}' not in {exp_csv}")
    merged = exp[["CAS", phase.exp_col]].merge(
        pred[["cas", phase.surrogate_col]], left_on="CAS", right_on="cas", how="inner"
    )
    sub = merged[[phase.exp_col, phase.surrogate_col]].apply(pd.to_numeric, errors="coerce").dropna()
    return sub[phase.exp_col].to_numpy(), sub[phase.surrogate_col].to_numpy()


def partition_experiment(
    exp_csv: Path,
    pred_csv: Path,
    repo_root: Optional[Path] = None,
) -> Callable[[Spec, Path], List[EvidenceItem]]:
    """Build the partition-validation ``ExperimentFn`` over the real UFZ data.

    Returns ``fn(spec, workspace_dir) -> [EvidenceItem]``: one measured EvidenceItem
    per phase carrying the coefficient of determination in ``Result.point`` (the
    statistic the frozen rule judges) plus the full per-phase stats as a JSON finding
    for the audit. The honest direction bearing (SUPPORTS iff R^2 >= threshold) is the
    experiment's own read; the binding verdict is still the engine's.

    Args:
        exp_csv: the UFZ input table with the measured ``*_exp`` columns.
        pred_csv: the xTB predictions table with the ``logK_*`` surrogate columns.
        repo_root: optional repo root for the provenance commit hash.
    """
    exp_csv = Path(exp_csv)
    pred_csv = Path(pred_csv)
    root = Path(repo_root) if repo_root else exp_csv.resolve().parent.parent

    def _run(spec: Spec, workspace_dir: Path) -> List[EvidenceItem]:
        commit = _git_commit(root)
        env = _environment()
        items: List[EvidenceItem] = []
        for phase in PHASES:
            measured, predicted = _load_paired(exp_csv, pred_csv, phase)
            n = int(len(measured))
            r2 = _coefficient_of_determination(predicted, measured)
            finding = {
                "phase": phase.key,
                "statistic": "coef_det_R2",
                "coef_det_R2": None if math.isnan(r2) else round(r2, 6),
                "pearson_r2": _round(_pearson_r2(predicted, measured)),
                "rmse_log_units": _round(_rmse(predicted, measured)),
                "n": n,
                "surrogate": phase.surrogate_col,
                "measured_column": phase.exp_col,
                "threshold": _R2_THRESHOLD,
                "metric_definition": "R2 = 1 - SS_res/SS_tot, surrogate used directly as prediction",
                "capability": PARTITION_CAPABILITY_ID,
            }
            # Honest read of the experiment's own statistic; engine remains authoritative.
            direction = (
                BearingDirection.SUPPORTS
                if (not math.isnan(r2) and r2 >= _R2_THRESHOLD)
                else BearingDirection.REFUTES
            )
            evidence = EvidenceItem(
                id=_evidence_id(phase.key),
                created_at=datetime.now(timezone.utc),
                spec_id=spec.id,
                kind=EvidenceKind.EXPERIMENT_RUN,
                provenance=Provenance(
                    code_ref=f"{commit}:rigor/partition_capability.py",
                    data_ref=(
                        f"UFZ-LSER PPLFER (measured={phase.exp_col}); "
                        f"xTB surrogate={phase.surrogate_col}; "
                        f"{exp_csv.name} x {pred_csv.name}"
                    ),
                    # The dependent variable is real measured UFZ experimental logK, so the
                    # R^2 is a statistic OF measured data -- the empirical adequacy gate is
                    # satisfied honestly (not a synthetic proxy).
                    data_source="measured",
                    environment=env,
                ),
                result=Result(
                    type="quantitative",
                    point=(None if math.isnan(r2) else float(r2)),
                    finding=json.dumps(finding, ensure_ascii=False),
                ),
                bears_on=[Bearing(target_id=phase.hyp_id, direction=direction)],
            )
            _save_evidence(evidence, spec, workspace_dir)
            items.append(evidence)
        return items

    return _run


# ---------------------------------------------------------------------------
# Helpers (local to the adapter; the kernel writes its own records).
# ---------------------------------------------------------------------------
def _round(x: float) -> Optional[float]:
    return None if (x is None or (isinstance(x, float) and math.isnan(x))) else round(float(x), 6)


def _evidence_id(phase_key: str) -> str:
    import uuid

    ts = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    return f"evi-partition-{phase_key}-{ts}-{uuid.uuid4().hex[:8]}"


def _git_commit(root: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _environment() -> str:
    try:
        import numpy
        import pandas

        return f"capability:{PARTITION_CAPABILITY_ID}, numpy:{numpy.__version__}, pandas:{pandas.__version__}"
    except Exception:
        return f"capability:{PARTITION_CAPABILITY_ID}"


def _save_evidence(evidence: EvidenceItem, spec: Spec, workspace_dir: Path) -> None:
    evidence_dir = Path(workspace_dir) / "runs" / spec.id / "evidence"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    (evidence_dir / f"{evidence.id}.json").write_text(
        json.dumps(evidence.model_dump(mode="json"), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


__all__ = [
    "PARTITION_CAPABILITY_ID",
    "Phase",
    "PHASES",
    "build_partition_spec",
    "partition_experiment",
]
