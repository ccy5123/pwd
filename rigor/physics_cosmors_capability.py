#!/usr/bin/env python3
"""
sci-adk capability: first-principles COSMO-RS partition prediction (muscle protein).

This is the PHYSICS run (proposal_physics.md), not the ML surrogate run. The partition
coefficient is COMPUTED from physical law -- quantum sigma-profiles -> COSMO-RS
infinite-dilution activity coefficients -> log K = (ln_gamma_water - ln_gamma_phase)/ln10
-- with the protein represented as an explicit amino-acid pseudo-solvent (actin+myosin
composition) and ZERO parameters fitted to the partition data. The prediction comes from
the repo's openCOSMO-RS pipeline (structural_protein/cosmo_kpw_openrs_v4.py); this
capability VALIDATES that a-priori physical prediction against the measured Endo set and
diagnoses where the physics systematically fails.

Scope (honest): only the muscle-protein phase has a committed first-principles result in
the repo (cosmo_work_v4/combined_results_46.json). The other four phases (membrane,
storage, albumin) need sigma-profile generation (ORCA QM), which is not available in
this sandbox -- their scripts exist but were never run/committed.

The engine judges three pre-registered, field-standard (NOT data-tuned) hypotheses:
  - useful        : RMSE <= 1.0 log units (useful for screening)
  - unbiased      : |mean signed error| <= 0.3 (no systematic offset)
  - missing-physics: the over-prediction is localized to polar/H-bond solutes --
                     mean(polar residual) - mean(nonpolar residual) >= 0.5 -- which
                     diagnoses the homogeneous-pseudo-solvent limitation (exposed polar
                     residues over-counted vs a folded protein).
"""

from __future__ import annotations

import ast
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, List, Optional

from sci_adk.core.evidence import (
    Bearing, BearingDirection, EvidenceItem, EvidenceKind, Provenance, Result,
)
from sci_adk.core.spec import (
    DecisionRule, DecisionRuleKind, Hypothesis, HypothesisMode, MethodPlan,
    RawProposal, Spec, TargetClaim, ToolRef,
)

CAPABILITY_ID = "cosmors-physics-muscle"
_SPEC_ID = "cosmors-physics-muscle"

# Pre-committed, field-standard bars (copied into rule params; judging path reads them
# only from there). 1.0 = useful screening; 0.3 = unbiased; 0.5 = a material localized gap.
_USEFUL_RMSE = 1.0
_UNBIASED = 0.3
_POLAR_GAP = 0.5

_PRIOR_WORK = {
    "searched_via": "agent web_search + repo data/README",
    "key_prior_art": [
        "Endo, Bauerfeind, Goss 2012, Environ. Sci. Technol. 46, 12697-12703, "
        "doi:10.1021/es303379y -- muscle protein-water partitioning (the measured set).",
        "Endo, Escher, Goss 2011, Environ. Sci. Technol. 45, 5912-5921, "
        "doi:10.1021/es200855w -- shows a homogeneous PC pseudo-solvent (COSMO-RS "
        "'Approach 2') is RMSE ~1.0 because a bilayer/protein is not a homogeneous "
        "solvent; depth/structure-resolved treatment (COSMOmic) is needed.",
        "openCOSMO-RS (open-source COSMO-RS, openCOSMORS24a parameterization) -- the "
        "thermodynamic engine used; parameters are universal, NOT fit to partition data.",
    ],
    "what_is_first_principles_here": (
        "log K is computed from quantum sigma-profiles + COSMO-RS activity coefficients "
        "with the protein as an amino-acid pseudo-solvent; NO parameter is fit to the "
        "measured partition coefficients. This is a-priori physical prediction, not a "
        "descriptor regression."
    ),
}


def _rule(statistic: str, op: str, value: float) -> DecisionRule:
    expr = {
        "rmse": f"parameter-free COSMO-RS RMSE vs measured logK <= {value} log units "
                f"=> support (useful screening); > {value} => refute",
        "abs_bias": f"|mean signed error (pred-exp)| <= {value} => support (no systematic "
                    f"offset); > {value} => refute (biased)",
        "polar_gap": f"mean(polar residual) - mean(nonpolar residual) >= {value} => support "
                     f"(over-prediction localized to polar/H-bond solutes -- diagnoses the "
                     f"homogeneous-pseudo-solvent missing physics); < {value} => refute",
    }[statistic]
    return DecisionRule(kind=DecisionRuleKind.THRESHOLD, expression=expr,
                        params={"statistic": statistic, "op": op, "value": value, "combine": "latest"})


def build_physics_spec(spec_id: str = _SPEC_ID) -> Spec:
    H = [
        Hypothesis(id="hyp-muscle-useful", referent="empirical", mode=HypothesisMode.EXPLORATORY,
                   statement="Parameter-free first-principles COSMO-RS (amino-acid pseudo-solvent) "
                             "predicts muscle protein-water logK with RMSE <= 1.0 log units",
                   decision_rule=_rule("rmse", "<=", _USEFUL_RMSE)),
        Hypothesis(id="hyp-muscle-unbiased", referent="empirical", mode=HypothesisMode.EXPLORATORY,
                   statement="That physical prediction has no systematic offset: "
                             "|mean signed error| <= 0.3 log units",
                   decision_rule=_rule("abs_bias", "<=", _UNBIASED)),
        Hypothesis(id="hyp-muscle-missing-physics", referent="empirical", mode=HypothesisMode.EXPLORATORY,
                   statement="The prediction error is a localized systematic over-prediction of "
                             "polar/H-bond solutes (mean polar residual exceeds mean nonpolar "
                             "residual by >= 0.5), diagnosing the homogeneous amino-acid "
                             "pseudo-solvent as the missing physics (exposed polar residues "
                             "over-counted vs a folded protein)",
                   decision_rule=_rule("polar_gap", ">=", _POLAR_GAP)),
    ]
    TC = [TargetClaim(id=f"tc-{h.id}", statement=h.statement[:60], answers=h.id) for h in H]
    return Spec(
        id=spec_id, version=1,
        raw_proposal=RawProposal(
            background=(
                "Tissue partition coefficients are needed a priori for PBPK/bioaccumulation. "
                "PP-LFER needs measured Abraham descriptors + coefficients FITTED to partition "
                "data. A first-principles route computes log K purely from quantum sigma-profiles "
                "+ COSMO-RS activity coefficients, with the protein as an explicit amino-acid "
                "pseudo-solvent and NO parameter fit to partition data. This run validates that "
                "a-priori physical prediction for muscle protein (the one phase with a committed "
                "COSMO-RS result; the others need ORCA sigma-profiles, unavailable here)."
            ),
            goal=("Quantify how well parameter-free COSMO-RS reproduces measured muscle logK, and "
                  "diagnose WHERE the physics systematically fails (the missing physics)."),
            method=("Protein = actin+myosin amino-acid composition pseudo-solvent (capped Ace-X-Nme "
                    "residues); solute sigma-profiles from quantum calc; openCOSMORS24a; "
                    "ln_gamma_inf in water and in the protein mix at 310.15 K; "
                    "log K = (ln_gamma_w - ln_gamma_p)/ln10, NO fit to partition data. Validate vs "
                    "the Endo measured set (n=46): RMSE, mean signed error (bias), R^2, slope, and "
                    "per-class residuals (polar O/N vs nonpolar). Bars are field-standard and "
                    "pre-committed: useful RMSE<=1.0, unbiased |bias|<=0.3, localized-gap>=0.5."),
            expected_output=("A SUPPORTED/REFUTED verdict per hypothesis, plus the systematic error "
                             "pattern that reveals the missing physics. Honest prior expectation "
                             "(Endo 2011): a homogeneous pseudo-solvent is ~RMSE 1.0 and over-counts "
                             "polar interactions. EXPLORATORY (results pre-exist in the repo); a "
                             "confirmatory test needs an untouched solute set."),
        ),
        hypotheses=H,
        method=MethodPlan(
            approaches=[
                "first-principles COSMO-RS: quantum sigma-profile -> activity coefficient -> log K",
                "protein as amino-acid pseudo-solvent at tissue composition (no partition-data fit)",
                "validation vs measured logK: RMSE / bias / R^2 / slope + per-class residual diagnosis",
            ],
            tools=[
                ToolRef(name="openCOSMO-RS (openCOSMORS24a)", kind="solver"),
                ToolRef(name="amino-acid pseudo-solvent (actin+myosin composition)", kind="model"),
                ToolRef(name="Endo muscle protein-water measured logK (n=46)", kind="dataset"),
            ],
        ),
        target_claims=TC,
    )


def _classify_polar(smiles: str) -> str:
    """polar = carries an O or N (H-bond donor/acceptor); else nonpolar (C/H/halogen/S)."""
    return "polar" if any(ch in smiles for ch in "ONon") else "nonpolar"


def physics_experiment(results_json: Path, source_py: Path,
                       repo_root: Optional[Path] = None) -> Callable[[Spec, Path], List[EvidenceItem]]:
    results_json = Path(results_json); source_py = Path(source_py)
    root = Path(repo_root) if repo_root else results_json.resolve().parents[2]

    def _run(spec: Spec, workspace_dir: Path) -> List[EvidenceItem]:
        import numpy as np

        commit = _git_commit(root)
        env = "capability:" + CAPABILITY_ID + "; physics=openCOSMO-RS (no partition-data fit)"
        items: List[EvidenceItem] = []

        # prior work (literature) recorded
        items.append(_save(_lit(spec, commit, env), spec, workspace_dir))

        # load the first-principles COSMO-RS predictions + the solute SMILES
        res = json.loads(results_json.read_text())
        rows = {r["name"]: (float(r["exp"]), float(r["pred"])) for r in res["rows"]}
        m = re.search(r"ENDO46\s*=\s*(\[.*?\])\n\n", source_py.read_text(), re.S)
        smiles = {n: s for (n, _e, s) in ast.literal_eval(m.group(1))}

        names = list(rows.keys())
        exp = np.array([rows[n][0] for n in names]); pred = np.array([rows[n][1] for n in names])
        resid = pred - exp
        rmse = float(np.sqrt(np.mean(resid ** 2)))
        bias = float(np.mean(resid))
        r2 = float(np.corrcoef(pred, exp)[0, 1] ** 2)
        slope, intercept = (float(v) for v in np.polyfit(exp, pred, 1))

        cls = {n: _classify_polar(smiles.get(n, "")) for n in names}
        polar_res = [resid[i] for i, n in enumerate(names) if cls[n] == "polar"]
        nonpolar_res = [resid[i] for i, n in enumerate(names) if cls[n] == "nonpolar"]
        polar_mean = float(np.mean(polar_res)); nonpolar_mean = float(np.mean(nonpolar_res))
        polar_gap = polar_mean - nonpolar_mean

        per_compound = sorted(
            [{"name": n, "exp": round(exp[i], 2), "pred": round(pred[i], 2),
              "resid": round(float(resid[i]), 2), "class": cls[n], "smiles": smiles.get(n, "")}
             for i, n in enumerate(names)], key=lambda d: -d["resid"])

        # full validation record (every compound + per-class -> nothing hidden)
        base_finding = {
            "phase": "muscle_protein", "method": "first-principles openCOSMO-RS (no fit)",
            "n": len(exp), "rmse": round(rmse, 4), "bias_mean_signed_error": round(bias, 4),
            "r2": round(r2, 4), "slope": round(slope, 4), "intercept": round(intercept, 4),
            "polar_mean_resid": round(polar_mean, 4), "nonpolar_mean_resid": round(nonpolar_mean, 4),
            "polar_minus_nonpolar": round(polar_gap, 4),
        }
        items.append(_save(EvidenceItem(
            id=_eid("validation-record"), spec_id=spec.id, kind=EvidenceKind.OBSERVATION,
            provenance=_prov(commit, env), bears_on=[],
            result=Result(type="qualitative", finding=json.dumps(
                {**base_finding, "per_compound": per_compound}, ensure_ascii=False)),
        ), spec, workspace_dir))

        # three judged hypotheses (engine reads Result.point)
        for hyp_id, point, op, thr, statname in [
            ("hyp-muscle-useful", rmse, "<=", _USEFUL_RMSE, "rmse"),
            ("hyp-muscle-unbiased", abs(bias), "<=", _UNBIASED, "abs_bias"),
            ("hyp-muscle-missing-physics", polar_gap, ">=", _POLAR_GAP, "polar_gap"),
        ]:
            met = (point <= thr) if op == "<=" else (point >= thr)
            d = BearingDirection.SUPPORTS if met else BearingDirection.REFUTES
            items.append(_save(EvidenceItem(
                id=_eid(hyp_id), spec_id=spec.id, kind=EvidenceKind.EXPERIMENT_RUN,
                provenance=_prov(commit, env),
                result=Result(type="quantitative", point=float(point),
                              finding=json.dumps({**base_finding, "statistic": statname,
                                                  "op": op, "threshold": thr}, ensure_ascii=False)),
                bears_on=[Bearing(target_id=hyp_id, direction=d)],
            ), spec, workspace_dir))
        return items

    return _run


def cosmors_physics_provider(results_json: Path, source_py: Path, repo_root: Path):
    from sci_adk.adapter.registry import CapabilityProvider

    def experiment_fn(**o):
        return physics_experiment(o["results_json"], o["source_py"], o.get("repo_root"))

    return CapabilityProvider(
        id=CAPABILITY_ID, experiment_fn=experiment_fn, demo_spec=build_physics_spec,
        demo_options=lambda: {"results_json": results_json, "source_py": source_py, "repo_root": repo_root},
    )


# ---- helpers ----
def _lit(spec, commit, env) -> EvidenceItem:
    return EvidenceItem(
        id=_eid("priorwork"), spec_id=spec.id, kind=EvidenceKind.LITERATURE,
        provenance=Provenance(code_ref=f"{commit}:rigor/physics_cosmors_capability.py",
                              data_ref="prior art (Endo 2011/2012; openCOSMO-RS)", environment=env),
        result=Result(type="qualitative", finding=json.dumps(_PRIOR_WORK, ensure_ascii=False)),
        bears_on=[])


def _prov(commit, env) -> Provenance:
    return Provenance(
        code_ref=f"{commit}:structural_protein/cosmo_kpw_openrs_v4.py (openCOSMO-RS, no fit)",
        data_ref="cosmo_work_v4/combined_results_46.json (computed) vs Endo measured logK",
        data_source="measured", environment=env)


def _eid(tag: str) -> str:
    import uuid
    return f"evi-cosmors-{tag}-{datetime.now(timezone.utc):%Y%m%d-%H%M%S}-{uuid.uuid4().hex[:8]}"


def _git_commit(root: Path) -> str:
    import subprocess
    try:
        return subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=10).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _save(ev: EvidenceItem, spec: Spec, workspace_dir: Path) -> EvidenceItem:
    d = Path(workspace_dir) / "runs" / spec.id / "evidence"; d.mkdir(parents=True, exist_ok=True)
    (d / f"{ev.id}.json").write_text(json.dumps(ev.model_dump(mode="json"), indent=2, ensure_ascii=False))
    return ev


__all__ = ["CAPABILITY_ID", "build_physics_spec", "physics_experiment", "cosmors_physics_provider"]
