#!/usr/bin/env python3
"""
Driver: compile the cross-validated multi-surrogate MODEL study into a sci-adk run.

This is the research-advance run (see partition_model_capability.py): it fits a
regularized multi-surrogate ridge model over the 72 xTB logK columns, evaluates it by
repeated cross-validation, and lets the sci-adk DecisionEngine judge each phase's
`useful` (mean CV R^2 >= 0.60) and `robust` (min CV R^2 >= 0.00) hypotheses.

Usage:
    python rigor/run_model.py
    sci-adk verify rigor/runs/partition-model-cv
"""

from __future__ import annotations

import sys
from pathlib import Path

RIGOR_DIR = Path(__file__).resolve().parent
REPO_ROOT = RIGOR_DIR.parent
sys.path.insert(0, str(RIGOR_DIR))

from partition_model_capability import (  # noqa: E402
    build_partition_model_spec,
    partition_model_experiment,
)

from sci_adk.core.validity import ValidityHalt  # noqa: E402
from sci_adk.loop.compiler import ResearchCompiler  # noqa: E402


def main() -> int:
    exp_csv = REPO_ROOT / "data" / "compounds_input_table1.csv"
    pred_csv = REPO_ROOT / "results" / "xtb_gfn12_predictions.csv"
    for p in (exp_csv, pred_csv):
        if not p.exists():
            print(f"error: required data file not found: {p}", file=sys.stderr)
            return 2

    spec = build_partition_model_spec()
    experiment = partition_model_experiment(exp_csv, pred_csv, repo_root=REPO_ROOT)

    compiler = ResearchCompiler(workspace_dir=RIGOR_DIR)
    try:
        result = compiler.compile("", spec=spec, experiment=experiment)
    except ValidityHalt as e:
        print(f"error: evidence-validity halt: {e.reason}", file=sys.stderr)
        return 2

    print(f"compiled Spec '{result.spec.id}' -> {result.run_dir}")
    print(f"  evidence: {len(result.evidence)} | claims: {len(result.claims)}\n")
    for claim in sorted(result.claims, key=lambda c: c.answers):
        print(f"    - {claim.answers:26s} {claim.status.value.upper()}")
    print(f"\n  paper draft: {result.paper_path}")
    print("\nRe-verify headlessly:")
    print(f"    sci-adk verify {result.run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
