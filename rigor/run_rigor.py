#!/usr/bin/env python3
"""
Driver: compile the partition-coefficient validation into a verifiable sci-adk run.

Runs the deterministic, zero-LLM path of the sci-adk research compiler over the real
partition-toolkit data:

    build_partition_spec()        -> frozen Spec (rigor/runs/<id>/spec.json)
    partition_experiment(...)     -> measured Evidence (rigor/runs/<id>/evidence/)
    ClaimUpdater + DecisionEngine -> autonomously judged Claims (rigor/runs/<id>/claims/)
    render_paper                  -> paper draft (rigor/runs/<id>/paper/draft.md)

Prereqs (see rigor/README.md):
    pip install "git+https://github.com/ccy5123/sci-adk@<commit>"   # sci_adk kernel
    pip install numpy pandas

Usage:
    python rigor/run_rigor.py
    sci-adk verify rigor/runs/partition-xtb-ufz     # third-party headless re-check
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make the capability importable whether run from the repo root or rigor/.
RIGOR_DIR = Path(__file__).resolve().parent
REPO_ROOT = RIGOR_DIR.parent
sys.path.insert(0, str(RIGOR_DIR))

from partition_capability import build_partition_spec, partition_experiment  # noqa: E402

from sci_adk.core.validity import ValidityHalt  # noqa: E402
from sci_adk.loop.compiler import ResearchCompiler  # noqa: E402


def main() -> int:
    exp_csv = REPO_ROOT / "data" / "compounds_input_table1.csv"
    pred_csv = REPO_ROOT / "results" / "xtb_gfn12_predictions.csv"
    for p in (exp_csv, pred_csv):
        if not p.exists():
            print(f"error: required data file not found: {p}", file=sys.stderr)
            return 2

    spec = build_partition_spec()
    experiment = partition_experiment(exp_csv, pred_csv, repo_root=REPO_ROOT)

    # Artifacts land under rigor/runs/<spec.id>/ (workspace_dir holds runs/).
    compiler = ResearchCompiler(workspace_dir=RIGOR_DIR)
    try:
        result = compiler.compile("", spec=spec, experiment=experiment)
    except ValidityHalt as e:
        print(f"error: evidence-validity halt: {e.reason}", file=sys.stderr)
        return 2

    print(f"compiled Spec '{result.spec.id}' -> {result.run_dir}")
    print(f"  evidence: {len(result.evidence)} | claims: {len(result.claims)}")
    for claim in sorted(result.claims, key=lambda c: c.answers):
        print(f"    - {claim.answers}: {claim.status.value.upper()}")
        print(f"        {claim.confidence.basis}")
    print(f"  paper draft: {result.paper_path}")
    print("\nRe-verify headlessly (no LLM, third-party reproducible):")
    print(f"    sci-adk verify {result.run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
