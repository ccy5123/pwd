#!/usr/bin/env python3
"""
Driver = the sci-adk composition root for the partition-model capability.

This drives sci-adk through its DESIGNED entry point (rigor-shell-architecture §3, F3):
it registers the partition-model capability in the adapter registry, then invokes the
real `sci-adk run --capability partition-model-cv` CLI -- the same path the kernel uses
for its own T-1 demo. The kernel resolves the ExperimentFn from the registry and renders
the verdicts itself; this script never judges and never calls an LLM (in-session-only).

So "is this actually run with sci-adk?" -> yes: the verdicts come out of
`sci_adk.cli.main(["run", "--capability", ...])`, not a hand-rolled compiler call.

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
    PARTITION_MODEL_CAPABILITY_ID,
    partition_model_provider,
)

import sci_adk.cli as sci_adk_cli  # noqa: E402
from sci_adk.adapter import registry  # noqa: E402


def main() -> int:
    exp_csv = REPO_ROOT / "data" / "compounds_input_table1.csv"
    pred_csv = REPO_ROOT / "results" / "xtb_gfn12_predictions.csv"
    for p in (exp_csv, pred_csv):
        if not p.exists():
            print(f"error: required data file not found: {p}", file=sys.stderr)
            return 2

    # Composition root: register the capability behind the adapter seam, then let the
    # real sci-adk CLI resolve + run it (idempotent if already registered).
    if PARTITION_MODEL_CAPABILITY_ID not in registry.available():
        registry.register(partition_model_provider(exp_csv, pred_csv, REPO_ROOT))

    print(f"registered capabilities: {registry.available()}")
    return sci_adk_cli.main([
        "run",
        "--capability", PARTITION_MODEL_CAPABILITY_ID,
        "-o", str(RIGOR_DIR),
        "--spec-id", PARTITION_MODEL_CAPABILITY_ID,
    ])


if __name__ == "__main__":
    raise SystemExit(main())
