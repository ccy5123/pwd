#!/usr/bin/env python3
"""
Composition root for the first-principles COSMO-RS physics run.

Registers the cosmors-physics-muscle capability in sci-adk's adapter registry and drives
the REAL `sci-adk run --capability cosmors-physics-muscle` CLI. The engine renders the
verdicts; this script never judges and never calls an LLM.

    python rigor/run_physics.py
    sci-adk verify rigor/runs/cosmors-physics-muscle
"""

from __future__ import annotations

import sys
from pathlib import Path

RIGOR_DIR = Path(__file__).resolve().parent
REPO_ROOT = RIGOR_DIR.parent
sys.path.insert(0, str(RIGOR_DIR))

from physics_cosmors_capability import CAPABILITY_ID, cosmors_physics_provider  # noqa: E402

import sci_adk.cli as sci_adk_cli  # noqa: E402
from sci_adk.adapter import registry  # noqa: E402


def main() -> int:
    results_json = REPO_ROOT / "structural_protein" / "cosmo_work_v4" / "combined_results_46.json"
    source_py = REPO_ROOT / "structural_protein" / "cosmo_kpw_openrs_v4.py"
    for p in (results_json, source_py):
        if not p.exists():
            print(f"error: required input not found: {p}", file=sys.stderr)
            return 2

    if CAPABILITY_ID not in registry.available():
        registry.register(cosmors_physics_provider(results_json, source_py, REPO_ROOT))
    print(f"registered capabilities: {registry.available()}")
    return sci_adk_cli.main([
        "run", "--capability", CAPABILITY_ID, "-o", str(RIGOR_DIR), "--spec-id", CAPABILITY_ID,
    ])


if __name__ == "__main__":
    raise SystemExit(main())
