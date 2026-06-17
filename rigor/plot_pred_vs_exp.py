#!/usr/bin/env python3
"""
Pred vs exp scatter for the first-principles COSMO-RS muscle-protein prediction.

Visualizes the run's key finding: polar/H-bond solutes sit ABOVE the 1:1 line
(systematic over-prediction), nonpolar solutes track it. Saves a PNG under the run dir.

    python rigor/plot_pred_vs_exp.py
"""
from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RIGOR = Path(__file__).resolve().parent
ROOT = RIGOR.parent
RES = ROOT / "structural_protein" / "cosmo_work_v4" / "combined_results_46.json"
SRC = ROOT / "structural_protein" / "cosmo_kpw_openrs_v4.py"
OUT = RIGOR / "runs" / "cosmors-physics-muscle" / "figures" / "pred_vs_exp.png"


def classify(smiles: str) -> str:
    return "polar (O/N, H-bond)" if any(c in smiles for c in "ONon") else "nonpolar (C/H/halogen)"


def main() -> int:
    rows = {r["name"]: (float(r["exp"]), float(r["pred"])) for r in json.loads(RES.read_text())["rows"]}
    smiles = {n: s for (n, _e, s) in ast.literal_eval(
        re.search(r"ENDO46\s*=\s*(\[.*?\])\n\n", SRC.read_text(), re.S).group(1))}
    names = list(rows)
    exp = np.array([rows[n][0] for n in names]); pred = np.array([rows[n][1] for n in names])
    cls = np.array([classify(smiles.get(n, "")) for n in names])

    rmse = float(np.sqrt(np.mean((pred - exp) ** 2)))
    bias = float(np.mean(pred - exp))
    r2 = float(np.corrcoef(pred, exp)[0, 1] ** 2)
    slope, intercept = (float(v) for v in np.polyfit(exp, pred, 1))

    lo, hi = min(exp.min(), pred.min()) - 0.4, max(exp.max(), pred.max()) + 0.4
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    ax.plot([lo, hi], [lo, hi], "k--", lw=1.2, label="1:1 (perfect)", zorder=1)
    xs = np.array([lo, hi])
    ax.plot(xs, slope * xs + intercept, color="0.45", lw=1.4, ls="-",
            label=f"OLS fit (slope={slope:.2f})", zorder=1)

    colors = {"polar (O/N, H-bond)": "#d6432f", "nonpolar (C/H/halogen)": "#2f6fd6"}
    for c in colors:
        mask = cls == c
        ax.scatter(exp[mask], pred[mask], s=46, c=colors[c], edgecolor="white",
                   linewidth=0.6, alpha=0.9, label=f"{c}  (n={int(mask.sum())})", zorder=3)

    # annotate the worst over-predictions (largest pred-exp)
    for i in np.argsort(-(pred - exp))[:5]:
        ax.annotate(names[i], (exp[i], pred[i]), textcoords="offset points",
                    xytext=(6, 3), fontsize=7.5, color="0.25")

    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("experimental  log K$_{muscle/water}$")
    ax.set_ylabel("COSMO-RS predicted  log K$_{muscle/water}$")
    ax.set_title("First-principles COSMO-RS vs experiment (muscle protein)\n"
                 "parameter-free; no fit to partition data", fontsize=11)
    ax.text(0.03, 0.97,
            f"n = {len(exp)}\nRMSE = {rmse:.2f}\nbias = {bias:+.2f}\n"
            f"R$^2$ = {r2:.2f}\nslope = {slope:.2f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.7"))
    ax.text(0.30, 0.88, "polar solutes sit ABOVE 1:1\n(systematic over-prediction)",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.5, color="#d6432f")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.grid(True, ls=":", alpha=0.4)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print(f"saved {OUT}  (RMSE={rmse:.2f}, bias={bias:+.2f}, R2={r2:.2f}, slope={slope:.2f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
