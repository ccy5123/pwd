#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scatter_all.py — 15-cell scatter grid (3 pose × 5 desolv)
predictions_all_cells.csv를 읽어서 전체 셀의 scatter plot 그리드 생성
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.metrics import r2_score

HERE = Path(__file__).resolve().parent
COL_ORDER = ['none', 'GFN1/gbsa', 'GFN1/alpb', 'GFN2/gbsa', 'GFN2/alpb']
POSE_ORDER = ['A', 'B1', 'B2']
NAVY, RED = '#004094', '#C00000'


def main():
    pred = pd.read_csv(HERE / "predictions_all_cells.csv")

    fig, axes = plt.subplots(3, 5, figsize=(20, 12))
    fig.suptitle('BSA Docking 15-Cell Factorial: Experimental vs Predicted log K_BSA/w',
                 fontsize=14, fontweight='bold')

    for i, pose in enumerate(POSE_ORDER):
        for j, desolv in enumerate(COL_ORDER):
            ax = axes[i, j]
            cell = pred[(pred.pose_method == pose) & (pred.desolv_method == desolv)]

            if len(cell) < 3:
                ax.text(0.5, 0.5, 'N/A', ha='center', va='center', transform=ax.transAxes)
                ax.set_xlim(0, 1); ax.set_ylim(0, 1)
                continue

            exp = cell.albumin_exp.values
            pr = cell.logK_pred.values

            r = pearsonr(exp, pr)[0]
            r2 = r2_score(exp, pr)

            ax.scatter(exp, pr, alpha=0.6, edgecolors='navy', facecolors=NAVY, s=30)

            # 1:1 line
            lims = [min(exp.min(), pr.min()) - 0.5, max(exp.max(), pr.max()) + 0.5]
            ax.plot(lims, lims, 'r--', lw=1, alpha=0.5)
            ax.set_xlim(lims); ax.set_ylim(lims)

            ax.set_xlabel(f'Experimental log K_BSA/w', fontsize=9)
            ax.set_ylabel(f'Predicted log K_BSA/w', fontsize=9)
            ax.set_title(f'{pose} × {desolv}\nR²={r2:+.3f}, r={r:+.3f}, n={len(exp)}',
                        fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outfile = HERE / "scatter_all_cells.png"
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"[saved] {outfile}")
    plt.close()


if __name__ == "__main__":
    main()
