#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reinterpret.py  —  BSA docking + xTB desolvation factorial 재해석 (docking 재실행 없음)
------------------------------------------------------------------------------------
읽기: predictions_all_cells.csv, results_factorial.csv, gsolv.csv, docking_poses.csv
출력: (1) 15셀 재해석 테이블  (2) 버그별 기여 분해  (3) scatter/heatmap PNG

핵심 원리 [Fact]
  logK_pred = -dG_corr/(2.302585*RT) - LOG_M_BSA
  -> offset(-1.823)·Gsolv 가산·pose 평균은 모두 pred를 '선형'으로만 바꿈.
  -> 따라서 pearson·calR²(=pearson², 선형보정 후 R² 상한)는 이 버그들에 '불변'.
  -> absR²만 왜곡됨. 즉 calR²가 각 셀의 '버그-무관 진짜 성능'.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.metrics import r2_score

HERE = Path(__file__).resolve().parent
RT = 0.0019872 * 310.15          # 0.6163 kcal/mol
LOG_M_BSA = np.log10(66.463)     # 1.8226
COL_ORDER = ['none', 'GFN1/gbsa', 'GFN1/alpb', 'GFN2/gbsa', 'GFN2/alpb']
POSE_ORDER = ['A', 'B1', 'B2']
NAVY, RED = '#004094', '#C00000'


def cell_metrics(e, p):
    """스케일-불변(calR², pearson) + 스케일-의존(absR², best-offset absR², bias, rmse)."""
    e = np.asarray(e, float); p = np.asarray(p, float)
    r = pearsonr(e, p)[0]
    bias = float(np.mean(e - p))                 # pred에 더할 최적 상수
    return dict(
        n=len(e),
        absR2=r2_score(e, p),                    # raw (offset -1.823 고정)
        absR2_off=r2_score(e, p + bias),         # offset만 자유(slope=1) — ① 상한
        calR2=r * r,                             # 완전 선형보정 R² — 버그-무관 상한
        pearson=r,
        bias=bias,
        rmse=float(np.sqrt(np.mean((e - p) ** 2))),
    )


def load():
    pred = pd.read_csv(HERE / "predictions_all_cells.csv")
    dock = pd.read_csv(HERE / "docking_poses.csv")
    return pred, dock


def table_15cell(pred):
    rows = []
    for pm in POSE_ORDER:
        for dm in COL_ORDER:
            c = pred[(pred.pose_method == pm) & (pred.desolv_method == dm)]
            if len(c) < 3:
                continue
            m = cell_metrics(c.albumin_exp.values, c.logK_pred.values)
            m.update(pose_method=pm, desolv_method=dm)
            rows.append(m)
    return pd.DataFrame(rows)


def report(df):
    print("=" * 92)
    print("15-CELL 재해석  (calR² = 버그-무관 진짜 성능 지표; absR²는 스케일 왜곡 포함)")
    print("=" * 92)
    for metric, title in [('calR2', 'calR²  (scale-invariant, TRUE signal)'),
                          ('absR2', 'absR²  (raw, offset+Gsolv 스케일 왜곡 포함)'),
                          ('absR2_off', 'absR²  (offset 자유화 후 — ① 기여 제거)')]:
        piv = df.pivot(index='pose_method', columns='desolv_method', values=metric)
        piv = piv.reindex(index=POSE_ORDER, columns=COL_ORDER)
        print(f"\n[{title}]")
        print(piv.to_string(float_format=lambda x: f"{x:+.3f}"))

    best = df.loc[df.calR2.idxmax()]
    print("\n" + "-" * 92)
    print(f"최고 calR² 셀 = {best.pose_method}×{best.desolv_method}: "
          f"calR²={best.calR2:.3f}  pearson={best.pearson:+.3f}  "
          f"absR²(raw)={best.absR2:+.3f}  absR²(offset자유)={best.absR2_off:+.3f}  "
          f"bias={best.bias:+.3f}  n={int(best.n)}")

    none_cal = df[(df.desolv_method == 'none')].set_index('pose_method').calR2
    print("\n[검정] desolvation이 calR²를 올리는 셀이 있는가?")
    beat = df[(df.desolv_method != 'none')]
    beat = beat[beat.apply(lambda r: r.calR2 > none_cal.get(r.pose_method, -9), axis=1)]
    if beat.empty:
        print("  없음. 모든 desolv 셀 calR² < 같은 pose의 none.  => desolvation 무익(버그-무관).")
    else:
        print("  있음:"); print(beat[['pose_method', 'desolv_method', 'calR2']].to_string(index=False))
    print("[검정] pose 순위(calR²): " +
          " > ".join(f"{p}({df[(df.pose_method==p)&(df.desolv_method=='none')].calR2.values[0]:.3f})"
                     for p in POSE_ORDER))
    print("-" * 92)


def figures(pred, df):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        rcParams['font.family'] = 'Arial'
        rcParams['axes.unicode_minus'] = False
    except Exception as ex:
        print(f"[figures skipped] matplotlib/font 문제: {ex}")
        return

    # (1) A×none scatter
    try:
        c = pred[(pred.pose_method == 'A') & (pred.desolv_method == 'none')]
        e, p = c.albumin_exp.values, c.logK_pred.values
        m = cell_metrics(e, p)
        fig, ax = plt.subplots(figsize=(5.2, 5.2))
        ax.scatter(e, p, s=28, c=NAVY, alpha=.75, edgecolors='white', linewidths=.5)
        lo, hi = min(e.min(), p.min()) - .5, max(e.max(), p.max()) + .5
        ax.plot([lo, hi], [lo, hi], '--', c=RED, lw=1.4, label='y = x')
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel('Experimental log $K_{BSA/w}$', fontsize=12)
        ax.set_ylabel('Predicted log $K_{BSA/w}$', fontsize=12)
        ax.set_title(f'A × none  (n={m["n"]})', fontsize=13)
        ax.text(.04, .96, f"pearson r = {m['pearson']:.2f}\n"
                          f"calR$^2$ = {m['calR2']:.2f}\n"
                          f"absR$^2$ = {m['absR2']:.2f}\n"
                          f"bias = {m['bias']:+.2f}",
                transform=ax.transAxes, va='top', ha='left', fontsize=10,
                bbox=dict(boxstyle='round', fc='white', ec=NAVY, alpha=.9))
        ax.legend(loc='lower right', fontsize=10)
        fig.tight_layout(); fig.savefig(HERE / "scatter_A_none.png", dpi=150); plt.close(fig)
        print("[saved] scatter_A_none.png")
    except Exception as ex:
        print(f"[scatter skipped] {ex}")

    # (2) calR² & absR² heatmaps
    for metric, fname, vlim in [('calR2', 'calR2_heatmap.png', (0, .3)),
                               ('absR2', 'absR2_heatmap.png', (-1, .3))]:
        try:
            piv = df.pivot(index='pose_method', columns='desolv_method', values=metric)
            piv = piv.reindex(index=POSE_ORDER, columns=COL_ORDER)
            fig, ax = plt.subplots(figsize=(9, 3.4))
            from matplotlib.colors import LinearSegmentedColormap
            cmap = LinearSegmentedColormap.from_list('nr', [RED, '#FFFFFF', NAVY])
            im = ax.imshow(piv.values, cmap=cmap, vmin=vlim[0], vmax=vlim[1], aspect='auto')
            ax.set_xticks(range(len(COL_ORDER))); ax.set_xticklabels(COL_ORDER, fontsize=10)
            ax.set_yticks(range(len(POSE_ORDER))); ax.set_yticklabels(POSE_ORDER, fontsize=11)
            for i in range(piv.shape[0]):
                for j in range(piv.shape[1]):
                    v = piv.values[i, j]
                    ax.text(j, i, 'N/A' if np.isnan(v) else f'{v:.2f}',
                            ha='center', va='center', fontsize=10,
                            color='black' if vlim[0] <= (v if not np.isnan(v) else 0) else 'white')
            ax.set_xlabel('Desolvation method', fontsize=11)
            ax.set_ylabel('Pose method', fontsize=11)
            ax.set_title(f'{metric}  (3 pose × 5 desolvation)', fontsize=12)
            fig.colorbar(im, ax=ax, fraction=.025)
            fig.tight_layout(); fig.savefig(HERE / fname, dpi=150); plt.close(fig)
            print(f"[saved] {fname}")
        except Exception as ex:
            print(f"[{fname} skipped] {ex}")


def main():
    pred, dock = load()
    # 정합성: box 좌표 상태 (② 실행가능성 판정)
    boxes = sorted(dock.box.unique())
    if boxes == [0]:
        print("[Fact] docking_poses.csv의 box 열이 전부 0 → 뜬 박스 필터(FIX②)·좌표 RMSD(FIX③) 불가")
        print("       (저장 시 box·coords 유실). ②③는 docking 재실행 필요.\n")
    df = table_15cell(pred)
    report(df)
    df.to_csv(HERE / "results_reinterpreted.csv", index=False)
    print(f"\n[saved] results_reinterpreted.csv")
    figures(pred, df)


if __name__ == "__main__":
    sys.exit(main())
