#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
box_filter_check.py — FIX② 목적 확인용 (docking 재실행 후 실행).
  1) docking_poses.csv의 box(0..8)별 pose 수·평균 dG 보고 -> '뜬 박스'는 평균 dG가 0 근처(약결합)로 드러남.
  2) A/B1/B2를 (전 박스) vs (EXCLUDE_BOXES 제외)로 재계산해 calR²/absR² 비교.
※ 뜬 박스를 코드에 하드코딩하지 않음. 아래 EXCLUDE_BOXES를 (1)의 데이터 보고 채워서 재실행.
※ 재실행 전(box 전부 0)이면 (1)만 의미. A/B1/B2는 현행과 동일하게 나옴(정상).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.metrics import r2_score

HERE = Path(__file__).resolve().parent
RT = 0.0019872 * 310.15
LOG_M_BSA = np.log10(66.463)

# (1) 보고를 보고 채운다. 예: 평균 dG가 0 근처인 box가 7,8이면 -> [7, 8]
EXCLUDE_BOXES: list = []


def conv(dg):
    return -np.asarray(dg, float) / (2.302585 * RT) - LOG_M_BSA


def metrics(e, p):
    e = np.asarray(e, float); p = np.asarray(p, float)
    r = pearsonr(e, p)[0]
    return r2_score(e, p), r * r, r  # absR2, calR2(=r^2), pearson


def pose_methods(dg):
    dg = np.asarray(dg, float)
    A = conv(dg.min())
    w = np.exp(-dg / RT)
    B1 = conv(np.average(dg, weights=w))
    B2 = conv(dg.mean())
    return A, B1, B2


def main():
    dock = pd.read_csv(HERE / "docking_poses.csv")
    pred = pd.read_csv(HERE / "predictions_all_cells.csv")
    exp_map = (pred[(pred.pose_method == "A") & (pred.desolv_method == "none")]
               .set_index("compound").albumin_exp.to_dict())

    boxes = sorted(dock.box.unique())
    print("=" * 78)
    print("(1) box별 pose 수 · 평균 dG  (평균 dG가 0 근처 = 뜬 박스 후보)")
    print("=" * 78)
    g = dock.groupby("box").agg(n=("dG", "size"), mean_dG=("dG", "mean"),
                                min_dG=("dG", "min"), max_dG=("dG", "max"))
    print(g.to_string(float_format=lambda x: f"{x:.3f}"))
    if boxes == [0]:
        print("\n[주의] box가 전부 0 → docking 재실행 전 상태. 재실행 후 다시 실행할 것.")

    print("\n" + "=" * 78)
    print(f"(2) A/B1/B2 재계산: 전 박스 vs 제외 {EXCLUDE_BOXES or '(없음)'}")
    print("=" * 78)
    for label, excl in [("all boxes", []), (f"exclude {EXCLUDE_BOXES}", EXCLUDE_BOXES)]:
        if excl == [] and label.startswith("exclude"):
            print(f"\n[{label}] EXCLUDE_BOXES 비어있음 → 생략(값 채우고 재실행).")
            continue
        rows = {"A": [], "B1": [], "B2": [], "exp": []}
        for cmpd, sub in dock.groupby("compound"):
            if cmpd not in exp_map:
                continue
            sub2 = sub[~sub.box.isin(excl)] if excl else sub
            if len(sub2) == 0:
                continue
            A, B1, B2 = pose_methods(sub2.dG.values)
            rows["A"].append(A); rows["B1"].append(B1); rows["B2"].append(B2)
            rows["exp"].append(exp_map[cmpd])
        e = rows["exp"]
        print(f"\n[{label}]  n={len(e)}")
        for m in ["A", "B1", "B2"]:
            a, c, r = metrics(e, rows[m])
            print(f"  {m:<3} absR²={a:+8.3f}  calR²={c:.3f}  pearson={r:+.3f}")
    print("\n해석: exclude 후 B1/B2 calR²가 A에 근접/역전하면 뜬 박스가 B1/B2를 오염시켰다는 [Fact].")
    print("      변화 없거나 여전히 A 최고면 -> 뜬 박스 필터해도 A>B1>B2 결론 유지.")


if __name__ == "__main__":
    main()
