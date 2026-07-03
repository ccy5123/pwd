#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
apply_fixes_23.py — bsa_docking_desolv_pipeline_v2.py 에 FIX②③ 실제 적용.
  ② box 식별자 보존: dock_compound가 (dG, coords, box_idx) 반환 -> dedup 통과 -> 저장 시 실제 box(0..8).
     (기존 L813 'box': pose_idx//len(=항상 0) 폐기)
  ③ RMSD 공식 교정: np.mean(diff**2)(3N으로 나눠 √3배 과소평가->과잉병합) -> np.mean(np.sum(diff**2,axis=1))
     = per-atom heavy-atom RMSD (2.0Å 컷오프 의도와 일치). ※기존 논리(heavy-atom RMSD<2.0Å) 유지.
  ④·① 은 건드리지 않음(기존 논리 유지). 이 스크립트 실행 후 docking 재실행 필요.
anchor 유일성 assert + [FIX2] 마커 멱등 가드 + ast.parse 게이트.
실행: (작업폴더) python apply_fixes_23.py
"""
import ast, sys
from pathlib import Path

F = Path("bsa_docking_desolv_pipeline_v2.py")
if not F.exists():
    print(f"[STOP] {F} 없음. 작업폴더 확인."); sys.exit(1)
src = F.read_text()

if "all_boxes = [b for _, _, b in dock_results]" in src or "box_idx" in src:
    print("[SKIP] '[FIX2]' 마커 존재 → 이미 적용됨. 아무것도 안 함."); sys.exit(0)


def repl(s, old, new, tag):
    n = s.count(old)
    assert n == 1, f"[STOP] [{tag}] anchor count={n} (기대 1). 파일이 이전과 다름 → 중단."
    return s.replace(old, new)


# ============ ③ + box 보존: deduplicate_poses 전체 교체 ============
old_dedup = '''def deduplicate_poses(coords: List[np.ndarray], energies: List[float]) -> Tuple[List[np.ndarray], List[float]]:
    """Remove duplicate poses based on heavy-atom RMSD < 2.0 Å."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    n = len(coords)
    keep = [True] * n
    unique_coords = []
    unique_energies = []

    for i in range(n):
        if not keep[i]:
            continue
        unique_coords.append(coords[i])
        unique_energies.append(energies[i])

        for j in range(i+1, n):
            if keep[j]:
                # Simple RMSD calculation
                diff = coords[i] - coords[j]
                rmsd = np.sqrt(np.mean(diff**2))
                if rmsd < RMSD_CUTOFF:
                    keep[j] = False  # Keep only the lowest energy (i comes first)

    return unique_coords, unique_energies'''

new_dedup = '''def deduplicate_poses(coords: List[np.ndarray], energies: List[float], boxes: List[int]) -> Tuple[List[np.ndarray], List[float], List[int]]:
    """Remove duplicate poses by heavy-atom RMSD < 2.0 A. Preserves box id. [FIX2/3]"""
    n = len(coords)
    keep = [True] * n
    unique_coords = []
    unique_energies = []
    unique_boxes = []

    for i in range(n):
        if not keep[i]:
            continue
        unique_coords.append(coords[i])
        unique_energies.append(energies[i])
        unique_boxes.append(boxes[i])  # [FIX2] box 식별자 보존

        for j in range(i+1, n):
            if keep[j]:
                diff = coords[i] - coords[j]
                # [FIX3] heavy-atom RMSD = sqrt(mean over atoms of per-atom squared distance).
                #        기존 np.mean(diff**2)는 3N으로 나눠 RMSD가 √3배 과소평가 -> 과잉 병합.
                #        sum(axis=1)로 per-atom 거리제곱을 먼저 합산해 교정(2.0A 컷오프 의도와 일치).
                rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
                if rmsd < RMSD_CUTOFF:
                    keep[j] = False  # Keep only the lowest energy (i comes first)

    return unique_coords, unique_energies, unique_boxes'''
src = repl(src, old_dedup, new_dedup, "dedup FIX3+box")

# ============ ② dock_compound: box 태그 ============
src = repl(src,
    '''                  single_box: bool = False) -> List[Tuple[float, np.ndarray]]:
    """Dock compound against all boxes. Returns list of (ΔG, coords)."""''',
    '''                  single_box: bool = False) -> List[Tuple[float, np.ndarray, int]]:
    """Dock compound against all boxes. Returns list of (ΔG, coords, box_idx). [FIX2]"""''',
    "dock_compound sig")

src = repl(src,
    "    for cx, cy, cz in centers_to_use:",
    "    for box_idx, (cx, cy, cz) in enumerate(centers_to_use):  # [FIX2] box 식별자",
    "dock_compound loop")

src = repl(src,
    "            results.append((float(dG), coords))",
    "            results.append((float(dG), coords, box_idx))  # [FIX2]",
    "dock_compound append")

# ============ ② smoke test 언패킹 ============
src = repl(src,
    "        best_dG = min(e for e, _ in dock_results)",
    "        best_dG = min(e for e, _, _ in dock_results)  # [FIX2] 3-tuple",
    "smoke unpack")

# ============ ② run_full_docking 언패킹 + 저장 ============
old_save = '''            all_coords = [coords for _, coords in dock_results]
            all_energies = [dG for dG, _ in dock_results]

            # Check atom count consistency
            n_atoms_list = [len(c) for c in all_coords]
            if len(set(n_atoms_list)) > 1:
                LOG.warning(f"Variable atom count in poses: {set(n_atoms_list)}, skipping dedup")
                # Still save, but without dedup
                unique_coords, unique_energies = all_coords, all_energies
            else:
                unique_coords, unique_energies = deduplicate_poses(all_coords, all_energies)

            LOG.info(f"  After RMSD deduplication: {len(unique_coords)} poses")

            # (c) Append deduplicated poses to CSV
            rows_to_append = []
            for pose_idx, (dG, coords) in enumerate(zip(unique_energies, unique_coords)):
                rows_to_append.append({
                    'compound': compound,
                    'SMILES': smiles,
                    'box': pose_idx // len(unique_coords),  # Re-index after dedup
                    'pose': pose_idx % len(unique_coords),'''

new_save = '''            all_coords = [c for _, c, _ in dock_results]        # [FIX2] 3-tuple
            all_energies = [dG for dG, _, _ in dock_results]    # [FIX2]
            all_boxes = [b for _, _, b in dock_results]         # [FIX2] 실제 box 식별자

            # Check atom count consistency
            n_atoms_list = [len(c) for c in all_coords]
            if len(set(n_atoms_list)) > 1:
                LOG.warning(f"Variable atom count in poses: {set(n_atoms_list)}, skipping dedup")
                # Still save, but without dedup
                unique_coords, unique_energies, unique_boxes = all_coords, all_energies, all_boxes
            else:
                unique_coords, unique_energies, unique_boxes = deduplicate_poses(all_coords, all_energies, all_boxes)

            LOG.info(f"  After RMSD deduplication: {len(unique_coords)} poses")

            # (c) Append deduplicated poses to CSV
            rows_to_append = []
            for pose_idx, (dG, coords, box_id) in enumerate(zip(unique_energies, unique_coords, unique_boxes)):
                rows_to_append.append({
                    'compound': compound,
                    'SMILES': smiles,
                    'box': box_id,              # [FIX2] 실제 box(0..8), 기존 pose_idx//len(=0) 폐기
                    'pose': pose_idx,'''
src = repl(src, old_save, new_save, "full_docking save")

ast.parse(src)
F.write_text(src)
print("[OK] FIX2/3 적용 완료 + ast.parse 통과.")
print("     다음: docking 재실행 필요(box·dedup 반영). 재실행 없이는 docking_poses.csv 안 바뀜.")
