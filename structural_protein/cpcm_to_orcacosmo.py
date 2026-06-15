#!/usr/bin/env python3
"""
.cpcm 파일을 .orcacosmo 형식으로 변환하는 스크립트
"""
import os
import re
from pathlib import Path
# from rdkit import Chem
# from rdkit.Chem import AllChem

def convert_cpcm_to_orcacosmo(cpcm_path, out_path, xyz_path, name, method="BP86def2-TZVPD"):
    """
    .cpcm 파일을 .orcacosmo 형식으로 변환

    Parameters
    ----------
    cpcm_path : Path
        .cpcm 파일 경로
    out_path : Path
        ORCA 출력 파일 경로 (에너지 추출용)
    xyz_path : Path
        xyz 구조 파일 경로
    name : str
        분자 이름
    method : str
        계산 방법
    """
    output_path = f"{cpcm_path.parent}/{name}.orcacosmo"

    with open(output_path, "w") as f:
        # 헤더
        f.write(f"{name} : {method}\n")

        # ENERGY 섹션 (ORCA 출력에서 추출)
        if out_path and out_path.exists():
            energy = None
            with open(out_path, "r") as out_file:
                for line in out_file:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        match = re.search(r"[-0-9.]+", line)
                        if match:
                            energy = match.group()
                            break
            if energy:
                f.write("\n" + "#" * 50 + "\n")
                f.write("#ENERGY\n")
                f.write(f"  {energy}\n")

        # XYZ_FILE 섹션 (표준 xyz 형식)
        if xyz_path and xyz_path.exists():
            f.write("\n" + "#" * 50 + "\n")
            f.write("#XYZ_FILE\n")
            with open(xyz_path, "r") as xyz_file:
                lines = xyz_file.readlines()
                # xyz 형식: 첫 줄은 원자 수, 두 번째 줄은 comment
                # 좌표는 Å 단위 그대로 사용 (변환 없음)
                for i, line in enumerate(lines):
                    # 두 번째 줄(comment)도 그대로 출력
                    f.write(line)

        # COSMO 섹션 (.cpcm 내용)
        f.write("\n" + "#" * 50 + "\n")
        f.write("#COSMO\n")
        with open(cpcm_path, "r") as cpcm_file:
            for line in cpcm_file:
                f.write(line)

        # ADJACENCY_MATRIX는 현재 비활성화 (필요시 추가)
        # ADJACENCY_MATRIX는 openCOSMO-RS에서 선택 사항일 수 있음

    # with 블록이 종료된 후 반환
    return output_path

def batch_convert(work_dir="./cosmo_work_v3/sigma_cache"):
    """
    캐시된 모든 분자를 .orcacosmo로 변환
    """
    work_path = Path(work_dir)
    converted = []

    for mol_dir in work_path.glob("MOL_*"):
        cpcm_file = mol_dir / "cosmo_sp.cpcm"
        out_file = mol_dir / "cosmo_sp.out"
        xyz_file = mol_dir / "xtbopt.xyz"  # xtb 결과 사용

        if not xyz_file.exists():
            xyz_file = mol_dir / "geom_opt.xyz"  # DFT 결과 사용

        if cpcm_file.exists() and xyz_file.exists():
            name = mol_dir.name
            try:
                result = convert_cpcm_to_orcacosmo(cpcm_file, out_file, xyz_file, name)
                converted.append(result)
                print(f"[OK] {name}")
            except Exception as e:
                print(f"[FAIL] {name}: {e}")

    return converted

def batch_convert_aa(work_dir="./cosmo_work_v3/sigma_cache"):
    """
    아미노산 캐시를 .orcacosmo로 변환
    """
    work_path = Path(work_dir)
    converted = []

    for aa_dir in work_path.glob("AA_*"):
        cpcm_file = aa_dir / "cosmo_sp.cpcm"
        out_file = aa_dir / "cosmo_sp.out"
        xyz_file = aa_dir / "xtbopt.xyz"

        if cpcm_file.exists() and xyz_file.exists():
            name = aa_dir.name
            try:
                result = convert_cpcm_to_orcacosmo(cpcm_file, out_file, xyz_file, name)
                converted.append(result)
                print(f"[OK] {name}")
            except Exception as e:
                print(f"[FAIL] {name}: {e}")

    return converted

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
    else:
        work_dir = "./cosmo_work_v3/sigma_cache"

    print("=== Converting molecules ===")
    batch_convert(work_dir)

    print("\n=== Converting amino acids ===")
    batch_convert_aa(work_dir)

    print("\n=== Conversion complete ===")