#!/usr/bin/env python3
"""
누락된 3개 분자 계산 스크립트 v4
xTB 구조 최적화 + ORCA COSMO single-point
"""
import subprocess
import hashlib
import re
from pathlib import Path

# 설정
WORK = Path("./cosmo_work_v4")
CACHE = WORK / "sigma_cache"
CACHE.mkdir(parents=True, exist_ok=True)

XTB = "/home1/s9383/Muscle_protein/xtb-6.5.1/bin/xtb"
ORCA = "/home1/s9383/Muscle_protein/orca_6_1_1_linux_x86-64_shared_openmpi418/orca"

# 3개 분자
MOLECULES = [
    ('phenanthrene', 'c1ccc2ccc3ccccc3c2c1', 3.14),
    ('chrysene', 'c1ccc2ccc3c4ccccc4ccc3c2c1', 4.53),
    ('BaP', 'c1cc2ccc3cccc4ccc5cccc(c12)c3c45', 4.88),
]


def generate_xyz_with_rdkit(smiles, xyz_path):
    """RDKit으로 초기 xyz 구조 생성"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(m)
    Chem.MolToXYZFile(m, str(xyz_path))
    return xyz_path


def run_xtb_optimization(xyz_path, work_dir, charge=0):
    """xTB로 구조 최적화"""
    with open(work_dir / "xtb_opt.err", "w") as f_err, \
         open(work_dir / "xtb_opt.out", "w") as f_out:
        result = subprocess.run(
            [XTB, str(xyz_path.resolve()), "--opt", "tight", "--gfn", "2", "-c", str(charge)],
            cwd=work_dir,
            stdout=f_out,
            stderr=f_err
        )
    return work_dir / "xtbopt.xyz"


def run_orca_cosmo_sp(xyz_path, work_dir, charge=0):
    """ORCA COSMO single-point 계산 (BP86/def2-TZVPD)"""
    # xyz 파일에서 좌표 추출 (처음 2줄 제외)
    lines = xyz_path.read_text().split("\n")
    coords = "\n".join(l for l in lines[2:] if l.strip())

    # ORCA 입력 파일 생성 (v3 방식: 좌표 직접 입력)
    inp = work_dir / "cosmo_sp.inp"
    inp.write_text(
        f"! BP86 def2-TZVPD SP CPCM\n"
        f"%pal nprocs 1 end\n"
        f"* xyz {charge} 1\n{coords}\n*\n"
    )

    # ORCA 실행
    with open(work_dir / "cosmo_sp.out", "w") as f_out, \
         open(work_dir / "cosmo_sp.err", "w") as f_err:
        subprocess.run(
            [ORCA, str(inp.resolve())],
            cwd=work_dir,
            stdout=f_out,
            stderr=f_err
        )

    # cosmo_sp.cpcm 파일 확인
    cpcm_file = work_dir / "cosmo_sp.cpcm"
    if not cpcm_file.exists():
        raise FileNotFoundError(f"ORCA COSMO output not found: {cpcm_file}")
    return cpcm_file


def convert_to_orcacosmo(mol_dir, name):
    """ORCA COSMO 결과를 .orcacosmo 형식으로 변환"""
    cpcm_file = mol_dir / "cosmo_sp.cpcm"
    out_file = mol_dir / "cosmo_sp.out"
    xyz_file = mol_dir / "xtbopt.xyz"

    if not cpcm_file.exists() or not xyz_file.exists():
        return None

    output = mol_dir / f"{mol_dir.name}.orcacosmo"

    # 에너지 추출 (ORCA 출력에서)
    energy = None
    if out_file.exists():
        with open(out_file) as f:
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:
                    match = re.search(r"[-0-9.]+", line)
                    if match:
                        energy = match.group()
                    break

    with open(output, "w") as f:
        # 헤더
        f.write(f"{name} : BP86def2-TZVPD\n")

        # ENERGY 섹션
        if energy:
            f.write("\n" + "#" * 50 + "\n")
            f.write("#ENERGY\n")
            f.write(f"  {energy}\n")

        # XYZ_FILE 섹션
        f.write("\n" + "#" * 50 + "\n")
        f.write("#XYZ_FILE\n")
        with open(xyz_file) as xyz_f:
            f.write(xyz_f.read())

        # COSMO 섹션
        f.write("\n" + "#" * 50 + "\n")
        f.write("#COSMO\n")
        with open(cpcm_file) as cpcm_f:
            f.write(cpcm_f.read())

    return output


# 메인 루프
for name, smiles, exp in MOLECULES:
    key = hashlib.md5(smiles.encode()).hexdigest()[:16]
    mol_dir = CACHE / f"MOL_{key}"
    mol_dir.mkdir(exist_ok=True)

    print(f"\n=== {name} ===")
    print(f"MD5 key: {key}")
    print(f"디렉토리: {mol_dir}")

    # 1. 초기 xyz 생성
    init_xyz = mol_dir / "init.xyz"
    if not init_xyz.exists():
        print("  [1/4] RDKit으로 초기 구조 생성...")
        generate_xyz_with_rdkit(smiles, init_xyz)
    else:
        print("  [1/4] 초기 구조 있음 (skip)")

    # 2. xTB 최적화
    opt_xyz = mol_dir / "xtbopt.xyz"
    if not opt_xyz.exists():
        print("  [2/4] xTB 최적화...")
        try:
            run_xtb_optimization(init_xyz, mol_dir)
        except Exception as e:
            print(f"    [ERROR] xTB 실패: {e}")
            continue
    else:
        print("  [2/4] xTB 결과 있음 (skip)")

    # 3. ORCA COSMO single-point
    cpcm_file = mol_dir / "cosmo_sp.cpcm"
    if not cpcm_file.exists():
        print("  [3/4] ORCA COSMO single-point...")
        try:
            run_orca_cosmo_sp(opt_xyz, mol_dir)
        except Exception as e:
            print(f"    [ERROR] ORCA COSMO 실패: {e}")
            continue
    else:
        print("  [3/4] ORCA COSMO 결과 있음 (skip)")

    # 4. .orcacosmo 변환
    orcacosmo = mol_dir / f"{mol_dir.name}.orcacosmo"
    if not orcacosmo.exists():
        print("  [4/4] .orcacosmo 변환...")
        convert_to_orcacosmo(mol_dir, name)
    else:
        print("  [4/4] .orcacosmo 있음 (skip)")

    print(f"  [완료] {orcacosmo}")

print("\n=== 전체 완료 ===")
