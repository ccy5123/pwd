#!/usr/bin/env python3
"""
benzo[a]pyrene 단일 분자 재생성 스크립트 (canonical-md5 통일)
- 기존 calc_missing_molecules_v4.py 파이프라인 재사용
- 올바른 SMILES + canonical-md5 + charge=0
"""
import subprocess
import re
from pathlib import Path

# ----------------------------------------------------------------------
# 설정 (절대경로)
# ----------------------------------------------------------------------
WORK = Path("./cosmo_work_v3")
CACHE = WORK / "sigma_cache"

XTB = "/home1/s9383/Muscle_protein/xtb-6.5.1/bin/xtb"
ORCA = "/home1/s9383/Muscle_protein/orca_6_1_1_linux_x86-64_shared_openmpi418/orca"

# ----------------------------------------------------------------------
# 재생성 대상 (단일, 확정값)
# ----------------------------------------------------------------------
BAP_NAME = "benzo[a]pyrene"
BAP_SMILES = "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34"  # 올바른 canonical SMILES
BAP_CHARGE = 0
BAP_EXP = 4.88  # chicken exp (참고용)

# 검증용 기대값
EXPECTED_INCHIKEY = "FMMWHPNWAFZXNH"  # 첫 블록
EXPECTED_DIR_HASH = "e2a7101ce1a41b29"  # canonical-md5

# ----------------------------------------------------------------------
# 캐시 키 (canonical-md5 통일)
# ----------------------------------------------------------------------
from cosmo_kpw_common_v3 import cache_key

# 사전검증: cache_key가 정말 e2a7101ce1a41b29 인지, SMILES가 진짜 BaP인지
print("=== 사전검증 ===")
computed_key = cache_key(BAP_SMILES)
print(f"cache_key  = {computed_key} (기대: {EXPECTED_DIR_HASH})")
print(f"일치       = {computed_key == EXPECTED_DIR_HASH}")

from rdkit import Chem
computed_inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(BAP_SMILES))
print(f"InChIKey   = {computed_inchikey} (기대: {EXPECTED_INCHIKEY}...)")
print(f"일치       = {computed_inchikey.split('-')[0] == EXPECTED_INCHIKEY}")

if computed_key != EXPECTED_DIR_HASH:
    raise ValueError(f"cache_key 불일치: {computed_key} != {EXPECTED_DIR_HASH}")
if computed_inchikey.split('-')[0] != EXPECTED_INCHIKEY:
    raise ValueError(f"InChIKey 불일치: {computed_inchikey} != {EXPECTED_INCHIKEY}...")

print("=== 사전검증 통과 ===\n")

# ----------------------------------------------------------------------
# 파이프라인 함수들 (calc_missing_molecules_v4.py 재사용)
# ----------------------------------------------------------------------
def generate_xyz_with_rdkit(smiles, xyz_path):
    """RDKit ETKDGv3 + MMFF로 초기 xyz 구조 생성"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(m)
    Chem.MolToXYZFile(m, str(xyz_path))
    return xyz_path


def run_xtb_optimization(xyz_path, work_dir, charge=0):
    """xTB로 구조 최적화 (--opt tight --gfn 2)"""
    with open(work_dir / "xtb_opt.err", "w") as f_err, \
         open(work_dir / "xtb_opt.out", "w") as f_out:
        result = subprocess.run(
            [XTB, str(xyz_path.resolve()), "--opt", "tight", "--gfn", "2", "-c", str(charge)],
            cwd=work_dir,
            stdout=f_out,
            stderr=f_err
        )
    if result.returncode != 0:
        raise RuntimeError(f"xTB optimization failed (code={result.returncode})")
    return work_dir / "xtbopt.xyz"


def run_orca_cosmo_sp(xyz_path, work_dir, charge=0):
    """ORCA COSMO single-point 계산 (BP86/def2-TZVPD, CPCM)"""
    lines = xyz_path.read_text().split("\n")
    coords = "\n".join(l for l in lines[2:] if l.strip())

    inp = work_dir / "cosmo_sp.inp"
    inp.write_text(
        f"! BP86 def2-TZVPD SP CPCM\n"
        f"%pal nprocs 1 end\n"
        f"* xyz {charge} 1\n{coords}\n*\n"
    )

    with open(work_dir / "cosmo_sp.out", "w") as f_out, \
         open(work_dir / "cosmo_sp.err", "w") as f_err:
        result = subprocess.run(
            [ORCA, str(inp.resolve())],
            cwd=work_dir,
            stdout=f_out,
            stderr=f_err
        )

    if result.returncode != 0:
        raise RuntimeError(f"ORCA COSMO failed (code={result.returncode})")

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

    # 에너지 추출
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


# ----------------------------------------------------------------------
# 메인 실행
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # 디렉토리명 = canonical-md5 (cache_key 사용)
    key = cache_key(BAP_SMILES)
    mol_dir = CACHE / f"MOL_{key}"
    mol_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== {BAP_NAME} 재생성 ===")
    print(f"SMILES    = {BAP_SMILES}")
    print(f"charge    = {BAP_CHARGE}")
    print(f"cache_key = {key}")
    print(f"디렉토리   = {mol_dir}")
    print(f"목표 파일  = {mol_dir / (mol_dir.name + '.orcacosmo')}")
    print()

    # 1. 초기 xyz 생성
    init_xyz = mol_dir / "init.xyz"
    if not init_xyz.exists():
        print("  [1/4] RDKit ETKDGv3 + MMFF로 초기 구조 생성...")
        generate_xyz_with_rdkit(BAP_SMILES, init_xyz)
    else:
        print("  [1/4] 초기 구조 있음 (skip)")

    # 2. xTB 최적화
    opt_xyz = mol_dir / "xtbopt.xyz"
    if not opt_xyz.exists():
        print("  [2/4] xTB 최적화 (--opt tight --gfn 2)...")
        run_xtb_optimization(init_xyz, mol_dir, charge=BAP_CHARGE)
    else:
        print("  [2/4] xTB 결과 있음 (skip)")

    # 3. ORCA COSMO single-point
    cpcm_file = mol_dir / "cosmo_sp.cpcm"
    if not cpcm_file.exists():
        print("  [3/4] ORCA COSMO single-point (BP86/def2-TZVPD, CPCM)...")
        run_orca_cosmo_sp(opt_xyz, mol_dir, charge=BAP_CHARGE)
    else:
        print("  [3/4] ORCA COSMO 결과 있음 (skip)")

    # 4. .orcacosmo 변환
    orcacosmo = mol_dir / f"{mol_dir.name}.orcacosmo"
    if not orcacosmo.exists():
        print("  [4/4] .orcacosmo 변환...")
        convert_to_orcacosmo(mol_dir, mol_dir.name)
    else:
        print("  [4/4] .orcacosmo 있음 (skip)")

    print(f"\n  [완료] {orcacosmo}")
    print(f"  크기   = {orcacosmo.stat().st_size} bytes")

    # 최종 검증
    print("\n=== 최종 검증 ===")
    print(f"dir 존재       = {mol_dir.exists()}")
    print(f".orcacosmo 존재 = {orcacosmo.exists()}")
    with open(orcacosmo, errors='ignore') as f:
        content = f.read()
    print(f"ENERGY 포함    = {'ENERGY' in content.upper()}")
    print(f"segment 포함   = {'segment' in content.lower()}")

    print("\n=== 전체 완료 ===")
