#!/usr/bin/env python3
"""
Muscle protein-water 분배계수 계산 (openCOSMO-RS_py 사용) v4
xTB만 사용하는 방식 - DFT 제거
"""
import os
import math
import numpy as np
from pathlib import Path

# openCOSMO-RS_py imports
from opencosmorspy.cosmors import COSMORS
from opencosmorspy.parameterization import openCOSMORS24a
from opencosmorspy.input_parsers import SigmaProfileParser

# ----------------------------------------------------------------------
# 설정
# ----------------------------------------------------------------------
WORK = Path("./cosmo_work_v4")
CACHE = WORK / "sigma_cache"
TEMP_K = 310.15  # 체온 37°C (K)

# 아미노산 조성 (muscle protein) - actin + myosin 서열 기준
MUSCLE_AA_COMPOSITION = {
    "ALA": 0.078, "ARG": 0.052, "ASN": 0.041, "ASP": 0.060, "CYS": 0.013,
    "GLN": 0.042, "GLU": 0.090, "GLY": 0.069, "HIS": 0.024, "ILE": 0.055,
    "LEU": 0.090, "LYS": 0.080, "MET": 0.025, "PHE": 0.038, "PRO": 0.040,
    "SER": 0.058, "THR": 0.051, "TRP": 0.012, "TYR": 0.032, "VAL": 0.060,
}

# 아미노산 20종 리스트 ( capped residue Ace-X-Nme )
AA_ORDER = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

# Endo et al. 2012 검증셋
ENDO46 = [
    ("n-nonane",3.86,"CCCCCCCCC"),("cyclooctane",2.94,"C1CCCCCCC1"),("cycloheptane",2.16,"C1CCCCCC1"),
    ("1-nonene",3.54,"C=CCCCCCCC"),("hept-1-yne",1.75,"C#CCCCCC"),("1-chlorooctane",2.93,"CCCCCCCCCl"),
    ("dibutylether",1.72,"CCCCOCCCC"),("2-octanone",1.57,"CCCCCCC(C)=O"),("2-nonanone",1.83,"CCCCCCCC(C)=O"),
    ("2-decanone",2.11,"CCCCCCCCC(C)=O"),("1-nitrohexane",1.34,"CCCCCC[N+](=O)[O-]"),
    ("TBP",1.69,"CCCCOP(=O)(OCCCC)OCCCC"),("1-octanol",1.42,"CCCCCCCCO"),("1-nonanol",1.97,"CCCCCCCCCO"),
    ("propylbenzene",1.98,"CCCc1ccccc1"),("pentylbenzene",2.71,"CCCCCc1ccccc1"),("hexylbenzene",3.19,"CCCCCCc1ccccc1"),
    ("124-TMB",1.90,"Mc1ccc(C)c(C)c1"),("12-DCB",2.02,"Clc1ccccc1Cl"),("1234-TeCB",2.80,"Clc1ccc(Cl)c(Cl)c1Cl"),
    ("1245-TeCB",2.94,"Clc1cc(Cl)c(Cl)cc1Cl"),("pentaCB",3.35,"Clc1cc(Cl)c(Cl)c(Cl)c1Cl"),
    ("14-DBB",2.18,"Brc1ccc(Br)cc1"),("dibenzothiophene",3.40,"c1ccc2c(c1)sc1ccccc12"),
    ("phenanthrene",3.14,"c1ccc2ccc3ccccc3c2c1"),("pyrene",3.67,"c1cc2ccc3cccc4ccc(c1)c2c34"),
    ("chrysene",4.53,"c1ccc2ccc3c4ccccc4ccc3c2c1"),("BaP",4.88,"c1cc2ccc3cccc4ccc5cccc(c12)c3c45"),
    ("4-nitroanisole",1.38,"COc1ccc([N+](=O)[O-])cc1"),("valerophenone",1.74,"CCCCC(=O)c1ccccc1"),
    ("benzophenone",2.01,"O=C(c1ccccc1)c1ccccc1"),("2-nitrotoluene",1.27,"Mc1ccccc1[N+](=O)[O-]"),
    ("24-DNT",1.22,"Mc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]"),("4-propylphenol",1.69,"CCCc1ccc(O)cc1"),
    ("4-fluorophenol",1.18,"Oc1ccc(F)cc1"),("3-chlorophenol",1.57,"Oc1cccc(Cl)c1"),("4-chlorophenol",1.55,"Oc1ccc(Cl)cc1"),
    ("4-bromophenol",1.70,"Oc1ccc(Br)cc1"),("4-iodophenol",1.93,"Oc1ccc(I)cc1"),("2-phenylphenol",2.00,"Oc1ccccc1-c1ccccc1"),
    ("4-chloroaniline",1.19,"Nc1ccc(Cl)cc1"),("4-iodoaniline",1.66,"Nc1ccc(I)cc1"),("4-nitroaniline",0.84,"Nc1ccc([N+](=O)[O-])cc1"),
    ("NN-diethylaniline",1.60,"CCN(CC)c1ccccc1"),("indole",1.40,"c1ccc2[nH]ccc2c1"),
    ("diazepam",2.01,"CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
]


def get_aa_orcacosmo_files():
    """단백질 pseudo-solvent를 구성하는 아미노산 .orcacosmo 파일들 반환"""
    aa_files = []
    for aa in AA_ORDER:
        # AA_XXX_hash 형식의 디렉토리 찾기
        for aa_dir in CACHE.glob(f"AA_{aa}_*"):
            orcacosmo = aa_dir / f"{aa_dir.name}.orcacosmo"
            if orcacosmo.exists():
                aa_files.append((str(orcacosmo), MUSCLE_AA_COMPOSITION[aa]))
                break
        else:
            print(f"[WARNING] AA_{aa} not found in cache")

    # 조성 정규화 (합=1.0)
    total = sum(frac for _, frac in aa_files)
    aa_files = [(path, frac/total) for path, frac in aa_files]

    return aa_files


def lng_inf(crs, solute_path, solvent_paths, solvent_x, T):
    """
    solute의 무한희석 lnγ∞ 계산

    Parameters
    ----------
    crs : COSMORS
        COSMORS 인스턴스
    solute_path : str
        용질의 .orcacosmo 파일 경로
    solvent_paths : list of str
        용매(들)의 .orcacosmo 파일 경로 리스트
    solvent_x : list of float
        용매(들)의 mole fraction (합=1)
    T : float
        온도 (K)

    Returns
    -------
    float
        lnγ∞ (용질의 무한희석 activity coefficient)
    """
    crs.clear_molecules()
    crs.add_molecule([solute_path])  # index 0 = solute (리스트 인자!)

    for p in solvent_paths:
        crs.add_molecule([p])

    crs.clear_jobs()
    x = np.array([0.0] + list(solvent_x))  # solute=0 → 무한희석, 나머지 합=1
    crs.add_job(x, T, refst='pure_component')

    result = crs.calculate()
    return result['tot']['lng'][0][0]  # [분자0=solute][job0]


def calculate_partition_coefficient(solute_orcacosmo, aa_files, water_orcacosmo, T=TEMP_K, c_unit=0.0):
    """
    분배계수 log K_pw 계산

    Parameters
    ----------
    solute_orcacosmo : str
        용질의 .orcacosmo 파일 경로
    aa_files : list of (str, float)
        아미노산 (파일경로, 몰분율) 리스트
    water_orcacosmo : str
        물의 .orcacosmo 파일 경로
    T : float
        온도 (K)
    c_unit : float
        단위 변환 상수 (log10(Vm_water/Vm_protein))

    Returns
    -------
    float or None
        log K_pw, 에러 시 None
    """
    try:
        crs = COSMORS(par=openCOSMORS24a())

        # --- (a) 물에서 ---
        aa_x = [frac for _, frac in aa_files]
        lng_w = lng_inf(crs, solute_orcacosmo, [water_orcacosmo], [1.0], T)

        # --- (b) 단백질 혼합상에서 ---
        aa_paths = [path for path, _ in aa_files]
        lng_p = lng_inf(crs, solute_orcacosmo, aa_paths, aa_x, T)

        # --- (c) 분배계수 ---
        ln_Kx = lng_w - lng_p
        log_Kpw = ln_Kx / math.log(10) + c_unit

        return float(log_Kpw)

    except Exception as e:
        import traceback
        print(f"  [COSMO-RS error: {e}]")
        traceback.print_exc()
        return None


def find_orcacosmo_file(cache_key_str):
    """캐시 키로 .orcacosmo 파일 찾기"""
    for mol_dir in CACHE.glob("MOL_*"):
        if cache_key_str in mol_dir.name:
            orcacosmo = mol_dir / f"{mol_dir.name}.orcacosmo"
            if orcacosmo.exists():
                return str(orcacosmo)
    return None


def validate_endo46(label="openrs_validation_v4"):
    """
    Endo 46종 검증 (WATER가 있을 때만 실행)
    v4: xTB만 사용
    """
    print("=== openCOSMO-RS 검증 시작 v4 (xTB only) ===")

    # 물의 .orcacosmo 파일 확인
    water_orcacosmo = str(CACHE / "WATER.orcacosmo")
    if not Path(water_orcacosmo).exists():
        print(f"[ERROR] {water_orcacosmo} not found.")
        return None

    print(f"Water file: {water_orcacosmo}")

    # 단백질 pseudo-solvent 파일들
    aa_files = get_aa_orcacosmo_files()
    print(f"단백질 아미노산 수: {len(aa_files)}/20")

    if len(aa_files) < 20:
        print("[WARNING] 일부 아미노산이 missing입니다")

    exp_vals, pred_vals, names = [], [], []
    results_list = []

    import hashlib

    for name, exp_k, smiles in ENDO46:
        # rdkit 없이 직접 SMILES 사용
        key = hashlib.md5(smiles.encode()).hexdigest()[:16]
        orcacosmo = find_orcacosmo_file(key)

        if orcacosmo is None:
            print(f"[skip] {name}: .orcacosmo 파일 없음")
            continue

        print(f"[*] {name}")
        log_k = calculate_partition_coefficient(
            orcacosmo, aa_files, water_orcacosmo, TEMP_K, c_unit=0.0
        )

        if log_k is not None:
            exp_vals.append(exp_k)
            pred_vals.append(log_k)
            names.append(name)
            results_list.append({"name": name, "exp": float(exp_k), "pred": float(log_k)})
            print(f"  -> exp={exp_k:.2f}, pred={log_k:.3f}")

    # 통계
    if len(exp_vals) >= 2:
        exp_arr = np.array(exp_vals)
        pred_arr = np.array(pred_vals)

        r2 = np.corrcoef(pred_arr, exp_arr)[0, 1]**2
        rmse = np.sqrt(np.mean((pred_arr - exp_arr)**2))

        import json
        output = {
            "label": label,
            "n": len(exp_arr),
            "r2": float(r2),
            "rmse": float(rmse),
            "rows": results_list
        }

        (WORK / f"results_{label}.json").write_text(json.dumps(output, indent=2))

        print(f"\n=== 결과 ===")
        print(f"n = {len(exp_arr)}")
        print(f"R² = {r2:.3f}")
        print(f"RMSE = {rmse:.3f} log units")

        return output
    else:
        print("데이터 부족으로 통계 계산 불가")
        return None


def test_sigma_profile():
    """Sigma profile 테스트"""
    # 첫 번째 아미노산으로 테스트
    aa_files = get_aa_orcacosmo_files()
    if not aa_files:
        print("[ERROR] 아미노산 .orcacosmo 파일 없음")
        return

    aa_path = aa_files[0][0]
    print(f"Testing sigma profile: {aa_path}")

    spp = SigmaProfileParser(aa_path)

    # 기본 정보 출력
    print(f"Name: {spp.get('name', 'N/A')}")
    print(f"Method: {spp.get('method', 'N/A')}")
    print(f"Area: {spp.get('area', 0):.2f} Å²")
    print(f"Volume: {spp.get('volume', 0):.2f} Å³")
    print(f"Atoms: {len(spp.get('atm_nr', []))}")
    print(f"Segments: {len(spp.get('seg_nr', []))}")

    # sigma moments (가능한 경우)
    try:
        spp.calculate_sigma_moments()
        print(f"sigma_moments: {spp['sigma_moments']}")
    except Exception as e:
        print(f"[WARNING] sigma_moments 계산 실패: {e}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        if sys.argv[1] == "test":
            test_sigma_profile()
        elif sys.argv[1] == "validate":
            validate_endo46()
    else:
        print("Usage:")
        print("  python cosmo_kpw_openrs_v4.py test      # sigma profile 테스트")
        print("  python cosmo_kpw_openrs_v4.py validate  # Endo 46종 검증")
