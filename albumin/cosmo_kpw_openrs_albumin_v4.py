#!/usr/bin/env python3
"""
Albumin(BSA)-water 분배계수 계산 (openCOSMO-RS_py 사용) v4-albumin
muscle 버전(cosmo_kpw_openrs_v4.py)에서 변경점:
  (1) 조성 MUSCLE_AA_COMPOSITION -> BSA_AA_COMPOSITION (UniProt P02769 mature 583)
  (2) 검증셋 ENDO46 -> ALBUMIN83 (Endo & Goss 2011, Table 1, 전부 37C)
나머지(openCOSMO-RS API, MD5 캐시키, TEMP_K=310.15, ORCA 설정)는 muscle 원본과 동일.

검증 해석 (Endo&Goss 2011 근거):
  - albumin은 용매가 아님(Kn/Kc~1, PP-LFER SD 0.41) -> pseudo-solvent가 muscle보다
    덜 맞는 게 정상.
  - McGowan V<1.2 (≈logKow<4) 작은 분자: 잘 맞을 것(논문 SD 0.28).
    V>1.2 큰 분자: pocket 크기 제한으로 어긋날 것(논문 rmse 0.71).
  -> 검증은 V(또는 logKow)로 stratify해서 볼 것.
"""
import os
import math
import numpy as np
from pathlib import Path

from opencosmorspy.cosmors import COSMORS
from opencosmorspy.parameterization import openCOSMORS24a
from opencosmorspy.input_parsers import SigmaProfileParser

# ----------------------------------------------------------------------
WORK = Path("./cosmo_work_v4_albumin")
CACHE = WORK / "sigma_cache"
TEMP_K = 310.15  # 37C — Endo&Goss 2011 측정온도 (muscle과 동일)

# 아미노산 조성 (BSA) - UniProt P02769 mature chain 25-607 (583 aa), 검증완료
BSA_AA_COMPOSITION = {
    "ALA": 0.0806, "ARG": 0.0395, "ASN": 0.0240, "ASP": 0.0686, "CYS": 0.0600,
    "GLN": 0.0343, "GLU": 0.1012, "GLY": 0.0274, "HIS": 0.0292, "ILE": 0.0240,
    "LEU": 0.1046, "LYS": 0.1012, "MET": 0.0069, "PHE": 0.0463, "PRO": 0.0480,
    "SER": 0.0480, "THR": 0.0566, "TRP": 0.0034, "TYR": 0.0343, "VAL": 0.0617,
}

AA_ORDER = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

# Endo & Goss 2011 검증셋 (BSA, 37C)
ALBUMIN83 = [
    ('1,2,4-Trimethylbenzene', 3.35, 'Cc1ccc(C)c(C)c1'),
    ('1,2-Dichlorobenzene', 3.03, 'Clc1ccccc1Cl'),
    ('1-Hexanol', 1.64, 'CCCCCCO'),
    ('Benzene', 1.58, 'c1ccccc1'),
    ('Carbon tetrachloride', 1.77, 'ClC(Cl)(Cl)Cl'),
    ('Halothane', 1.62, 'FC(F)(F)C(Cl)Br'),
    ('Enflurane', 1.59, 'FC(F)OC(F)(F)C(F)Cl'),
    ('Isoflurane', 1.58, 'FC(F)OC(Cl)C(F)(F)F'),
    ('Chlorobenzene', 2.32, 'Clc1ccccc1'),
    ('Cycloheptane', 2.52, 'C1CCCCCC1'),
    ('Cyclohexane', 2.01, 'C1CCCCC1'),
    ('Ethylbenzene', 2.7, 'CCc1ccccc1'),
    ('Heptane', 3.59, 'CCCCCCC'),
    ('Hexafluorobenzene', 1.55, 'Fc1c(F)c(F)c(F)c(F)c1F'),
    ('Hexane', 3.09, 'CCCCCC'),
    ('Methoxyflurane', 1.77, 'COC(F)(F)C(Cl)Cl'),
    ('Methylpentafluorobenzene', 2.32, 'Cc1c(F)c(F)c(F)c(F)c1F'),
    ('Nonane', 4.45, 'CCCCCCCCC'),
    ('Propylbenzene', 2.95, 'CCCc1ccccc1'),
    ('Tetrachloroethene', 2.4, 'ClC(Cl)=C(Cl)Cl'),
    ('Toluene', 2.26, 'Cc1ccccc1'),
    ('Trichloroethene', 1.88, 'ClC=C(Cl)Cl'),
    ('1-heptanol', 2.18, 'CCCCCCCO'),
    ('di-n-butyl ether', 2.01, 'CCCCOCCCC'),
    ('Cyclooctane', 2.98, 'C1CCCCCCC1'),
    ('Octan-1-ol', 2.74, 'CCCCCCCCO'),
    ('Nonan-1-ol', 3.1, 'CCCCCCCCCO'),
    ('1-Chlorooctane', 3.85, 'CCCCCCCCCl'),
    ('1-Nonene', 4.22, 'C=CCCCCCCC'),
    ('1,2,4-Trichlorobenzene', 3.6, 'Clc1ccc(Cl)c(Cl)c1'),
    ('Dipentyl ether', 3.0, 'CCCCCOCCCCC'),
    ('2-octanone', 2.09, 'CCCCCCC(C)=O'),
    ('2-nonanone', 2.48, 'CCCCCCCC(C)=O'),
    ('2-Nitrotoluene', 2.12, 'Cc1ccccc1[N+](=O)[O-]'),
    ('4-Nitroanisole', 2.48, 'COc1ccc([N+](=O)[O-])cc1'),
    ('4-Chlorophenol', 2.43, 'Oc1ccc(Cl)cc1'),
    ('4-Ethyl-3-hexanol', 1.48, 'CCC(O)C(CC)CC'),
    ('2,4-Dinitrotoluene', 1.73, 'Cc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]'),
    ('Styrene', 2.76, 'C=Cc1ccccc1'),
    ('Octane', 4.01, 'CCCCCCCC'),
    ('3-Chlorophenol', 2.35, 'Oc1cccc(Cl)c1'),
    ('4-bromophenol', 2.81, 'Oc1ccc(Br)cc1'),
    ('4-Iodoaniline', 2.95, 'Nc1ccc(I)cc1'),
    ('N,N-Diethylaniline', 2.27, 'CCN(CC)c1ccccc1'),
    ('4-n-Propylphenol', 2.59, 'CCCc1ccc(O)cc1'),
    ('4-iodophenol', 3.41, 'Oc1ccc(I)cc1'),
    ('Phenanthrene', 4.15, 'c1ccc2c(c1)ccc1ccccc12'),
    ('Fluoranthene', 4.28, 'c1ccc2c(c1)-c1cccc3cccc-2c13'),
    ('Pyrene', 4.76, 'c1cc2ccc3cccc4ccc(c1)c2c34'),
    ('1-heptyne', 2.49, 'C#CCCCCC'),
    ('tribromomethane', 1.95, 'BrC(Br)Br'),
    ('γ-hexachlorocyclohexane', 2.46, 'ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl'),
    ('2-decanone', 2.88, 'CCCCCCCCC(C)=O'),
    ('1-nitrooctane', 3.38, 'CCCCCCCC[N+](=O)[O-]'),
    ('tri-n-butyl phosphate', 2.47, 'CCCCOP(=O)(OCCCC)OCCCC'),
    ('1,4-dibromobenzene', 3.97, 'Brc1ccc(Br)cc1'),
    ('1,2,3,4-tetrachlorobenzene', 4.21, 'Clc1ccc(Cl)c(Cl)c1Cl'),
    ('indene', 2.92, 'C1=Cc2ccccc2C1'),
    ('naphthalene', 3.56, 'c1ccc2ccccc2c1'),
    ('dibenzofuran', 3.79, 'c1ccc2c(c1)oc1ccccc12'),
    ('dibenzothiophene', 4.16, 'c1ccc2c(c1)sc1ccccc12'),
    ('chrysene', 4.46, 'c1ccc2c(c1)ccc1c3ccccc3ccc21'),
    ('benzo[b]fluoranthene', 4.42, 'c1ccc2c(c1)-c1cccc3c1c-2cc1ccccc13'),
    ('benzo[ghi]perylene', 4.76, 'c1cc2ccc3ccc4ccc5cccc6c(c1)c2c3c4c56'),
    ('anisole', 2.16, 'COc1ccccc1'),
    ('valerophenone', 2.7, 'CCCCC(=O)c1ccccc1'),
    ('benzophenone', 2.62, 'O=C(c1ccccc1)c1ccccc1'),
    ('di-n-propyl phthalate', 2.84, 'CCCOC(=O)c1ccccc1C(=O)OCCC'),
    ('1-nitronaphthalene', 3.17, 'O=[N+]([O-])c1cccc2ccccc12'),
    ('4-chlorobenzyl alcohol', 2.1, 'OCc1ccc(Cl)cc1'),
    ('2-phenylphenol', 2.62, 'Oc1ccccc1-c1ccccc1'),
    ('4-fluorophenol', 1.57, 'Oc1ccc(F)cc1'),
    ('bisphenol A', 2.88, 'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1'),
    ('4-nitroaniline', 1.69, 'Nc1ccc([N+](=O)[O-])cc1'),
    ('2-chloroaniline', 1.95, 'Nc1ccccc1Cl'),
    ('4-aminobiphenyl', 2.55, 'Nc1ccc(-c2ccccc2)cc1'),
    ('indole', 2.25, 'c1ccc2[nH]ccc2c1'),
    ('carbazole', 3.52, 'c1ccc2c(c1)[nH]c1ccccc12'),
    ('metolachlor', 1.74, 'CCc1cccc(C)c1N(C(=O)CCl)C(C)COC'),
    ('atrazine', 1.77, 'CCNc1nc(Cl)nc(NC(C)C)n1'),
    ('diazepam', 2.68, 'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21'),
    ('estrone', 2.69, 'CC12CCC3c4ccc(O)cc4CCC3C1CCC2=O'),
    ('endosulfan α', 3.24, 'O=S1OCC2C(CO1)C1(Cl)C(Cl)=C(Cl)C2(Cl)C1(Cl)Cl'),
]


def get_aa_orcacosmo_files():
    """단백질 pseudo-solvent를 구성하는 아미노산 .orcacosmo 파일들 반환"""
    aa_files = []
    for aa in AA_ORDER:
        for aa_dir in CACHE.glob(f"AA_{aa}_*"):
            orcacosmo = aa_dir / f"{aa_dir.name}.orcacosmo"
            if orcacosmo.exists():
                aa_files.append((str(orcacosmo), BSA_AA_COMPOSITION[aa]))
                break
        else:
            print(f"[WARNING] AA_{aa} not found in cache")
    total = sum(frac for _, frac in aa_files)
    aa_files = [(path, frac/total) for path, frac in aa_files]
    return aa_files


def lng_inf(crs, solute_path, solvent_paths, solvent_x, T):
    """solute의 무한희석 lnγ∞ 계산 (muscle 원본과 동일)"""
    crs.clear_molecules()
    crs.add_molecule([solute_path])
    for p in solvent_paths:
        crs.add_molecule([p])
    crs.clear_jobs()
    x = np.array([0.0] + list(solvent_x))
    crs.add_job(x, T, refst='pure_component')
    result = crs.calculate()
    return result['tot']['lng'][0][0]


def calculate_partition_coefficient(solute_orcacosmo, aa_files, water_orcacosmo, T=TEMP_K, c_unit=0.0):
    """log K_pw 계산 (muscle 원본과 동일)"""
    try:
        crs = COSMORS(par=openCOSMORS24a())
        aa_x = [frac for _, frac in aa_files]
        lng_w = lng_inf(crs, solute_orcacosmo, [water_orcacosmo], [1.0], T)
        aa_paths = [path for path, _ in aa_files]
        lng_p = lng_inf(crs, solute_orcacosmo, aa_paths, aa_x, T)
        ln_Kx = lng_w - lng_p
        log_Kpw = ln_Kx / math.log(10) + c_unit
        return float(log_Kpw)
    except Exception as e:
        import traceback
        print(f"  [COSMO-RS error: {e}]")
        traceback.print_exc()
        return None


def find_orcacosmo_file(cache_key_str):
    """캐시 키로 .orcacosmo 파일 찾기 (muscle 원본과 동일)"""
    for mol_dir in CACHE.glob("MOL_*"):
        if cache_key_str in mol_dir.name:
            orcacosmo = mol_dir / f"{mol_dir.name}.orcacosmo"
            if orcacosmo.exists():
                return str(orcacosmo)
    return None


def validate_albumin83(label="openrs_albumin_v4"):
    """Albumin 83종 검증 (WATER 있을 때만)"""
    print("=== openCOSMO-RS Albumin 검증 시작 v4 (37C) ===")
    water_orcacosmo = str(CACHE / "WATER.orcacosmo")
    if not Path(water_orcacosmo).exists():
        print(f"[ERROR] {water_orcacosmo} not found."); return None
    print(f"Water file: {water_orcacosmo}")

    aa_files = get_aa_orcacosmo_files()
    print(f"단백질 아미노산 수: {len(aa_files)}/20")
    if len(aa_files) < 20:
        print("[WARNING] 일부 아미노산이 missing입니다")

    exp_vals, pred_vals, names = [], [], []
    results_list = []
    import hashlib

    for name, exp_k, smiles in ALBUMIN83:
        key = hashlib.md5(smiles.encode()).hexdigest()[:16]
        orcacosmo = find_orcacosmo_file(key)
        if orcacosmo is None:
            print(f"[skip] {name}: .orcacosmo 파일 없음"); continue
        print(f"[*] {name}")
        log_k = calculate_partition_coefficient(orcacosmo, aa_files, water_orcacosmo, TEMP_K, c_unit=0.0)
        if log_k is not None:
            exp_vals.append(exp_k); pred_vals.append(log_k); names.append(name)
            results_list.append({"name": name, "exp": float(exp_k), "pred": float(log_k)})
            print(f"  -> exp={exp_k:.2f}, pred={log_k:.3f}")

    if len(exp_vals) >= 2:
        exp_arr = np.array(exp_vals); pred_arr = np.array(pred_vals)
        r2 = np.corrcoef(pred_arr, exp_arr)[0, 1]**2
        rmse = np.sqrt(np.mean((pred_arr - exp_arr)**2))
        bias = float(np.mean(pred_arr - exp_arr))
        import json
        output = {"label": label, "n": len(exp_arr), "r2": float(r2),
                  "rmse": float(rmse), "bias": bias, "rows": results_list}
        (WORK / f"results_{label}.json").write_text(json.dumps(output, indent=2))
        print(f"\n=== 결과 ===")
        print(f"n = {len(exp_arr)}")
        print(f"R2 = {r2:.3f}")
        print(f"RMSE = {rmse:.3f} log units  (bias={bias:+.3f})")
        print(f"제안 c_unit(전체 평균 offset) = {-bias:+.3f}")
        return output
    else:
        print("데이터 부족"); return None


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "validate":
        validate_albumin83()
    else:
        print("Usage: python cosmo_kpw_openrs_albumin_v4.py validate")
