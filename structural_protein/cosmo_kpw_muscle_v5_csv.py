#!/usr/bin/env python3
"""
Muscle protein-water 분배계수 계산 (openCOSMO-RS_py 사용) v5 CSV
- CSV 단일소스 리팩터
- canonical-md5 캐시 조회 통일
"""
import os
import math
import csv
import numpy as np
from pathlib import Path

# openCOSMO-RS_py imports
from opencosmorspy.cosmors import COSMORS
from opencosmorspy.parameterization import openCOSMORS24a
from opencosmorspy.input_parsers import SigmaProfileParser

# ----------------------------------------------------------------------
# 설정
# ----------------------------------------------------------------------
WORK = Path("./cosmo_work_v3")
CACHE = WORK / "sigma_cache"
TEMP_K = 310.15  # 체온 37°C (K)

# 아미노산 조성 (muscle protein) - actin + myosin 서열 기준
MUSCLE_AA_COMPOSITION = {
    # Chicken (Gallus gallus) myofibrillar core-6, sequence-based, stoichiometry actin:myosin:Tm:TnT:TnI:TnC = 7:1:1:1:1:1
    # Sequences (UniProt): actin P68139, myosin P13538, Tm P04268, TnT P12620, TnI P68246, TnC P02588 (all OX=9031)
    # Stoichiometry: Potter 1974 (Arch.Biochem.Biophys.162:436); Parvatiyar&Pinto 2022 (10.1016/j.abb.2022.109301); Yates 1983
    "ALA": 0.0841, "ARG": 0.0526, "ASN": 0.0329, "ASP": 0.0587, "CYS": 0.0115,
    "GLN": 0.0439, "GLU": 0.1128, "GLY": 0.0556, "HIS": 0.0219, "ILE": 0.0623,
    "LEU": 0.0841, "LYS": 0.0837, "MET": 0.0367, "PHE": 0.0291, "PRO": 0.0329,
    "SER": 0.0541, "THR": 0.0563, "TRP": 0.0075, "TYR": 0.0307, "VAL": 0.0486,
}

# 아미노산 20종 리스트 ( capped residue Ace-X-Nme )
AA_ORDER = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

# CSV 경로 (__file__ 기준 후보탐색, cwd-비의존)
_HERE = Path(__file__).resolve().parent
_CSV_CANDIDATES = [
    _HERE / "3D" / "compounds_input_table1.csv",          # Muscle_protein/3D/
    _HERE.parent / "3D" / "compounds_input_table1.csv",   # /gpfs/home1/s9383/3D/
]
CSV_PATH = next((str(p) for p in _CSV_CANDIDATES if p.exists()), None)
if CSV_PATH is None:
    raise FileNotFoundError(f"compounds_input_table1.csv not found in {_CSV_CANDIDATES}")


def load_validation_set_from_csv():
    """
    CSV에서 muscle_chicken_exp notna 행을 (compound, exp, SMILES)로 로드
    """
    rows = []
    with open(CSV_PATH) as f:
        reader = csv.DictReader(f)
        for r in reader:
            exp = r['muscle_chicken_exp']
            if exp not in ('', 'None', None):
                rows.append((r['compound'], float(exp), r['SMILES']))
    return rows


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
    """캐시 키로 .orcacosmo 파일 찾기 (canonical-md5 기준)"""
    for mol_dir in CACHE.glob("MOL_*"):
        if cache_key_str in mol_dir.name:
            orcacosmo = mol_dir / f"{mol_dir.name}.orcacosmo"
            if orcacosmo.exists():
                return str(orcacosmo)
    return None


def validate(label="chicken_core6_seq_v5"):
    """
    CSV 검증셋 실행 (WATER가 있을 때만 실행)
    v5: CSV 단일소스, canonical-md5 캐시 조회
    """
    print("=== openCOSMO-RS 검증 시작 v5 (CSV single-source, canonical-md5) ===")

    # CSV에서 검증셋 로드
    val_set = load_validation_set_from_csv()
    print(f"CSV에서 {len(val_set)}종 로드 (46이어야)")
    assert len(val_set) == 46, f"CSV 46종 아님: {len(val_set)}"

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

    # canonical-md5 통일: cosmo_kpw_common_v3의 cache_key 사용
    from cosmo_kpw_common_v3 import cache_key

    for name, exp_k, smiles in val_set:
        # canonical-md5로 캐시 조회
        key = cache_key(smiles)
        orcacosmo = find_orcacosmo_file(key)

        if orcacosmo is None:
            print(f"[skip] {name}: .orcacosmo 파일 없음 (key={key})")
            continue

        print(f"[*] {name}")
        log_k = calculate_partition_coefficient(
            orcacosmo, aa_files, water_orcacosmo, TEMP_K, c_unit=-0.80  # = log10(Vm_water / M_residue), residue-MW basis (rho=1); recomputed for this composition = -0.799
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

        # Pearson R² (offset 불변)
        r2_pearson = np.corrcoef(pred_arr, exp_arr)[0, 1]**2
        # CoD R² (1:1 기준, c_unit 효과 반영됨)
        ss_res = np.sum((pred_arr - exp_arr)**2)
        ss_tot = np.sum((exp_arr - exp_arr.mean())**2)
        r2_cod = 1 - ss_res / ss_tot
        # RMSE
        rmse = np.sqrt(np.mean((pred_arr - exp_arr)**2))
        # Within 10× (log 단위 1 = 10배)
        within_10x = np.mean(np.abs(pred_arr - exp_arr) <= 1.0) * 100

        import json
        output = {
            "label": label,
            "n": len(exp_arr),
            "r2_pearson": float(r2_pearson),
            "r2_cod": float(r2_cod),
            "rmse": float(rmse),
            "within_10x_pct": float(within_10x),
            "c_unit": -0.80,
            "c_unit_basis": "log10(Vm_water / M_residue), residue-MW basis (rho=1); recomputed for this composition = -0.799",
            "rows": results_list
        }

        (WORK / f"results_{label}.json").write_text(json.dumps(output, indent=2))

        print(f"\n=== 결과 ===")
        print(f"n = {len(exp_arr)}")
        print(f"R²_cod = {r2_cod:.3f} (1:1 기준)")
        print(f"R²_pearson = {r2_pearson:.3f} (offset 불변)")
        print(f"RMSE = {rmse:.3f} log units")
        print(f"Within 10× = {within_10x:.1f}%")

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
            validate()
    else:
        print("Usage:")
        print("  python cosmo_kpw_muscle_v5_csv.py test      # sigma profile 테스트")
        print("  python cosmo_kpw_muscle_v5_csv.py validate  # CSV 검증셋 검증")
