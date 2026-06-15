# Albumin (BSA) K_pw 검증 — openCOSMO-RS, muscle v4 포팅

muscle 파이프라인(v4)을 albumin으로 포팅. **변경점은 조성·검증셋·온도라벨뿐**,
openCOSMO-RS API·MD5 캐시키·TEMP_K(310.15)·ORCA 설정은 muscle 원본과 동일.

## 파일
- `calc_sigma_albumin_v4.py` — σ-profile 생성 (muscle calc_missing_molecules_v4.py 포팅)
  - 대상: ALBUMIN83 solute + capped residue 20종 + water
  - 캐시키: solute=MOL_<md5>, 아미노산=AA_<3letter>_<md5>, water=WATER.orcacosmo
- `cosmo_kpw_openrs_albumin_v4.py` — K_pw 계산·검증 (muscle cosmo_kpw_openrs_v4.py 포팅)
  - 조성 = BSA_AA_COMPOSITION (UniProt P02769 mature 583, 검증완료)
  - 검증셋 = ALBUMIN83 (Endo & Goss 2011 Table 1, 전부 37°C, CSV와 14/14 대조 일치)
- `run_calc_albumin_v4.sh` — SLURM (muscle run_calc_3mol_v4.sh 포팅)
- `_albumin83_list.py` — 검증 데이터 원본 리스트

## 실행
```bash
# XTB/ORCA 경로를 calc_sigma_albumin_v4.py 상단에서 서버 실제 경로로 확인
sbatch run_calc_albumin_v4.sh
# 또는 수동:
python calc_sigma_albumin_v4.py residues   # 잔기20+water
python calc_sigma_albumin_v4.py solutes    # 83종
python cosmo_kpw_openrs_albumin_v4.py validate
```

## muscle 잔기 재사용 (재계산 절약)
muscle에서 이미 만든 AA_*/WATER σ-profile을 그대로 복사하면 잔기 재계산 0:
```bash
mkdir -p cosmo_work_v4_albumin/sigma_cache
cp -r <muscle>/cosmo_work_v4/sigma_cache/AA_*      cosmo_work_v4_albumin/sigma_cache/
cp    <muscle>/cosmo_work_v4/sigma_cache/WATER.orcacosmo cosmo_work_v4_albumin/sigma_cache/
# 이러면 solutes만 새로 계산됨
python calc_sigma_albumin_v4.py solutes
python cosmo_kpw_openrs_albumin_v4.py validate
```
※ 단 capped residue가 **중성 protonation form**으로 muscle과 동일해야 재사용 가능
  (이 코드의 CAPPED_RESIDUE_SMILES가 muscle과 같은지 확인). 다르면 잔기도 새로 계산.

## 검증 해석 (Endo & Goss 2011 근거) — ★ 중요
- albumin은 **용매가 아님**: 저자가 Kn/Kc≈1, PP-LFER SD=0.41로 입증.
  → pseudo-solvent가 muscle(R²0.66)보다 **덜 맞는 게 물리적으로 정상**.
- **크기로 stratify 필수** (McGowan V≈1.2 ≈ logKow 4):
  - 작은 분자(V<1.2): 논문 PP-LFER SD=0.28 → pseudo-solvent도 잘 맞을 것 = "비특이 dissolution 재현"
  - 큰 분자(V>1.2): 논문 rmse=0.71 → pocket 크기제한으로 어긋날 것 = "특이결합 한계"
  → RMSE 하나로 보지 말고 V(또는 logKow) 경계로 나눠 해석.
- c_unit: 1차 c_unit=0.0으로 돌려 bias 확인 후, 전체 평균 또는 (더 깨끗하게) 작은분자(V<1.2)만으로 offset 결정.

## 한계
- openCOSMO-RS 24a neutral 전용: BSA charged residue(~34%) 중성처리. 단 solute 83종 전부 중성이라 solute측 무관.
- CYS 35개 free thiol 노출(실제 17 SS+1 free) 근사.
