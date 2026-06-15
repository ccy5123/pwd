# Structural Protein (Muscle) K_pw 검증 — openCOSMO-RS

## 데이터셋
- **MUSCLE46**: muscle protein-water partition coefficient
- 46종 화합물, n = 46
- 온도: 310.15 K (37°C)

## 결과
- **R² = 0.66**
- RMSE: (결과 파일 참조)

## 파일
- `calc_missing_molecules_v4.py` — σ-profile 생성
- `cosmo_kpw_openrs_v4.py` — K_pw 계산·검증
- `cpcm_to_orcacosmo.py` — CPQM → ORCA COSMO 변환
- `cosmo_work_v4/combined_results_46.json` — 최종 결과 (R²=0.66)
- `cosmo_work_v3/results_openrs_validation.json` — 검증 결과 (v3)

## 실행 방법
```bash
python calc_missing_molecules_v4.py residues   # 잔기 계산
python calc_missing_molecules_v4.py solutes    # 용질 계산
python cosmo_kpw_openrs_v4.py validate
```

## 한계
- openCOSMO-RS 24a neutral 전용
- Muscle protein은 복잡한 구조로 예측이 어려움
