# BSA Docking + xTB Desolvation Factorial Experiment

AutoDock Vina docking에 xTB desolvation 보정을 더한 BSA(Albumin)-water 분배계수 예측 실험.

## 연구 목적

중성 유기물의 albumin-water 분배계수(log K_BSA/w) 예측에서, docking 친화도(ΔG)에 xTB 수상 desolvation(Gsolv)을 더한 보정이 예측 정확도를 개선하는지 15셀 factorial로 검증.

## 실험 설계

### 변수

| 변수 | 수준 | 설명 |
|------|------|------|
| **포즈 방식** | 3 | A(최저ΔG 1개), B1(Boltzmann 가중평균), B2(산술평균) |
| **Desolvation** | 5 | none, GFN1/gbsa, GFN1/alpb, GFN2/gbsa, GFN2/alpb |
| **총 셀** | 15 | 3 × 5 = 15 combinations |

### 보정식

```
ΔG_corrected = ΔG_contact − Gsolv
log K_pw = logKA − 1.823
logKA = −ΔG / (2.302585 × RT)
RT = 0.616 kcal/mol (310.15 K)
```

## 실험 프로토콜

### Docking

- 단백질: BSA 4F5S (PDB ID: 4F5S), chain A
- 박스: 3×3×1 = 9 centers, 26×26×26 Å, spacing 1.0 Å
- Vina: exhaustiveness=32, n_poses=20, seed=42
- Ligand: RDKit ETKDGv3 (seed=42) → MMFF → Meeko PDBQT, 중성 유지

### Desolvation

- xTB GFN1-xTB / GFN2-xTB
- Solvation: GBSA(water) / ALPB(water)
- Conformer: ETKDGv3 seed 고정 (재현성)

## 실험 결과

### 15셀 Factorial 결과

| Pose | Desolv | absR² | calR² | Pearson | Spearman | RMSE |
|------|--------|--------|-------|---------|----------|------|
| A | none | **-0.676** | 0.219 | 0.468 | 0.321 | 1.116 |
| A | GFN1/gbsa | -24.818 | 0.125 | 0.353 | 0.439 | 4.380 |
| A | GFN1/alpb | -24.280 | 0.147 | 0.383 | 0.413 | 4.334 |
| A | GFN2/gbsa | -25.484 | 0.083 | 0.288 | 0.284 | 4.436 |
| A | GFN2/alpb | -30.304 | 0.128 | 0.358 | 0.384 | 4.823 |
| B1 | none | -1.553 | 0.208 | 0.457 | 0.286 | 1.377 |
| B1 | GFN1/gbsa | -29.407 | 0.116 | 0.341 | 0.421 | 4.753 |
| B1 | GFN1/alpb | -28.813 | 0.138 | 0.371 | 0.399 | 4.706 |
| B1 | GFN2/gbsa | -30.097 | 0.077 | 0.277 | 0.279 | 4.807 |
| B1 | GFN2/alpb | -35.126 | 0.121 | 0.348 | 0.374 | 5.181 |
| B2 | none | -5.244 | 0.154 | 0.393 | 0.187 | 2.154 |
| B2 | GFN1/gbsa | -41.938 | 0.069 | 0.263 | 0.307 | 5.648 |
| B2 | GFN1/alpb | -41.053 | 0.088 | 0.297 | 0.292 | 5.590 |
| B2 | GFN2/gbsa | -42.354 | 0.041 | 0.203 | 0.165 | 5.676 |
| B2 | GFN2/alpb | -47.636 | 0.083 | 0.289 | 0.327 | 6.011 |

### 핵심 발견

**Desolvation 보정은 absR²를 모든 셀에서 크게 하락시킴**

| Pose | none (absR²) | 최고 desolv (absR²) | 변화 |
|------|--------------|---------------------|------|
| A | -0.676 | -24.28 | △ -23.60 |
| B1 | -1.553 | -28.81 | △ -27.26 |
| B2 | -5.244 | -41.05 | △ -35.81 |

**결론**: Desolvation 미적용(A×none)이 최고 성능 (absR² = -0.676). 보정식 ΔG_corr = ΔG − Gsolv이 상관구조를 망가뜨림.

## 파일 설명

### 코드

| 파일 | 설명 |
|------|------|
| `bsa_docking_desolv_pipeline_v2.py` | 메인 파이프라인 (docking + 분석) |
| `gsolv_fixed.py` | xTB Gsolv 계산 모듈 |
| `run_gsolv.py` | Gsolv 일괄 계산 스크립트 |

### 데이터

| 파일 | 설명 |
|------|------|
| `docking_poses.csv` | Docking 결과 (2867 poses, 83 compounds) |
| `gsolv.csv` | Gsolv 값 (332 rows = 83 × 4 combo) |
| `predictions_all_cells.csv` | 15셀 예측 값 (1245 rows = 83 × 15) |
| `results_factorial.csv` | 15셀 요약 통계 |

### 시각화

| 파일 | 설명 |
|------|------|
| `absR2_heatmap.png` | 3×5 absR² 히트맵 |
| `calR2_heatmap.png` | 3×5 calR² 히트맵 |

## 실행 방법

### Docking만 실행

```bash
conda activate rapids-cuml
python bsa_docking_desolv_pipeline_v2.py --docking-only
```

### Smoke test

```bash
python bsa_docking_desolv_pipeline_v2.py --smoke --docking-only
```

### 분석 (결과 있을 때)

```bash
python bsa_docking_desolv_pipeline_v2.py --analyze
```

## 의존성

- Python 3.x
- vina (AutoDock Vina Python wrapper)
- meeko (PDBQT 변환)
- rdkit (화학 구조 처리)
- xtb (GFN-xTB, desolvation 계산)
- pdb2pqr (단백질 protonation)
- openbabel (파일 변환)
- scikit-learn (지표 계산)

## 참고 문헌

- Chen, R., Muensterman, D., Field, M., Ng, C. (2025). "Deriving Membrane–Water and Protein–Water Partition Coefficients from In Vitro Experiments for Per- and Polyfluoroalkyl Substances (PFAS)." *Environmental Science & Technology*, 59(1), 82-91. (본문·SI read ✓)

## 실험 일지

- 실험 완료일: 2026-07-02
- Docking job: 718126 (COMPLETED, 10분)
- Gsolv job: 718148 (COMPLETED, 3분)
- 총 화합물: 83 (중성 유기물)
- 총 포즈: 2867 (RMSD dedup 후)
