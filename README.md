# Partition Coefficient Toolkit

SMILES → tissue/blood partition coefficient 예측 파이프라인.

## 모듈

- `shared/` — 공통 유틸 (SMILES 파서, descriptor, 단위 변환)
- `membrane_lipid/` — k_M (membrane lipid–water)
- `storage_lipid/` — P_SL (storage lipid–water)
- `albumin/` — albumin partition / binding
- `structural_protein/` — structural protein partition / binding
- `pipelines/` — 모듈 통합

## xTB 계산 (GFN1-xTB, GFN2-xTB)

### 스크립트
- `pipelines/xtb_gfn12_batch.py` — GFN1-xTB, GFN2-xTB로 분배 계수 계산

### 실행 방법
```bash
cd /home1/s9383/partition coefficient/pipelines/
python xtb_gfn12_batch.py --input ../data/compounds_input_table1.csv \
                           --output ../results/xtb_gfn12_predictions.csv \
                           --ncpu 8
```

### 결과 요약 (440 화합물)

| 상 (State) | 최고 R² | 방법 | 샘플 수 |
|------------|---------|------|--------|
| **Phospholipid (membrane)** | **0.81** | GFN2-GBSA-acetonitrile | 207 |
| **Storage lipid** | **0.84** | GFN1-ALPB-hexane | 247 |
| **Albumin** | **0.31** | GFN1-ALPB-woctanol | 83 |
| **Muscle (chicken)** | **-0.03** | GFN1-ALPB-octanol | 46 |
| **Muscle (fish)** | **-0.01** | GFN1-ALPB-octanol | 45 |

### 계산 방법
- **GFN1-xTB:** 1 vac + 25 ALPB + 12 GBSA = 38 SP/분자
- **GFN2-xTB:** 1 vac + 25 ALPB + 14 GBSA = 40 SP/분자
- 총 78 SP/분자 × 440 화합물 = 34,320 SP
- 실행 시간: 약 12시간 (8 workers)

### 데이터 파일
- 입력: `data/compounds_input_table1.csv` (440 화합물)
- 결과: `results/xtb_gfn12_predictions.csv` (159 컬럼)

---

## 시작하기

```bash
pip install -e ".[dev]"
pytest
ruff format .
ruff check .
```
