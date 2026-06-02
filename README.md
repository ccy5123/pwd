# Partition Coefficient Toolkit

SMILES → tissue/blood partition coefficient 예측 파이프라인.

## 모듈

- `shared/` — 공통 유틸 (SMILES 파서, descriptor, 단위 변환)
- `membrane_lipid/` — k_M (membrane lipid–water)
- `storage_lipid/` — P_SL (storage lipid–water)
- `albumin/` — albumin partition / binding
- `structural_protein/` — structural protein partition / binding
- `pipelines/` — 모듈 통합

## 시작하기

```bash
pip install -e ".[dev]"
pytest
ruff format .
ruff check .
```
