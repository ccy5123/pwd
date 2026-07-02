# Muscle Protein K_pw v5 - CSV 리팩터 & 캐시 재생성 Handoff

> **마지막 업데이트**: 2026-07-02
> **현재 상태**: Phase 0 실행 대기
> **진행 상황**: 읽기전용 검증 완료, Phase 0 보고 완료

---

## 작업 개요

muscle protein-water 분배계수(log K_pw) openCOSMO-RS 검증 파이프라인을 CSV 단일 소스로 리팩터하고, 오염/누락 캐시 4-6종을 재생성하는 작업.

**목표**: chicken 46종 전체 검증 (현재 n=43 → 46)

---

## 문제 정리

### 현재 이슈 (서버 실측)

| 문제 | 원인 | 영향 |
|------|------|------|
| **124-TMB, 2-nitrotoluene, 24-DNT** | SMILES 메틸기 `M` 오타 → RDKit 파싱 불가 | 캐시 미생성, 3종 miss |
| **benzo[a]pyrene (BaP)** | 틀린 구조 캐시 존재 (InChIKey QXXOVHAZDCQCPN ≠ 진짜 FMMWHPNWAFZXNH) | 잘못된 구조로 계산됨 |
| **phenanthrene, chrysene** | 캐시는 있으나 옛 SMILES 문자열 해시로 명명됨 | CSV canonical SMILES와 md5 불일치 가능성 |

### 근본 원인

- SMILES를 스크립트에 하드코딩 → CSV와 이중 관리
- `md5(raw SMILES)` 사용 → canonical/raw 드리프트

---

## 서버 상태 (Fact)

### 파일 구조

| 항목 | 경로/값 |
|------|---------|
| **작업 디렉토리** | `/gpfs/home1/s9383/Muscle_protein` |
| **WORK** | `./cosmo_work_v3` |
| **CACHE** | `./cosmo_work_v3/sigma_cache` |
| **입력 CSV** | `../3D/compounds_input_table1.csv` |
| **진단 번들** | `/gpfs/home1/s9383/Muscle_protein/diag_bundle.tgz` (4.9MB) |

### 캐시 현황 (cosmo_work_v3/sigma_cache)

| 항목 | 개수/상태 |
|------|----------|
| `AA_*_*` 디렉토리 | 20개 ✓ |
| `WATER.orcacosmo` | 존재 ✓ |
| `MOL_*` 디렉토리 | 49개 |
| **최근 검증 결과** | n=43, 실패 3종 (TMB, nitrotoluene, DNT) |

### CSV 현황

| 항목 | 값 |
|------|-----|
| 경로 | `../3D/compounds_input_table1.csv` |
| 총 행 | 440 |
| muscle_chicken_exp notna | **46** ✓ |

---

## 작업 규칙 (엄수)

```bash
# 모든 bash는 새 셸이므로 매 호출을 이렇게 묶을 것
cd /gpfs/home1/s9383/Muscle_protein && <명령>
```

| 제약 | 설명 |
|------|------|
| 기존 결과 덮어쓰기 금지 | `*_corrected.py`, v4 JSON 보존 |
| 라벨 유지 | `chicken_core6_seq_v5` |
| c_unit 고정 | `-0.80` (데이터 피팅값으로 변경 금지) |
| 조성/온도/검증셋 변경 금지 | core-6 서열, 310.15K, chicken 46종 |
| 재계산은 sbatch만 | 포그라운드 금지, 48h 제한 |
| 커밋/푸시 | repo `tlagustn123-source`, CLAUDE.md 규약, 승인 후 |
| 보고 구분 | [Fact/서버확인] vs [추론] 분리 |

---

## Phase별 진행 계획

### ✅ Phase 0 — 읽기전용 검증 (완료)

| 항목 | 확인 사항 | 결과 |
|------|----------|------|
| WORK/CACHE 경로 | `cosmo_kpw_muscle_v5.py` | ✓ v3 사용 |
| CSV 확인 | `3D/compounds_input_table1.csv` | ✓ 46종 |
| 캐시 개수 | AA=20, WATER=1, MOL=49 | ✓ |
| **4종 올바른 SMILES 캐시** | canonical md5로 존재 여부 | 보고 완료 |

**올바른 SMILES** (CSV 기준):
```
124-TMB:           Cc1ccc(C)c(C)c1
2-nitrotoluene:    Cc1ccccc1[N+](=O)[O-]
24-DNT:            Cc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]
benzo[a]pyrene:    c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34
```

### Phase 1 — CSV 리팩터 (대기)

**새 파일**: `cosmo_kpw_muscle_v5_csv.py`

- 하드코딩 ENDO46 삭제 → CSV에서 muscle_chicken_exp notna 46종 로드
- 캐시 키: `cosmo_kpw_common_v3.py`의 `cache_key()` import (= md5(canonical)[:16])
- WORK=./cosmo_work_v3, CACHE=WORK/sigma_cache
- c_unit=-0.80 고정
- JSON/콘솔 출력: r2_pearson, r2_cod, rmse, within_10x_pct, c_unit, c_unit_basis

### Phase 2 — 캐시 드라이런 (대기)

- CSV 46종 canonical md5 → MOL_* 존재 + .orcacosmo 완전성 검사
- n_hit/n_miss + miss 목록·사유 보고

**예상 miss**: {124-TMB, 2-nitrotoluene, 24-DNT, BaP, phenanthrene, chrysene}

### Phase 2.5 — 오염/누락 재생성 (대기)

**생성기**: `gen_cache_from_csv.py` (신규) 또는 `calc_missing_molecules_v4.py` 수정본

- 대상 6종의 CSV SMILES 사용
- 파이프라인: RDKit ETKDGv3+MMFF → xTB(--opt tight --gfn 2) → ORCA BP86/def2-TZVPD CPCM → .orcacosmo
- 디렉토리명 = `cache_key()`(canonical md5) 통일
- 기존 오염 디렉토리 → `sigma_cache/_deprecated/` 이동
- sbatch (48h)

### Phase 3 — 검증 (대기)

- `run_cosmo_muscle_v5_csv.sh` 수정/생성
- sbatch (48h)
- 결과 n==46 확인

### Phase 4 — 비교 보고 (대기)

| metric | v5 현재(n=43, BaP오염) | v5 수정(n=46) |
|---|---|:---|
| R²_cod | 0.646 | ? |
| Pearson r² | 0.650 | ? |
| RMSE | 0.547 | ? |
| within10× | 93.0% | ? |

---

## 다음 작업 (새 대화방)

### 실행할 명령어 복사

```bash
# Phase 0 재검증 필요시
cd /gpfs/home1/s9383/Muscle_protein && grep -E "WORK|CACHE|ENDO46|md5" cosmo_kpw_muscle_v5.py | head -20

# 4종 올바른 SMILES 캐시 존재 확인 (Phase 0 핵심)
cd /gpfs/home1/s9383/Muscle_protein && python3 -c "
import hashlib

smiles_list = [
    ('124-TMB', 'Cc1ccc(C)c(C)c1'),
    ('2-nitrotoluene', 'Cc1ccccc1[N+](=O)[O-]'),
    ('24-DNT', 'Cc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]'),
    ('benzo[a]pyrene', 'c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34')
]

for name, smiles in smiles_list:
    key = hashlib.md5(smiles.encode()).hexdigest()[:16]
    mol_dir = f'cosmo_work_v3/sigma_cache/MOL_{key}'
    print(f'{name}: MOL_{key}')
    
    import os
    if os.path.isdir(mol_dir):
        orcacosmo = f'{mol_dir}/{os.path.basename(mol_dir)}.orcacosmo'
        if os.path.exists(orcacosmo):
            print(f'  ✓ 캐시 존재')
        else:
            print(f'  ✗ .orcacosmo 없음')
    else:
        print(f'  ✗ 디렉토리 없음')
"
```

### 첫 프롬프트 (새 대화방 시작 시 붙여넣기)

```
# 작업: muscle K_pw v5 — CSV 단일소스 리팩터 + 오염 캐시 4~6종 재생성

/home1/s9383/Muscle_protein/HANDOFF.md 를 읽고 작업을 계속해주세요.
현재 Phase 0 완료 상태이고, Phase 1부터 진행해야 합니다.
모든 bash는 cd /gpfs/home1/s9383/Muscle_protein && ... 로 묶어주세요.
각 Phase 끝에서 멈추고 보고해주세요.
```

---

## 참고 파일

| 파일 | 경로 |
|------|------|
| HANDOFF 본문 | `/home1/s9383/Muscle_protein/HANDOFF.md` |
| 진단 번들 | `/gpfs/home1/s9383/Muscle_protein/diag_bundle.tgz` |
| 입력 CSV | `/gpfs/home1/s9383/3D/compounds_input_table1.csv` |
| 현재 스크립트 | `/home1/s9383/Muscle_protein/cosmo_kpw_muscle_v5.py` |
| 캐시 함수 참고 | `/home1/s9383/Muscle_protein/cosmo_kpw_common_v3.py` |

---

*이 문서는 새 대화방에서 작업을 계속할 때 사용합니다.*
