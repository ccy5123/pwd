# Phospholipid membrane K_lipw 검증 — openCOSMO-RS (Approach 2: DMPC pseudo-solvent)

## 데이터셋 (만들어 둠)
- `_membrane207_list.py` — **MEMBRANE207**: Endo, Escher, Goss 2011 (ES&T 45, 5912)
  liposome(PC)-water logK_lipw, **207종**. CSV의 membrane_exp 컬럼에서 추출.
  품질: SMILES 'M'오타 0, 중복 0, 빈값 0, logK 범위 −1.24~7.86 (논문 −0.8~7.9 범위 내).
  논문 Approach 2/3 테스트셋(n=207)과 일치.

## ★ 핵심 — 이건 성능 테스트(실패 재현)이지 정답 방법 아님 ★
storage/albumin/muscle는 등방성 bulk라 pseudo-solvent가 맞지만,
**막은 이방성 이중층**이라 단일 용매로 모델링하면 부정확.
Endo 2011이 직접 입증:
- Approach 2 (DMPC를 균질 용매로): COSMOthermX RMSE 1.01, SPARC RMSE 1.07
- Approach 3 (COSMOmic, 깊이별 적분): RMSE 0.79 ← 정답 방법
지금 openCOSMO-RS/ORCA로 가능한 건 Approach 2뿐이라 이걸로 감.
→ 기대: RMSE~1.0, 화합물별 오차(PAH 과소예측 / OH 2개 과대예측) 재현.
   이게 나오면 (a) 파이프라인 정상, (b) "membrane엔 COSMOmic 필요"의 정량 근거.

## 설정 (storage/albumin/muscle와 다른 점)
- **용매 = DMPC 단일 분자** (1,2-dimyristoyl-PC, 14:0/14:0, zwitterion 중성).
  근거: Endo 2011 Approach 2 — POPC/DMPC/DPPC 시험, PC 종류 둔감 → DMPC 채택.
  RDKit 검증: C36H72NO8P, MW 677.9, 순전하 0.
- **온도 = 298.15 K (25°C).** 근거: Endo 2011 "All calculations performed for 25 C".
  ※ storage/albumin/muscle(37°C)과 다름 — 논문 기준 따름.

## 파일
- `cosmo_kpw_openrs_membrane_v4.py` — K_lipw 계산·검증 (용매 DMPC 단일)
- `calc_sigma_membrane_v4.py` — σ-profile 생성 (DMPC+water+207)
- `run_calc_membrane_v4.sh` — SLURM
- `_membrane207_list.py` — 검증 데이터

## 실행
```bash
sbatch run_calc_membrane_v4.sh
# 또는:
python calc_sigma_membrane_v4.py solvent   # DMPC + water (DMPC는 큰 분자라 느림)
python calc_sigma_membrane_v4.py solutes   # 207종
python cosmo_kpw_openrs_membrane_v4.py validate
```

## 용질 σ-profile 재사용
캐시 키가 SMILES MD5라, storage/albumin/muscle의 MOL_<md5>를 복사하면 겹치는 용질 재사용:
```bash
cp -r <other>/sigma_cache/MOL_*  cosmo_work_v4_membrane/sigma_cache/
```
(DMPC는 membrane 전용이라 새로 계산. 단 온도 25C라 add_job에서 처리되므로
 σ-profile 자체는 온도무관 → 다른 상 용질 캐시 그대로 재사용 가능)

## 한계 / 다음
- Approach 2는 원리적 한계(이방성 무시). 정답은 COSMOmic.
- 다음: openCOSMO-RS의 COSMOmic 지원여부 조사 / TURBOMOLE BP-TZVP+COSMOmic 트랙.
- 막상태: 논문 데이터는 Tc 위 액정상(DMPC≥28°C 등)만. pseudo-solvent엔 막상태 개념 없음.
- openCOSMO-RS 24a neutral 전용: 이온화 화합물 중성종 기준. PFAS 음이온 범위 밖.
