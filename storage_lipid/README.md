# Storage lipid-water K_pw 검증 — openCOSMO-RS (muscle/albumin v4 포팅)

## 데이터셋 (만들어 둠)
- `_storage247_list.py` — **STORAGE247**: Geisler, Endo, Goss 2012 (ES&T 46, 9519)
  storage lipid-water logK, **247종, 전부 37°C**. CSV의 storage_exp 컬럼에서 추출.
  품질: SMILES 'M'오타 0, 중복 0, 빈값 0, logK 범위 −2.66~9.88 (논문과 일치).

## 핵심 차이 (muscle/albumin 대비)
**용매 = NTG(trinonanoylglycerol) 단일 분자.** 단백질처럼 잔기 혼합·조성이 없음.
- 근거: Geisler 2015 — COSMOtherm 계산에 NTG/triolein 사용, 둘 거의 동일(평균차 0.25 log),
  NTG가 작아 계산 빠름. Geisler 2012 — 지방산 조성 무관, 모든 storage lipid 단일상.
- 코드: `lng_inf(solute, [NTG], x=[1.0])` — 순수 NTG에 무한희석.
- σ-profile 만들 분자 = NTG + water 단 2개 (+ 용질 247).

## 파일
- `cosmo_kpw_openrs_storage_v4.py` — K 계산·검증 (용매 NTG 단일)
- `calc_sigma_storage_v4.py` — σ-profile 생성 (NTG+water+247)
- `run_calc_storage_v4.sh` — SLURM
- `_storage247_list.py` — 검증 데이터

## 실행
```bash
# XTB/ORCA 경로 확인 (calc_sigma_storage_v4.py 상단)
sbatch run_calc_storage_v4.sh
# 또는:
python calc_sigma_storage_v4.py solvent   # NTG + water
python calc_sigma_storage_v4.py solutes   # 247종
python cosmo_kpw_openrs_storage_v4.py validate
```

## 용질 σ-profile 재사용 (muscle/albumin과 겹치는 것)
캐시 키가 SMILES MD5라, muscle/albumin에서 만든 MOL_<md5>를 이 폴더로 복사하면
겹치는 용질은 재계산 0:
```bash
cp -r <muscle/albumin>/sigma_cache/MOL_*  cosmo_work_v4_storage/sigma_cache/
```
(NTG·water는 storage 전용이라 새로 계산)

## 성능 기대 & 해석 (Geisler 2015 COSMOtherm 기준)
- 단순 화합물 rmse~0.45, **H결합 donor rmse~0.35**(octanol 대용보다 훨씬 좋음), 복잡~0.71.
- storage lipid는 비구조적(유기용매類) 상 → **4상 중 COSMO가 가장 잘 맞는 상.**
  albumin(pocket 한계)·membrane(이방성) 같은 구조적 문제 없음.
- 고불소화(PFAS 중성종)에도 storage pp-LFER가 olive oil식보다 우수 → 유리.

## 한계
- openCOSMO-RS 24a neutral 전용. 이온화 화합물은 중성종 기준(논문도 동일).
- PFAS 음이온은 범위 밖.
- (선택) 복잡·다관능 화합물 외삽 테스트: Geisler 2015 Table 2의 24종(호르몬·농약·
  마이코톡신)을 SI에서 추출해 보조 테스트셋으로 추가하면 PFAS 대비 외삽 성능 평가 가능.
