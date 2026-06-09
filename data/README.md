# data/

입력 데이터 보관 폴더.

작은 파일만 git에 commit. 큰 파일은 외부 클라우드 + 여기에 링크.

---

## UFZ_PPLFER_partition_data_combined_v2.xlsx

- **출처**:
  - UFZ-LSER database (Helmholtz Centre for Environmental Research – UFZ, Leipzig, Germany)
  - https://www.ufz.de/lserd
  - Maintainers: N. Ulrich, S. Endo, T. N. Brown, N. Watanabe, G. Bronner, M. H. Abraham, K.-U. Goss
  - 원본 논문 4개:
    - Endo, S.; Escher, B. I.; Goss, K.-U. *Environ. Sci. Technol.* **2011**, *45*, 5912–5921. DOI: 10.1021/es200855w (membrane lipid, 207 + 57 excluded)
    - Geisler, A.; Endo, S.; Goss, K.-U. *Environ. Sci. Technol.* **2012**, *46*, 9519–9524. DOI: 10.1021/es301921w (storage lipid, 250)
    - Endo, S.; Bauerfeind, J.; Goss, K.-U. *Environ. Sci. Technol.* **2012**, *46*, 12697–12703. DOI: 10.1021/es303379y (muscle protein, 67)
    - Endo, S.; Brown, T. N.; Goss, K.-U. *Environ. Sci. Technol.* **2013**, *47*, 6630–6639. DOI: 10.1021/es401772m (albumin + general tissue model)

- **내용**:
  - 네 phase에 대한 평형 분배계수: phospholipid–water (membrane), storage lipid–water, BSA/albumin–water, muscle protein–water
  - 440종 유기화합물 × 4 phase
  - 화합물 식별자: CAS, IUPAC 이름, canonical SMILES (RDKit)
  - PP-LFER descriptor (Abraham): E, S, A, B, V, L (250종)
  - 화학분류: PFAS / Synthetic Musk / Siloxane / Other

- **라이선스**:
  - **원본 UFZ-LSER database**: free for academic use after registration on UFZ website. 상용 사용은 별도 문의 필요.
  - **원본 논문 데이터**: American Chemical Society (ACS) 저작권. 학술 연구 목적 fair use 허용. 재배포·상용 시 ACS 허가 필요.
  - **combined v2 워크북 (이 파일)**: 위 데이터를 통합·가공한 derivative work. 학술 연구용 내부 사용 목적. 외부 공개·재배포 시 원본 출처를 모두 명시해야 함.

- **파일 정보**:
  - 시트 수: 7개 (README, combined, membrane_lipid, storage_lipid, albumin_measured, albumin_literature, muscle_protein, phase_classification)
  - 화합물 수: 440종 (CAS 유일, canonical SMILES 유일)
  - phase 커버리지: membrane 207, storage 247, albumin 83, muscle chicken 63 / fish 67 / collagen 15
  - 검증: 원본 PDF와 cell-by-cell 0 불일치, 시트 간 cross-check 0 불일치
  - 버전: v2 (membrane lipid 통합 추가 + SMILES 생성 + 화학분류)
