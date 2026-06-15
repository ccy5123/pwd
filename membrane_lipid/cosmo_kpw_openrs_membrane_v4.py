#!/usr/bin/env python3
"""
Phospholipid membrane(liposome)-water 분배계수 K_lipw 계산 (openCOSMO-RS_py) v4-membrane
=== Approach 2 (Endo, Escher, Goss 2011): 막을 균질 PC 용매로 모델링 ===

storage/albumin/muscle 버전에서 변경점:
  (1) 용매: DMPC 단일 분자 (pseudo-solvent). 순수 DMPC(x=1)에 무한희석.
      근거: Endo 2011 Approach 2 — POPC/DMPC/DPPC 시험, PC 종류에 둔감 → DMPC 채택.
  (2) 검증셋: MEMBRANE207 (Endo 2011 Approach 2/3 test set)
  (3) 온도: TEMP_K=298.15 (25C). 근거: Endo 2011 "All calculations were performed
      for 25 C." (Approach 1/2/3 모두). ※ storage/albumin/muscle(37C)과 다름 — 논문 기준.

★★ 중요 — 이건 "성능 테스트(실패 재현)"이지 정답 방법이 아님 ★★
  Endo 2011이 Approach 2(DMPC pseudo-solvent)는 RMSE~1.0으로 부정확하다고 입증:
    SPARC RMSE 1.07, COSMOthermX RMSE 1.01 (vs COSMOmic RMSE 0.79).
  막은 이방성 이중층이라 단일 용매로 깊이 의존성을 못 담음:
    - 비극성(특히 PAH): 막 코어에 위치 못해 과소예측
    - OH 2개 화합물: 균질상서 양쪽 OH가 PC와 동시 최적결합 → 과대예측
  → 기대: RMSE~1.0, 위 방향 오차 재현. 이게 나오면 (a)파이프라인 정상,
     (b)"membrane은 COSMOmic 필요"의 정량 근거가 됨.
  정답 방법은 COSMOmic(깊이별 적분) — openCOSMO-RS 지원여부/TURBOMOLE 트랙 별도.
"""
import os, math
import numpy as np
from pathlib import Path

from opencosmorspy.cosmors import COSMORS
from opencosmorspy.parameterization import openCOSMORS24a
from opencosmorspy.input_parsers import SigmaProfileParser

WORK = Path("./cosmo_work_v4_membrane")
CACHE = WORK / "sigma_cache"
TEMP_K = 298.15  # 25C — Endo 2011 예측계산 기준 (Approach 1/2/3 모두 25C)

# 용매 모델: DMPC = 1,2-dimyristoyl-sn-glycero-3-phosphocholine (14:0/14:0)
# 양쪽성 이온(zwitterion): 포스페이트 음이온 + 콜린 4급암모늄 양이온 → 전체 중성
DMPC_SMILES = "CCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCC"

MEMBRANE207 = [('1,2-Dichlorobenzene', 3.64, 'Clc1ccccc1Cl'), ('1,3-Dichlorobenzene', 3.71, 'Clc1cccc(Cl)c1'), ('1-Butanol', 0.51, 'CCCCO'), ('1-Hexanol', 1.88, 'CCCCCCO'), ('1-Pentanol', 1.08, 'CCCCCO'), ('1-Propanol', 0.17, 'CCCO'), ('2,2,4-Trimethylpentane', 4.61, 'CC(C)CC(C)(C)C'), ('2-Butoxyethanol', 0.6, 'CCCCOCCO'), ('2-Methyl-2-propanol', 0.16, 'CC(C)(C)O'), ('2-Propanol', -0.04, 'CC(C)O'), ('Acetone', 0.06, 'CC(C)=O'), ('Carbon tetrachloride', 2.61, 'ClC(Cl)(Cl)Cl'), ('Chlorobenzene', 2.91, 'Clc1ccccc1'), ('Cyclohexane', 3.27, 'C1CCCCC1'), ('Ethanol', -0.26, 'CCO'), ('Ethyl acetate', 0.46, 'CCOC(C)=O'), ('Heptane', 4.55, 'CCCCCCC'), ('Hexane', 3.91, 'CCCCCC'), ('Methanol', -0.53, 'CO'), ('p-Xylene', 2.98, 'Cc1ccc(C)cc1'), ('Propyl acetate', 1.01, 'CCCOC(C)=O'), ('Tetrachloroethene', 3.08, 'ClC(Cl)=C(Cl)Cl'), ('Trichloroethene', 2.43, 'ClC=C(Cl)Cl'), ('1-heptanol', 2.38, 'CCCCCCCO'), ('benzyl alcohol', 1.14, 'OCc1ccccc1'), ('cyclopentanone', 0.3, 'O=C1CCCC1'), ('di-n-butyl ether', 2.78, 'CCCCOCCCC'), ('N,N-dimethylaniline', 2.33, 'CN(C)c1ccccc1'), ('Octan-1-ol', 2.66, 'CCCCCCCCO'), ('1,2,4-Trichlorobenzene', 4.2, 'Clc1ccc(Cl)c(Cl)c1'), ('Dipentyl ether', 3.77, 'CCCCCOCCCCC'), ('2-octanone', 2.42, 'CCCCCCC(C)=O'), ('2-nonanone', 2.83, 'CCCCCCCC(C)=O'), ('Nitrobenzene', 2.01, 'O=[N+]([O-])c1ccccc1'), ('2-Nitrotoluene', 2.41, 'Cc1ccccc1[N+](=O)[O-]'), ('4-Chlorophenol', 2.73, 'Oc1ccc(Cl)cc1'), ('Octane', 4.67, 'CCCCCCCC'), ('3-Chlorophenol', 2.78, 'Oc1cccc(Cl)c1'), ('4-bromophenol', 2.4, 'Oc1ccc(Br)cc1'), ('4-n-Propylphenol', 2.92, 'CCCc1ccc(O)cc1'), ('4-iodophenol', 2.55, 'Oc1ccc(I)cc1'), ('Anthracene', 5.21, 'c1ccc2cc3ccccc3cc2c1'), ('Phenanthrene', 4.95, 'c1ccc2c(c1)ccc1ccccc12'), ('Fluoranthene', 5.58, 'c1ccc2c(c1)-c1cccc3cccc-2c13'), ('Pyrene', 5.71, 'c1cc2ccc3cccc4ccc(c1)c2c34'), ('tribromomethane', 2.33, 'BrC(Br)Br'), ('2-decanone', 3.16, 'CCCCCCCCC(C)=O'), ('1,4-dibromobenzene', 4.3, 'Brc1ccc(Br)cc1'), ('chrysene', 6.4, 'c1ccc2c(c1)ccc1c3ccccc3ccc21'), ('benzo[b]fluoranthene', 7.11, 'c1ccc2c(c1)-c1cccc3c1c-2cc1ccccc13'), ('benzo[ghi]perylene', 7.78, 'c1cc2ccc3ccc4ccc5cccc6c(c1)c2c3c4c56'), ('2-phenylphenol', 3.43, 'Oc1ccccc1-c1ccccc1'), ('4-fluorophenol', 2.19, 'Oc1ccc(F)cc1'), ('bisphenol A', 3.92, 'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1'), ('diazepam', 2.99, 'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21'), ('estrone', 3.59, 'CC12CCC3c4ccc(O)cc4CCC3C1CCC2=O'), ('1,2,4,5-tetrachlorobenzene', 4.73, 'Clc1cc(Cl)c(Cl)cc1Cl'), ('pentachlorobenzene', 5.18, 'Clc1cc(Cl)c(Cl)c(Cl)c1Cl'), ('benzo[a]pyrene', 7.19, 'c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34'), ('di-n-butyl phthalate', 3.87, 'CCCCOC(=O)c1ccccc1C(=O)OCCCC'), ('estradiol', 3.33, 'CC12CCC3c4ccc(O)cc4CCC3C1CCC2O'), ('4-nitrophenol', 2.72, 'O=[N+]([O-])c1ccc(O)cc1'), ('2,6-diethylphenol', 2.73, 'CCc1cccc(CC)c1O'), ('4-n-nonylphenol', 5.61, 'CCCCCCCCCc1ccc(O)cc1'), ('diethyl malonate', 0.5, 'CCOC(=O)CC(=O)OCC'), ('4-methylphenol', 2.35, 'Cc1ccc(O)cc1'), ('1,4-dichlorobenzene', 3.57, 'Clc1ccc(Cl)cc1'), ('ethylene glycol', -0.79, 'OCCO'), ('3-methylphenol', 2.34, 'Cc1cccc(O)c1'), ('1,3,5-trichlorobenzene', 4.16, 'Clc1cc(Cl)cc(Cl)c1'), ('cyclohexanol', 1.01, 'OC1CCCCC1'), ('cyclohexanone', 0.54, 'O=C1CCCCC1'), ('phenol', 1.96, 'Oc1ccccc1'), ('1,5-pentanediol', -0.7, 'OCCCCCO'), ('1,2-octanediol', 1.93, 'CCCCCCC(O)CO'), ('hexachlorobenzene', 5.64, 'Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl'), ('anthranilic acid', 2.08, 'Nc1ccccc1C(=O)O'), ('2-isopropyl-4,6-dinitrophenol', 3.14, 'CC(C)c1cc([N+](=O)[O-])cc([N+](=O)[O-])c1O'), ('tetrachlorocatechol', 4.41, 'Oc1c(O)c(Cl)c(Cl)c(Cl)c1Cl'), ('2,4-dichlorophenol', 3.57, 'Oc1ccc(Cl)cc1Cl'), ('4-ethylphenol', 2.78, 'CCc1ccc(O)cc1'), ('norfluoxetine', 4.49, 'NCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1'), ('4-phenylbutylamine', 2.41, 'NCCCCc1ccccc1'), ('lidocaine', 2.15, 'CCN(CC)CC(=O)Nc1c(C)cccc1C'), ('3,5-dibromo-4-methylphenol', 4.51, 'Cc1c(Br)cc(O)cc1Br'), ('4-tert-octylphenol', 5.39, 'CC(C)(C)CC(C)(C)c1ccc(O)cc1'), ('dibutyl succinate', 2.4, 'CCCCOC(=O)CCC(=O)OCCCC'), ('diethyl adipate', 1.66, 'CCOC(=O)CCCCC(=O)OCC'), ('2-tert-butyl-4,6-dinitrophenol', 4.1, 'CC(C)(C)c1cc([N+](=O)[O-])cc([N+](=O)[O-])c1O'), ('trans-1,2-cyclohexanediol', 0.23, 'O[C@@H]1CCCC[C@H]1O'), ('4-(methylsulfonyl)phenol', 1.27, 'CS(=O)(=O)c1ccc(O)cc1'), ('erythritol', -1.24, 'OC[C@@H](O)[C@@H](O)CO'), ('4-n-pentylphenol', 4.31, 'CCCCCc1ccc(O)cc1'), ('diclofenac', 4.45, 'O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl'), ('Ibuprofen', 3.8, 'CC(C)Cc1ccc(C(C)C(=O)O)cc1'), ('S-13b', 6.44, 'CC(C)(C)c1cc(Cl)c(O)c(C(=O)Nc2ccc(Cl)cc2[N+](=O)[O-])c1'), ('4-n-butylphenol', 3.13, 'CCCCc1ccc(O)cc1'), ('3,5-dibromo-4-hydroxybenzonitrile', 3.16, 'N#Cc1cc(Br)c(O)c(Br)c1'), ('2-allylphenol', 3.06, 'C=CCc1ccccc1O'), ('benzo[e]pyrene', 7.08, 'c1ccc2c(c1)c1cccc3ccc4cccc2c4c31'), ('indeno[1,2,3-cd]pyrene', 7.86, 'c1ccc2c(c1)-c1ccc3ccc4cccc5cc-2c1c3c45'), ('benzo[k]fluoranthene', 7.13, 'c1ccc2cc3c(cc2c1)-c1cccc2cccc-3c12'), ('dibenz[a,c]anthracene', 7.49, 'c1ccc2cc3c4ccccc4c4ccccc4c3cc2c1'), ('5-phenylvaleric acid', 3.06, 'O=C(O)CCCCc1ccccc1'), ('2-ethyl-1-pentanol', 1.81, 'CCCC(CC)CO'), ('3-hexyn-2,5-diol', -0.22, 'CC(O)C#CC(C)O'), ('2,3,5,6-tetrachloroaniline', 4.6, 'Nc1c(Cl)c(Cl)cc(Cl)c1Cl'), ("2,2',5,5'-tetrachlorobiphenyl (PCB 52)", 5.94, 'Clc1ccc(Cl)c(-c2cc(Cl)ccc2Cl)c1'), ("2,2',3,3',6,6'-hexachlorobiphenyl (PCB 136)", 6.5, 'Clc1ccc(Cl)c(-c2c(Cl)ccc(Cl)c2Cl)c1Cl'), ('4-tert-butyl-2,6-dinitrophenol', 3.81, 'CC(C)(C)c1cc([N+](=O)[O-])c(O)c([N+](=O)[O-])c1'), ('2-ethyl-4,6-dinitrophenol', 3.02, 'CCc1cc([N+](=O)[O-])cc([N+](=O)[O-])c1O'), ('methyl-4-chloro-2-nitrobenzoate', 2.45, 'COC(=O)c1ccc(Cl)cc1[N+](=O)[O-]'), ('trans-1,2-cyclooctanediol', 1.11, 'O[C@@H]1CCCCCC[C@H]1O'), ('2,3,4,5-tetrachlorophenol', 4.76, 'Oc1cc(Cl)c(Cl)c(Cl)c1Cl'), ('1-hepten-3-ol', 1.61, 'C=CC(O)CCCC'), ('estriol', 1.96, 'CC12CCC3c4ccc(O)cc4CCC3C1CC(O)C2O'), ('cycloheptanol', 1.51, 'OC1CCCCCC1'), ('cycloheptanone', 0.98, 'O=C1CCCCCC1'), ('cyclooctanone', 1.4, 'O=C1CCCCCCC1'), ('1,3-cyclohexanedione', 0.49, 'O=C1CCCC(=O)C1'), ('2,4-dinitrophenol', 2.67, 'O=[N+]([O-])c1ccc(O)c([N+](=O)[O-])c1'), ('propranolol', 3.32, 'CC(C)NCC(O)COc1cccc2ccccc12'), ('pentachloroaniline', 4.69, 'Nc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl'), ('3,4,5-trimethylphenol', 2.66, 'Cc1cc(O)cc(C)c1C'), ('dibenz[a,h]anthracene', 7.72, 'c1ccc2c(c1)ccc1cc3c(ccc4ccccc43)cc12'), ('2-methyl-4,6-dinitrophenol', 2.71, 'Cc1cc([N+](=O)[O-])cc([N+](=O)[O-])c1O'), ('1,2-pentanediol', 0.26, 'CCCC(O)CO'), ('dimethyl-2-amino-p-phthalate', 2.53, 'COC(=O)c1ccc(C(=O)OC)c(N)c1'), ('butyramide', -0.2, 'CCCC(N)=O'), ('fluoxetine', 4.35, 'CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1'), ('3-nitrophenol', 2.56, 'O=[N+]([O-])c1cccc(O)c1'), ('diethylstilbestrol', 4.76, 'CC/C(=C(/CC)c1ccc(O)cc1)c1ccc(O)cc1'), ('benzo[a]anthracene', 6.44, 'c1ccc2cc3c(ccc4ccccc43)cc2c1'), ('glycerol', -1.04, 'OCC(O)CO'), ("2,2',4,6,6'-pentachlorobiphenyl (PCB 104)", 6.13, 'Clc1cc(Cl)c(-c2c(Cl)cccc2Cl)c(Cl)c1'), ('17α-ethynylestradiol', 3.74, 'C#CC1(O)CCC2C3CCc4cc(O)ccc4C3CCC21C'), ('progesterone', 3.28, 'CC(=O)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C'), ('2,6-dinitrophenol', 2.03, 'O=[N+]([O-])c1cccc([N+](=O)[O-])c1O'), ('2,6-dimethylphenol', 2.47, 'Cc1cccc(C)c1O'), ('3,4-dinitrophenol', 3.17, 'O=[N+]([O-])c1ccc(O)cc1[N+](=O)[O-]'), ('2,3,4,6-tetrachlorophenol', 4.46, 'Oc1c(Cl)cc(Cl)c(Cl)c1Cl'), ('3-pentanol', 1.0, 'CCC(O)CC'), ('3-tert-butylphenol', 3.25, 'CC(C)(C)c1cccc(O)c1'), ('4-heptanol', 1.7, 'CCCC(O)CCC'), ('procaine', 2.38, 'CCN(CC)CCOC(=O)c1ccc(N)cc1'), ('4-chloro-3-methylphenol', 3.32, 'Cc1cc(O)ccc1Cl'), ('2,4-dimethyl-3-pentanol', 1.71, 'CC(C)C(O)C(C)C'), ("2,2',4,5',6-pentachlorobiphenyl (PCB 103)", 6.32, 'Clc1cc(Cl)c(-c2cc(Cl)ccc2Cl)c(Cl)c1'), ("2,2',3,4,4',5,6'-heptachlorobiphenyl (PCB 182)", 6.83, 'Clc1cc(Cl)c(-c2cc(Cl)c(Cl)c(Cl)c2Cl)c(Cl)c1'), ('3,4,5-trichlorophenol', 4.71, 'Oc1cc(Cl)c(Cl)c(Cl)c1'), ('4-methyl-2,6-dinitrophenol', 2.34, 'Cc1cc([N+](=O)[O-])c(O)c([N+](=O)[O-])c1'), ("4,4'-dihydroxybenzophenone", 2.7, 'O=C(c1ccc(O)cc1)c1ccc(O)cc1'), ('aniline', 1.63, 'Nc1ccccc1'), ('2,4-pentanediol', -0.3, 'CC(O)CC(C)O'), ("2,2',4,6-tetrachlorobiphenyl (PCB 50)", 5.92, 'Clc1cc(Cl)c(-c2ccccc2Cl)c(Cl)c1'), ('1,6-hexanediol', -0.1, 'OCCCCCCO'), ('1,7-hexanediol', 0.32, 'OCCCCCCCO'), ('1,8-octanediol', 0.81, 'OCCCCCCCCO'), ('1,2,3,5-tetrachlorobenzene', 4.77, 'Clc1cc(Cl)c(Cl)c(Cl)c1'), ('2,4,5-trichloroaniline', 4.16, 'Nc1cc(Cl)c(Cl)cc1Cl'), ('2-n-propylphenol', 3.13, 'CCCc1ccccc1O'), ('2,4,5-trichlorotoluene', 4.72, 'Cc1cc(Cl)c(Cl)cc1Cl'), ('salicylic acid', 2.55, 'O=C(O)c1ccccc1O'), ('1,2-hexanediol', 0.81, 'CCCCC(O)CO'), ('cyclooctanol', 2.09, 'OC1CCCCCCC1'), ("3,3',4,5-tetrachlorobiphenyl (PCB 78)", 6.53, 'Clc1cccc(-c2cc(Cl)c(Cl)c(Cl)c2)c1'), ('trimipramine', 4.74, 'CC(CN(C)C)CN1c2ccccc2CCc2ccccc21'), ('4-cyanophenol', 2.11, 'N#Cc1ccc(O)cc1'), ('4-tert-amylphenol', 3.54, 'CCC(C)(C)c1ccc(O)cc1'), ('warfarin', 3.42, 'CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O'), ('5-hexen-1-ol', 1.2, 'C=CCCCCO'), ('meso-hexestrol', 4.37, 'CCC(c1ccc(O)cc1)C(CC)c1ccc(O)cc1'), ('dienestrol', 5.08, 'CC=C(C(=CC)c1ccc(O)cc1)c1ccc(O)cc1'), ('diethylphthalate', 1.77, 'CCOC(=O)c1ccccc1C(=O)OCC'), ('butylbenzylphthalate', 4.3, 'CCCCc1ccc(C(=O)[O-])c(C(=O)[O-])c1Cc1ccccc1'), ('1,2,3-trichlorobenzene', 4.19, 'Clc1cccc(Cl)c1Cl'), ('2,6-dichlorophenol', 2.86, 'Oc1c(Cl)cccc1Cl'), ('pentachlorophenol', 5.1, 'Oc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl'), ('2,4,6-trimethylaniline', 2.38, 'Cc1cc(C)c(N)c(C)c1'), ('2,4,6-trichlorophenol', 3.8, 'Oc1c(Cl)cc(Cl)cc1Cl'), ('2-tert-butylphenol', 3.51, 'CC(C)(C)c1ccccc1O'), ('2-nitrophenol', 1.89, 'O=[N+]([O-])c1ccccc1O'), ('2-sec-butyl-4,6-dinitrophenol', 3.73, 'CCC(C)c1cc([N+](=O)[O-])cc([N+](=O)[O-])c1O'), ('amlodipine', 3.75, 'CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl'), ('2-sec-butylphenol', 3.47, 'CCC(C)c1ccccc1O'), ('2-ethylphenol', 2.81, 'CCc1ccccc1O'), ('quinoline', 1.67, 'c1ccc2ncccc2c1'), ('4-phenylphenol', 3.52, 'Oc1ccc(-c2ccccc2)cc1'), ('1,5-hexanediol', -0.22, 'CC(O)CCCCO'), ('5-hexyn-1-ol', 0.78, 'C#CCCCCO'), ('trans-2-hexen-1-ol', 1.36, 'CCC/C=C/CO'), ('cis-3-hexen-1-ol', 1.2, 'CC/C=C\\CCO'), ('trans-3-hexen-1-ol', 1.23, 'CC/C=C/CCO'), ('benzyl-4-hydroxybenzoate', 3.75, 'O=C(OCc1ccccc1)c1ccc(O)cc1'), ('tetracaine', 3.23, 'CCCCNc1ccc(C(=O)OCCN(C)C)cc1'), ('butyl-4-hydroxybenzoate', 3.4, 'CCCCOC(=O)c1ccc(O)cc1'), ('2-methylphenol', 2.45, 'Cc1ccccc1O'), ('2-chlorophenol', 2.76, 'Oc1ccccc1Cl'), ('3,4-dimethylaniline', 2.11, 'Cc1ccc(N)cc1C'), ('3,4-dichlorophenol', 3.76, 'Oc1ccc(Cl)c(Cl)c1'), ('2,4,5-trichlorophenol', 4.46, 'Oc1cc(Cl)c(Cl)cc1Cl'), ('cyclopentanol', 0.52, 'OC1CCCC1'), ('3-trifluoromethylphenol', 3.25, 'Oc1cccc(C(F)(F)F)c1'), ('4-tert-butylphenol', 3.48, 'CC(C)(C)c1ccc(O)cc1'), ('3-nitroaniline', 2.17, 'Nc1cccc([N+](=O)[O-])c1'), ('4-sec-butylphenol', 3.26, 'CCC(C)c1ccc(O)cc1'), ('4-isopropylphenol', 3.25, 'CC(C)c1ccc(O)cc1')]


def lng_inf(crs, solute_path, solvent_paths, solvent_x, T):
    """solute의 무한희석 lnγ∞ (muscle 원본과 동일)"""
    crs.clear_molecules()
    crs.add_molecule([solute_path])
    for p in solvent_paths:
        crs.add_molecule([p])
    crs.clear_jobs()
    x = np.array([0.0] + list(solvent_x))
    crs.add_job(x, T, refst='pure_component')
    result = crs.calculate()
    return result['tot']['lng'][0][0]


def calculate_partition_coefficient(solute_orcacosmo, dmpc_orcacosmo, water_orcacosmo, T=TEMP_K, c_unit=0.0):
    """log K_lipw 계산 (Approach 2). 용매=순수 DMPC(x=1), 물=순수(x=1)."""
    try:
        crs = COSMORS(par=openCOSMORS24a())
        lng_w = lng_inf(crs, solute_orcacosmo, [water_orcacosmo], [1.0], T)   # 물에서
        lng_m = lng_inf(crs, solute_orcacosmo, [dmpc_orcacosmo],  [1.0], T)   # DMPC(막)에서
        ln_Kx = lng_w - lng_m
        log_K = ln_Kx / math.log(10) + c_unit
        return float(log_K)
    except Exception as e:
        import traceback; print(f"  [COSMO-RS error: {e}]"); traceback.print_exc()
        return None


def find_orcacosmo_file(cache_key_str):
    for mol_dir in CACHE.glob("MOL_*"):
        if cache_key_str in mol_dir.name:
            oc = mol_dir / f"{mol_dir.name}.orcacosmo"
            if oc.exists():
                return str(oc)
    return None


def validate_membrane207(label="openrs_membrane_v4"):
    print("=== openCOSMO-RS Membrane 검증 v4 (25C, Approach2: DMPC pseudo-solvent) ===")
    print("    ※ Endo 2011 기준 RMSE~1.0 예상(성능테스트). COSMOmic 아님.")
    water = str(CACHE / "WATER.orcacosmo")
    dmpc  = str(CACHE / "DMPC.orcacosmo")
    for f,nm in [(water,"WATER"),(dmpc,"DMPC")]:
        if not Path(f).exists():
            print(f"[ERROR] {nm}.orcacosmo 없음 → calc_sigma_membrane_v4.py 먼저 실행"); return None
    print(f"Solvent: DMPC={dmpc}\n         WATER={water}")

    import hashlib
    exp_vals, pred_vals, names, rows = [], [], [], []
    for name, exp_k, smiles in MEMBRANE207:
        key = hashlib.md5(smiles.encode()).hexdigest()[:16]
        oc = find_orcacosmo_file(key)
        if oc is None:
            print(f"[skip] {name}: .orcacosmo 없음"); continue
        log_k = calculate_partition_coefficient(oc, dmpc, water, TEMP_K, c_unit=0.0)
        if log_k is not None:
            exp_vals.append(exp_k); pred_vals.append(log_k); names.append(name)
            rows.append({"name":name,"exp":float(exp_k),"pred":float(log_k)})
            print(f"[*] {name}: exp={exp_k:.2f}, pred={log_k:.3f}")

    if len(exp_vals) >= 2:
        e=np.array(exp_vals); p=np.array(pred_vals)
        r2=np.corrcoef(p,e)[0,1]**2; rmse=float(np.sqrt(np.mean((p-e)**2)))
        bias=float(np.mean(p-e))
        import json
        out={"label":label,"n":len(e),"r2":float(r2),"rmse":rmse,"bias":bias,"rows":rows}
        (WORK/f"results_{label}.json").write_text(json.dumps(out,indent=2))
        print(f"\n=== 결과 ===\nn={len(e)}\nR2={r2:.3f}\nRMSE={rmse:.3f} (bias={bias:+.3f})")
        print(f"제안 c_unit = {-bias:+.3f}")
        print("해석: RMSE~1.0이면 논문 Approach2 재현(성공). 화합물별 오차방향")
        print("      (PAH 과소/OH2개 과대)도 확인 → COSMOmic 필요성 근거.")
        return out
    print("데이터 부족"); return None


if __name__ == "__main__":
    import sys
    if len(sys.argv)>1 and sys.argv[1]=="validate":
        validate_membrane207()
    else:
        print("Usage: python cosmo_kpw_openrs_membrane_v4.py validate")
