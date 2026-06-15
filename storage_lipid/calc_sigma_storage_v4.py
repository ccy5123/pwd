#!/usr/bin/env python3
"""
Storage lipid σ-profile 계산 v4
muscle calc_missing_molecules_v4.py와 동일 워크플로:
  RDKit ETKDGv3+MMFF -> xTB(GFN2) opt -> ORCA BP86/def2-TZVPD CPCM -> .orcacosmo
대상: STORAGE247 solute + NTG(용매) + water
캐시 규약: solute=MOL_<md5>, NTG=NTG.orcacosmo(고정명), water=WATER.orcacosmo

실행:
  python calc_sigma_storage_v4.py solvent   # NTG + water (2개)
  python calc_sigma_storage_v4.py solutes   # 247종
  python calc_sigma_storage_v4.py all
"""
import subprocess, hashlib, re, sys, shutil
from pathlib import Path

WORK = Path("./cosmo_work_v4_storage")
CACHE = WORK / "sigma_cache"
CACHE.mkdir(parents=True, exist_ok=True)

# ★ 서버 실제 경로로 수정 (muscle 원본 경로 그대로) ★
XTB  = "/home1/s9383/Muscle_protein/xtb-6.5.1/bin/xtb"
ORCA = "/home1/s9383/Muscle_protein/orca_6_1_1_linux_x86-64_shared_openmpi418/orca"

NTG_SMILES = "CCCCCCCCC(=O)OCC(COC(=O)CCCCCCCC)OC(=O)CCCCCCCC"  # trinonanoylglycerol

STORAGE247 = [('1,1,1,2-Tetrachloroethane', 2.85, 'ClCC(Cl)(Cl)Cl'), ('1,1,1-Trichloroethane', 2.66, 'CC(Cl)(Cl)Cl'), ('1,1,2,2-Tetrachloroethane', 2.46, 'ClC(Cl)C(Cl)Cl'), ('1,1,2-Trichloroethane', 2.09, 'ClCC(Cl)Cl'), ('1,1-Dichloro-1-fluoroethane', 1.34, 'CC(F)(Cl)Cl'), ('1,1-Dichloroethane', 1.74, 'CC(Cl)Cl'), ('1,1-Dichloroethene', 2.19, 'C=C(Cl)Cl'), ('1,2,3-Trimethylbenzene', 3.32, 'Cc1cccc(C)c1C'), ('1,2,4-Trifluorobenzene', 2.67, 'Fc1ccc(F)c(F)c1'), ('1,2,4-Trimethylbenzene', 3.43, 'Cc1ccc(C)c(C)c1'), ('1,2-Dibromoethane', 1.67, 'BrCCBr'), ('1,2-Dichlorobenzene', 3.71, 'Clc1ccccc1Cl'), ('1,2-Dichloroethane', 1.6, 'ClCCCl'), ('1,2-Dichloropropane', 1.99, 'CC(Cl)CCl'), ('1,2-Difluorobenzene', 2.62, 'Fc1ccccc1F'), ('1,3,5-Trifluorobenzene', 2.93, 'Fc1cc(F)cc(F)c1'), ('1,3,5-Trimethylbenzene', 3.49, 'Cc1cc(C)cc(C)c1'), ('1,3-Dichlorobenzene', 3.84, 'Clc1cccc(Cl)c1'), ('1,4-Difluorobenzene', 2.58, 'Fc1ccc(F)cc1'), ('1-Bromo-2-chloroethane', 2.07, 'ClCCBr'), ('1-Butanol', 0.02, 'CCCCO'), ('1-Chlorobutane', 2.78, 'CCCCCl'), ('1-Chloropentane', 3.48, 'CCCCCCl'), ('1-Chloropropane', 2.22, 'CCCCl'), ('1-Hexanol', 1.36, 'CCCCCCO'), ('1-Methoxy-2-propanol', -1.53, 'COCC(C)O'), ('1-Nitropropane', 0.97, 'CCC[N+](=O)[O-]'), ('1-Pentanol', 0.54, 'CCCCCO'), ('1-Propanol', -0.48, 'CCCO'), ('2,2,4-Trimethylpentane', 4.64, 'CC(C)CC(C)(C)C'), ('2,2-Dichloro-1,1,1-\ntrifluoroethane', 1.81, 'FC(F)(F)C(Cl)Cl'), ('2,2-Dimethylbutane', 3.79, 'CCC(C)(C)C'), ('2,3,4-Trimethylpentane', 5.05, 'CC(C)C(C)C(C)C'), ('2-Butoxyethanol', -0.21, 'CCCCOCCO'), ('2-Chloropropane', 2.03, 'CC(C)Cl'), ('2-Ethoxyethanol', -1.27, 'CCOCCO'), ('2-Fluoropropane', 1.38, 'CC(C)F'), ('2-Heptanone', 1.81, 'CCCCCC(C)=O'), ('2-Hexanone', 1.19, 'CCCCC(C)=O'), ('2-Isopropoxyethanol', -1.16, 'CC(C)OCCO'), ('2-Methoxyethanol', -1.67, 'COCCO'), ('2-Methyl-1-propanol', -0.06, 'CC(C)CO'), ('2-Methyl-2-propanol', -0.73, 'CC(C)(C)O'), ('2-Methylpentane', 3.99, 'CCCC(C)C'), ('2-Nitropropane', 0.61, 'CC(C)[N+](=O)[O-]'), ('2-Pentanone', 0.65, 'CCCC(C)=O'), ('2-Propanol', -0.82, 'CC(C)O'), ('3-Methyl-1-butanol', 0.41, 'CC(C)CCO'), ('3-Methylhexane', 4.61, 'CCCC(C)CC'), ('3-Methylpentane', 4.07, 'CCC(C)CC'), ('3-Pentanone', 0.56, 'CCC(=O)CC'), ('4-Methyl-2-pentanone', 1.05, 'CC(=O)CC(C)C'), ('Acetone', -0.32, 'CC(C)=O'), ('Allylbenzene', 2.96, 'C=CCc1ccccc1'), ('Benzene', 2.12, 'c1ccccc1'), ('Bromochloromethane', 1.49, 'ClCBr'), ('2-Butanone', 0.18, 'CCC(C)=O'), ('Butane', 3.03, 'CCCC'), ('Butyl acetate', 1.59, 'CCCCOC(C)=O'), ('Carbon tetrachloride', 3.18, 'ClC(Cl)(Cl)Cl'), ('1,2-\nDichlorotetrafluoroethane', 3.0, 'FC(F)(Cl)C(F)(F)Cl'), ('1,1,2,2,3,3,4,4-\nOctafluorobutane', 2.5, 'FC(F)C(F)(F)C(F)(F)C(F)F'), ('1,1,2,2,3,3-\nHexafluoropropane', 1.49, 'FC(F)C(F)(F)C(F)F'), ('1,1,2,2-Tetrafluoroethane', 0.92, 'FC(F)C(F)F'), ('1,1-Difluoroethane', 0.91, 'CC(F)F'), ('Halothane', 2.28, 'FC(F)(F)C(Cl)Br'), ('Teflurane', 1.8, 'FC(F)C(F)(F)Br'), ('1,1,1,2-Tetrafluoroethane', 1.03, 'FCC(F)(F)F'), ('Fluroxene', 1.26, 'C=COCC(F)(F)F'), ('Carbon tetrafluoride', 1.12, 'FC(F)(F)F'), ('1,3-Difluoropropane', 0.66, 'FCCCF'), ('Enflurane', 2.03, 'FC(F)OC(F)(F)C(F)Cl'), ('Isoflurane', 1.99, 'FC(F)OC(Cl)C(F)(F)F'), ('Desflurane', 1.88, 'FC(F)OC(F)C(F)(F)F'), ('Sevoflurane', 1.87, 'FCOC(C(F)(F)F)C(F)(F)F'), ('Chlorobenzene', 2.73, 'Clc1ccccc1'), ('Chlorodibromomethane', 2.12, 'ClC(Br)Br'), ('Chloroethane', 1.58, 'CCCl'), ('Chloroform', 2.02, 'ClC(Cl)Cl'), ('cis-1,2-Dichloroethene', 1.59, 'Cl/C=C\\Cl'), ('Cycloheptane', 4.82, 'C1CCCCCC1'), ('Cyclohexane', 4.05, 'C1CCCCC1'), ('Cyclopentane', 3.5, 'C1CCCC1'), ('Cyclopropane', 1.73, 'C1CC1'), ('Decane', 6.39, 'CCCCCCCCCC'), ('Dibromethane', 1.76, 'BrCBr'), ('Dichloromethane', 1.28, 'ClCCl'), ('Diethyl ether', 0.79, 'CCOCC'), ('Difluoromethane', 0.53, 'FCF'), ('Divinyl ether', 1.69, 'C=COC=C'), ('Ethane', 1.73, 'CC'), ('Ethanol', -1.14, 'CCO'), ('Ethene', 0.84, 'C=C'), ('Ethyl acetate', 0.46, 'CCOC(C)=O'), ('Ethyl tert-butyl ether', 1.37, 'CCOC(C)(C)C'), ('Ethyl tert-pentyl ether', 1.93, 'CCOC(C)(C)CC'), ('Ethylbenzene', 3.14, 'CCc1ccccc1'), ('Fluorobenzene', 2.41, 'Fc1ccccc1'), ('Fluoroethane', 0.62, 'CCF'), ('Fluorochloromethane', 0.77, 'FCCl'), ('Heptane', 4.71, 'CCCCCCC'), ('Hexafluorobenzene', 2.45, 'Fc1c(F)c(F)c(F)c(F)c1F'), ('Hexane', 4.11, 'CCCCCC'), ('Isobutyl acetate', 1.6, 'CC(=O)OCC(C)C'), ('Isopentyl acetate', 2.11, 'CC(=O)OCCC(C)C'), ('Isopropyl acetate', 1.02, 'CC(=O)OC(C)C'), ('Isopropylbenzene', 3.52, 'CC(C)c1ccccc1'), ('Methoxyflurane', 2.11, 'COC(F)(F)C(Cl)Cl'), ('Methane', 0.75, 'C'), ('Methanol', -1.95, 'CO'), ('Methyl acetate', -0.02, 'COC(C)=O'), ('Methyl chloride', 0.8, 'CCl'), ('Methylcyclopentane', 4.0, 'CC1CCCC1'), ('Methylpentafluorobenzene', 3.27, 'Cc1c(F)c(F)c(F)c(F)c1F'), ('m-Methylstyrene', 3.23, 'C=Cc1cccc(C)c1'), ('m-Xylene', 3.16, 'Cc1cccc(C)c1'), ('Nonane', 5.82, 'CCCCCCCCC'), ('o-Xylene', 3.12, 'Cc1ccccc1C'), ('p-Xylene', 3.16, 'Cc1ccc(C)cc1'), ('Pentachloroethane', 2.93, 'ClC(Cl)C(Cl)(Cl)Cl'), ('Pentafluorobenzene', 2.31, 'Fc1cc(F)c(F)c(F)c1F'), ('Pentane', 3.49, 'CCCCC'), ('Pentyl acetate', 2.11, 'CCCCCOC(C)=O'), ('p-Methylstyrene', 3.21, 'C=Cc1ccc(C)cc1'), ('Propane', 2.35, 'CCC'), ('Propyl acetate', 1.02, 'CCCOC(C)=O'), ('Propylbenzene', 3.65, 'CCCc1ccccc1'), ('Tetrachloroethene', 3.57, 'ClC(Cl)=C(Cl)Cl'), ('Toluene', 2.67, 'Cc1ccccc1'), ('trans-1,2-Dichloroethene', 2.06, 'Cl/C=C/Cl'), ('Trichloroethene', 2.8, 'ClC=C(Cl)Cl'), ('tridecane', 8.16, 'CCCCCCCCCCCCC'), ('triethylamine', 1.04, 'CCN(CC)CC'), ('undecane', 7.03, 'CCCCCCCCCCC'), ('alpha-pinene', 4.58, 'CC1=CCC2CC1C2(C)C'), ('1,2-dimethoxyethane', -0.33, 'COCCOC'), ('1,4-dioxane', -0.24, 'C1COCCO1'), ('1-heptanol', 1.84, 'CCCCCCCO'), ('2-methylpyridine', 0.61, 'Cc1ccccn1'), ('3-methylpyridine', 0.84, 'Cc1cccnc1'), ('4-methylpyridine', 0.82, 'Cc1ccncc1'), ('benzyl alcohol', 0.12, 'OCc1ccccc1'), ('bromobenzene', 3.27, 'Brc1ccccc1'), ('butyl formate', 1.46, 'CCCCOC=O'), ('butyl propanoate', 2.35, 'CCCCOC(=O)CC'), ('butylbenzene', 4.13, 'CCCCc1ccccc1'), ('1,1-difluoro-2-chloroethene', 1.63, 'FC(F)=CCl'), ('1-chloro-2,2,2-\ntrifluoroethane', 1.34, 'FC(F)(F)CCl'), ('bis-(2,2,2-\ntrifluoroethyl)ether', 1.92, 'FC(F)(F)COCC(F)(F)F'), ('cyclohexene', 3.33, 'C1=CCCCC1'), ('cyclopentanone', 0.39, 'O=C1CCCC1'), ('difluorochloromethane', 0.72, 'FC(F)Cl'), ('diisopropyl ether', 1.29, 'CC(C)OC(C)C'), ('dimethoxymethane', 0.18, 'COCOC'), ('dimethyl ether', 0.13, 'COC'), ('dimethylacetamide', -1.25, 'CC(=O)N(C)C'), ('dimethylformamide', -1.57, 'CN(C)C=O'), ('di-n-butyl ether', 3.14, 'CCCCOCCCC'), ('dodecane', 7.59, 'CCCCCCCCCCCC'), ('ethyl formate', 0.16, 'CCOC=O'), ('ethyl propanoate', 1.13, 'CCOC(=O)CC'), ('iodoethane', 1.85, 'CCI'), ('fluorotrichloromethane', 2.28, 'FC(Cl)(Cl)Cl'), ('methyl formate', -0.48, 'COC=O'), ('methylcyclohexane', 4.54, 'CC1CCCCC1'), ('nitroethane', 0.45, 'CC[N+](=O)[O-]'), ('nitromethane', -0.02, 'C[N+](=O)[O-]'), ('N,N-dimethylaniline', 2.29, 'CN(C)c1ccccc1'), ('pentadecane', 9.31, 'CCCCCCCCCCCCCCC'), ('piperidine', 1.01, 'C1CCNCC1'), ('propyl bromide', 2.47, 'CCCBr'), ('propyl formate', 0.87, 'CCCOC=O'), ('pyridine', 0.15, 'c1ccncc1'), ('tetradecane', 8.74, 'CCCCCCCCCCCCCC'), ('tetrahydrofuran', 0.62, 'C1CCOC1'), ('halopropane', 2.91, 'FC(F)(F)C(F)(Cl)Br'), ('CF3CHFOCHFCl', 2.15, 'FC(Cl)OC(F)C(F)(F)F'), ('Cyclooctane', 5.34, 'C1CCCCCCC1'), ('Octan-1-ol', 2.46, 'CCCCCCCCO'), ('Nonan-1-ol', 3.01, 'CCCCCCCCCO'), ('Hexanal', 1.92, 'CCCCCC=O'), ('Heptanal', 2.46, 'CCCCCCC=O'), ('Octanal', 3.01, 'CCCCCCCC=O'), ('Nonanal / n-Nonyl\nAldehyde', 3.42, 'CCCCCCCCC=O'), ('1-Chloroheptane', 4.37, 'CCCCCCCCl'), ('1-Chlorooctane', 4.95, 'CCCCCCCCCl'), ('1-Hexene', 3.51, 'C=CCCCC'), ('1-Heptene', 4.04, 'C=CCCCCC'), ('1-Octene', 4.6, 'C=CCCCCCC'), ('1-Nonene', 5.15, 'C=CCCCCCCC'), ('1-Decene', 5.64, 'C=CCCCCCCCC'), ('1,2,4-Trichlorobenzene', 4.17, 'Clc1ccc(Cl)c(Cl)c1'), ('Di-n-propyl ether', 2.06, 'CCCOCCC'), ('Dipentyl ether', 4.16, 'CCCCCOCCCCC'), ('2-octanone', 2.31, 'CCCCCCC(C)=O'), ('2-nonanone', 2.79, 'CCCCCCCC(C)=O'), ('1-Nitrobutane', 1.53, 'CCCC[N+](=O)[O-]'), ('1-Nitrohexane', 2.58, 'CCCCCC[N+](=O)[O-]'), ('4-Ethylpyridine', 1.14, 'CCc1ccncc1'), ('1-Chloro-4-nitrobenzene', 2.38, 'O=[N+]([O-])c1ccc(Cl)cc1'), ('Nitrobenzene', 1.92, 'O=[N+]([O-])c1ccccc1'), ('2-Nitrotoluene', 2.4, 'Cc1ccccc1[N+](=O)[O-]'), ('2,6-Dinitrotoluene', 1.94, 'Cc1c([N+](=O)[O-])cccc1[N+](=O)[O-]'), ('4-Nitroanisole', 2.32, 'COc1ccc([N+](=O)[O-])cc1'), ('1,4-Dimethoxybenzene', 2.12, 'COc1ccc(OC)cc1'), ('4-Chlorophenol', 1.51, 'Oc1ccc(Cl)cc1'), ('Ethyl benzoate', 2.59, 'CCOC(=O)c1ccccc1'), ('2-Ethyl-1-hexanol', 2.04, 'CCCCC(CC)CO'), ('3-Ethyl-3-hexanol', 1.7, 'CCCC(O)(CC)CC'), ('4-Ethyl-3-hexanol', 1.92, 'CCC(O)C(CC)CC'), ('3-Ethyl-3-pentanol', 1.13, 'CCC(O)(CC)CC'), ('2,4-Dinitrotoluene', 2.34, 'Cc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]'), ('CF3CH2OCF2Cl', 2.73, 'FC(F)(F)COC(F)(F)Cl'), ('1-fluropropane', 0.98, 'CCCF'), ('Hexachloroethane', 3.71, 'ClC(Cl)(Cl)C(Cl)(Cl)Cl'), ('Biphenyl', 4.14, 'c1ccc(-c2ccccc2)cc1'), ('1,1,1,2,3,4,4,4-\nOctafluorobutane', 2.44, 'FC(C(F)C(F)(F)F)C(F)(F)F'), ('hexadecane', 9.88, 'CCCCCCCCCCCCCCCC'), ('isopropyl bromide', 2.28, 'CC(C)Br'), ('1,1,1,2,2,3,3,4,4-\nnonafluorobutane', 2.82, 'FC(F)C(F)(F)C(F)(F)C(F)(F)F'), ('beta-pinene', 4.11, 'C=C1CCC2CC1C2(C)C'), ('limonene', 4.17, 'C=C(C)C1CC=C(C)CC1'), ('Fluoromethane', 0.0, 'CF'), ('Tricyclo[5.2.1.0(2,6)]decane', 4.62, 'C1CC2C3CCC(C3)C2C1'), ('Methyl tert-butyl ether', 0.89, 'COC(C)(C)C'), ('Sulfur hexafluoride', 1.85, 'FS(F)(F)(F)(F)F'), ('vinyl chloride', 1.58, 'C=CCl'), ('dimethyl sulfoxide', -2.66, 'CS(C)=O'), ('formic acid', -1.69, 'O=CO'), ('3-carene', 4.71, 'CC1=CCC2C(C1)C2(C)C'), ('vinyl bromide', 1.34, 'C=CBr'), ('allyl chloride', 1.84, 'C=CCCl'), ('Styrene', 2.68, 'C=Cc1ccccc1'), ('Octane', 5.27, 'CCCCCCCC'), ('3-Chlorophenol', 1.66, 'Oc1cccc(Cl)c1'), ('Benzyl acetate', 1.7, 'CC(=O)OCc1ccccc1'), ('1-Naphthol', 2.19, 'Oc1cccc2ccccc12'), ('4-bromophenol', 1.61, 'Oc1ccc(Br)cc1'), ('N,N-Diethylaniline', 3.17, 'CCN(CC)c1ccccc1'), ('4-n-Propylphenol', 2.19, 'CCCc1ccc(O)cc1'), ('1,3-Dinitrobenzene', 1.42, 'O=[N+]([O-])c1cccc([N+](=O)[O-])c1'), ('Anthracene', 4.83, 'c1ccc2cc3ccccc3cc2c1'), ('Phenanthrene', 4.8, 'c1ccc2c(c1)ccc1ccccc12'), ('Fluoranthene', 5.18, 'c1ccc2c(c1)-c1cccc3cccc-2c13'), ('Pyrene', 5.26, 'c1cc2ccc3cccc4ccc(c1)c2c34'), ('Fluorene', 4.39, 'c1ccc2c(c1)Cc1ccccc1-2'), ('Acenaphthene', 3.97, 'c1cc2c3c(cccc3c1)CC2')]


def generate_xyz_with_rdkit(smiles, xyz_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(m)
    Chem.MolToXYZFile(m, str(xyz_path))
    return xyz_path

def run_xtb_optimization(xyz_path, work_dir, charge=0):
    with open(work_dir/"xtb_opt.err","w") as fe, open(work_dir/"xtb_opt.out","w") as fo:
        subprocess.run([XTB, str(xyz_path.resolve()), "--opt","tight","--gfn","2","-c",str(charge)],
                       cwd=work_dir, stdout=fo, stderr=fe)
    return work_dir/"xtbopt.xyz"

def run_orca_cosmo_sp(xyz_path, work_dir, charge=0):
    lines = xyz_path.read_text().split("\n")
    coords = "\n".join(l for l in lines[2:] if l.strip())
    inp = work_dir/"cosmo_sp.inp"
    inp.write_text(f"! BP86 def2-TZVPD SP CPCM\n%pal nprocs 1 end\n* xyz {charge} 1\n{coords}\n*\n")
    with open(work_dir/"cosmo_sp.out","w") as fo, open(work_dir/"cosmo_sp.err","w") as fe:
        subprocess.run([ORCA, str(inp.resolve())], cwd=work_dir, stdout=fo, stderr=fe)
    cpcm = work_dir/"cosmo_sp.cpcm"
    if not cpcm.exists(): raise FileNotFoundError(f"ORCA COSMO output not found: {cpcm}")
    return cpcm

def convert_to_orcacosmo(mol_dir, name):
    cpcm=mol_dir/"cosmo_sp.cpcm"; outf=mol_dir/"cosmo_sp.out"; xyz=mol_dir/"xtbopt.xyz"
    if not cpcm.exists() or not xyz.exists(): return None
    output=mol_dir/f"{mol_dir.name}.orcacosmo"
    energy=None
    if outf.exists():
        for line in open(outf):
            if "FINAL SINGLE POINT ENERGY" in line:
                mt=re.search(r"[-0-9.]+",line); energy=mt.group() if mt else None; break
    with open(output,"w") as f:
        f.write(f"{name} : BP86def2-TZVPD\n")
        if energy: f.write("\n"+"#"*50+"\n#ENERGY\n"+f"  {energy}\n")
        f.write("\n"+"#"*50+"\n#XYZ_FILE\n"); f.write(open(xyz).read())
        f.write("\n"+"#"*50+"\n#COSMO\n");    f.write(open(cpcm).read())
    return output

def process(name, smiles, dirname):
    mol_dir = CACHE/dirname
    mol_dir.mkdir(exist_ok=True)
    print(f"\n=== {name} ({dirname}) ===")
    init=mol_dir/"init.xyz"
    if not init.exists(): generate_xyz_with_rdkit(smiles, init)
    opt=mol_dir/"xtbopt.xyz"
    if not opt.exists():
        try: run_xtb_optimization(init, mol_dir)
        except Exception as e: print(f"  [ERR xtb] {e}"); return
    if not (mol_dir/"cosmo_sp.cpcm").exists():
        try: run_orca_cosmo_sp(opt, mol_dir)
        except Exception as e: print(f"  [ERR orca] {e}"); return
    oc=mol_dir/f"{mol_dir.name}.orcacosmo"
    if not oc.exists(): convert_to_orcacosmo(mol_dir, name)
    print(f"  [done] {oc}")
    return mol_dir

def do_solvent():
    # NTG -> NTG.orcacosmo (고정명)
    wkey="NTG_TMP"
    d=process("NTG", NTG_SMILES, wkey)
    if d:
        src=d/f"{d.name}.orcacosmo"
        if src.exists(): shutil.copy(src, CACHE/"NTG.orcacosmo"); print("  NTG.orcacosmo 생성")
    # water -> WATER.orcacosmo (고정명)
    d=process("WATER","O","WATER_TMP")
    if d:
        src=d/f"{d.name}.orcacosmo"
        if src.exists(): shutil.copy(src, CACHE/"WATER.orcacosmo"); print("  WATER.orcacosmo 생성")

def do_solutes():
    for name,_exp,smi in STORAGE247:
        key=hashlib.md5(smi.encode()).hexdigest()[:16]
        process(name, smi, f"MOL_{key}")

if __name__=="__main__":
    mode=sys.argv[1] if len(sys.argv)>1 else "all"
    if mode in ("solvent","all"): do_solvent()
    if mode in ("solutes","all"): do_solutes()
    print("\n=== 완료 ===")
