#!/usr/bin/env python3
"""
Albumin σ-profile 계산 스크립트 v4
muscle calc_missing_molecules_v4.py와 동일 워크플로:
  RDKit ETKDGv3+MMFF -> xTB(GFN2) opt -> ORCA BP86/def2-TZVPD CPCM -> .orcacosmo
변경점: 대상 분자 = ALBUMIN83 solute + capped residue 20종 + water
캐시 규약(muscle 동일): solute=MOL_<md5>, 아미노산=AA_<3letter>_<md5>, water=WATER.orcacosmo
"""
import subprocess, hashlib, re, sys
from pathlib import Path

WORK = Path("./cosmo_work_v4_albumin")
CACHE = WORK / "sigma_cache"
CACHE.mkdir(parents=True, exist_ok=True)

# ★ 서버 실제 경로로 수정 (muscle 원본 경로 그대로 둠) ★
XTB  = "/home1/s9383/Muscle_protein/xtb-6.5.1/bin/xtb"
ORCA = "/home1/s9383/Muscle_protein/orca_6_1_1_linux_x86-64_shared_openmpi418/orca"

ALBUMIN83 = [('1,2,4-Trimethylbenzene', 3.35, 'Cc1ccc(C)c(C)c1'), ('1,2-Dichlorobenzene', 3.03, 'Clc1ccccc1Cl'), ('1-Hexanol', 1.64, 'CCCCCCO'), ('Benzene', 1.58, 'c1ccccc1'), ('Carbon tetrachloride', 1.77, 'ClC(Cl)(Cl)Cl'), ('Halothane', 1.62, 'FC(F)(F)C(Cl)Br'), ('Enflurane', 1.59, 'FC(F)OC(F)(F)C(F)Cl'), ('Isoflurane', 1.58, 'FC(F)OC(Cl)C(F)(F)F'), ('Chlorobenzene', 2.32, 'Clc1ccccc1'), ('Cycloheptane', 2.52, 'C1CCCCCC1'), ('Cyclohexane', 2.01, 'C1CCCCC1'), ('Ethylbenzene', 2.7, 'CCc1ccccc1'), ('Heptane', 3.59, 'CCCCCCC'), ('Hexafluorobenzene', 1.55, 'Fc1c(F)c(F)c(F)c(F)c1F'), ('Hexane', 3.09, 'CCCCCC'), ('Methoxyflurane', 1.77, 'COC(F)(F)C(Cl)Cl'), ('Methylpentafluorobenzene', 2.32, 'Cc1c(F)c(F)c(F)c(F)c1F'), ('Nonane', 4.45, 'CCCCCCCCC'), ('Propylbenzene', 2.95, 'CCCc1ccccc1'), ('Tetrachloroethene', 2.4, 'ClC(Cl)=C(Cl)Cl'), ('Toluene', 2.26, 'Cc1ccccc1'), ('Trichloroethene', 1.88, 'ClC=C(Cl)Cl'), ('1-heptanol', 2.18, 'CCCCCCCO'), ('di-n-butyl ether', 2.01, 'CCCCOCCCC'), ('Cyclooctane', 2.98, 'C1CCCCCCC1'), ('Octan-1-ol', 2.74, 'CCCCCCCCO'), ('Nonan-1-ol', 3.1, 'CCCCCCCCCO'), ('1-Chlorooctane', 3.85, 'CCCCCCCCCl'), ('1-Nonene', 4.22, 'C=CCCCCCCC'), ('1,2,4-Trichlorobenzene', 3.6, 'Clc1ccc(Cl)c(Cl)c1'), ('Dipentyl ether', 3.0, 'CCCCCOCCCCC'), ('2-octanone', 2.09, 'CCCCCCC(C)=O'), ('2-nonanone', 2.48, 'CCCCCCCC(C)=O'), ('2-Nitrotoluene', 2.12, 'Cc1ccccc1[N+](=O)[O-]'), ('4-Nitroanisole', 2.48, 'COc1ccc([N+](=O)[O-])cc1'), ('4-Chlorophenol', 2.43, 'Oc1ccc(Cl)cc1'), ('4-Ethyl-3-hexanol', 1.48, 'CCC(O)C(CC)CC'), ('2,4-Dinitrotoluene', 1.73, 'Cc1ccc([N+](=O)[O-])cc1[N+](=O)[O-]'), ('Styrene', 2.76, 'C=Cc1ccccc1'), ('Octane', 4.01, 'CCCCCCCC'), ('3-Chlorophenol', 2.35, 'Oc1cccc(Cl)c1'), ('4-bromophenol', 2.81, 'Oc1ccc(Br)cc1'), ('4-Iodoaniline', 2.95, 'Nc1ccc(I)cc1'), ('N,N-Diethylaniline', 2.27, 'CCN(CC)c1ccccc1'), ('4-n-Propylphenol', 2.59, 'CCCc1ccc(O)cc1'), ('4-iodophenol', 3.41, 'Oc1ccc(I)cc1'), ('Phenanthrene', 4.15, 'c1ccc2c(c1)ccc1ccccc12'), ('Fluoranthene', 4.28, 'c1ccc2c(c1)-c1cccc3cccc-2c13'), ('Pyrene', 4.76, 'c1cc2ccc3cccc4ccc(c1)c2c34'), ('1-heptyne', 2.49, 'C#CCCCCC'), ('tribromomethane', 1.95, 'BrC(Br)Br'), ('γ-hexachlorocyclohexane', 2.46, 'ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl'), ('2-decanone', 2.88, 'CCCCCCCCC(C)=O'), ('1-nitrooctane', 3.38, 'CCCCCCCC[N+](=O)[O-]'), ('tri-n-butyl phosphate', 2.47, 'CCCCOP(=O)(OCCCC)OCCCC'), ('1,4-dibromobenzene', 3.97, 'Brc1ccc(Br)cc1'), ('1,2,3,4-tetrachlorobenzene', 4.21, 'Clc1ccc(Cl)c(Cl)c1Cl'), ('indene', 2.92, 'C1=Cc2ccccc2C1'), ('naphthalene', 3.56, 'c1ccc2ccccc2c1'), ('dibenzofuran', 3.79, 'c1ccc2c(c1)oc1ccccc12'), ('dibenzothiophene', 4.16, 'c1ccc2c(c1)sc1ccccc12'), ('chrysene', 4.46, 'c1ccc2c(c1)ccc1c3ccccc3ccc21'), ('benzo[b]fluoranthene', 4.42, 'c1ccc2c(c1)-c1cccc3c1c-2cc1ccccc13'), ('benzo[ghi]perylene', 4.76, 'c1cc2ccc3ccc4ccc5cccc6c(c1)c2c3c4c56'), ('anisole', 2.16, 'COc1ccccc1'), ('valerophenone', 2.7, 'CCCCC(=O)c1ccccc1'), ('benzophenone', 2.62, 'O=C(c1ccccc1)c1ccccc1'), ('di-n-propyl phthalate', 2.84, 'CCCOC(=O)c1ccccc1C(=O)OCCC'), ('1-nitronaphthalene', 3.17, 'O=[N+]([O-])c1cccc2ccccc12'), ('4-chlorobenzyl alcohol', 2.1, 'OCc1ccc(Cl)cc1'), ('2-phenylphenol', 2.62, 'Oc1ccccc1-c1ccccc1'), ('4-fluorophenol', 1.57, 'Oc1ccc(F)cc1'), ('bisphenol A', 2.88, 'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1'), ('4-nitroaniline', 1.69, 'Nc1ccc([N+](=O)[O-])cc1'), ('2-chloroaniline', 1.95, 'Nc1ccccc1Cl'), ('4-aminobiphenyl', 2.55, 'Nc1ccc(-c2ccccc2)cc1'), ('indole', 2.25, 'c1ccc2[nH]ccc2c1'), ('carbazole', 3.52, 'c1ccc2c(c1)[nH]c1ccccc12'), ('metolachlor', 1.74, 'CCc1cccc(C)c1N(C(=O)CCl)C(C)COC'), ('atrazine', 1.77, 'CCNc1nc(Cl)nc(NC(C)C)n1'), ('diazepam', 2.68, 'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21'), ('estrone', 2.69, 'CC12CCC3c4ccc(O)cc4CCC3C1CCC2=O'), ('endosulfan α', 3.24, 'O=S1OCC2C(CO1)C1(Cl)C(Cl)=C(Cl)C2(Cl)C1(Cl)Cl')]

CAPPED_RESIDUE_SMILES = {'GLY': 'CC(=O)NCC(=O)NC', 'ALA': 'CC(=O)N[C@@H](C)C(=O)NC', 'VAL': 'CC(=O)N[C@@H](C(C)C)C(=O)NC', 'LEU': 'CC(=O)N[C@@H](CC(C)C)C(=O)NC', 'ILE': 'CC(=O)N[C@@H]([C@@H](C)CC)C(=O)NC', 'PRO': 'CC(=O)N1CCC[C@H]1C(=O)NC', 'PHE': 'CC(=O)N[C@@H](Cc1ccccc1)C(=O)NC', 'TRP': 'CC(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)NC', 'MET': 'CC(=O)N[C@@H](CCSC)C(=O)NC', 'SER': 'CC(=O)N[C@@H](CO)C(=O)NC', 'THR': 'CC(=O)N[C@@H]([C@@H](C)O)C(=O)NC', 'CYS': 'CC(=O)N[C@@H](CS)C(=O)NC', 'TYR': 'CC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NC', 'ASN': 'CC(=O)N[C@@H](CC(N)=O)C(=O)NC', 'GLN': 'CC(=O)N[C@@H](CCC(N)=O)C(=O)NC', 'ASP': 'CC(=O)N[C@@H](CC(O)=O)C(=O)NC', 'GLU': 'CC(=O)N[C@@H](CCC(O)=O)C(=O)NC', 'LYS': 'CC(=O)N[C@@H](CCCCN)C(=O)NC', 'ARG': 'CC(=O)N[C@@H](CCCNC(N)=N)C(=O)NC', 'HIS': 'CC(=O)N[C@@H](Cc1c[nH]cn1)C(=O)NC'}


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
    if not cpcm.exists():
        raise FileNotFoundError(f"ORCA COSMO output not found: {cpcm}")
    return cpcm

def convert_to_orcacosmo(mol_dir, name):
    cpcm = mol_dir/"cosmo_sp.cpcm"; outf = mol_dir/"cosmo_sp.out"; xyz = mol_dir/"xtbopt.xyz"
    if not cpcm.exists() or not xyz.exists(): return None
    output = mol_dir/f"{mol_dir.name}.orcacosmo"
    energy=None
    if outf.exists():
        for line in open(outf):
            if "FINAL SINGLE POINT ENERGY" in line:
                mt=re.search(r"[-0-9.]+", line);  energy=mt.group() if mt else None; break
    with open(output,"w") as f:
        f.write(f"{name} : BP86def2-TZVPD\n")
        if energy:
            f.write("\n"+"#"*50+"\n#ENERGY\n"+f"  {energy}\n")
        f.write("\n"+"#"*50+"\n#XYZ_FILE\n"); f.write(open(xyz).read())
        f.write("\n"+"#"*50+"\n#COSMO\n");    f.write(open(cpcm).read())
    return output

def process(name, smiles, prefix):
    key = hashlib.md5(smiles.encode()).hexdigest()[:16]
    mol_dir = CACHE / (f"{prefix}_{key}" if prefix=="MOL" else f"{prefix}_{key}")
    mol_dir.mkdir(exist_ok=True)
    print(f"\n=== {name} ({mol_dir.name}) ===")
    init = mol_dir/"init.xyz"
    if not init.exists(): generate_xyz_with_rdkit(smiles, init)
    opt = mol_dir/"xtbopt.xyz"
    if not opt.exists():
        try: run_xtb_optimization(init, mol_dir)
        except Exception as e: print(f"  [ERR xtb] {e}"); return
    if not (mol_dir/"cosmo_sp.cpcm").exists():
        try: run_orca_cosmo_sp(opt, mol_dir)
        except Exception as e: print(f"  [ERR orca] {e}"); return
    oc = mol_dir/f"{mol_dir.name}.orcacosmo"
    if not oc.exists(): convert_to_orcacosmo(mol_dir, name)
    print(f"  [done] {oc}")

def do_residues():
    for aa, smi in CAPPED_RESIDUE_SMILES.items():
        process(aa, smi, f"AA_{aa}")  # -> AA_<3letter>_<md5>
    # water
    process("WATER", "O", "WATER_TMP")
    # water는 고정 이름으로 복사 (muscle 규약: WATER.orcacosmo)
    import hashlib, shutil
    wkey = hashlib.md5("O".encode()).hexdigest()[:16]
    wdir = CACHE/f"WATER_TMP_{wkey}"
    src = wdir/f"{wdir.name}.orcacosmo"
    if src.exists(): shutil.copy(src, CACHE/"WATER.orcacosmo"); print("  WATER.orcacosmo 생성")

def do_solutes():
    for name, _exp, smi in ALBUMIN83:
        process(name, smi, "MOL")

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv)>1 else "all"
    if mode in ("residues","all"): do_residues()
    if mode in ("solutes","all"):  do_solutes()
    print("\n=== 완료 ===")
