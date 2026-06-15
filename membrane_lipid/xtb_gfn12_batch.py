#!/usr/bin/env python3
"""
xTB Multi-Method Batch Calculator for Partition Coefficients
============================================================

Computes log K (solvent/water) for a set of compounds using
GFN1-xTB and GFN2-xTB, each with ALPB and GBSA solvation models.

Per molecule (78 SP total):
  GFN1-xTB:  1 vac + 25 ALPB + 12 GBSA = 38 SP
  GFN2-xTB:  1 vac + 25 ALPB + 14 GBSA = 40 SP

Solvent lists verified against:
  - tblite/src/tblite/solvation/data/shift/  (live library)
  - github.com/grimme-lab/gbsa-parameters    (Ehlert official repo)

References:
  - Ehlert, Stahn, Spicher, Grimme. JCTC 2021, 17, 4250-4261.
    (ALPB; recommended: GFN2-xTB(ALPB), GFN1-xTB(GBSA))
  - Bannwarth, Ehlert, Grimme. JCTC 2019, 15, 1652-1671.
    (GFN2-xTB original — GBSA parameterization)
  - Grimme, Bannwarth, Shushkov. JCTC 2017, 13, 1989-2009.
    (GFN1-xTB original — GBSA first coupling)

Usage:
  python xtb_gfn12_batch.py --input compounds_input_table1.csv \\
                            --output xtb_gfn12_predictions.csv \\
                            --ncpu 16

  # Sanity test (10 compounds):
  python xtb_gfn12_batch.py --input compounds_input_table1.csv \\
                            --output test.csv --ncpu 4 --limit 10
"""
import os
import sys
import argparse
import pickle
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tblite.interface import Calculator

# ----------------------------------------------------------------------------
# Solvent definitions (verified 2026-06-12 against tblite library)
# ----------------------------------------------------------------------------
ALPB_SOLVENTS = [
    'acetone', 'acetonitrile', 'aniline', 'benzaldehyde', 'benzene',
    'ch2cl2', 'chcl3', 'cs2', 'dioxane', 'dmf', 'dmso', 'ethanol',
    'ether', 'ethyl acetate', 'furane', 'hexadecane', 'hexane',
    'methanol', 'nitromethane', 'octanol', 'phenol', 'thf',
    'toluene', 'water', 'woctanol',
]  # 25 solvents — Ehlert 2021 Table 1 (24) + ethanol (PR #387 by M. Stahn)

# GBSA solvent availability differs by GFN method
GBSA_SOLVENTS_GFN1 = [
    'acetone', 'acetonitrile', 'benzene', 'ch2cl2', 'chcl3', 'cs2',
    'dmso', 'ether', 'methanol', 'thf', 'toluene', 'water',
]  # 12 — verified via gbsa-parameters/gfn1-0-1/

GBSA_SOLVENTS_GFN2 = [
    'acetone', 'acetonitrile', 'benzene', 'ch2cl2', 'chcl3', 'cs2',
    'dmf', 'dmso', 'ether', 'hexane', 'methanol', 'thf', 'toluene', 'water',
]  # 14 — verified via gbsa-parameters/gfn2-0-1/  (GFN1 + dmf + hexane)

METHODS = ['GFN1-xTB', 'GFN2-xTB']

T_KELVIN = 310.15  # body temperature for biological partition
RT_kcal = 1.987204e-3 * T_KELVIN  # kcal/mol
LN10 = np.log(10.0)


# ----------------------------------------------------------------------------
# 3D structure generation
# ----------------------------------------------------------------------------
def smiles_to_3d(smiles, seed=42):
    """Generate one 3D conformer with ETKDG + MMFF94.

    Returns (atomic_numbers, coords_bohr) or (None, None) on failure.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        if AllChem.EmbedMolecule(mol, params) != 0:
            # Try without random seed as a fallback
            params.randomSeed = -1
            if AllChem.EmbedMolecule(mol, params) != 0:
                return None, None
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            pass  # MMFF may fail on exotic structures; use unoptimized geom
        conf = mol.GetConformer()
        atoms = np.array([a.GetAtomicNum() for a in mol.GetAtoms()])
        coords_ang = np.array([[conf.GetAtomPosition(i).x,
                                conf.GetAtomPosition(i).y,
                                conf.GetAtomPosition(i).z]
                               for i in range(mol.GetNumAtoms())])
        coords_bohr = coords_ang / 0.529177  # Å → Bohr
        return atoms, coords_bohr
    except Exception:
        return None, None


# ----------------------------------------------------------------------------
# xTB single-point energy
# ----------------------------------------------------------------------------
def xtb_energy(atoms, coords, method, model=None, solvent=None,
               charge=0, uhf=0):
    """Run a single-point xTB calculation.

    Args:
        atoms: array of atomic numbers
        coords: array (N, 3) in Bohr
        method: 'GFN1-xTB' or 'GFN2-xTB'
        model: None | 'alpb' | 'gbsa'  (None = vacuum)
        solvent: solvent name (required if model is not None)

    Returns:
        Total energy in Hartree, or None on failure.
    """
    try:
        c = Calculator(method, atoms, coords, charge=float(charge), uhf=int(uhf))
        if model is not None and solvent is not None:
            c.add(f'{model}-solvation', solvent)
        c.set('verbosity', 0)
        c.set('max-iter', 250)
        res = c.singlepoint()
        return float(res['energy'])
    except Exception:
        return None


# ----------------------------------------------------------------------------
# Per-compound worker (multiprocessing-safe)
# ----------------------------------------------------------------------------
def process_compound(args):
    """Compute all required SP energies for one compound.

    Returns dict with keys:
      cas, compound, smiles, status,
      E_<method>_vac,
      E_<method>_alpb_<solvent>,
      E_<method>_gbsa_<solvent>,
      logK_<method>_<model>_<solvent>  (= log10(K_solvent/water))
    """
    idx, row = args
    out = {
        'cas': row.get('CAS', ''),
        'compound': row.get('compound', ''),
        'smiles': row.get('SMILES', ''),
        'status': 'pending',
    }

    smi = row.get('SMILES', '')
    if not smi or pd.isna(smi):
        out['status'] = 'missing_smiles'
        return out

    atoms, coords = smiles_to_3d(smi)
    if atoms is None:
        out['status'] = 'geometry_failed'
        return out

    # Energy dictionary
    energies = {}  # {(method, model, solvent): E_hartree}, model='vac' for vacuum

    # 1) Vacuum (per method)
    for method in METHODS:
        E = xtb_energy(atoms, coords, method, model=None, solvent=None)
        energies[(method, 'vac', None)] = E
        out[f'E_{method.replace("-xTB","")}_vac'] = E

    # 2) ALPB (25 solvents, both methods)
    for method in METHODS:
        for sv in ALPB_SOLVENTS:
            E = xtb_energy(atoms, coords, method, model='alpb', solvent=sv)
            energies[(method, 'alpb', sv)] = E
            key = f'E_{method.replace("-xTB","")}_alpb_{sv.replace(" ","_")}'
            out[key] = E

    # 3) GBSA (method-specific solvent list)
    gbsa_map = {'GFN1-xTB': GBSA_SOLVENTS_GFN1, 'GFN2-xTB': GBSA_SOLVENTS_GFN2}
    for method, sv_list in gbsa_map.items():
        for sv in sv_list:
            E = xtb_energy(atoms, coords, method, model='gbsa', solvent=sv)
            energies[(method, 'gbsa', sv)] = E
            key = f'E_{method.replace("-xTB","")}_gbsa_{sv.replace(" ","_")}'
            out[key] = E

    # 4) Compute log K for each (method, model, solvent) using water reference
    # log K_solvent/water = (E_water - E_solvent) / (RT * ln 10)
    # Energy in Hartree → kcal/mol: 1 Eh = 627.5094740631 kcal/mol
    EH_to_kcal = 627.5094740631
    for method in METHODS:
        # ALPB water reference
        E_water_alpb = energies.get((method, 'alpb', 'water'))
        for sv in ALPB_SOLVENTS:
            if sv == 'water':
                continue
            E_solv = energies.get((method, 'alpb', sv))
            if E_water_alpb is not None and E_solv is not None:
                dG_kcal = (E_water_alpb - E_solv) * EH_to_kcal
                logK = dG_kcal / (RT_kcal * LN10)
                key = f'logK_{method.replace("-xTB","")}_alpb_{sv.replace(" ","_")}'
                out[key] = logK
        # GBSA water reference
        E_water_gbsa = energies.get((method, 'gbsa', 'water'))
        sv_list = gbsa_map[method]
        for sv in sv_list:
            if sv == 'water':
                continue
            E_solv = energies.get((method, 'gbsa', sv))
            if E_water_gbsa is not None and E_solv is not None:
                dG_kcal = (E_water_gbsa - E_solv) * EH_to_kcal
                logK = dG_kcal / (RT_kcal * LN10)
                key = f'logK_{method.replace("-xTB","")}_gbsa_{sv.replace(" ","_")}'
                out[key] = logK

    # Status
    n_total = len(energies)
    n_ok = sum(1 for v in energies.values() if v is not None)
    if n_ok == n_total:
        out['status'] = 'ok'
    elif n_ok > n_total // 2:
        out['status'] = f'partial_{n_ok}/{n_total}'
    else:
        out['status'] = f'mostly_failed_{n_ok}/{n_total}'

    return out


# ----------------------------------------------------------------------------
# Main pipeline
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input', required=True,
                    help='Input CSV with columns: CAS, compound, SMILES, ...')
    ap.add_argument('--output', required=True,
                    help='Output CSV path for predictions')
    ap.add_argument('--energies-pkl', default=None,
                    help='Optional: save raw energy dict as pickle')
    ap.add_argument('--ncpu', type=int, default=4,
                    help='Number of parallel processes (default 4)')
    ap.add_argument('--limit', type=int, default=None,
                    help='Limit to first N compounds (for testing)')
    args = ap.parse_args()

    print(f'Loading {args.input} ...')
    df = pd.read_csv(args.input)
    if args.limit:
        df = df.head(args.limit)
        print(f'  Limited to {len(df)} compounds')
    print(f'  Loaded {len(df)} compounds')

    # SP count per molecule
    sp_per_mol = (2  # 2 methods × 1 vac
                  + 2 * len(ALPB_SOLVENTS)
                  + len(GBSA_SOLVENTS_GFN1) + len(GBSA_SOLVENTS_GFN2))
    print(f'  SP per molecule: {sp_per_mol}')
    print(f'  Total SP: {sp_per_mol * len(df):,}')
    print(f'  Workers: {args.ncpu}')
    print()

    work = list(df.iterrows())
    t0 = time.time()

    if args.ncpu > 1:
        with Pool(args.ncpu) as pool:
            results = []
            for i, r in enumerate(pool.imap_unordered(process_compound, work, chunksize=1)):
                results.append(r)
                if (i + 1) % max(1, len(work)//20) == 0 or i + 1 == len(work):
                    el = time.time() - t0
                    rate = (i + 1) / el if el > 0 else 0
                    eta = (len(work) - i - 1) / rate if rate > 0 else float('inf')
                    print(f'  [{i+1:4d}/{len(work)}] {el:.0f}s elapsed, '
                          f'{rate:.2f} mol/s, ETA {eta:.0f}s')
    else:
        results = [process_compound(w) for w in work]

    elapsed = time.time() - t0
    print(f'\nCompleted in {elapsed:.1f}s ({elapsed/60:.1f} min)')

    # Preserve input column order
    result_df = pd.DataFrame(results)
    # Merge input metadata back (in case extra columns existed)
    extra_cols = [c for c in df.columns if c not in result_df.columns
                  and c not in ['CAS', 'compound', 'SMILES']]
    if extra_cols:
        df_meta = df[['CAS'] + extra_cols].rename(columns={'CAS': 'cas'})
        result_df = result_df.merge(df_meta, on='cas', how='left')

    result_df.to_csv(args.output, index=False)
    print(f'Saved predictions: {args.output}  ({len(result_df.columns)} columns)')

    # Status summary
    print('\nStatus breakdown:')
    print(result_df['status'].value_counts().to_string())

    # Optional: save raw energy structure
    if args.energies_pkl:
        with open(args.energies_pkl, 'wb') as f:
            pickle.dump(results, f)
        print(f'Saved raw energies pickle: {args.energies_pkl}')


if __name__ == '__main__':
    main()
