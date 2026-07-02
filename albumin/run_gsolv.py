#!/usr/bin/env python3
"""
Gsolv calculation for BSA factorial pipeline.

Calculates Gsolv for 83 compounds × 4 combinations:
  GFN1/gbsa, GFN1/alpb, GFN2/gbsa, GFN2/alpb

Output: gsolv.csv with columns [compound, method, gsolv]
"""

import pandas as pd
import logging
from pathlib import Path
from tqdm import tqdm
import sys

# Import gsolv_fixed module
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gsolv_fixed import calculate_gsolv, SEED

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler('gsulf_calculation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
LOG = logging.getLogger(__name__)

# Constants
HERE = Path(__file__).resolve().parent
INPUT_CSV = HERE / "docking_poses.csv"
OUTPUT_CSV = HERE / "gsolv.csv"

# 4 combinations
COMBINATIONS = [
    (1, "gbsa"),
    (1, "alpb"),
    (2, "gbsa"),
    (2, "alpb"),
]

def main():
    LOG.info("=" * 60)
    LOG.info("GSOLV CALCULATION: 83 compounds × 4 combinations")
    LOG.info("=" * 60)

    # Load unique compounds from docking_poses.csv
    LOG.info(f"Loading compounds from {INPUT_CSV}")
    df_docking = pd.read_csv(INPUT_CSV)
    compounds = df_docking['compound'].unique().tolist()
    n_compounds = len(compounds)

    LOG.info(f"Unique compounds: {n_compounds}")
    if n_compounds != 83:
        LOG.warning(f"Expected 83 compounds, got {n_compounds}")

    # Get SMILES for each compound
    compound_smiles = df_docking[['compound', 'SMILES']].drop_duplicates('compound')
    compound_smiles = compound_smiles.set_index('compound')['SMILES'].to_dict()

    LOG.info(f"Combinations to calculate: {COMBINATIONS}")
    total_calculations = n_compounds * len(COMBINATIONS)
    LOG.info(f"Total calculations: {total_calculations}")

    # Calculate Gsolv for all combinations
    results = []

    for gfn, solv in COMBINATIONS:
        method = f"GFN{gfn}/{solv}"
        LOG.info(f"\n--- Calculating {method} ---")

        for compound in compounds:
            smiles = compound_smiles.get(compound)
            if smiles is None:
                LOG.warning(f"No SMILES found for {compound}, skipping")
                continue

            LOG.info(f"  {compound} ({method})...")
            gsolv = calculate_gsolv(smiles, gfn=gfn, solvation=solv)

            if gsolv is not None:
                results.append({
                    'compound': compound,
                    'method': method,
                    'gsolv': gsolv
                })
                LOG.info(f"    → Gsolv = {gsolv:.3f} kcal/mol")
            else:
                LOG.error(f"    → FAILED for {compound}")

    # Save results
    LOG.info(f"\nSaving results to {OUTPUT_CSV}")
    df_results = pd.DataFrame(results)
    df_results.to_csv(OUTPUT_CSV, index=False)

    # Summary
    LOG.info("\n" + "=" * 60)
    LOG.info("GSOLV CALCULATION COMPLETE")
    LOG.info("=" * 60)
    LOG.info(f"Total rows in gsolv.csv: {len(df_results)}")
    LOG.info(f"Expected: {total_calculations}")

    # Check for duplicates (different compounds with same Gsolv)
    LOG.info("\n--- Checking for duplicate Gsolv values ---")
    for method in df_results['method'].unique():
        df_method = df_results[df_results['method'] == method]
        n_unique_gsolv = df_method['gsolv'].nunique()
        n_total = len(df_method)
        dup_ratio = 100.0 * (n_total - n_unique_gsolv) / n_total if n_total > 0 else 0
        LOG.info(f"{method}: {n_unique_gsolv} unique / {n_total} total ({dup_ratio:.1f}% duplicates)")

        if dup_ratio > 5.0:
            LOG.error(f"FAIL: {method} has {dup_ratio:.1f}% duplicates > 5%")
            LOG.error("STOPPING - Please check xTB output")
            sys.exit(1)

    # Sample values
    LOG.info("\n--- Sample Gsolv values (first 6) ---")
    print(df_results.head(6).to_string(index=False))

    LOG.info("\n✓ Gsolv calculation complete - all checks passed")

if __name__ == "__main__":
    main()
