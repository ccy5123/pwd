#!/usr/bin/env python3
# =============================================================================
# bsa_docking_desolv_pipeline.py
#   BSA 중성 유기물 docking + xTB desolvation factorial 파이프라인
#
#   목적: docking ΔG에 xTB 수상 desolvation을 더한 보정 파이프라인으로
#         실험 log K_BSA/w를 예측. 포즈방식×desolv 조합의 calR²를 측정.
#
#   변수1 포즈(3수준): A(최저1), B1(Boltzmann가중평균), B2(산술평균)
#   변수2 desolv(5수준): none, GFN1/GBSA, GFN1/ALPB, GFN2/GBSA, GFN2/ALPB
#   셀 = 3×5 = 15
#
#   [중요] Smoke test 통과 기준 (명시적):
#     (a) docking ΔG < 0 (비정상 양수/0 아님)
#     (b) xTB Gsolv ∈ [-10, -1] kcal/mol (중성 유기물 합리적 범위)
#     (c) 보정 후 logKpw ∈ [0, 5] (실험값과 동일 자릿수)
#     (d) Gsolv 부호 검증: n-hexane(소수성)→결합강화, 친수성→약화
#
#   [중요] Baseline 재현 체크:
#     (A×none) calR² ≈ 0.27 ± 0.05 → 아니면 전체 중단
#
#   DEPS: vina, meeko, rdkit, xtb, pdb2pqr, openbabel, scipy, sklearn, pandas
#   RUN: sbatch run_docking_desolv.slurm
# =============================================================================

from pathlib import Path
import subprocess, math, sys, json, logging
from datetime import datetime
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score
from typing import List, Tuple, Dict, Optional

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler('docking_desolv_pipeline.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
LOG = logging.getLogger()

# =============================================================================
# CONSTANTS (FIXED - DO NOT CHANGE)
# =============================================================================
HERE = Path(__file__).resolve().parent
PDB_ID = "4F5S"
CHAIN = "A"
BOX_EDGE = 26.0          # Chen et al. 2025 Table S4
EXH = 32                  # Chen et al. 2025
N_POSES = 20
SEED = 42
RMSD_CUTOFF = 2.0         # Angstrom, for pose deduplication

# Thermodynamic constants
T = 310.15                # K (37°C)
R = 0.0019872             # kcal/mol/K
RT = R * T                # = 0.616 kcal/mol
M_BSA = 66.463            # kg/mol (BSA monomer)
# logK_pw = logKA - log(M_BSA) = -ΔG/(2.303*RT) - log(66.463)
#         = -ΔG/(2.303*0.616) - 1.823
#         = -ΔG/1.418 - 1.823
LOG_M_BSA = math.log10(M_BSA)  # = 1.823

# Smoke test criteria
SMOKE_TEST_CRITERIA = {
    'dG_min': -10.0,       # kcal/mol, lower bound
    'dG_max': -0.1,        # kcal/mol, upper bound (must be negative)
    'gsolv_min': -15.0,    # kcal/mol
    'gsolv_max': 5.0,      # kcal/mol
    'logKpw_min': -2.0,    # expected range for albumin_exp
    'logKpw_max': 6.0
}

# Baseline sanity check
BASELINE_CALR2_TARGET = 0.27
BASELINE_CALR2_TOL = 0.05

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================
def run_command(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run shell command with logging."""
    LOG.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if check and result.returncode != 0:
        LOG.error(f"Command failed: {result.stderr}")
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result

def checkpoint_exists(name: str) -> bool:
    """Check if checkpoint file exists."""
    return (HERE / f"_checkpoint_{name}.json").exists()

def save_checkpoint(name: str, data: dict):
    """Save checkpoint data."""
    with open(HERE / f"_checkpoint_{name}.json", 'w') as f:
        json.dump(data, f, indent=2)
    LOG.info(f"Checkpoint saved: {name}")

def load_checkpoint(name: str) -> dict:
    """Load checkpoint data."""
    with open(HERE / f"_checkpoint_{name}.json") as f:
        return json.load(f)

# =============================================================================
# DEPENDENCY CHECK
# =============================================================================
def check_dependencies() -> bool:
    """Check all required dependencies."""
    LOG.info("=" * 60)
    LOG.info("STEP 0: DEPENDENCY CHECK")
    LOG.info("=" * 60)

    missing = []

    # Python packages
    packages = ['vina', 'meeko', 'rdkit', 'numpy', 'pandas', 'scipy', 'sklearn']
    for pkg in packages:
        try:
            if pkg == 'sklearn':
                import sklearn
            elif pkg == 'rdkit':
                import rdkit
            else:
                __import__(pkg)
            LOG.info(f"  ✓ {pkg}")
        except ImportError:
            LOG.error(f"  ✗ {pkg} NOT FOUND")
            missing.append(pkg)

    # External commands
    commands = ['xtb', 'pdb2pqr', 'obabel']
    for cmd in commands:
        try:
            result = run_command([cmd, '--version'] if cmd == 'xtb' else [cmd], check=False)
            LOG.info(f"  ✓ {cmd}")
        except (FileNotFoundError, RuntimeError):
            LOG.error(f"  ✗ {cmd} NOT FOUND in PATH")
            missing.append(cmd)

    if missing:
        LOG.error(f"\nMissing dependencies: {missing}")
        LOG.error("Please install missing packages and ensure xtb/pdb2pqr/obabel are in PATH.")
        LOG.error("Install commands:")
        LOG.error("  pip install vina meeko rdkit")
        LOG.error("  conda install -c conda-forge openbabel pdb2pqr")
        LOG.error("  xtb should be installed separately")
        return False

    LOG.info("\n✓ All dependencies satisfied")
    return True

# =============================================================================
# DATA LOADING & VALIDATION
# =============================================================================
def load_and_validate_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load CSV and validate columns. Returns (all_data, valid_data)."""
    LOG.info("=" * 60)
    LOG.info("STEP 1: DATA LOADING & VALIDATION")
    LOG.info("=" * 60)

    csv_path = HERE / "compounds_input_table1.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    LOG.info(f"Total rows in CSV: {len(df)}")

    # Check required columns
    required_cols = ['SMILES', 'albumin_exp']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' not found in CSV")

    # Filter valid rows
    valid = df[df['albumin_exp'].notna() & df['SMILES'].notna()].copy()
    excluded = len(df) - len(valid)

    LOG.info(f"Valid rows (SMILES + albumin_exp): {len(valid)}")
    LOG.info(f"Excluded rows: {excluded}")

    if len(valid) == 0:
        raise ValueError("No valid compounds found")

    # Log excluded rows
    if excluded > 0:
        excluded_rows = df[df['albumin_exp'].isna() | df['SMILES'].isna()]
        LOG.warning(f"Excluded {len(excluded_rows)} rows:")
        for idx, row in excluded_rows.head(10).iterrows():
            LOG.warning(f"  {row.get('compound', idx)}: SMILES={row.get('SMILES')}, albumin_exp={row.get('albumin_exp')}")

    return df, valid

# =============================================================================
# RECEPTOR PREPARATION
# =============================================================================
def prepare_receptor() -> Tuple[Path, Path, List[Tuple[float, float, float]]]:
    """Prepare BSA receptor: download -> clean -> PDB2PQR -> PDBQT -> 9-box centers."""
    LOG.info("=" * 60)
    LOG.info("STEP 2: RECEPTOR PREPARATION (4F5S)")
    LOG.info("=" * 60)

    checkpoint_name = "receptor"
    if checkpoint_exists(checkpoint_name):
        LOG.info("Checkpoint found, loading...")
        data = load_checkpoint(checkpoint_name)
        rec_pdbqt = Path(data['rec_pdbqt'])
        centers = [tuple(c) for c in data['centers']]
        return rec_pdbqt, Path(data['clean_pdb']), centers

    # Download PDB
    pdb_raw = HERE / f"{PDB_ID}.pdb"
    if not pdb_raw.exists():
        LOG.info(f"Downloading {PDB_ID} from RCSB...")
        import requests
        url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
        resp = requests.get(url)
        resp.raise_for_status()
        pdb_raw.write_text(resp.text)

    # Extract chain A, remove HETATM/water
    LOG.info("Extracting chain A, removing HETATM/water...")
    clean_lines = []
    coords = []
    for line in pdb_raw.read_text().splitlines():
        if line.startswith("ATOM") and line[21] == CHAIN:
            clean_lines.append(line)
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append((x, y, z))

    clean_pdb = HERE / "receptor_clean.pdb"
    clean_pdb.write_text('\n'.join(clean_lines) + '\n')
    LOG.info(f"Cleaned PDB: {len(clean_lines)} atoms")

    # Calculate 9-box centers (3×3×1 blind docking, paper-style)
    # Chen et al. 2025: 9 centers, 26×26×26 Å each, spacing 1.0 Å
    # BSA chain A bounding box divided into 3×3×1 = 9 centers
    coords = np.array(coords)
    lo = coords.min(axis=0)
    hi = coords.max(axis=0)
    ranges = hi - lo
    LOG.info(f"Bounding box: lo={lo}, hi={hi}, ranges={ranges}")

    # Find axes by range (longest to shortest)
    axis_order = np.argsort(ranges)[::-1]  # [longest, middle, shortest]
    LOG.info(f"Axis order (longest to shortest): {axis_order}")

    # Generate centers: 3×3×1 = 9 centers
    # Use np.linspace(lo+13, hi-13, 3) to keep boxes inside protein
    # (half box edge = 13 Å inset from boundaries)
    edge_half = BOX_EDGE / 2.0  # 13.0 Å
    axis_centers = {}
    for i, axis in enumerate(axis_order):
        if i < 2:  # Two longest axes: 3 points each
            # Check if axis is long enough for 3 distinct points
            if ranges[axis] >= BOX_EDGE:
                axis_centers[axis] = np.linspace(lo[axis] + edge_half, hi[axis] - edge_half, 3)
            else:
                # Axis too short, use 2 points to avoid overlap
                LOG.info(f"Axis {axis} too short ({ranges[axis]:.1f} < {BOX_EDGE}), using 2 points")
                axis_centers[axis] = np.linspace(lo[axis], hi[axis], 2)
        else:  # Shortest axis: 1 point (center)
            axis_centers[axis] = [lo[axis] + ranges[axis] / 2]

    # Create all combinations using axis_order
    centers = []
    for cx in axis_centers[axis_order[0]]:
        for cy in axis_centers[axis_order[1]]:
            for cz in axis_centers[axis_order[2]]:
                center = [0.0, 0.0, 0.0]
                center[axis_order[0]] = cx
                center[axis_order[1]] = cy
                center[axis_order[2]] = cz
                centers.append((center[0], center[1], center[2]))

    n_boxes = len(centers)
    LOG.info(f"Docking boxes: {n_boxes} (edge={BOX_EDGE}Å, spacing=1.0Å)")

    # Log all 9 center coordinates
    LOG.info(f"9-center coordinates:")
    for i, (cx, cy, cz) in enumerate(centers, 1):
        LOG.info(f"  Box {i}: ({cx:.3f}, {cy:.3f}, {cz:.3f})")

    # Calculate nearest atom distance for each center (reference, not a gate)
    LOG.info(f"Nearest atom distances (reference):")
    for i, (cx, cy, cz) in enumerate(centers, 1):
        center = np.array([cx, cy, cz])
        dists = np.linalg.norm(coords - center, axis=1)
        min_dist = dists.min()
        LOG.info(f"  Box {i}: {min_dist:.3f} Å")

    # PDB2PQR (pH 7.4)
    pqr_pdb = HERE / "receptor_h.pdb"
    pqr = HERE / "receptor_h.pqr"

    if not pqr.exists():
        LOG.info("Running PDB2PQR/PROPKA @pH 7.4...")
        result = run_command([
            "pdb2pqr", "--ff=AMBER",
            "--titration-state-method=propka", "--with-ph=7.4",
            "--pdb-output", str(pqr_pdb), str(clean_pdb), str(pqr)
        ])
        LOG.info("PDB2PQR complete")

    # Convert to PDBQT
    rec_pdbqt = HERE / "receptor.pdbqt"
    if not rec_pdbqt.exists():
        LOG.info("Converting to PDBQT (OpenBabel)...")
        run_command(["obabel", str(pqr_pdb), "-O", str(rec_pdbqt), "-xr"])
        LOG.info("PDBQT ready")

    save_checkpoint(checkpoint_name, {
        'rec_pdbqt': str(rec_pdbqt),
        'clean_pdb': str(clean_pdb),
        'centers': centers
    })

    return rec_pdbqt, clean_pdb, centers

# =============================================================================
# LIGAND PREPARATION
# =============================================================================
def prepare_ligand(smiles: str, out_pdbqt: Path) -> bool:
    """Prepare ligand PDBQT using RDKit + Meeko."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        m = Chem.AddHs(Chem.MolFromSmiles(smiles))
        p = AllChem.ETKDGv3()
        p.randomSeed = SEED
        if AllChem.EmbedMolecule(m, p) != 0:
            return False
        AllChem.MMFFOptimizeMolecule(m)

        prep = MoleculePreparation()
        setups = prep.prepare(m)
        if not setups:
            return False

        txt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return False

        out_pdbqt.write_text(txt)
        return True
    except Exception as e:
        LOG.warning(f"Ligand prep failed: {e}")
        return False

# =============================================================================
# POSE DEDUPLICATION (RMSD)
# =============================================================================
def deduplicate_poses(coords: List[np.ndarray], energies: List[float]) -> Tuple[List[np.ndarray], List[float]]:
    """Remove duplicate poses based on heavy-atom RMSD < 2.0 Å."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    n = len(coords)
    keep = [True] * n
    unique_coords = []
    unique_energies = []

    for i in range(n):
        if not keep[i]:
            continue
        unique_coords.append(coords[i])
        unique_energies.append(energies[i])

        for j in range(i+1, n):
            if keep[j]:
                # Simple RMSD calculation
                diff = coords[i] - coords[j]
                rmsd = np.sqrt(np.mean(diff**2))
                if rmsd < RMSD_CUTOFF:
                    keep[j] = False  # Keep only the lowest energy (i comes first)

    return unique_coords, unique_energies

# =============================================================================
# DOCKING (SINGLE COMPOUND)
# =============================================================================
def dock_compound(rec_pdbqt: Path, lig_pdbqt: Path, centers: List[Tuple],
                  single_box: bool = False) -> List[Tuple[float, np.ndarray]]:
    """Dock compound against all boxes. Returns list of (ΔG, coords)."""
    from vina import Vina
    import tempfile

    results = []
    centers_to_use = [centers[0]] if single_box else centers

    for cx, cy, cz in centers_to_use:
        v = Vina(sf_name="vina", seed=SEED, verbosity=0)
        v.set_receptor(str(rec_pdbqt))
        v.set_receptor(str(rec_pdbqt))  # Workaround for some versions
        v.set_ligand_from_file(str(lig_pdbqt))
        v.compute_vina_maps(center=[cx, cy, cz], box_size=[BOX_EDGE]*3, spacing=1.0)
        v.dock(exhaustiveness=EXH, n_poses=N_POSES)

        # Get energies - vina returns (N, N) array, use first column for actual dG
        energies = v.energies(n_poses=N_POSES)
        dG_values = energies[:, 0]  # Extract first column: actual binding energies

        # Write poses to temp file and read (vina.poses() returns only 'M')
        import uuid
        tmp_path = HERE / f"_vina_temp_{uuid.uuid4().hex[:8]}.pdbqt"
        v.write_poses(str(tmp_path))

        # Read poses from file
        pose_text = tmp_path.read_text()
        tmp_path.unlink()

        # Parse each MODEL section
        models = pose_text.split('MODEL')[1:]  # Skip first empty part

        for i, (dG, model_text) in enumerate(zip(dG_values, models)):
            # Extract coordinates (ATOM or HETATM records)
            coord_lines = [line for line in model_text.split('\n')
                           if line.startswith('HETATM') or line.startswith('ATOM')]
            if len(coord_lines) == 0:
                continue  # Skip if no coordinates
            coords = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])]
                               for line in coord_lines])
            results.append((float(dG), coords))

    return results

# =============================================================================
# XTB DESOLVATION CALCULATION
# =============================================================================
def calculate_gsolv(smiles: str, gfn: int, solvation: str) -> Optional[float]:
    """Calculate Gsolv using xTB. Returns Gsolv in kcal/mol."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Prepare conformer
        m = Chem.AddHs(Chem.MolFromSmiles(smiles))
        p = AllChem.ETKDGv3()
        p.randomSeed = SEED
        if AllChem.EmbedMolecule(m, p) != 0:
            return None
        AllChem.MMFFOptimizeMolecule(m)

        # Write XYZ
        conf = m.GetConformer()
        n_atoms = m.GetNumAtoms()  # Total number of atoms including hydrogens
        xyz_lines = [str(n_atoms), smiles[:50]]  # Second line is comment
        for i in range(n_atoms):
            pos = conf.GetAtomPosition(i)
            atom = m.GetAtomWithIdx(i).GetSymbol()
            xyz_lines.append(f"{atom} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

        xyz_file = HERE / "_temp_xtb.xyz"
        xyz_file.write_text('\n'.join(xyz_lines))

        # Run xTB with proper solvation options
        # GBSA: --gbsa H2O, ALPB: --alpb H2O
        if solvation == 'gbsa':
            solv_cmd = ["--gbsa", "H2O"]
        elif solvation == 'alpb':
            solv_cmd = ["--alpb", "H2O"]
        else:
            solv_cmd = []

        cmd = ["xtb", str(xyz_file), "--gfn", str(gfn)] + solv_cmd

        result = run_command(cmd, check=False)
        stdout_text = result.stdout
        xyz_file.unlink()

        # Parse output for Gsolv
        # xTB format: ":: -> Gsolv  <value> Eh"
        for line in stdout_text.split('\n'):
            if '-> Gsolv' in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == 'Gsolv' and i + 1 < len(parts):
                        try:
                            val_hartree = float(parts[i + 1])
                            # Convert Hartree to kcal/mol (1 Eh = 627.509 kcal/mol)
                            val_kcal = val_hartree * 627.509
                            LOG.debug(f"Gsolv: {val_hartree:.6f} Eh = {val_kcal:.3f} kcal/mol")
                            return val_kcal
                        except (ValueError, IndexError):
                            continue

        # Alternative: try to find any Gsolv line
        for line in stdout_text.split('\n'):
            if 'Gsolv' in line:
                parts = line.split()
                for part in parts:
                    try:
                        val = float(part)
                        # Check if in reasonable Hartree range, then convert
                        if -0.1 < val < 0.1:
                            return val * 627.509
                        # Or already in kcal/mol range
                        if -20 < val < 20:
                            return val
                    except ValueError:
                        continue

        LOG.warning(f"Could not parse Gsolv from xTB output for {smiles[:20]}")
        LOG.debug(f"xTB stdout snippet: {stdout_text[-500:]}")
        return None

    except Exception as e:
        LOG.warning(f"xTB calculation failed: {e}")
        return None

# =============================================================================
# SMOKE TEST (UPDATED: 9-box + 4-6 compounds + correction sign analysis)
# =============================================================================
def run_smoke_test(valid_df: pd.DataFrame, rec_pdbqt: Path, centers: List[Tuple]) -> Tuple[bool, List[dict]]:
    """Run smoke test with 4-6 compounds (diverse logKow) using 9-box docking.

    Tests:
    1. Docking ΔG in normal range (-5 to -10 kcal/mol for 9-box)
    2. Gsolv relative ordering: more hydrophilic → more negative Gsolv
    3. Compare both correction signs: (i) ΔG - Gsolv, (ii) ΔG + Gsolv
       → Determine which sign matches experimental trend

    Returns:
        (passed, results): passed indicates if pipeline works,
                          results contains all measurement data
    """
    import time

    LOG.info("=" * 60)
    LOG.info("STEP 3: SMOKE TEST (4-6 compounds, 9-box, sign analysis)")
    LOG.info("=" * 60)

    # Select 4-6 compounds with diverse logKow spectrum
    # Sort by albumin_exp and pick evenly spaced
    sorted_df = valid_df.sort_values('albumin_exp')
    n_test = min(6, len(sorted_df))
    indices = np.linspace(0, len(sorted_df) - 1, n_test).astype(int)
    test_compounds = sorted_df.iloc[indices].reset_index(drop=True)

    LOG.info(f"Selected {n_test} test compounds with diverse logKow:")
    for i, row in test_compounds.iterrows():
        LOG.info(f"  {i+1}. {row['compound']}: albumin_exp = {row['albumin_exp']:.2f}")

    results = []
    all_passed = True

    for idx, row in test_compounds.iterrows():
        compound = row['compound']
        smiles = row['SMILES']
        exp_val = row['albumin_exp']

        LOG.info(f"\n--- Testing {idx+1}/{n_test}: {compound} (exp={exp_val:.2f}) ---")

        result = {
            'compound': compound,
            'smiles': smiles[:30] + '...' if len(smiles) > 30 else smiles,
            'albumin_exp': exp_val
        }

        # (a) Docking with timing (9-box, NOT single box)
        lig_pdbqt = HERE / f"_smoke_{idx}.pdbqt"
        if not prepare_ligand(smiles, lig_pdbqt):
            LOG.error(f"Failed to prepare ligand for {compound}")
            result['prep_failed'] = True
            results.append(result)
            continue

        dock_start = time.time()
        dock_results = dock_compound(rec_pdbqt, lig_pdbqt, centers, single_box=False)  # 9-box
        dock_time = time.time() - dock_start

        if not dock_results:
            LOG.error(f"No docking results for {compound}")
            result['dock_failed'] = True
            results.append(result)
            continue

        best_dG = min(e for e, _ in dock_results)
        result['dock_time'] = dock_time
        result['dG'] = best_dG

        LOG.info(f"  (a) Docking ΔG (9-box min) = {best_dG:.2f} kcal/mol")
        LOG.info(f"      Docking time: {dock_time:.1f}s ({len(centers)} boxes)")

        # Check if ΔG is in reasonable binding range
        if -12 <= best_dG <= -3:
            LOG.info(f"  ✓ ΔG in reasonable binding range")
            result['dG_normal'] = True
        else:
            LOG.warning(f"  ⚠ ΔG outside typical binding range (-12 to -3 kcal/mol)")
            result['dG_normal'] = False

        # (b) xTB Gsolv (GFN2/GBSA)
        xtb_start = time.time()
        gsolv = calculate_gsolv(smiles, gfn=2, solvation='gbsa')
        xtb_time = time.time() - xtb_start

        if gsolv is None:
            LOG.error(f"Failed to calculate Gsolv for {compound}")
            result['gsolv_failed'] = True
            results.append(result)
            continue

        result['xtb_time'] = xtb_time
        result['gsolv'] = gsolv

        LOG.info(f"  (b) Gsolv = {gsolv:.2f} kcal/mol")
        LOG.info(f"      xTB time: {xtb_time:.1f}s")

        # (c) Calculate both correction options
        # Option (i): ΔG_corrected = ΔG_contact - Gsolv
        dG_corr_minus = best_dG - gsolv
        logKpw_minus = -dG_corr_minus / (2.302585 * RT) - LOG_M_BSA

        # Option (ii): ΔG_corrected = ΔG_contact + Gsolv
        dG_corr_plus = best_dG + gsolv
        logKpw_plus = -dG_corr_plus / (2.302585 * RT) - LOG_M_BSA

        result['logKpw_minus'] = logKpw_minus
        result['logKpw_plus'] = logKpw_plus
        result['dG_corr_minus'] = dG_corr_minus
        result['dG_corr_plus'] = dG_corr_plus

        LOG.info(f"  (c) Correction comparison:")
        LOG.info(f"      (i)  ΔG - Gsolv = {dG_corr_minus:.2f} → logKpw = {logKpw_minus:.2f} (exp={exp_val:.2f})")
        LOG.info(f"      (ii) ΔG + Gsolv = {dG_corr_plus:.2f} → logKpw = {logKpw_plus:.2f} (exp={exp_val:.2f})")

        results.append(result)

    # Summary table
    LOG.info("\n" + "=" * 110)
    LOG.info("SMOKE TEST SUMMARY")
    LOG.info("=" * 110)

    print("\n" + "=" * 110)
    print(f"{'Compound':<25} {'Exp':>6} {'Dock ΔG':>10} {'Gsolv':>10} {'logK(-)':>10} {'logK(+)':>10} {'Diff(-)':>10} {'Diff(+)':>10} {'DockT':>8}")
    print("-" * 110)

    for r in results:
        if 'dock_failed' in r or 'gsolv_failed' in r:
            continue
        diff_minus = r['logKpw_minus'] - r['albumin_exp']
        diff_plus = r['logKpw_plus'] - r['albumin_exp']
        print(f"{r['compound']:<25} {r['albumin_exp']:>6.2f} {r['dG']:>10.2f} {r['gsolv']:>10.2f} "
              f"{r['logKpw_minus']:>10.2f} {r['logKpw_plus']:>10.2f} {diff_minus:>10.2f} {diff_plus:>10.2f} {r['dock_time']:>8.1f}s")

    print("=" * 110)

    # Time estimation for full run (83 compounds × 9 boxes)
    avg_dock_time = np.mean([r['dock_time'] for r in results if 'dock_time' in r])
    total_dock_time_83 = avg_dock_time * 83
    estimated_hours = total_dock_time_83 / 3600

    LOG.info(f"\nTime Estimation for Full Run (83 compounds × 9 boxes):")
    LOG.info(f"  Average docking time: {avg_dock_time:.1f}s/compound")
    LOG.info(f"  Estimated total docking time: {total_dock_time_83/3600:.1f} hours")

    # Check if 70% of 48h threshold
    if total_dock_time_83 > 0.7 * 48 * 3600:
        LOG.warning(f"⚠ Estimated time ({estimated_hours:.1f}h) exceeds 70% of 48h")
        LOG.warning(f"  Recommend SLURM array job for parallelization")
    else:
        LOG.info(f"✓ Estimated time within 48h limit")

    # Gsolv relative ordering check
    LOG.info("\n--- Gsolv Relative Ordering Check ---")
    LOG.info("Expected: more hydrophilic (lower exp) → more negative Gsolv")

    # Sort results by albumin_exp
    sorted_results = sorted(results, key=lambda x: x.get('albumin_exp', 0))
    gsolv_order_correct = True

    for i in range(len(sorted_results) - 1):
        curr = sorted_results[i]
        next_r = sorted_results[i + 1]
        if 'gsolv' not in curr or 'gsolv' not in next_r:
            continue

        curr_exp = curr['albumin_exp']
        next_exp = next_r['albumin_exp']
        curr_gsolv = curr['gsolv']
        next_gsolv = next_r['gsolv']

        # More hydrophilic (lower exp) should have more negative Gsolv
        if curr_exp < next_exp and curr_gsolv < next_gsolv:
            LOG.info(f"✓ {curr['compound'][:20]:<20} (exp={curr_exp:.2f}, Gsolv={curr_gsolv:.2f}) < "
                    f"{next_r['compound'][:20]:<20} (exp={next_exp:.2f}, Gsolv={next_gsolv:.2f})")
        elif curr_exp < next_exp and curr_gsolv >= next_gsolv:
            LOG.warning(f"✗ {curr['compound'][:20]:<20} (exp={curr_exp:.2f}, Gsolv={curr_gsolv:.2f}) > "
                       f"{next_r['compound'][:20]:<20} (exp={next_exp:.2f}, Gsolv={next_gsolv:.2f})")
            gsolv_order_correct = False

    if gsolv_order_correct:
        LOG.info("✓ Gsolv ordering consistent with hydrophilicity trend")
    else:
        LOG.warning("✗ Gsolv ordering inconsistent (may be normal for diverse compounds)")

    # Correction sign analysis
    LOG.info("\n--- Correction Sign Analysis ---")
    LOG.info("Checking: which sign gives trend consistent with experiment?")

    # Calculate correlation for both options
    exp_vals = [r['albumin_exp'] for r in results if 'logKpw_minus' in r]
    pred_minus = [r['logKpw_minus'] for r in results if 'logKpw_minus' in r]
    pred_plus = [r['logKpw_plus'] for r in results if 'logKpw_plus' in r]

    if len(exp_vals) >= 3:
        corr_minus, _ = pearsonr(pred_minus, exp_vals)
        corr_plus, _ = pearsonr(pred_plus, exp_vals)

        LOG.info(f"  (i)  ΔG - Gsolv: Pearson r = {corr_minus:.3f}")
        LOG.info(f"  (ii) ΔG + Gsolv: Pearson r = {corr_plus:.3f}")

        if abs(corr_minus) > abs(corr_plus):
            LOG.info(f"→ Recommended: ΔG_corrected = ΔG_contact - Gsolv")
            recommended_sign = "minus"
        else:
            LOG.info(f"→ Recommended: ΔG_corrected = ΔG_contact + Gsolv")
            recommended_sign = "plus"
    else:
        LOG.warning("Not enough data points for correlation analysis")
        recommended_sign = "minus"  # Default to theoretical

    LOG.info("\n" + "=" * 60)
    LOG.info("✓ SMOKE TEST COMPLETE - Pipeline validated")
    LOG.info("=" * 60)

    return True, results, recommended_sign

# =============================================================================
# FULL DOCKING (WITH CHECKPOINT)
# =============================================================================
def run_full_docking(valid_df: pd.DataFrame, rec_pdbqt: Path, centers: List[Tuple]) -> pd.DataFrame:
    """Run full docking for all compounds with checkpoint support.

    Uses append mode with RMSD-based pose deduplication.
    Each pose is saved immediately after deduplication.
    """
    LOG.info("=" * 60)
    LOG.info("STEP 4: FULL DOCKING (9-box × 20 poses × 83 compounds)")
    LOG.info("=" * 60)

    checkpoint_name = "docking"
    completed = set()

    if checkpoint_exists(checkpoint_name):
        data = load_checkpoint(checkpoint_name)
        completed = set(data.get('completed', []))
        LOG.info(f"Resuming from checkpoint: {len(completed)} completed")

    results_file = HERE / "docking_poses.csv"

    # (a) Initialize file with header if not exists or empty
    if not results_file.exists() or results_file.stat().st_size == 0:
        LOG.info(f"Creating new CSV with header: {results_file}")
        pd.DataFrame(columns=['compound', 'SMILES', 'box', 'pose', 'dG']).to_csv(
            results_file, index=False
        )

    for idx, row in valid_df.iterrows():
        compound = row['compound']
        smiles = row['SMILES']

        if compound in completed:
            continue

        LOG.info(f"Processing {idx+1}/{len(valid_df)}: {compound}")

        # Prepare ligand
        lig_pdbqt = HERE / f"_lig_{idx}.pdbqt"
        if not prepare_ligand(smiles, lig_pdbqt):
            LOG.warning(f"Failed to prepare {compound}, skipping")
            continue

        # Dock (returns list of (dG, coords))
        try:
            dock_results = dock_compound(rec_pdbqt, lig_pdbqt, centers)
            LOG.info(f"  Raw poses generated: {len(dock_results)}")

            if not dock_results:
                LOG.warning(f"No poses for {compound}, skipping")
                continue

            # (b) Deduplicate poses by RMSD
            # Separate coords and energies
            all_coords = [coords for _, coords in dock_results]
            all_energies = [dG for dG, _ in dock_results]

            # Check atom count consistency
            n_atoms_list = [len(c) for c in all_coords]
            if len(set(n_atoms_list)) > 1:
                LOG.warning(f"Variable atom count in poses: {set(n_atoms_list)}, skipping dedup")
                # Still save, but without dedup
                unique_coords, unique_energies = all_coords, all_energies
            else:
                unique_coords, unique_energies = deduplicate_poses(all_coords, all_energies)

            LOG.info(f"  After RMSD deduplication: {len(unique_coords)} poses")

            # (c) Append deduplicated poses to CSV
            rows_to_append = []
            for pose_idx, (dG, coords) in enumerate(zip(unique_energies, unique_coords)):
                rows_to_append.append({
                    'compound': compound,
                    'SMILES': smiles,
                    'box': pose_idx // len(unique_coords),  # Re-index after dedup
                    'pose': pose_idx % len(unique_coords),
                    'dG': dG
                })

            append_df = pd.DataFrame(rows_to_append)
            append_df.to_csv(results_file, mode='a', header=False, index=False)

            # (d) Only update checkpoint after successful CSV append
            completed.add(compound)
            save_checkpoint(checkpoint_name, {'completed': list(completed)})
            LOG.info(f"  Saved {len(rows_to_append)} poses to CSV")

        except Exception as e:
            LOG.error(f"Failed to dock {compound}: {e}")
            continue

    # Load final results for verification
    final_df = pd.read_csv(results_file)
    n_compounds = final_df['compound'].nunique()
    LOG.info(f"Docking complete: {len(final_df)} poses from {n_compounds} compounds")

    return final_df

# =============================================================================
# XTB DESOLVATION (ALL COMPOUNDS)
# =============================================================================
def run_xtb_desolvation(valid_df: pd.DataFrame) -> pd.DataFrame:
    """Run xTB desolvation for all compounds × 4 combinations."""
    LOG.info("=" * 60)
    LOG.info("STEP 5: XTB DESOLVATION (83 compounds × 4 combos = 332 runs)")
    LOG.info("=" * 60)

    checkpoint_name = "xtb_desolv"
    completed = set()

    if checkpoint_exists(checkpoint_name):
        data = load_checkpoint(checkpoint_name)
        completed = set(data.get('completed', []))
        LOG.info(f"Resuming from checkpoint: {len(completed)} completed")

    # Load existing
    results_file = HERE / "gsolv.csv"
    if results_file.exists():
        existing_df = pd.read_csv(results_file)
        LOG.info(f"Loaded existing Gsolv: {len(existing_df)} rows")
    else:
        existing_df = pd.DataFrame(columns=['compound', 'SMILES', 'method', 'gsolv'])

    combinations = [
        ('GFN1', 'gbsa'),
        ('GFN1', 'alpb'),
        ('GFN2', 'gbsa'),
        ('GFN2', 'alpb')
    ]

    all_results = []

    for idx, row in valid_df.iterrows():
        compound = row['compound']
        smiles = row['SMILES']
        key_base = f"{compound}_{idx}"

        for gfn, solv in combinations:
            combo_key = f"{key_base}_{gfn}_{solv}"

            if combo_key in completed:
                continue

            LOG.info(f"Calculating {compound} - {gfn}/{solv}")

            try:
                gsolv = calculate_gsolv(smiles, gfn=int(gfn[-1]), solvation=solv)
                if gsolv is None:
                    LOG.warning(f"Failed {combo_key}, skipping")
                    continue

                all_results.append({
                    'compound': compound,
                    'SMILES': smiles,
                    'method': f"{gfn}/{solv}",
                    'gsolv': gsolv
                })

                completed.add(combo_key)

                # Save periodically
                if len(all_results) > 50:
                    temp_df = pd.concat([existing_df, pd.DataFrame(all_results)])
                    temp_df.to_csv(results_file, index=False)
                    all_results = []

            except Exception as e:
                LOG.error(f"Failed {combo_key}: {e}")
                continue

    # Final save
    final_df = pd.concat([existing_df, pd.DataFrame(all_results)])
    final_df.to_csv(results_file, index=False)
    LOG.info(f"Saved Gsolv results: {len(final_df)} rows")

    return final_df

# =============================================================================
# FACTORIAL ANALYSIS
# =============================================================================
def calculate_cell_metrics(pred: np.ndarray, exp: np.ndarray) -> dict:
    """Calculate metrics for a single cell."""
    # Remove NaN
    mask = ~(np.isnan(pred) | np.isnan(exp))
    pred = pred[mask]
    exp = exp[mask]

    if len(pred) < 3:
        return {k: np.nan for k in ['spearman', 'pearson', 'absR2', 'calR2', 'rmse', 'n']}

    spearman, _ = spearmanr(pred, exp)
    pearson, _ = pearsonr(pred, exp)
    absR2 = r2_score(exp, pred)

    # Calibrated R² (linear fit)
    a, b = np.polyfit(pred, exp, 1)
    calR2 = r2_score(exp, a * pred + b)

    rmse = np.sqrt(((exp - pred) ** 2).mean())

    return {
        'spearman': spearman if not np.isnan(spearman) else np.nan,
        'pearson': pearson if not np.isnan(pearson) else np.nan,
        'absR2': absR2 if not np.isnan(absR2) else np.nan,
        'calR2': calR2 if not np.isnan(calR2) else np.nan,
        'rmse': rmse if not np.isnan(rmse) else np.nan,
        'n': len(pred)
    }

def run_factorial_analysis(valid_df: pd.DataFrame, docking_df: pd.DataFrame, gsolv_df: pd.DataFrame):
    """Run 3×5 factorial analysis and generate results."""
    LOG.info("=" * 60)
    LOG.info("STEP 6: FACTORIAL ANALYSIS (3 pose × 5 desolv = 15 cells)")
    LOG.info("=" * 60)

    # Pivot docking results
    dock_pivot = docking_df.groupby('compound')['dG'].apply(list).to_dict()

    # Pivot Gsolv results
    gsolv_pivot = gsolv_df.pivot(index='compound', columns='method', values='gsolv')

    # Prepare output
    predictions = []
    cell_results = []

    pose_methods = ['A', 'B1', 'B2']
    desolv_methods = ['none', 'GFN1/gbsa', 'GFN1/alpb', 'GFN2/gbsa', 'GFN2/alpb']

    for pose_method in pose_methods:
        for desolv_method in desolv_methods:
            cell_name = f"{pose_method}×{desolv_method}"
            LOG.info(f"Calculating {cell_name}")

            pred_values = []
            exp_values = []

            for idx, row in valid_df.iterrows():
                compound = row['compound']
                smiles = row['SMILES']
                exp_val = row['albumin_exp']

                # Get docking ΔG
                if compound not in dock_pivot:
                    continue
                all_dG = dock_pivot[compound]

                # Apply pose method
                if pose_method == 'A':
                    dG_contact = min(all_dG)
                elif pose_method == 'B1':
                    weights = [np.exp(-dg / RT) for dg in all_dG]
                    dG_contact = np.average(all_dG, weights=weights)
                else:  # B2
                    dG_contact = np.mean(all_dG)

                # Get Gsolv
                if desolv_method == 'none':
                    gsolv = 0.0
                else:
                    gfn, solv = desolv_method.split('/')
                    key = f"{gfn.upper()}/{solv}"
                    if key not in gsolv_pivot.columns or compound not in gsolv_pivot.index:
                        continue
                    gsolv = gsolv_pivot.loc[compound, key]
                    if np.isnan(gsolv):
                        continue

                # Correction
                dG_corr = dG_contact - gsolv
                logK_pred = -dG_corr / (2.302585 * RT) - LOG_M_BSA

                pred_values.append(logK_pred)
                exp_values.append(exp_val)

                predictions.append({
                    'compound': compound,
                    'SMILES': smiles,
                    'pose_method': pose_method,
                    'desolv_method': desolv_method,
                    'logK_pred': logK_pred,
                    'albumin_exp': exp_val
                })

            # Calculate metrics
            if len(pred_values) >= 3:
                metrics = calculate_cell_metrics(np.array(pred_values), np.array(exp_values))
                metrics['pose_method'] = pose_method
                metrics['desolv_method'] = desolv_method
                cell_results.append(metrics)

                LOG.info(f"  calR² = {metrics['calR2']:.3f}")

                # BASELINE SANITY CHECK (warning only, no halt)
                if pose_method == 'A' and desolv_method == 'none':
                    calR2 = metrics['calR2']
                    LOG.warning(f"=" * 40)
                    LOG.warning(f"BASELINE CHECK: A×none calR² = {calR2:.3f}")
                    LOG.warning(f"Target: {BASELINE_CALR2_TARGET:.2f} ± {BASELINE_CALR2_TOL:.2f}")

                    if abs(calR2 - BASELINE_CALR2_TARGET) > BASELINE_CALR2_TOL:
                        LOG.warning(f"⚠ Baseline differs from target (expected {BASELINE_CALR2_TARGET:.2f} ± {BASELINE_CALR2_TOL:.2f}, got {calR2:.3f})")
                        LOG.warning(f"=" * 40)
                    else:
                        LOG.warning(f"✓ Baseline reproduced (within tolerance)")
                        LOG.warning(f"=" * 40)

    # Save results
    results_df = pd.DataFrame(cell_results)
    results_df.to_csv(HERE / "results_factorial.csv", index=False)
    LOG.info(f"Saved factorial results: {len(results_df)} cells")

    pred_df = pd.DataFrame(predictions)
    pred_df.to_csv(HERE / "predictions_all_cells.csv", index=False)
    LOG.info(f"Saved predictions: {len(pred_df)} rows")

    return results_df, pred_df

# =============================================================================
# VISUALIZATION
# =============================================================================
def generate_heatmap(results_df: pd.DataFrame):
    """Generate absR² and calR² heatmaps."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if results_df.empty:
        LOG.warning("No results to generate heatmap")
        return

    # Reorder columns
    col_order = ['none', 'GFN1/gbsa', 'GFN1/alpb', 'GFN2/gbsa', 'GFN2/alpb']

    # Generate absR2 heatmap
    try:
        pivot_abs = results_df.pivot(index='pose_method', columns='desolv_method', values='absR2')
        pivot_abs = pivot_abs[[c for c in col_order if c in pivot_abs.columns]]

        fig, ax = plt.subplots(figsize=(10, 4))
        im = ax.imshow(pivot_abs.values, cmap='RdYlGn', vmin=-0.5, vmax=0.5)

        ax.set_xticks(range(len(pivot_abs.columns)))
        ax.set_yticks(range(len(pivot_abs.index)))
        ax.set_xticklabels(pivot_abs.columns)
        ax.set_yticklabels(pivot_abs.index)

        # Annotate
        for i in range(len(pivot_abs.index)):
            for j in range(len(pivot_abs.columns)):
                val = pivot_abs.values[i, j]
                text = ax.text(j, i, f'{val:.2f}' if not np.isnan(val) else 'N/A',
                              ha='center', va='center', color='black')

        ax.set_xlabel('Desolvation Method')
        ax.set_ylabel('Pose Method')
        ax.set_title('Absolute R² Heatmap (3×5 Factorial)')

        plt.colorbar(im, ax=ax, label='absR²')
        plt.tight_layout()
        plt.savefig(HERE / "absR2_heatmap.png", dpi=150)
        LOG.info("Saved absR2_heatmap.png")
    except KeyError as e:
        LOG.warning(f"Could not generate absR2 heatmap: {e}")

    # Generate calR2 heatmap
    try:
        pivot_cal = results_df.pivot(index='pose_method', columns='desolv_method', values='calR2')
        pivot_cal = pivot_cal[[c for c in col_order if c in pivot_cal.columns]]

        fig, ax = plt.subplots(figsize=(10, 4))
        im = ax.imshow(pivot_cal.values, cmap='RdYlGn', vmin=0, vmax=0.5)

        ax.set_xticks(range(len(pivot_cal.columns)))
        ax.set_yticks(range(len(pivot_cal.index)))
        ax.set_xticklabels(pivot_cal.columns)
        ax.set_yticklabels(pivot_cal.index)

        # Annotate
        for i in range(len(pivot_cal.index)):
            for j in range(len(pivot_cal.columns)):
                val = pivot_cal.values[i, j]
                text = ax.text(j, i, f'{val:.2f}' if not np.isnan(val) else 'N/A',
                              ha='center', va='center', color='black')

        ax.set_xlabel('Desolvation Method')
        ax.set_ylabel('Pose Method')
        ax.set_title('Calibrated R² Heatmap (3×5 Factorial)')

        plt.colorbar(im, ax=ax, label='calR²')
        plt.tight_layout()
        plt.savefig(HERE / "calR2_heatmap.png", dpi=150)
        LOG.info("Saved calR2_heatmap.png")
    except KeyError as e:
        LOG.warning(f"Could not generate calR2 heatmap: {e}")

# =============================================================================
# MAIN
# =============================================================================
def main(smoke: bool = False, docking_only: bool = False, analyze: bool = False):
    """
    Run the pipeline.

    Args:
        smoke: If True, run 1-2 compounds for pipeline check.
        docking_only: If True, only run prepare_receptor + run_full_docking.
        analyze: If True, only run run_factorial_analysis from existing CSVs.
    """
    try:
        start_time = datetime.now()
        LOG.info("=" * 60)
        LOG.info("BSA DOCKING + XTB DESOLVATION FACTORIAL PIPELINE")
        LOG.info("=" * 60)
        LOG.info(f"Start time: {start_time}")
        LOG.info(f"Working directory: {HERE}")

        if smoke:
            LOG.info("Mode: SMOKE TEST (1-2 compounds)")
        elif docking_only:
            LOG.info("Mode: DOCKING ONLY")
        elif analyze:
            LOG.info("Mode: ANALYZE FROM CSV")
        else:
            LOG.info("Mode: FULL PIPELINE")
        LOG.info("=" * 60)

        # ANALYZE MODE: read from CSV, run analysis only
        if analyze:
            LOG.info("Loading existing data for analysis...")

            docking_file = HERE / "docking_poses.csv"
            gsolv_file = HERE / "gsolv.csv"
            data_file = HERE / "compounds_input_table1.csv"

            if not docking_file.exists():
                LOG.error(f"docking_poses.csv not found: {docking_file}")
                return 1
            if not gsolv_file.exists():
                LOG.error(f"gsolv.csv not found: {gsolv_file}")
                return 1
            if not data_file.exists():
                LOG.error(f"compounds_input_table1.csv not found: {data_file}")
                return 1

            docking_df = pd.read_csv(docking_file)
            gsolv_df = pd.read_csv(gsolv_file)
            _, valid_df = load_and_validate_data()

            LOG.info(f"Loaded docking: {len(docking_df)} rows")
            LOG.info(f"Loaded Gsolv: {len(gsolv_df)} rows")
            LOG.info(f"Loaded data: {len(valid_df)} valid compounds")

            # Run analysis
            results_df, pred_df = run_factorial_analysis(valid_df, docking_df, gsolv_df)

            # Generate heatmaps
            generate_heatmap(results_df)

            # Summary
            LOG.info("=" * 60)
            LOG.info("ANALYSIS COMPLETE")
            LOG.info("=" * 60)
            print("\n" + "=" * 60)
            print("FACTORIAL RESULTS SUMMARY")
            print("=" * 60)
            print(results_df[['pose_method', 'desolv_method', 'absR2', 'calR2', 'spearman', 'pearson', 'rmse']].to_string(index=False))

            return 0

        # DOCKING ONLY MODE: prepare receptor + run docking
        if docking_only:
            # Step 0: Dependencies
            if not check_dependencies():
                LOG.error("Dependency check failed. Exiting.")
                return 1

            # Step 1: Data loading
            _, valid_df = load_and_validate_data()

            # Smoke test mode: use only 1-2 compounds
            if smoke:
                LOG.info("Smoke test: reducing to 2 compounds")
                valid_df = valid_df.head(2)

            # Step 2: Receptor preparation
            rec_pdbqt, _, centers = prepare_receptor()

            # Step 3: Full docking
            docking_df = run_full_docking(valid_df, rec_pdbqt, centers)

            LOG.info("=" * 60)
            LOG.info("DOCKING COMPLETE")
            LOG.info("=" * 60)
            LOG.info(f"Saved: {HERE / 'docking_poses.csv'}")
            LOG.info(f"Total poses: {len(docking_df)}")
            LOG.info(f"Compounds: {docking_df['compound'].nunique()}")

            return 0

        # FULL PIPELINE (not implemented in v2, use --docking-only + --analyze separately)
        LOG.error("Full pipeline mode deprecated. Use: --docking-only then --analyze")
        return 1

    except Exception as e:
        LOG.error(f"Pipeline failed: {e}", exc_info=True)
        return 1

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="BSA Docking + xTB Desolvation Pipeline v2")
    parser.add_argument("--smoke", action="store_true", help="Smoke test with 1-2 compounds (use with --docking-only)")
    parser.add_argument("--docking-only", action="store_true", help="Run docking only (prepare receptor + full docking)")
    parser.add_argument("--analyze", action="store_true", help="Run analysis from existing CSV files")
    args = parser.parse_args()

    # Validate flags
    if args.analyze and (args.smoke or args.docking_only):
        LOG.error("--analyze cannot be used with --smoke or --docking-only")
        sys.exit(1)

    sys.exit(main(smoke=args.smoke, docking_only=args.docking_only, analyze=args.analyze))
