#!/usr/bin/env python3
"""
Fixed Gsolv calculation module.

Key fixes:
1. Unique temporary directory per compound (no _temp_xtb.xyz race condition)
2. Parse only "-> Gsolv" line (single format, Hartree→kcal conversion)
3. No fallback parsing (strict parsing)
4. ETKDGv3 seed fixed for reproducibility

Author: Claude Code
Date: 2025-07-01
"""

import tempfile
import logging
from pathlib import Path
from typing import Optional
from subprocess import run, CalledProcessError

LOG = logging.getLogger(__name__)

# Constants
SEED = 42
HARTREE_TO_KCAL = 627.509  # 1 Eh = 627.509 kcal/mol


def calculate_gsolv(smiles: str, gfn: int, solvation: str) -> Optional[float]:
    """
    Calculate Gsolv using xTB with strict parsing.

    Args:
        smiles: SMILES string
        gfn: GFN level (1 or 2)
        solvation: 'gbsa' or 'alpb'

    Returns:
        Gsolv in kcal/mol, or None if failed
    """
    tmpdir = None
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # === Step 1: Generate conformer with FIXED seed ===
        m = Chem.AddHs(Chem.MolFromSmiles(smiles))
        p = AllChem.ETKDGv3()
        p.randomSeed = SEED  # Fixed seed for reproducibility
        if AllChem.EmbedMolecule(m, p) != 0:
            LOG.warning(f"Conformer generation failed for {smiles[:30]}")
            return None
        AllChem.MMFFOptimizeMolecule(m)

        # === Step 2: Write XYZ in UNIQUE temporary directory ===
        tmpdir = tempfile.mkdtemp(prefix=f"gsolv_{gfn}_{solvation}_")
        conf = m.GetConformer()
        n_atoms = m.GetNumAtoms()
        xyz_lines = [str(n_atoms), smiles[:50]]  # Second line is comment
        for i in range(n_atoms):
            pos = conf.GetAtomPosition(i)
            atom = m.GetAtomWithIdx(i).GetSymbol()
            xyz_lines.append(f"{atom} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

        xyz_file = Path(tmpdir) / "input.xyz"
        xyz_file.write_text('\n'.join(xyz_lines))

        # === Step 3: Run xTB ===
        if solvation == 'gbsa':
            solv_cmd = ["--gbsa", "H2O"]
        elif solvation == 'alpb':
            solv_cmd = ["--alpb", "H2O"]
        else:
            raise ValueError(f"Invalid solvation: {solvation}")

        cmd = ["xtb", str(xyz_file), "--gfn", str(gfn)] + solv_cmd

        LOG.debug(f"Running: {' '.join(cmd)}")
        result = run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode != 0:
            LOG.warning(f"xtb failed (exit {result.returncode}) for {smiles[:30]}")
            LOG.debug(f"stderr: {result.stderr[-500:]}")
            return None

        stdout_text = result.stdout

        # === Step 4: Parse ONLY "-> Gsolv" line (strict, no fallback) ===
        # xTB format: ":: -> Gsolv  <value> Eh"
        for line in stdout_text.split('\n'):
            if '-> Gsolv' in line:
                parts = line.split()
                for i, part in enumerate(parts):
                    if part == 'Gsolv' and i + 1 < len(parts):
                        try:
                            val_hartree = float(parts[i + 1])
                            val_kcal = val_hartree * HARTREE_TO_KCAL
                            LOG.debug(f"Gsolv parsed: {val_hartree:.6f} Eh = {val_kcal:.3f} kcal/mol")
                            return val_kcal
                        except (ValueError, IndexError):
                            LOG.warning(f"Failed to parse Gsolv from line: {line}")
                            return None

        # No "-> Gsolv" line found = failure
        LOG.warning(f"No '-> Gsolv' line in xTB output for {smiles[:30]}")
        LOG.debug(f"stdout snippet: {stdout_text[-500:]}")
        return None

    except CalledProcessError as e:
        LOG.warning(f"xtb crashed: {e}")
        return None
    except Exception as e:
        LOG.warning(f"Gsolv calculation error: {e}")
        return None
    finally:
        # Clean up temporary directory
        if tmpdir:
            import shutil
            try:
                shutil.rmtree(tmpdir, ignore_errors=True)
            except Exception:
                pass


def validate_gsolv_reproducibility(smiles_list, n_trials=2):
    """
    Validate Gsolv reproducibility: same compound should give same value.

    Returns: dict with {compound: [trial1, trial2, ...]}
    """
    results = {}
    for smiles in smiles_list:
        key = smiles[:30]  # Short key
        results[key] = []
        for _ in range(n_trials):
            val = calculate_gsolv(smiles, gfn=2, solvation='gbsa')
            if val is not None:
                results[key].append(val)
            else:
                results[key].append(None)
    return results


if __name__ == "__main__":
    # Quick test
    logging.basicConfig(level=logging.DEBUG)
    test_smiles = ["CCCCCC", "c1ccccc1", "CC(=O)O"]
    for smi in test_smiles:
        result = calculate_gsolv(smi, gfn=2, solvation='gbsa')
        print(f"{smi[:20]:20} -> Gsolv = {result}")
