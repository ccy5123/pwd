#!/bin/bash
#==============================================================================
# run_cosmo_muscle_v5.sh
# openCOSMO-RS muscle K_pw validation (chicken core-6 sequence-based composition)
#==============================================================================
#SBATCH --job-name=cosmo_muscle_v5
#SBATCH --output=logs/cosmo_muscle_v5_%j.out
#SBATCH --error=logs/cosmo_muscle_v5_%j.err
#SBATCH --time=48:00:00              # 48h wall-time limit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=gpu1,gpu2,gpu3,gpu4,gpu5,gpu6
#SBATCH --cpus-per-task=16           # CPU; xTB/ORCA/openCOSMO-RS는 CPU 위주
#SBATCH --mem=64G

set -euo pipefail

# ---- 0. 환경 ----
ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT"
mkdir -p logs

source /gpfs/home1/s9383/.bashrc
conda activate rapids-cuml

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="/home1/s9383/Muscle_protein/openCOSMO-RS_conformer_pipeline${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

echo "[INFO] host=$(hostname)  cpus=${OMP_NUM_THREADS}  start=$(date)"
echo "[INFO] CACHE=./cosmo_work_v3"

# ---- 1. 사전 점검 (cosmo_work_v3/sigma_cache에 AA 20종 + WATER 있어야 함) ----
WATER="cosmo_work_v3/sigma_cache/WATER.orcacosmo"
if [[ ! -f "$WATER" ]]; then
  echo "[ERROR] $WATER 없음. 아미노산/물 .orcacosmo 캐시를 먼저 생성하세요." >&2
  exit 1
fi
naa=$(ls -d cosmo_work_v3/sigma_cache/AA_*_* 2>/dev/null | wc -l)
echo "[INFO] AA pseudo-solvent dirs found: ${naa}/20"
nmol=$(ls -d cosmo_work_v3/sigma_cache/MOL_* 2>/dev/null | wc -l)
echo "[INFO] MOL dirs found: ${nmol}"

# ---- 2. 실행 (48h 안에서 완료) ----
#   - CSV 단일소스 리팩터 버전 (v5_csv) 사용
timeout 47h python cosmo_kpw_muscle_v5_csv.py validate

# ---- 3. 결과 위치 ----
echo "[INFO] 결과: cosmo_work_v3/results_chicken_core6_seq_v5.json"
echo "[INFO] done=$(date)"
