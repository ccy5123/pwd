#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=gpu1,gpu2,gpu3,gpu4,gpu5,gpu6
#SBATCH --cpus-per-task=4
#SBATCH --job-name=gen_bap
#SBATCH --time=48:00:00
#SBATCH -o ./logs/gen_bap_%j.out
#SBATCH -e ./logs/gen_bap_%j.err

echo "start at: $(date)"
echo "node: $HOSTNAME"
echo "jobid: $SLURM_JOB_ID"

source /gpfs/home1/s9383/.bashrc
conda activate rapids-cuml

# ORCA용 LD_LIBRARY_PATH (unbound 오류 방지)
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

echo "=== benzo[a]pyrene 단일 분자 재생성 시작 (canonical-md5) ==="
cd /home1/s9383/Muscle_protein
timeout 47h $CONDA_PREFIX/bin/python gen_bap.py
echo "=== 완료 ==="
