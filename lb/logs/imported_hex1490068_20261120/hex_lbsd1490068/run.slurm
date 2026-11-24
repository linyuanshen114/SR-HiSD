#!/bin/bash
#SBATCH --job-name=BCCHEX_LBSD_HISTV4
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --nodelist=c02
#SBATCH --output=/bicmr/home/yslin114/LB/experiments/bcc2hex_lbsd_historical_v4_rebenchmark_20261117/logs/%x_%j.out
#SBATCH --error=/bicmr/home/yslin114/LB/experiments/bcc2hex_lbsd_historical_v4_rebenchmark_20261117/logs/%x_%j.err
set -euo pipefail
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
TASK_ROOT=/bicmr/home/yslin114/LB/experiments/bcc2hex_lbsd_historical_v4_rebenchmark_20261117
cd "$TASK_ROOT/code"
module add matlab/r2023a
echo "job=$SLURM_JOB_ID host=$(hostname) date=$(date -Is)"
awk '/Cpus_allowed_list/ {print}' /proc/self/status
srun --exact --cpus-per-task=4 --hint=nomultithread --cpu-bind=cores matlab -batch "maxNumCompThreads(4); test_n8; test_new_protocol; test_original_hvp_n8"
srun --exact --cpus-per-task=4 --hint=nomultithread --cpu-bind=cores matlab -batch "maxNumCompThreads(4); run_historical_route('$TASK_ROOT/results/search')"
