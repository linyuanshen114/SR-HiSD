#!/bin/bash
#SBATCH --job-name=BCCHEX_GF_SEQUENTIAL_BOTH
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --nodelist=c03
#SBATCH --output=/bicmr/home/yslin114/LB/experiments/bcc2hex_lbsd_sequential_bilateral_20261119/logs/%x_%j.out
#SBATCH --error=/bicmr/home/yslin114/LB/experiments/bcc2hex_lbsd_sequential_bilateral_20261119/logs/%x_%j.err
set -euo pipefail
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
TASK_ROOT=/bicmr/home/yslin114/LB/experiments/bcc2hex_lbsd_sequential_bilateral_20261119
SOURCE_STATUS=/bicmr/home/yslin114/LB/experiments/bcc2hex_lbsd_historical_v4_rebenchmark_20261117/results/search/status.json
SELECTION="$TASK_ROOT/results/source_selection.json"
cd "$TASK_ROOT"
echo "job=$SLURM_JOB_ID host=$(hostname) date=$(date -Is)"
awk '/Cpus_allowed_list/ {print}' /proc/self/status
python3 code/test_select_sequential_source.py
python3 code/select_sequential_source.py "$SOURCE_STATUS" "$SELECTION"
python3 - "$SELECTION" <<'PY'
import hashlib, json, sys
from pathlib import Path
s = json.loads(Path(sys.argv[1]).read_text())
for key in ("candidate", "spectrum", "gate"):
    h = hashlib.sha256()
    with open(s[key], "rb") as f:
        for b in iter(lambda: f.read(8*1024*1024), b""):
            h.update(b)
    if h.hexdigest() != s[key + "_sha256"]:
        raise SystemExit("SOURCE_HASH_MISMATCH: " + key)
print("SEQUENTIAL_SOURCE_HASH_PASS", flush=True)
PY
SOURCE_CANDIDATE=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["candidate"])' "$SELECTION")
SOURCE_SPECTRUM=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["spectrum"])' "$SELECTION")
module add matlab/r2023a
cd "$TASK_ROOT/gf_code"
srun --exact --cpus-per-task=4 --hint=nomultithread --cpu-bind=cores matlab -batch "maxNumCompThreads(4); test_bilateral_gf_n8; test_bilateral_failover"
srun --exact --cpus-per-task=4 --hint=nomultithread --cpu-bind=cores matlab -batch "maxNumCompThreads(4); run_lbsd_bilateral('$SOURCE_CANDIDATE','$SOURCE_SPECTRUM','$TASK_ROOT/results/flow')"
