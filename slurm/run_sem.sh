#!/bin/bash
# Prepare SEM inputs (Python) and fit the parental-age SEM (R / lavaan) as one SLURM job.
#   bash slurm/run_sem.sh                    # submits; logs in logs/sem/
#   DEPEND=123:456 bash slurm/run_sem.sh     # start after those jobs succeed
set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$REPO/config/slurm.env"
mkdir -p "$REPO/logs/sem"

sbatch --parsable --job-name=gxe_sem --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" \
    --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=32G --time=02:00:00 \
    --output="$REPO/logs/sem/sem_%j.out" --error="$REPO/logs/sem/sem_%j.err" \
    --export=ALL,GXE_REPO="$REPO" ${DEPEND:+--dependency=afterok:$DEPEND} --wrap="
set -euo pipefail
cd '$REPO'
export PYTHONPATH='$REPO' OPENBLAS_NUM_THREADS=1
echo \"Start: \$(date)  Host: \$(hostname)\"
'$GXE_PYTHON' -u 02_analysis/sem/01_prepare_sem_data.py
'$GXE_RSCRIPT' 02_analysis/sem/02_sem_lavaan.R '$REPO'
'$GXE_PYTHON' -u 02_analysis/sem/03_sem_report.py
'$GXE_PYTHON' -u 02_analysis/sem/04_sem_figures.py
echo \"End: \$(date)\""
