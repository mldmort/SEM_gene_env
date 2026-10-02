#!/bin/bash
# Run one pipeline script inside a SLURM job: run_step.sh <script.py> [args...]
# Resources and log paths are set by the caller (see run_data_assembly.sh).
set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
# Inside a job, BASH_SOURCE points at SLURM's spool copy; GXE_REPO is passed in by the submitter.
REPO=${GXE_REPO:-$REPO}
source "$REPO/config/slurm.env"
export BEDTOOLS PYTHONPATH="$REPO${PYTHONPATH:+:$PYTHONPATH}"

echo "Host:  $(hostname)"
echo "Start: $(date)"
echo "Run:   $*"
echo
"$GXE_PYTHON" -u "$REPO/$1" "${@:2}"
echo
echo "End:   $(date)"
