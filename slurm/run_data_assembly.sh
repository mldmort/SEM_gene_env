#!/bin/bash
# Submit the whole data-assembly stage (01_data_assembly) as SLURM jobs with dependencies.
#
#   bash slurm/run_data_assembly.sh                 # all cohorts, all steps
#   bash slurm/run_data_assembly.sh spark_wgs ssc   # selected cohorts
#   STEPS="07 08" bash slurm/run_data_assembly.sh   # selected steps (earlier outputs must exist)
#
# Per cohort: 01 burden, 02 annotate -> 03 filter -> 04 counts, 07 QC; 05/06 phenotypes once;
# 08 assemble waits for the cohort's 01, 04, 07 and its study's phenotype step.
set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$REPO/config/slurm.env"
if (( $# )); then COHORTS=("$@"); else COHORTS=(spark_wes spark_wgs ssc); fi
STEPS=${STEPS:-"01 02 03 04 05 06 07 08"}
LOG=$REPO/logs/01_data_assembly
mkdir -p "$LOG"

want() { [[ " $STEPS " == *" $1 "* ]]; }

# submit <name> <cpus> <mem> <time> <dependency job ids, space-separated or ""> <script> [args]
submit() {
    local name=$1 cpus=$2 mem=$3 time=$4 deps=$5; shift 5
    local dep_arg=()
    deps=$(echo "$deps" | xargs)
    [[ -n $deps ]] && dep_arg=(--dependency="afterok:${deps// /:}")
    sbatch --parsable --job-name="gxe_$name" --account="$SLURM_ACCOUNT" \
        --partition="$SLURM_PARTITION" --nodes=1 --ntasks=1 \
        --cpus-per-task="$cpus" --mem="$mem" --time="$time" \
        --output="$LOG/${name}_%j.out" --error="$LOG/${name}_%j.err" \
        --export=ALL,GXE_REPO="$REPO" "${dep_arg[@]}" \
        "$REPO/slurm/run_step.sh" "$@"
}

S=01_data_assembly
declare -A PHENO
want 05 && PHENO[spark]=$(submit 05_spark_pheno 4 32G 01:00:00 "" $S/05_spark_phenotypes.py) && echo "05 spark phenotypes: ${PHENO[spark]}"
want 06 && PHENO[ssc]=$(submit 06_ssc_pheno 2 16G 00:30:00 "" $S/06_ssc_phenotypes.py) && echo "06 ssc phenotypes:   ${PHENO[ssc]}"

for c in "${COHORTS[@]}"; do
    j01= j02= j03= j04= j07= j08=
    want 01 && j01=$(submit 01_burden_$c 8 64G 08:00:00 "" $S/01_inherited_burden.py --cohort "$c")
    want 02 && j02=$(submit 02_annotate_$c 8 96G 04:00:00 "" $S/02_annotate_dnms.py --cohort "$c")
    want 03 && j03=$(submit 03_filter_$c 4 48G 01:00:00 "$j02" $S/03_filter_dnms.py --cohort "$c")
    want 04 && j04=$(submit 04_counts_$c 4 32G 01:00:00 "$j03" $S/04_dnm_counts.py --cohort "$c")
    want 07 && j07=$(submit 07_qc_$c 2 16G 00:30:00 "" $S/07_sequencing_qc.py --cohort "$c")
    # Must match `study` in config/config.yaml.
    case $c in ssc) study=ssc ;; *) study=spark ;; esac
    want 08 && j08=$(submit 08_assemble_$c 4 64G 01:00:00 "$j01 $j04 $j07 ${PHENO[$study]:-}" \
                         $S/08_assemble.py --cohort "$c")
    echo "$c: 01=$j01 02=$j02 03=$j03 04=$j04 07=$j07 08=$j08"
done
