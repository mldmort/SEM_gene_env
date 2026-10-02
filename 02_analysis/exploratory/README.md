# Exploratory plots

Quick-look figures from the assembled tables. Outputs go to `results/exploratory/`.

| Script | Figure |
|---|---|
| `01_dnm_vs_parental_age_edu_pgs.py` | Passing DNMs per child vs paternal / maternal age and father's / mother's EA PGS, SPARK WGS and SSC |

DNM children are selected by `gxe.analysis.dnm_children` (in the DNM callset,
zero-count filter and a median + 5 MAD outlier cutoff from `config.yaml: analysis.dnm_children`).

Run on a compute node:

```bash
sbatch --account=ddp195 --partition=ind-shared --nodes=1 --ntasks=1 -c 2 --mem=16G -t 00:20:00 \
  --export=ALL,GXE_REPO=$PWD -o logs/exploratory_%j.out -e logs/exploratory_%j.err \
  slurm/run_step.sh 02_analysis/exploratory/01_dnm_vs_parental_age_edu_pgs.py
```
