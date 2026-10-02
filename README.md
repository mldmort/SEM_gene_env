# GxE_analysis

Gene × environment analyses in SPARK (WES, WGS) and SSC: rare inherited variant
burden, de novo mutations, polygenic scores, parental age and environmental
exposures (air pollution, area deprivation, income, education, parental smoking).

The pipeline has two stages:

| Stage | Folder | What it does |
|---|---|---|
| 1. Data assembly | `01_data_assembly/` | Builds one table per cohort, one row per sequenced individual, holding every burden count and phenotype |
| 2. Analysis | `02_analysis/` | Models and figures, reading only the stage 1 tables |

## Layout

```
config/config.yaml      all input paths, cohort definitions and thresholds
config/slurm.env        cluster account, partition, python and bedtools paths
gxe/                    shared Python helpers (config, I/O, phenotype codings)
01_data_assembly/       steps 01-08 (see its README)
02_analysis/            analysis scripts
slurm/                  launchers
docs/                   coding notes
data/, logs/            outputs, not tracked by git
```

## Running

Requirements: Python 3.11+ with pandas, polars, pyarrow, pyyaml (see
`environment.yml`) and bedtools. On Expanse the `bioinf2` micromamba env has them.

```bash
bash slurm/run_data_assembly.sh
```

This submits every step for all three cohorts with SLURM dependencies. Pass
cohort names to limit it, or set `STEPS` to rerun some steps:

```bash
STEPS="05 08" bash slurm/run_data_assembly.sh spark_wgs
```

Final tables land in `data/assembled/gxe_<cohort>.{parquet,tsv.gz}` with a
column inventory in `gxe_<cohort>.columns.tsv`.

## Moving to another machine

Edit `config/config.yaml` (input paths) and `config/slurm.env` (cluster). No
script holds a path of its own. `GXE_CONFIG=/path/to/other.yaml` points the
scripts at a different config.
