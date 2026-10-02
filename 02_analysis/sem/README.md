# Parental age → ASD structural equation model

Separates the parental age → ASD association into a direct path, an indirect causal path
through damaging de novo mutations, and non-causal paths through shared genetics
(inherited rare variants, EA PGS, ASD PGS).

```
parent_age ~ r_p·RV + e_p·EDU + p_p·PGS + birth order
EDU        ~ r_e·RV + p_e·PGS
DNM        ~ γ·parent_age
ASD        ~ c·parent_age + δ·DNM + r_a·RV + e_a·EDU + p_a·PGS + birth order   (probit)
```

Every equation also has sex, ancestry PCs 1–10, child / mother / father coverage,
mother / father CRAM mapped %, and (SPARK WES) genotype-source indicators. Defined effects:
indirect `γ·δ`, confounding through each genetic node, and their totals.

| Script | Does |
|---|---|
| `01_prepare_sem_data.py` | Analysis-ready table per cohort (EUR cases and control siblings in the DNM callset) |
| `02_sem_lavaan.R` | Fits each cohort × parent (father, mother) × RV definition, all children and controls only; WLSMV, family-clustered SEs; fixed-effect meta |
| `03_sem_report.py` | `results/sem/sem_report.md`: per-cohort and meta tables, RV definitions compared |
| `04_sem_figures.py` | Path diagrams and forest plot in `results/sem/figures/` |

Settings are under `analysis.sem` in `config/config.yaml`: cohorts, ancestry, DNM and RV
definitions, sequencing covariates, birth-order coding (`linear` or `categorical`).

Run everything (needs the `sem` R environment, `environment_sem.yml`):

```bash
bash slurm/run_sem.sh
```
