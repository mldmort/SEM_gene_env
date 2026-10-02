# 01_data_assembly

Builds `data/assembled/gxe_<cohort>.parquet` for `spark_wes`, `spark_wgs` and `ssc`:
one row per sequenced individual (children and parents) with genetics, de novo
counts, sequencing QC and phenotypes.

| Step | Script | Per | Output (under `data/intermediate/`) |
|---|---|---|---|
| 01 | `01_inherited_burden.py` | cohort | `inherited_burden/inherited_burden_<cohort>.parquet` |
| 02 | `02_annotate_dnms.py` | cohort | `dnm/dnm_annotated_<cohort>.parquet` |
| 03 | `03_filter_dnms.py` | cohort | `dnm/dnm_filtered_<cohort>.parquet` |
| 04 | `04_dnm_counts.py` | cohort | `dnm/dnm_counts_<cohort>.parquet` |
| 05 | `05_spark_phenotypes.py` | once | `phenotypes/spark_phenotypes.parquet` |
| 06 | `06_ssc_phenotypes.py` | once | `phenotypes/ssc_phenotypes.parquet` |
| 07 | `07_sequencing_qc.py` | cohort | `qc/sequencing_qc_<cohort>.parquet` |
| 08 | `08_assemble.py` | cohort | `data/assembled/gxe_<cohort>.{parquet,tsv.gz,columns.tsv}` |

## Definitions

**Inherited burden (01).** Rare inherited variants per sample (`carrier == "1"`), after
removing the child's de novo candidate sites. Two sets:
- DNM-consistent tiers, computed with the same code and sources as the de novo tiers
  (`gxe/annotation.py`): `rv_lof_t1-t3` (LOFTEE HC x GeneBayes s_het) and `rv_miss_t1-t4`
  (number of ClinPred / AlphaMissense / popEVE / MPC dbNSFP MANE rankscores passing).
- The original metrics: `lof_tier1` (LoF HC, LOEUF < 0.25), `lof_tier2` (0.25 <= LOEUF < 0.5),
  `mis_both` (missense, MPC > 1.75 and AlphaMissense > 0.997), `mis_am_only`, `mis_mpc_only`.

**De novo mutations (02-04).** All synthdnm candidates are annotated (VEP,
GeneBayes s_het, dbNSFP scores, gnomAD 4.1 frequencies), then kept when
`synthdnm_class == 1`, cohort AF < 0.001, gnomAD POPMAX AF < 0.001 or missing,
and outside RepeatMasker, simple repeats and segmental duplications.

- LoF tiers: HC with s_het ≥ 0.18 (t1), 0.03–0.18 (t2), < 0.03 or unknown (t3).
- Missense tiers: number of tools passing (ClinPred, AlphaMissense, popEVE, MPC); t1 = all 4, t4 = 1.

Counts are `dnm_<category>`, `dnm_snv_<category>`, `dnm_indel_<category>` for
`all, exonic, lof_t1-3, miss_t1-4, synonymous`, over autosomes only: chrX and
chrY calls are left out of the counts (`dnm.exclude_chroms_from_counts`) but
kept in the annotated and filtered parquets. `in_dnm_callset` marks children
that were in the DNM call set, so a 0 means no DNM rather than not called.

**Phenotypes (05, 06).** Per individual: birth year, parental ages at birth,
parental education (raw and harmonised 1–10, see `docs/education_coding.md`),
household income (SPARK, 1–9), ADI (SPARK; 2015, 2019, 2020, 2022), air
pollution (prenatal and at enrolment), parental smoking, core descriptive
measures (IQ, SCQ, RBS-R, Vineland, ...), and for SPARK the Adult Self-Report
substance-use scores (on the adult's row and copied to their children as
`mother_asr_*` / `father_asr_*`).

**Assembly (08).** Joins everything onto the manifest, z-scores the 19 PGS
within the cohort (the manifest values are on a raw scale, SD ~ 2e-7), adds `role`
(child / mother / father / other), drops SPARK withdrawn individuals, and adds
`mother_<col>` / `father_<col>` copies of the columns listed in
`config.yaml: assembly.parent_lookup_columns`.

## Inputs

All paths are in `config/config.yaml`. SPARK phenotypes come from release
2026-09-28, except air pollution, which is not in that release and still
comes from the 2026-03-23 folder.
