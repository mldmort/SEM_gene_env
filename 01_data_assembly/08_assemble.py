#!/usr/bin/env python3
"""
Assemble the analysis table for one cohort: one row per sequenced individual.

Columns, in order:
  manifest      sample_id, source, FID, IID, PAT, MAT, SEX, case_status, ancestry,
                PGS (z-scored within the cohort), proj_PCs (raw PC1..PC10 dropped,
                see config `assembly`); extra PGS from config `inputs.extra_pgs`
                (asd_loSSC), each sample scored from its own genotype source
  role          child (PAT or MAT set), father / mother (listed as a child's PAT / MAT), other
  burden        rv_lof_t1..t3, rv_miss_t1..t4 (DNM-consistent tiers) and the original
                lof_tier1, lof_tier2, mis_both, mis_am_only, mis_mpc_only
                (rare inherited variants; 0 when the sample carries none)
  de novo       in_dnm_callset, then 30 dnm_* counts (missing when not in the callset)
  QC            coverage_mean, cram_*
  phenotypes    SPARK (step 05) or SSC (step 06) table
  parents       mother_<col> / father_<col> for config `assembly.parent_lookup_columns`,
                taken from the parent's own row in this cohort

SPARK individuals on the withdrawn list are dropped.

Outputs: data/assembled/gxe_<cohort>.parquet, gxe_<cohort>.tsv.gz,
         gxe_<cohort>.columns.tsv (column, dtype, non-missing count)

Usage: python 08_assemble.py --cohort spark_wgs
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from gxe import (assembled_path, cohort_manifest, cohort_names, intermediate_path,
                 left_merge, load_config, log, write_table)
from gxe.pgs import add_extra_pgs

BURDEN = ["rv_lof_t1", "rv_lof_t2", "rv_lof_t3", "rv_miss_t1", "rv_miss_t2", "rv_miss_t3", "rv_miss_t4",
          "lof_tier1", "lof_tier2", "mis_both", "mis_am_only", "mis_mpc_only"]


def assign_role(df):
    pat = df["PAT"].fillna("0")
    mat = df["MAT"].fillna("0")
    role = pd.Series("other", index=df.index)
    role[df["sample_id"].isin(set(mat))] = "mother"
    role[df["sample_id"].isin(set(pat))] = "father"
    role[(pat != "0") | (mat != "0")] = "child"
    return role


def add_parent_columns(df, cols):
    keyed = df.set_index("sample_id")
    for side, id_col in (("mother", "MAT"), ("father", "PAT")):
        ids = df[id_col].fillna("0")
        for c in cols:
            df[f"{side}_{c}"] = ids.map(keyed[c])
    return df


def main():
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--cohort", required=True, choices=cohort_names(cfg))
    args = ap.parse_args()
    cohort = args.cohort
    study = cfg["cohorts"][cohort]["study"]

    df = cohort_manifest(cohort, cfg=cfg)
    df = df.drop(columns=[c for c in cfg["assembly"]["drop_manifest_columns"] if c in df.columns])
    df.insert(df.columns.get_loc("SEX") + 1, "role", assign_role(df))
    df = add_extra_pgs(df, cfg["inputs"].get("extra_pgs", {}))
    for c in cfg["assembly"]["pgs_columns"]:
        x = pd.to_numeric(df[c], errors="coerce")
        df[c] = (x - x.mean()) / x.std()
    log(f"{cohort}: {len(df):,} manifest samples; roles: {df['role'].value_counts().to_dict()}")

    burden = pd.read_parquet(intermediate_path("inherited_burden", f"inherited_burden_{cohort}.parquet"))
    df = left_merge(df, burden)
    df[BURDEN] = df[BURDEN].fillna(0).astype(int)

    dnm = pd.read_parquet(intermediate_path("dnm", f"dnm_counts_{cohort}.parquet"))
    df["in_dnm_callset"] = df["sample_id"].isin(set(dnm["sample_id"]))
    df = left_merge(df, dnm)
    log(f"  in DNM callset: {df['in_dnm_callset'].sum():,} "
        f"(children: {(df['in_dnm_callset'] & (df['role'] == 'child')).sum():,})")

    df = left_merge(df, pd.read_parquet(intermediate_path("qc", f"sequencing_qc_{cohort}.parquet")))
    log(f"  with coverage: {df['coverage_mean'].notna().sum():,}")

    pheno = pd.read_parquet(intermediate_path("phenotypes", f"{study}_phenotypes.parquet"))
    df = left_merge(df, pheno)
    log(f"  with phenotype record: {df['sample_id'].isin(set(pheno['sample_id'])).sum():,}")

    if "spark_withdrawn" in df.columns:
        withdrawn = df["spark_withdrawn"].astype("boolean").fillna(False).astype(bool)
        log(f"  dropping {withdrawn.sum():,} withdrawn individuals")
        df = df[~withdrawn].reset_index(drop=True)

    df = add_parent_columns(df, cfg["assembly"]["parent_lookup_columns"])

    out = assembled_path(f"gxe_{cohort}.parquet")
    write_table(df, out, tsv=True)
    cols = pd.DataFrame({"column": df.columns, "dtype": df.dtypes.astype(str).values,
                         "n_non_missing": df.notna().sum().values})
    cols.to_csv(out.with_suffix(".columns.tsv"), sep="\t", index=False)
    log(f"  wrote {out.with_suffix('.columns.tsv')}")


if __name__ == "__main__":
    main()
