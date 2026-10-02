#!/usr/bin/env python3
"""
Per-sample sequencing QC for one cohort.

  coverage_mean                 mosdepth mean over the target ("total_region" row)
  cram_total                    samtools flagstat total reads
  cram_mapped_pct, cram_primary_mapped_pct, cram_properly_paired_pct
  cram_primary_duplicates_pct   100 * primary_duplicates / primary

Output: data/intermediate/qc/sequencing_qc_<cohort>.parquet

Usage: python 07_sequencing_qc.py --cohort spark_wgs
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from gxe import cohort_names, intermediate_path, load_config, log, write_table


def load_coverage(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["chrom"] == "total_region"]
    df = df.assign(coverage_mean=pd.to_numeric(df["mean"], errors="coerce"))
    return df.rename(columns={"SAMPLE": "sample_id"})[["sample_id", "coverage_mean"]].drop_duplicates("sample_id")


def load_cram_stats(path):
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).rename(columns={
        "spid": "sample_id", "total": "cram_total", "mapped_pct": "cram_mapped_pct",
        "primary_mapped_pct": "cram_primary_mapped_pct",
        "properly_paired_pct": "cram_properly_paired_pct",
    })
    cols = ["cram_total", "cram_mapped_pct", "cram_primary_mapped_pct", "cram_properly_paired_pct"]
    for c in cols + ["primary", "primary_duplicates"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["cram_primary_duplicates_pct"] = 100.0 * df["primary_duplicates"] / df["primary"]
    return df[["sample_id"] + cols + ["cram_primary_duplicates_pct"]].drop_duplicates("sample_id")


def main():
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--cohort", required=True, choices=cohort_names(cfg))
    args = ap.parse_args()

    cov = load_coverage(cfg["inputs"]["coverage"][args.cohort])
    cram = load_cram_stats(cfg["inputs"]["cram_stats"][args.cohort])
    log(f"{args.cohort}: coverage for {len(cov):,} samples, cram stats for {len(cram):,}")
    df = cov.merge(cram, on="sample_id", how="outer")
    write_table(df, intermediate_path("qc", f"sequencing_qc_{args.cohort}.parquet"))


if __name__ == "__main__":
    main()
