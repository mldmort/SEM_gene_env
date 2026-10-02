#!/usr/bin/env python3
"""
Per-child de novo counts for one cohort: 10 categories x 3 variant classes = 30 columns.

  categories: all, exonic, lof_t1..t3, miss_t1..t4, synonymous
    exonic = lof t1-3 + miss t1-4 + synonymous (mutually exclusive, so a count of the union)
  classes:    dnm_<cat> (all variants), dnm_snv_<cat>, dnm_indel_<cat>
    SNV = single-base REF and ALT_specific; everything else is an indel

DNMs on the chromosomes in config `dnm.exclude_chroms_from_counts` (chrX, chrY)
are not counted.

Rows: every sample that appears as a child in the synthdnm scored callset (any
candidate, before filtering), so a child with no surviving DNM gets 0 rather
than being confused with a sample that was never called. Step 08 leaves
dnm_* missing for samples outside the callset.

Output: data/intermediate/dnm/dnm_counts_<cohort>.parquet

Usage: python 04_dnm_counts.py --cohort ssc
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import polars as pl

from gxe import cohort_names, intermediate_path, load_config, log

CATEGORIES = ["all", "exonic", "lof_t1", "lof_t2", "lof_t3",
              "miss_t1", "miss_t2", "miss_t3", "miss_t4", "synonymous"]
TIERS = ["lof_t1", "lof_t2", "lof_t3", "miss_t1", "miss_t2", "miss_t3", "miss_t4"]


def main():
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--cohort", required=True, choices=cohort_names(cfg))
    args = ap.parse_args()

    callset = (pl.scan_parquet(cfg["inputs"]["dnm_scored"][args.cohort])
                 .select(pl.col("SAMPLE").unique()).collect())
    log(f"{args.cohort}: {len(callset):,} children in the scored callset")

    df = pl.read_parquet(intermediate_path("dnm", f"dnm_filtered_{args.cohort}.parquet"))
    log(f"  filtered DNMs: {len(df):,}")
    excluded = cfg["dnm"]["exclude_chroms_from_counts"]
    on_excluded = pl.col("#CHROM").is_in(excluded)
    n_excluded = df.select(on_excluded.sum()).item()
    df = df.filter(~on_excluded)
    log(f"  dropped {n_excluded:,} on {', '.join(excluded)}; counting {len(df):,}")

    tier = {t: pl.col(t).fill_null(False).cast(pl.Boolean) for t in TIERS}
    # A bare pl.lit(True) sums to 1, not to the row count, so "all" uses a real column.
    every_row = pl.col("SAMPLE").is_not_null()
    cat = {
        "all": every_row,
        **tier,
        "synonymous": pl.col("Consequence").fill_null("").str.contains("synonymous_variant"),
    }
    cat["exonic"] = pl.any_horizontal(*tier.values(), cat["synonymous"])
    is_snv = (pl.col("REF").str.len_chars() == 1) & (pl.col("ALT_specific").str.len_chars() == 1)
    classes = {"dnm_": every_row, "dnm_snv_": is_snv, "dnm_indel_": ~is_snv}

    aggs = [(cls & cat[c]).sum().alias(f"{prefix}{c}")
            for prefix, cls in classes.items() for c in CATEGORIES]
    counts = df.group_by("SAMPLE").agg(aggs)

    out = (callset.join(counts, on="SAMPLE", how="left")
                  .with_columns(pl.exclude("SAMPLE").fill_null(0).cast(pl.Int32))
                  .rename({"SAMPLE": "sample_id"})
                  .sort("sample_id"))
    for c in CATEGORIES:
        s = out[f"dnm_{c}"]
        log(f"  dnm_{c}: mean {s.mean():.2f}, median {s.median():.0f}, max {s.max()}")
    path = intermediate_path("dnm", f"dnm_counts_{args.cohort}.parquet")
    out.write_parquet(path)
    log(f"  wrote {path} ({len(out):,} rows)")


if __name__ == "__main__":
    main()
