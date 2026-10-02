#!/usr/bin/env python3
"""
Per-sample counts of rare inherited variants (carrier == "1"), for one cohort.

Two sets of metrics.

DNM-consistent tiers (same sources and rules as the de novo tiers, gxe/annotation.py):
  rv_lof_t1..t3   LOFTEE HC x GeneBayes s_het (t1 >= 0.18, t2 0.03-0.18, t3 < 0.03 or missing)
  rv_miss_t1..t4  missense with 4 / 3 / 2 / 1 of ClinPred, AlphaMissense, popEVE, MPC
                  rankscores (dbNSFP MANE select) at or above their thresholds

Original metrics (thresholds in config `inherited_burden`), kept for comparison:
  lof_tier1     LoF == "HC" & LOEUF < 0.25
  lof_tier2     LoF == "HC" & 0.25 <= LOEUF < 0.5
  mis_both      missense & MPC > 1.75 & AM > 0.997
  mis_am_only   missense & AM > 0.997 & not MPC > 1.75
  mis_mpc_only  missense & MPC > 1.75 & not AM > 0.997

A child's de novo candidate sites (every variant in the synthdnm scored callset for
that child) are removed first when config `inherited_burden.exclude_dnm_candidates`
is true, so the counts are inherited variants only and do not overlap the DNM counts.

Output: data/intermediate/inherited_burden/inherited_burden_<cohort>.parquet
  sample_id + all metrics, one row per sample with >= 1 carrier row.
  Step 08 zero-fills the remaining cohort samples.

Usage: python 01_inherited_burden.py --cohort spark_wes
"""
import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import polars as pl

from gxe import cohort_names, intermediate_path, load_config, log
from gxe.annotation import TIER_COLS, add_dbnsfp_scores, add_genebayes, add_tiers

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
SITE = ["SAMPLE", "POS", "REF", "ALT_specific"]
READ_COLS = ["SAMPLE", "POS", "REF", "ALT_specific", "carrier", "LoF", "Consequence",
             "Gene", "MPC_max", "AM_max", "gnomad_lof_oe_ci_upper_min"]
OLD_METRICS = ["lof_tier1", "lof_tier2", "mis_both", "mis_am_only", "mis_mpc_only"]
NEW_METRICS = [f"rv_{t}" for t in TIER_COLS]
METRICS = NEW_METRICS + OLD_METRICS


def old_metric_exprs(t):
    loeuf = pl.col("gnomad_lof_oe_ci_upper_min").cast(pl.Float64, strict=False)
    mpc = (pl.col("MPC_max").cast(pl.Float64, strict=False) > t["mpc_min"]).fill_null(False)
    am = (pl.col("AM_max").cast(pl.Float64, strict=False) > t["am_min"]).fill_null(False)
    hc = (pl.col("LoF") == "HC").fill_null(False)
    miss = pl.col("Consequence").str.contains("missense_variant").fill_null(False)
    return {
        "lof_tier1": hc & (loeuf < t["loeuf_tier1_max"]).fill_null(False),
        "lof_tier2": hc & ((loeuf >= t["loeuf_tier1_max"]) & (loeuf < t["loeuf_tier2_max"])).fill_null(False),
        "mis_both": miss & mpc & am,
        "mis_am_only": miss & am & ~mpc,
        "mis_mpc_only": miss & mpc & ~am,
    }


def chrom_counts(path, chrom, cfg, dnm_scored):
    """Per-sample counts for one chromosome, plus carrier rows read and removed as DNM sites."""
    ann, dnm_cfg = cfg["inputs"]["annotation"], cfg["dnm"]
    lf = (pl.scan_parquet(path).select(READ_COLS)
            .filter(pl.col("carrier") == "1")
            .with_columns(pl.col("POS").cast(pl.Int64), pl.lit(chrom).alias("#CHROM")))
    n_rows = lf.select(pl.len()).collect().item()
    n_removed = 0
    if dnm_scored is not None:
        cand = (pl.scan_parquet(dnm_scored).filter(pl.col("#CHROM") == chrom)
                  .select(SITE).with_columns(pl.col("POS").cast(pl.Int64)).unique())
        kept = lf.join(cand, on=SITE, how="anti")
        n_removed = n_rows - kept.select(pl.len()).collect().item()
        lf = kept
    score_cols = list(dnm_cfg["missense_rankscore_thresholds"])
    lf = (lf.pipe(add_genebayes, ann["genebayes"])
            .pipe(add_dbnsfp_scores, ann["dbnsfp_scores_dir"], chrom, score_cols)
            .pipe(add_tiers, dnm_cfg))
    old = old_metric_exprs(cfg["inherited_burden"])
    counts = (lf.group_by("SAMPLE")
                .agg(*(pl.col(t).sum().alias(f"rv_{t}") for t in TIER_COLS),
                     *(e.sum().alias(k) for k, e in old.items()))
                .collect())
    return counts, n_rows, n_removed


def main():
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--cohort", required=True, choices=cohort_names(cfg))
    args = ap.parse_args()

    template = cfg["inputs"]["rare_variants"]
    exclude = cfg["inherited_burden"].get("exclude_dnm_candidates", False)
    dnm_scored = cfg["inputs"]["dnm_scored"][args.cohort] if exclude else None
    log(f"Excluding de novo candidate sites: {exclude}")
    parts, removed = [], 0
    for chrom in CHROMS:
        path = Path(template.format(cohort=args.cohort, chrom=chrom))
        if not path.exists():
            log(f"  {chrom}: not found, skipped ({path})")
            continue
        t0 = time.time()
        counts, n_rows, n_rm = chrom_counts(path, chrom, cfg, dnm_scored)
        parts.append(counts)
        removed += n_rm
        log(f"  {chrom}: {n_rows:,} carrier rows, {n_rm:,} at DNM candidate sites removed, "
            f"{counts.height:,} samples ({time.time() - t0:.0f}s)")

    out = (pl.concat(parts).group_by("SAMPLE").agg(pl.all().sum())
             .with_columns(pl.exclude("SAMPLE").cast(pl.Int64))
             .rename({"SAMPLE": "sample_id"}).select("sample_id", *METRICS).sort("sample_id"))
    log(f"\n{args.cohort}: {out.height:,} samples; {removed:,} carrier rows at DNM candidate sites removed")
    for m in METRICS:
        log(f"  {m:14s} total {out[m].sum():>10,}  mean per sample {out[m].mean():.4f}")
    path = intermediate_path("inherited_burden", f"inherited_burden_{args.cohort}.parquet")
    out.write_parquet(path)
    log(f"  wrote {path}")


if __name__ == "__main__":
    main()
