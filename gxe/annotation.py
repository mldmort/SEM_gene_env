"""Variant annotation shared by de novo (step 02) and inherited rare-variant (step 01) counts.

Both use the same sources and the same tier rules, so the LoF / missense tiers mean the
same thing for de novo and inherited variants:
  GeneBayes s_het      per-gene posterior mean, joined on Ensembl gene id
  dbNSFP MANE-select   ClinPred / AlphaMissense / popEVE / MPC rankscores, joined on
                       (#CHROM, POS, REF, ALT_specific)
  lof_t1..t3           LOFTEE HC x s_het (t1 >= shet_hi, t2 in [shet_lo, shet_hi), t3 below or missing)
  miss_t1..t4          missense with 4 / 3 / 2 / 1 of the rankscores at or above their thresholds
"""
from pathlib import Path

import polars as pl

KEY = ["#CHROM", "POS", "REF", "ALT_specific"]
TIER_COLS = ["lof_t1", "lof_t2", "lof_t3", "miss_t1", "miss_t2", "miss_t3", "miss_t4"]


def add_genebayes(lf, path):
    gb = pl.scan_csv(path, separator="\t").select(
        pl.col("ensg").alias("Gene"),
        pl.col("post_mean").cast(pl.Float64, strict=False).alias("genebayes_post_mean"),
    )
    return lf.join(gb, on="Gene", how="left")


def add_dbnsfp_scores(lf, scores_dir, chrom, score_cols):
    src = Path(scores_dir) / f"{chrom}.parquet"
    if not src.exists():
        return lf.with_columns(pl.lit(None, dtype=pl.Float64).alias(c) for c in score_cols)
    db = pl.scan_parquet(src).select(
        pl.lit(chrom).alias("#CHROM"),
        pl.col("pos(1-based)").cast(pl.Int64, strict=False).alias("POS"),
        pl.col("ref").cast(pl.Utf8).alias("REF"),
        pl.col("alt").cast(pl.Utf8).alias("ALT_specific"),
        *(pl.col(c).cast(pl.Float64, strict=False) for c in score_cols),
    ).unique(subset=KEY, keep="first")
    return lf.join(db, on=KEY, how="left")


def add_tiers(lf, dnm_cfg):
    """Needs columns LoF, Consequence, genebayes_post_mean and the rankscores."""
    is_hc = pl.col("LoF") == "HC"
    shet = pl.col("genebayes_post_mean")
    hi, lo = dnm_cfg["shet_hi"], dnm_cfg["shet_lo"]
    is_miss = pl.col("Consequence").str.contains("missense_variant")
    n_pass = sum(
        (pl.col(c).fill_null(-1.0) >= t).cast(pl.Int32)
        for c, t in dnm_cfg["missense_rankscore_thresholds"].items()
    )
    return lf.with_columns(
        (is_hc & shet.is_not_null() & (shet >= hi)).fill_null(False).alias("lof_t1"),
        (is_hc & shet.is_not_null() & (shet >= lo) & (shet < hi)).fill_null(False).alias("lof_t2"),
        (is_hc & (shet.is_null() | (shet < lo))).fill_null(False).alias("lof_t3"),
        (is_miss & (n_pass == 4)).fill_null(False).alias("miss_t1"),
        (is_miss & (n_pass == 3)).fill_null(False).alias("miss_t2"),
        (is_miss & (n_pass == 2)).fill_null(False).alias("miss_t3"),
        (is_miss & (n_pass == 1)).fill_null(False).alias("miss_t4"),
    )
