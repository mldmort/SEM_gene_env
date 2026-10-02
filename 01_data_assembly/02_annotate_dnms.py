#!/usr/bin/env python3
"""
Annotate all synthdnm-scored de novo candidates for one cohort (no filtering).

Per chromosome:
  1. cohort_af = AC / (2 * n_samples_vcf)
  2. VEP: Consequence, SYMBOL, IMPACT, Gene, LoF
  3. GeneBayes s_het posterior mean (per gene)
  4. dbNSFP missense rankscores (ClinPred, AlphaMissense, popEVE, MPC)
  5. dbNSFP gnomAD v4.1 joint AF / POPMAX AF / nhomalt / flag
  6. lof_t1..t3 (LOFTEE HC x s_het) and miss_t1..t4 (number of tools passing; t1 = all 4)

Ported from spark_gene_env/03_denovos/annotate_denovos/17_annotate_dnms.py.

Output: data/intermediate/dnm/dnm_annotated_<cohort>.parquet

Usage: python 02_annotate_dnms.py --cohort ssc
"""
import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import polars as pl

from gxe import cohort_names, intermediate_path, load_config, log
from gxe.annotation import KEY, TIER_COLS, add_dbnsfp_scores, add_genebayes, add_tiers

AF_COLS = {
    "gnomAD4.1_joint_AF": pl.Float64,
    "gnomAD4.1_joint_POPMAX_AF": pl.Float64,
    "gnomAD4.1_joint_nhomalt": pl.Int64,
    "gnomAD4.1_joint_flag": pl.Utf8,
}


def add_vep(lf, vep_dir, cohort):
    vep = pl.scan_parquet(Path(vep_dir) / cohort / "*.parquet").rename({
        "ALT": "ALT_specific", "CSQConsequence": "Consequence", "CSQSYMBOL": "SYMBOL",
        "CSQIMPACT": "IMPACT", "CSQGene": "Gene", "CSQLoF": "LoF",
    })
    return lf.join(vep, on=KEY, how="left")


def add_dbnsfp_af(lf, af_dir, chrom):
    src = Path(af_dir) / f"{chrom}.parquet"
    if not src.exists():
        return lf.with_columns(pl.lit(None, dtype=t).alias(c) for c, t in AF_COLS.items())
    af = pl.scan_parquet(src).select(
        pl.lit(chrom).alias("#CHROM"),
        pl.col("pos").cast(pl.Int64, strict=False).alias("POS"),
        pl.col("ref").cast(pl.Utf8).alias("REF"),
        pl.col("alt").cast(pl.Utf8).alias("ALT_specific"),
        *(pl.col(c).cast(t, strict=False) for c, t in AF_COLS.items()),
    )
    return lf.join(af, on=KEY, how="left")


def main():
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--cohort", required=True, choices=cohort_names(cfg))
    args = ap.parse_args()
    cohort = args.cohort

    ann = cfg["inputs"]["annotation"]
    dnm_cfg = cfg["dnm"]
    score_cols = list(dnm_cfg["missense_rankscore_thresholds"])
    src = cfg["inputs"]["dnm_scored"][cohort]
    n = cfg["cohorts"][cohort]["n_samples_vcf"]
    out = intermediate_path("dnm", f"dnm_annotated_{cohort}.parquet")
    partials_dir = out.parent / f"_partials_{cohort}"
    partials_dir.mkdir(exist_ok=True)

    log(f"Cohort {cohort}: {src}\n  n_samples_vcf = {n:,}")
    chroms = pl.scan_parquet(src).select(pl.col("#CHROM").unique()).collect()["#CHROM"].to_list()
    main_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chroms = [c for c in main_chroms if c in chroms] + sorted(c for c in chroms if c not in main_chroms)

    partials = []
    for chrom in chroms:
        part = partials_dir / f"{chrom}.parquet"
        log(f"  {chrom}")
        (pl.scan_parquet(src)
           .filter(pl.col("#CHROM") == chrom)
           .with_columns(pl.lit(cohort).alias("cohort"),
                         (pl.col("AC").cast(pl.Float64) / (2.0 * n)).alias("cohort_af"))
           .pipe(add_vep, ann["vep_dir"], cohort)
           .pipe(add_genebayes, ann["genebayes"])
           .pipe(add_dbnsfp_scores, ann["dbnsfp_scores_dir"], chrom, score_cols)
           .pipe(add_dbnsfp_af, ann["dbnsfp_af_dir"], chrom)
           .pipe(add_tiers, dnm_cfg)
           .sink_parquet(part, engine="streaming"))
        partials.append(part)

    pl.concat([pl.scan_parquet(p) for p in partials]).sink_parquet(out, engine="streaming")
    shutil.rmtree(partials_dir)

    summary = pl.scan_parquet(out).select(
        pl.len().alias("rows"),
        *(pl.col(c).cast(pl.Int64).sum() for c in TIER_COLS),
        pl.col("gnomAD4.1_joint_POPMAX_AF").is_not_null().sum().alias("with_gnomad_popmax"),
    ).collect()
    log(f"\nWrote {out}")
    for c, v in summary.row(0, named=True).items():
        log(f"  {c}: {v:,}")


if __name__ == "__main__":
    main()
