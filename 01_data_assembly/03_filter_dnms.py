#!/usr/bin/env python3
"""
Filter annotated DNMs for one cohort.

Keeps a candidate when all of these hold:
  1. synthdnm_class == 1
  2. cohort_af < max_cohort_af
  3. gnomAD4.1_joint_POPMAX_AF is missing or < max_gnomad_popmax_af
  4. its position (POS-1, POS) does not overlap RepeatMasker, simple repeats or segdups

Step 4 runs bedtools on a 4-column BED of positions only, so the annotation
columns never go through a TSV round trip.

Output: data/intermediate/dnm/dnm_filtered_<cohort>.parquet

Usage: python 03_filter_dnms.py --cohort ssc     (needs bedtools on PATH or $BEDTOOLS)
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import polars as pl

from gxe import cohort_names, intermediate_path, load_config, log


def mask_bed(spec, workdir, name):
    """Write chrom/start/end of a mask file (header skipped) to a plain 3-column BED."""
    cols = ",".join(str(c) for c in spec["cols"])
    out = Path(workdir) / f"{name}.bed"
    subprocess.run(f"tail -n +2 '{spec['path']}' | cut -f{cols} > '{out}'",
                   shell=True, check=True, executable="/bin/bash")
    return out


def outside_masks(df, masks, bedtools, workdir):
    """Row ids of `df` whose position overlaps none of the mask BEDs."""
    bed = Path(workdir) / "dnms.bed"
    (df.select(pl.col("#CHROM"), (pl.col("POS") - 1).alias("start"), pl.col("POS").alias("end"),
               pl.col("_rid"))
       .write_csv(bed, separator="\t", include_header=False))
    cmd = f"'{bedtools}' intersect -v -a '{bed}' -b '{masks[0]}'"
    for m in masks[1:]:
        cmd += f" | '{bedtools}' intersect -v -a - -b '{m}'"
    res = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True,
                         executable="/bin/bash")
    return [int(line.rsplit("\t", 1)[1]) for line in res.stdout.splitlines() if line]


def main():
    cfg = load_config()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--cohort", required=True, choices=cohort_names(cfg))
    args = ap.parse_args()

    bedtools = os.environ.get("BEDTOOLS") or shutil.which("bedtools")
    if not bedtools:
        sys.exit("bedtools not found: put it on PATH or set $BEDTOOLS")

    d = cfg["dnm"]
    src = intermediate_path("dnm", f"dnm_annotated_{args.cohort}.parquet")
    out = intermediate_path("dnm", f"dnm_filtered_{args.cohort}.parquet")

    n_in = pl.scan_parquet(src).select(pl.len()).collect().item()
    popmax = pl.col("gnomAD4.1_joint_POPMAX_AF")
    df = (pl.scan_parquet(src)
            .filter((pl.col("synthdnm_class") == 1)
                    & (pl.col("cohort_af") < d["max_cohort_af"])
                    & (popmax.is_null() | (popmax < d["max_gnomad_popmax_af"])))
            .collect()
            .with_row_index("_rid"))
    log(f"{args.cohort}: {n_in:,} candidates -> {len(df):,} after class / AF filters")

    with tempfile.TemporaryDirectory(dir=out.parent) as tmp:
        masks = [mask_bed(spec, tmp, name) for name, spec in cfg["inputs"]["region_masks"].items()]
        keep = outside_masks(df, masks, bedtools, tmp)
    df = df.filter(pl.col("_rid").is_in(keep)).drop("_rid")
    log(f"  -> {len(df):,} outside repeat / segdup regions")

    df.write_parquet(out)
    log(f"  wrote {out}")


if __name__ == "__main__":
    main()
