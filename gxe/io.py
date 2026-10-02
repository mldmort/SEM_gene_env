"""Logging, manifest loading and table writing."""
import sys

import pandas as pd

from gxe.config import load_config


def log(msg=""):
    print(msg, file=sys.stderr, flush=True)


def read_manifest(cfg=None):
    """Full sample manifest as strings, one row per sample_id.

    The manifest can list a sample under more than one source; the first row is kept.
    """
    cfg = cfg or load_config()
    df = pd.read_csv(cfg["inputs"]["manifest"], sep="\t", dtype=str, low_memory=False)
    n = len(df)
    df = df.drop_duplicates(subset="sample_id", keep="first").reset_index(drop=True)
    if len(df) < n:
        log(f"  manifest: dropped {n - len(df):,} duplicate sample_id rows")
    return df


def cohort_manifest(cohort, manifest=None, cfg=None):
    """Manifest rows whose `source` belongs to `cohort`."""
    cfg = cfg or load_config()
    manifest = read_manifest(cfg) if manifest is None else manifest
    sources = cfg["cohorts"][cohort]["sources"]
    return manifest[manifest["source"].isin(sources)].reset_index(drop=True)


def left_merge(left, right, on="sample_id"):
    """Left join that refuses overlapping column names or a non-unique right key."""
    overlap = (set(left.columns) & set(right.columns)) - {on}
    if overlap:
        raise ValueError(f"columns present on both sides of the merge: {sorted(overlap)}")
    return left.merge(right, on=on, how="left", validate="many_to_one")


def write_table(df, path, tsv=False):
    """Write parquet (and optionally a gzipped TSV next to it)."""
    df.to_parquet(path, index=False)
    log(f"  wrote {path} ({len(df):,} rows x {df.shape[1]} cols)")
    if tsv:
        tsv_path = path.with_suffix(".tsv.gz")
        df.to_csv(tsv_path, sep="\t", index=False, compression="gzip")
        log(f"  wrote {tsv_path}")
