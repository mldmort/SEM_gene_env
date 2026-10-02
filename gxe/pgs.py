"""Polygenic scores that are not in the manifest, read from per-source PLINK .sscore files."""
import pandas as pd

from gxe.io import log


def read_sscore(path, strip_suffix=False):
    """sample_id -> score from a PLINK2 .sscore (IID column, *_AVG score column).

    strip_suffix drops a trailing _SF###### / _SP###### from array IIDs, matching
    assemble_pgs_results_v2.py and build_asd_loSSC_table.py.
    """
    d = pd.read_csv(path, sep="\t")
    iid = next(c for c in d.columns if c.lstrip("#").upper() == "IID")
    score = next(c for c in d.columns if c.upper().endswith("_AVG"))
    ids = d[iid].astype(str)
    if strip_suffix:
        ids = ids.str.replace(r"_(SF|SP)\d+$", "", regex=True)
    out = pd.Series(pd.to_numeric(d[score], errors="coerce").values, index=ids)
    return out[~out.index.duplicated(keep="first")]


def add_extra_pgs(df, extra_cfg):
    """Add each extra PGS as a column, taking every sample's score from its own `source`."""
    for name, by_source in extra_cfg.items():
        col = pd.Series(float("nan"), index=df.index)
        for source, spec in by_source.items():
            rows = df["source"] == source
            if not rows.any():
                continue
            scores = read_sscore(spec["path"], spec.get("strip_suffix", False))
            col[rows] = df.loc[rows, "sample_id"].map(scores)
            log(f"  {name} [{source}]: {col[rows].notna().sum():,} / {rows.sum():,} scored")
        df[name] = col
    return df
