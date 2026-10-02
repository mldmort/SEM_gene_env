"""Loading the assembled tables for analysis."""
import pandas as pd

from gxe.config import REPO_ROOT, assembled_path, load_config
from gxe.io import log


def results_path(*parts, cfg=None):
    """Path under the results dir; parent directories are created."""
    cfg = cfg or load_config()
    out = (REPO_ROOT / cfg["analysis"]["results_dir"]).joinpath(*parts)
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def load_cohort(cohort, columns=None):
    return pd.read_parquet(assembled_path(f"gxe_{cohort}.parquet"), columns=columns)


def dnm_children(cohort, columns=None, cfg=None):
    """Children in the DNM callset after the zero-count and outlier filters in config.

    Returns (DataFrame, dict of sample counts at each filter step).
    """
    cfg = cfg or load_config()
    rules = cfg["analysis"]["dnm_children"]
    need = None if columns is None else sorted({"role", "in_dnm_callset", "dnm_all", *columns})
    df = load_cohort(cohort, need)
    flow = {"children": int((df["role"] == "child").sum())}

    df = df[(df["role"] == "child") & df["in_dnm_callset"]]
    flow["in_dnm_callset"] = len(df)

    if rules["drop_zero_dnm_all"][cohort]:
        df = df[df["dnm_all"] > 0]
        flow["dnm_all_gt_0"] = len(df)

    x, k = df["dnm_all"], rules["outlier_k"]
    if rules["outlier_method"] == "median_mad":
        cutoff = x.median() + k * 1.4826 * (x - x.median()).abs().median()
    elif rules["outlier_method"] == "mean_sd":
        cutoff = x.mean() + k * x.std()
    else:
        raise ValueError(f"unknown outlier_method {rules['outlier_method']!r}")
    df = df[df["dnm_all"] <= cutoff]
    flow[f"dnm_all_le_{cutoff:.0f}"] = len(df)

    log(f"{cohort}: " + " -> ".join(f"{k} {v:,}" for k, v in flow.items()))
    return df.reset_index(drop=True), flow
