#!/usr/bin/env python3
"""
Passing DNMs per child (dnm_all) vs paternal / maternal age at birth and the
father's / mother's educational-attainment PGS (edu_2022), for SPARK WGS and SSC.

Children: in the DNM callset, after the zero-count and outlier filters in config
`analysis.dnm_children`. Each panel uses the children with that predictor present.

Per panel: one grey dot per child (counts jittered vertically by +/-0.3 for display only;
x-axis shows the 0.5-99.5th percentiles, the fit uses all children),
means within deciles of the predictor (95% CI), and an OLS line, with the slope,
its 95% CI, p and n. Slopes are DNMs per year of age or per 1 SD of PGS.

Outputs (results/exploratory/):
  dnm_vs_parental_age_edu_pgs.{png,pdf}
  dnm_vs_parental_age_edu_pgs.tsv      slope, CI, p, R^2, n per cohort x predictor

Usage: python 01_dnm_vs_parental_age_edu_pgs.py [--ancestry EUR]
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import statsmodels.api as sm

from gxe import log
from gxe.analysis import dnm_children, results_path

COHORTS = {"spark_wgs": "SPARK WGS", "ssc": "SSC (WGS)"}
PREDICTORS = {
    "father_age_at_child_birth": ("Paternal age at birth (years)", "per year"),
    "mother_age_at_child_birth": ("Maternal age at birth (years)", "per year"),
    "father_edu_2022": ("Father's EA PGS (SD)", "per SD"),
    "mother_edu_2022": ("Mother's EA PGS (SD)", "per SD"),
}
Y = "dnm_all"

# Reference palette (dataviz skill): slot 1 blue, slot 2 orange; ink for text.
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK_2, GRID, DOTS = "#0b0b0b", "#52514e", "#e4e3df", "#a3a29c"


def ols(x, y):
    fit = sm.OLS(y, sm.add_constant(x)).fit()
    lo, hi = fit.conf_int()[1]
    return {"slope": fit.params[1], "ci_low": lo, "ci_high": hi,
            "p": fit.pvalues[1], "r2": fit.rsquared, "n": int(fit.nobs),
            "intercept": fit.params[0]}


def decile_means(x, y):
    bins = pd.qcut(x, 10, duplicates="drop")
    g = pd.DataFrame({"x": x, "y": y}).groupby(bins, observed=True)
    out = g.agg(x=("x", "mean"), y=("y", "mean"), sd=("y", "std"), n=("y", "size"))
    out["ci"] = 1.96 * out["sd"] / np.sqrt(out["n"])
    return out


def fmt_p(p):
    return f"p = {p:.2g}" if p >= 1e-300 else "p < 1e-300"


def panel(ax, x, y, label, unit, rng):
    st = ols(x, y)
    jitter = rng.uniform(-0.3, 0.3, len(y))
    alpha = min(0.35, 1500 / len(x))   # fainter as the panel gets denser
    ax.scatter(x, y + jitter, s=3, color=DOTS, alpha=alpha, linewidths=0, rasterized=True)
    dm = decile_means(x, y)
    ax.errorbar(dm["x"], dm["y"], yerr=dm["ci"], fmt="o", ms=5, color=BLUE,
                mec="white", mew=1, elinewidth=1.5, capsize=0, zorder=3, label="Decile mean (95% CI)")
    # Show the central 99% of x; a few extreme values (e.g. paternal age > 70) would
    # otherwise stretch the axis. The fit still uses every child.
    x_lo, x_hi = np.percentile(x, [0.5, 99.5])
    pad = 0.03 * (x_hi - x_lo)
    ax.set_xlim(x_lo - pad, x_hi + pad)
    xs = np.linspace(x_lo - pad, x_hi + pad, 50)
    ax.plot(xs, st["intercept"] + st["slope"] * xs, color=ORANGE, lw=2, zorder=4, label="OLS fit")
    ax.text(0.03, 0.97,
            f"{st['slope']:+.3f} DNMs {unit}\n95% CI {st['ci_low']:+.3f} to {st['ci_high']:+.3f}\n"
            f"{fmt_p(st['p'])}, R² = {st['r2']:.3f}\nn = {st['n']:,}",
            transform=ax.transAxes, va="top", ha="left", fontsize=8, color=INK,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GRID, alpha=0.9))
    ax.set_xlabel(label, color=INK_2, fontsize=9)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    st["n_outside_shown_x"] = int(((x < x_lo - pad) | (x > x_hi + pad)).sum())
    return st


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--ancestry", default=None,
                    help="Restrict to one predicted_ancestry (e.g. EUR); default all")
    args = ap.parse_args()

    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(len(COHORTS), len(PREDICTORS), figsize=(16, 8.2),
                             constrained_layout=True)
    rows = []
    for r, (cohort, title) in enumerate(COHORTS.items()):
        df, flow = dnm_children(cohort, columns=["predicted_ancestry", *PREDICTORS])
        if args.ancestry:
            df = df[df["predicted_ancestry"] == args.ancestry]
            log(f"  {args.ancestry} only: {len(df):,}")
        # Keep the y-axis on the bulk of the data; the fit uses every child.
        y_top = np.percentile(df[Y], 99.5) * 1.1
        for c, (pred, (label, unit)) in enumerate(PREDICTORS.items()):
            ax = axes[r, c]
            d = df[[pred, Y]].dropna()
            st = panel(ax, d[pred].to_numpy(float), d[Y].to_numpy(float), label, unit, rng)
            ax.set_ylim(0, y_top)
            rows.append({"cohort": cohort, "predictor": pred,
                         "ancestry": args.ancestry or "all", **st})
            log(f"  {pred}: slope {st['slope']:+.4f} (p {st['p']:.2g}, n {st['n']:,})")
        axes[r, 0].set_ylabel(f"{title}\npassing DNMs per child", color=INK, fontsize=10)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", frameon=False, fontsize=9, ncol=2)
    anc = f", {args.ancestry} ancestry" if args.ancestry else ", all ancestries"
    fig.suptitle(f"Passing de novo mutations per child vs parental age and parental EA PGS{anc}",
                 x=0.01, ha="left", fontsize=12, color=INK)

    stem = "dnm_vs_parental_age_edu_pgs" + (f"_{args.ancestry}" if args.ancestry else "")
    for ext in ("png", "pdf"):
        out = results_path("exploratory", f"{stem}.{ext}")
        fig.savefig(out, dpi=200)
        log(f"wrote {out}")
    tsv = results_path("exploratory", f"{stem}.tsv")
    pd.DataFrame(rows).drop(columns="intercept").to_csv(tsv, sep="\t", index=False)
    log(f"wrote {tsv}")


if __name__ == "__main__":
    main()
